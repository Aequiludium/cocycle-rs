//! Finite nested chain windows, coordinate inclusions and rank readouts.
use super::matrix::{ColumnBasis, Decomposition, invalid, zeros};
use super::{ClassQuery, F2Matrix, HomologyOperator, MinimumChain, minimum_affine};
use crate::Result;
use crate::complex::FilteredComplex;
use crate::execution::{Execution, WorkBudget};
use std::collections::HashMap;
use std::hash::Hash;
use std::sync::Arc;

/// One deterministic interval generator, bound to the actual family run.
///
/// Repeated endpoints have distinct occurrences and complete corrected histories.
/// A fresh family rejects this reference even when all its input entries agree.
/// Cloning the family preserves the binding. No cross-process serialization or
/// authentication of caller-declared source labels is implied.
#[derive(Clone, Debug)]
pub struct Feature {
    run: Arc<()>,
    birth: usize,
    death: Option<usize>,
    occurrence: usize,
    history: Vec<Vec<u8>>,
}
impl Feature {
    /// First live stage position.
    pub fn birth(&self) -> usize {
        self.birth
    }
    /// First dead stage, or unknown continuation beyond the observed family.
    pub fn death(&self) -> Option<usize> {
        self.death
    }
    /// Zero-based occurrence among equal endpoint pairs.
    pub fn occurrence(&self) -> usize {
        self.occurrence
    }
}

/// Cost compatibility policy for a finite operator family.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WeightPolicy {
    /// Shared coordinate IDs retain identical costs and units.
    Inherited,
    /// Costs may change; geometry is always read in the target window's costs.
    /// Chain compatibility alone establishes no isometry or mass bound.
    Variable,
}

/// A half-open interval in a finite sequence of observed stages.
///
/// Stage positions, including repeated scales, are distinct. `death=None`
/// means alive at the last observed stage; it is not proof of essentiality
/// beyond the finite sequence. Multiplicity is retained explicitly.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StageInterval {
    /// First live stage position.
    pub birth: usize,
    /// First dead stage position, or an unknown end beyond the observations.
    pub death: Option<usize>,
    /// Number of copies with these same stage endpoints.
    pub multiplicity: usize,
}

/// An owned finite family of feasible F2 operators on nested named bases.
///
/// Adjacent coordinate inclusions are checked in all three degrees, with both
/// chain equations. Degrees, declared source snapshots and cost units must
/// agree. Reordered bases are supported; arbitrary linear maps are excluded.
/// Empty families are rejected; empty individual stages are valid.
#[derive(Clone, Debug)]
pub struct OperatorFamily<Id = usize> {
    scales: Vec<f64>,
    operators: Vec<HomologyOperator<Id>>,
    policy: WeightPolicy,
    run: Arc<()>,
    simplex_source: Option<Arc<crate::complex::SimplicialComplex>>,
}

impl<Id: Clone + Eq + Hash> OperatorFamily<Id> {
    /// Validate a finite sequence using inherited positive integer costs.
    pub fn new(scales: Vec<f64>, operators: Vec<HomologyOperator<Id>>) -> Result<Self> {
        Self::with_weight_policy(scales, operators, WeightPolicy::Inherited)
    }
    /// Validate a finite sequence with an explicit cost policy. Nondecreasing
    /// finite scales are required; repeated scales retain separate stage IDs.
    pub fn with_weight_policy(
        scales: Vec<f64>,
        operators: Vec<HomologyOperator<Id>>,
        policy: WeightPolicy,
    ) -> Result<Self> {
        Self::new_work(scales, operators, policy, &mut WorkBudget::unlimited())
    }
    fn new_work<const C: bool>(
        scales: Vec<f64>,
        operators: Vec<HomologyOperator<Id>>,
        policy: WeightPolicy,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Self> {
        budget.step_by(scales.len())?;
        if operators.is_empty()
            || scales.len() != operators.len()
            || scales.iter().any(|x| !x.is_finite())
            || scales.windows(2).any(|w| w[0] > w[1])
        {
            return Err(invalid(
                "family needs nonempty stages and finite nondecreasing scales",
            ));
        }
        let family = Self {
            scales,
            operators,
            policy,
            run: Arc::new(()),
            simplex_source: None,
        };
        for i in 0..family.operators.len() - 1 {
            budget.step()?;
            let s = family.operators[i].window();
            let t = family.operators[i + 1].window();
            if s.degree() != t.degree() || s.source() != t.source() || s.unit() != t.unit() {
                return Err(invalid(
                    "family degree, source snapshot or cost units disagree",
                ));
            }
            let maps: Vec<_> = (0..3)
                .map(|d| {
                    budget.step_by(s.bases()[d].len())?;
                    budget.step_by(t.bases()[d].len())?;
                    inclusion(&s.bases()[d], &t.bases()[d])
                })
                .collect::<Result<_>>()?;
            for j in 0..s.weights().len() {
                budget.step()?;
                if t.boundaries()[0].column(maps[1][j])?
                    != embed(
                        &s.boundaries()[0].column(j)?,
                        &maps[0],
                        t.bases()[0].len(),
                        budget,
                    )?
                {
                    return Err(invalid("current inclusion does not commute with A"));
                }
                if policy == WeightPolicy::Inherited {
                    let target = maps[1][j];
                    if s.weights()[j] != t.weights()[target] {
                        return Err(invalid("inherited coordinate costs changed"));
                    }
                }
            }
            for j in 0..s.bases()[2].len() {
                budget.step()?;
                if embed(
                    &s.boundaries()[1].column(j)?,
                    &maps[1],
                    t.weights().len(),
                    budget,
                )? != t.boundaries()[1].column(maps[2][j])?
                {
                    return Err(invalid("next inclusion does not commute with D"));
                }
            }
        }
        budget.check()?;
        Ok(family)
    }
    /// Observed scales, retaining repeated values and stage order.
    pub fn scales(&self) -> &[f64] {
        &self.scales
    }
    /// Accepted cost policy; variable costs change target geometry only.
    pub fn weight_policy(&self) -> WeightPolicy {
        self.policy
    }
    /// Read one accepted operator without recomputation or copying buffers.
    pub fn stage(&self, i: usize) -> Result<&HomologyOperator<Id>> {
        self.operators
            .get(i)
            .ok_or(invalid("stage index out of range"))
    }
    fn pair(&self, i: usize, j: usize) -> Result<()> {
        self.stage(i)?;
        self.stage(j)?;
        if i > j {
            return Err(invalid("transport requires source stage <= target stage"));
        }
        Ok(())
    }
    /// Exact induced transport matrix in deterministic source/target homology
    /// bases, not an identification of same-shaped chain coordinates.
    pub fn transport(&self, i: usize, j: usize) -> Result<F2Matrix> {
        self.transport_work(i, j, &mut WorkBudget::unlimited())
    }
    /// As [`Self::transport`], sharing one work budget across source/target
    /// basis extraction, projection actions and coordinate solves.
    pub fn transport_with(
        &self,
        i: usize,
        j: usize,
        execution: &Execution<'_>,
    ) -> Result<F2Matrix> {
        self.transport_work(i, j, &mut WorkBudget::new(execution)?)
    }
    fn transport_work<const C: bool>(
        &self,
        i: usize,
        j: usize,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<F2Matrix> {
        budget.check()?;
        self.pair(i, j)?;
        let s = self.stage(i)?;
        let t = self.stage(j)?;
        if i == j {
            return identity(s.betti(), budget);
        }
        let source =
            F2Matrix::from_columns(&s.homology_basis_work(budget)?, s.window().weights().len())?;
        let target =
            F2Matrix::from_columns(&t.homology_basis_work(budget)?, t.window().weights().len())?
                .prepare(true, budget)?;
        self.transport_bases(i, j, &source, &target, budget)
    }
    fn transport_bases<const C: bool>(
        &self,
        i: usize,
        j: usize,
        source: &F2Matrix,
        target: &Decomposition,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<F2Matrix> {
        let s = self.stage(i)?;
        let t = self.stage(j)?;
        let indices = inclusion(&s.window().bases()[1], &t.window().bases()[1])?;
        let mut columns = Vec::new();
        for k in 0..source.ncols() {
            let z = source.column(k)?;
            let p = t.project_work(
                &embed(&z, &indices, t.window().weights().len(), budget)?,
                budget,
            )?;
            columns.push(
                target
                    .solve(&p, budget)?
                    .ok_or(invalid("transport escaped target homology image"))?,
            );
        }
        F2Matrix::from_columns(&columns, t.betti())
    }
    /// Rank of the induced finite-stage class map over F2.
    pub fn transport_rank(&self, i: usize, j: usize) -> Result<usize> {
        self.transport(i, j)?.rank()
    }
    /// Track a source cycle through the actual coordinate inclusion, reading
    /// its representative and checked mass under the target projection/costs.
    pub fn track_class(&self, z: &[u8], i: usize, j: usize) -> Result<ClassQuery<Id>> {
        self.track_class_with(z, i, j, &Execution::default())
    }
    /// As [`Self::track_class`], budgeting source-cycle validation, inclusion,
    /// target action and mass together. Raw coordinates use the source basis.
    pub fn track_class_with(
        &self,
        z: &[u8],
        i: usize,
        j: usize,
        execution: &Execution<'_>,
    ) -> Result<ClassQuery<Id>> {
        self.track_work(z, i, j, &mut WorkBudget::new(execution)?)
    }
    /// Transport an owned query from the actual source operator run. Equal
    /// shapes, source labels or projection entries do not identify foreign runs.
    pub fn track_query(
        &self,
        query: &ClassQuery<Id>,
        i: usize,
        j: usize,
    ) -> Result<ClassQuery<Id>> {
        self.track_query_with(query, i, j, &Execution::default())
    }
    /// As [`Self::track_query`], with one fresh cooperative transport budget.
    pub fn track_query_with(
        &self,
        query: &ClassQuery<Id>,
        i: usize,
        j: usize,
        execution: &Execution<'_>,
    ) -> Result<ClassQuery<Id>> {
        let mut budget = WorkBudget::new(execution)?;
        self.pair(i, j)?;
        if !std::sync::Arc::ptr_eq(&query.run, &self.stage(i)?.run) {
            return Err(invalid("query does not belong to the source stage run"));
        }
        self.track_work(query.representative(), i, j, &mut budget)
    }
    fn track_work<const C: bool>(
        &self,
        z: &[u8],
        i: usize,
        j: usize,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<ClassQuery<Id>> {
        budget.check()?;
        self.pair(i, j)?;
        let s = self.stage(i)?;
        let t = self.stage(j)?;
        if s.window().boundaries()[0]
            .apply_with(z, budget)?
            .iter()
            .any(|&x| x != 0)
        {
            return Err(invalid("tracking requires a source cycle"));
        }
        budget.step_by(s.window().bases()[1].len())?;
        budget.step_by(t.window().bases()[1].len())?;
        let e = inclusion(&s.window().bases()[1], &t.window().bases()[1])?;
        t.query_work(&embed(z, &e, t.window().weights().len(), budget)?, budget)
    }
    /// Read multiplicities from the exact rank invariant of this finite family.
    /// This returns stage intervals rather than silently creating a complete PH
    /// diagram from observations with unknown continuation.
    pub fn barcode(&self) -> Result<Vec<StageInterval>> {
        self.barcode_with(&Execution::default())
    }
    /// As [`Self::barcode`], sharing one budget across all induced transports,
    /// their rank reductions and the finite multiplicity readout.
    pub fn barcode_with(&self, execution: &Execution<'_>) -> Result<Vec<StageInterval>> {
        let mut budget = WorkBudget::new(execution)?;
        let n = self.operators.len();
        let mut ranks = F2RankTable::new(n)?;
        // Keep only two stage bases during preparation. All retained maps
        // below are homology-coordinate maps; no persistent family cache grows.
        let mut source = F2Matrix::from_columns(
            &self.operators[0].homology_basis_work(&mut budget)?,
            self.operators[0].window().weights().len(),
        )?;
        let mut adjacent = Vec::new();
        for j in 1..n {
            let target = F2Matrix::from_columns(
                &self.operators[j].homology_basis_work(&mut budget)?,
                self.operators[j].window().weights().len(),
            )?;
            let prepared = target.prepare(true, &mut budget)?;
            adjacent.push(self.transport_bases(j - 1, j, &source, &prepared, &mut budget)?);
            source = target;
        }
        for i in 0..n {
            budget.step()?;
            ranks.set(i, i, self.operators[i].betti());
            let mut action = identity(self.operators[i].betti(), &mut budget)?;
            for j in i + 1..n {
                action = adjacent[j - 1].product(&action, &mut budget)?;
                ranks.set(i, j, action.reduce(&mut budget)?.1.len());
            }
        }
        let mut intervals = Vec::new();
        for b in 0..n {
            for d in b + 1..=n {
                budget.step()?;
                let prev = if b == 0 { 0 } else { ranks.get(b - 1, d - 1) };
                let here = ranks
                    .get(b, d - 1)
                    .checked_sub(prev)
                    .ok_or(invalid("nonmonotone rank invariant"))?;
                let after = if d == n {
                    0
                } else {
                    ranks
                        .get(b, d)
                        .checked_sub(if b == 0 { 0 } else { ranks.get(b - 1, d) })
                        .ok_or(invalid("nonmonotone rank invariant"))?
                };
                let multiplicity = here
                    .checked_sub(after)
                    .ok_or(invalid("negative interval multiplicity"))?;
                if multiplicity != 0 {
                    intervals.push(StageInterval {
                        birth: b,
                        death: (d < n).then_some(d),
                        multiplicity,
                    });
                }
            }
        }
        budget.check()?;
        Ok(intervals)
    }
    /// Read interval generators from this family's adjacent induced maps.
    /// Complete histories are corrected backwards at deaths, rather than
    /// interpreting an endpoint pair as a unique cycle. Preparation, histories
    /// and independent raw-boundary verification share one execution budget.
    /// One incremental column basis is maintained per adjacent map and only two
    /// temporary homology bases are retained. Full chain histories remain owned
    /// output; this is a bounded exact query, not a streaming or memory-capped API.
    pub fn features(&self, execution: &Execution<'_>) -> Result<Vec<Feature>> {
        self.features_work(&mut WorkBudget::new(execution)?)
    }
    fn features_work<const C: bool>(&self, budget: &mut WorkBudget<'_, C>) -> Result<Vec<Feature>> {
        let mut source = F2Matrix::from_columns(
            &self.operators[0].homology_basis_work(budget)?,
            self.operators[0].window().weights().len(),
        )?;
        let unit = |n: usize, i: usize| {
            let mut v = vec![0; n];
            v[i] = 1;
            v
        };
        let mut alive: Vec<_> = (0..source.ncols())
            .map(|i| Ok((0, unit(source.ncols(), i), vec![source.column(i)?])))
            .collect::<Result<_>>()?;
        let mut ended = Vec::new();
        for stage in 1..self.operators.len() {
            let target = F2Matrix::from_columns(
                &self.operators[stage].homology_basis_work(budget)?,
                self.operators[stage].window().weights().len(),
            )?;
            let prepared = target.prepare(true, budget)?;
            let map = self.transport_bases(stage - 1, stage, &source, &prepared, budget)?;
            let mut span = ColumnBasis::new(map.nrows())?;
            let mut following: Vec<(usize, Vec<u8>, Vec<Vec<u8>>)> = Vec::new();
            for (birth, vector, mut history) in alive {
                let image = map.apply_with(&vector, budget)?;
                if let Some(coefficients) = span.insert(&image, following.len(), budget)? {
                    for (bit, (older_birth, _, older_history)) in
                        coefficients.iter().zip(&following)
                    {
                        if *bit != 0 {
                            for s in birth..stage {
                                budget.step_by(history[s - birth].len())?;
                                for (x, y) in history[s - birth]
                                    .iter_mut()
                                    .zip(&older_history[s - older_birth])
                                {
                                    *x ^= y;
                                }
                            }
                        }
                    }
                    ended.push((birth, Some(stage), history));
                } else {
                    history.push(target.apply_with(&image, budget)?);
                    following.push((birth, image, history));
                }
            }
            for i in 0..map.nrows() {
                let v = unit(map.nrows(), i);
                if span.insert(&v, following.len(), budget)?.is_none() {
                    following.push((stage, v, vec![target.column(i)?]));
                }
            }
            alive = following;
            source = target;
        }
        ended.extend(
            alive
                .into_iter()
                .map(|(birth, _, history)| (birth, None, history)),
        );
        ended.sort_by_key(|(b, d, _)| (*b, d.unwrap_or(self.operators.len())));
        let mut counts = HashMap::new();
        let mut features = Vec::new();
        for (birth, death, history) in ended {
            let occurrence = counts.entry((birth, death)).or_insert(0);
            features.push(Feature {
                run: self.run.clone(),
                birth,
                death,
                occurrence: *occurrence,
                history,
            });
            *occurrence += 1;
        }
        self.validate_features(&features, budget)?;
        Ok(features)
    }
    fn validate_features<const C: bool>(
        &self,
        features: &[Feature],
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<()> {
        // Verify against raw boundaries, rather than trusting coordinate solves.
        for (i, op) in self.operators.iter().enumerate() {
            let active = features
                .iter()
                .filter(|f| f.birth <= i && f.death.is_none_or(|d| i < d))
                .collect::<Vec<_>>();
            let d = &op.window().boundaries()[1];
            let mut columns = (0..d.ncols())
                .map(|j| d.column(j))
                .collect::<Result<Vec<_>>>()?;
            for f in &active {
                let z = &f.history[i - f.birth];
                if op.window().boundaries()[0]
                    .apply_with(z, budget)?
                    .iter()
                    .any(|&x| x != 0)
                    || op.project_work(z, budget)? != *z
                {
                    return Err(invalid(
                        "feature history is not in the accepted homology image",
                    ));
                }
                if i > f.birth {
                    let previous = &self.operators[i - 1];
                    let indices =
                        inclusion(&previous.window().bases()[1], &op.window().bases()[1])?;
                    let image = embed(
                        &f.history[i - f.birth - 1],
                        &indices,
                        op.window().weights().len(),
                        budget,
                    )?;
                    if op.project_work(&image, budget)? != *z {
                        return Err(invalid("feature history does not commute"));
                    }
                }
                columns.push(z.clone());
            }
            // Only dying histories need a raw-D solve. Rank-only stages must
            // not allocate a full row transform merely to certify independence.
            let boundary_space = d.prepare(features.iter().any(|f| f.death == Some(i)), budget)?;
            let rank_d = boundary_space.rank();
            let rank_joined = F2Matrix::from_columns(&columns, d.nrows())?
                .reduce(budget)?
                .1
                .len();
            if active.len() != op.betti() || rank_joined != rank_d + op.betti() {
                return Err(invalid("feature history is not a full quotient basis"));
            }
            for f in features.iter().filter(|f| f.death == Some(i)) {
                let previous = &self.operators[i - 1];
                let indices = inclusion(&previous.window().bases()[1], &op.window().bases()[1])?;
                let image = embed(
                    f.history.last().ok_or(invalid("missing feature history"))?,
                    &indices,
                    d.nrows(),
                    budget,
                )?;
                if boundary_space.solve(&image, budget)?.is_none() {
                    return Err(invalid("feature does not die in the raw boundary space"));
                }
            }
        }
        budget.check()
    }
    fn check_feature(&self, feature: &Feature) -> Result<()> {
        if !Arc::ptr_eq(&feature.run, &self.run) {
            return Err(invalid("feature belongs to another family run"));
        }
        Ok(())
    }
    /// Resolve a bound generator at a stage. `None` means not born or already
    /// dead; a foreign reference or invalid stage is an error. The returned
    /// query belongs to that stage's actual operator run.
    pub fn resolve_feature(
        &self,
        feature: &Feature,
        stage: usize,
        execution: &Execution<'_>,
    ) -> Result<Option<ClassQuery<Id>>> {
        let mut budget = WorkBudget::new(execution)?;
        self.check_feature(feature)?;
        let op = self.stage(stage)?;
        if stage < feature.birth || feature.death.is_some_and(|d| stage >= d) {
            return Ok(None);
        }
        Ok(Some(op.query_work(
            &feature.history[stage - feature.birth],
            &mut budget,
        )?))
    }
}

impl OperatorFamily<crate::complex::SimplexId> {
    /// Extract one stage per simplex insertion in the frozen source order.
    /// Tied scales remain distinct face-first events. Costs follow the full
    /// degree-k basis and are inherited. The owned source preserves original
    /// vertex IDs and certifies event provenance; no PH reducer supplies results.
    /// Unlike `from_stages`, this constructor can certify simplex pair queries.
    /// This eager small-scale reference extracts each insertion window, then
    /// independently constructs and retains one operator per insertion. Total
    /// storage is the sum of all stage windows/actions and can grow quadratically
    /// or worse. Work controls are cooperative, not a memory quota; callers must
    /// bound input size. Use the persistence APIs for scalable diagram-only work.
    pub fn from_simplices(
        complex: &crate::complex::SimplicialComplex,
        degree: usize,
        weights: Vec<u64>,
        source: &str,
        unit: &str,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        let mut budget = WorkBudget::new(execution)?;
        let ids = complex
            .cells()
            .filter(|&id| FilteredComplex::dimension(complex, id) == degree)
            .collect::<Vec<_>>();
        if weights.len() != ids.len() || complex.is_empty() {
            return Err(invalid(
                "nonempty simplex source and one cost per degree-k simplex required",
            ));
        }
        let mut operators = Vec::new();
        for i in 0..complex.len() {
            budget.step()?;
            let costs = ids
                .iter()
                .zip(&weights)
                .filter(|(id, _)| id.index() <= i)
                .map(|(_, &w)| w)
                .collect();
            let window = simplex_window(complex, degree, i + 1, costs, source, unit, &mut budget)?;
            operators.push(HomologyOperator::feasible_work(window, &mut budget)?);
        }
        let mut family = Self::new_work(
            complex.simplices().iter().map(|s| s.value()).collect(),
            operators,
            WeightPolicy::Inherited,
            &mut budget,
        )?;
        budget.step_by(complex.len())?;
        family.simplex_source = Some(Arc::new(complex.clone()));
        budget.check()?;
        Ok(family)
    }
    /// Certified birth/death simplex IDs for a bound feature. Algebra-only or
    /// scale-grouped families have no one-simplex provenance and are rejected.
    /// The IDs refer to this family's retained frozen source.
    pub fn feature_pairing(
        &self,
        feature: &Feature,
    ) -> Result<(crate::complex::SimplexId, Option<crate::complex::SimplexId>)> {
        self.check_feature(feature)?;
        let source = self
            .simplex_source
            .as_ref()
            .ok_or(invalid("one-simplex source required"))?;
        let degree = self.operators[0].window().degree();
        let birth = crate::complex::SimplexId(feature.birth);
        if FilteredComplex::dimension(source.as_ref(), birth) != degree {
            return Err(invalid("birth event has the wrong degree"));
        }
        let at_birth = self.stage(feature.birth)?.window();
        let index = at_birth.bases()[1]
            .iter()
            .position(|&id| id == birth)
            .ok_or(invalid("missing birth simplex"))?;
        if feature.history[0][index] == 0 {
            return Err(invalid("history does not contain its birth simplex"));
        }
        let death = if let Some(stage) = feature.death {
            let death = crate::complex::SimplexId(stage);
            if FilteredComplex::dimension(source.as_ref(), death) != degree + 1 {
                return Err(invalid("death event has the wrong degree"));
            }
            Some(death)
        } else {
            None
        };
        Ok((birth, death))
    }
    /// Original vertex lists for a raw current-degree chain at a stage.
    /// A raw chain need not be a cycle. Geometry uses original vertex IDs;
    /// Euclidean/periodic coordinates and their units remain caller data.
    pub fn locate_chain(&self, stage: usize, chain: &[u8]) -> Result<Vec<Vec<usize>>> {
        let op = self.stage(stage)?;
        super::matrix::check_vector(chain, op.window().weights().len())?;
        let source = self
            .simplex_source
            .as_ref()
            .ok_or(invalid("one-simplex source required"))?;
        chain
            .iter()
            .zip(&op.window().bases()[1])
            .filter(|(b, _)| **b != 0)
            .map(|(_, &id)| {
                Ok(source
                    .simplex(id)
                    .ok_or(invalid("missing source simplex"))?
                    .vertices()
                    .to_vec())
            })
            .collect()
    }
    /// Construct the explicitly costed upper-degree window at a finite death.
    /// Weights follow its ordered upper-simplex basis; the cost unit is separate
    /// from both lower-chain costs and filtration scales.
    pub fn filling_window(
        &self,
        feature: &Feature,
        weights: Vec<u64>,
        unit: &str,
        execution: &Execution<'_>,
    ) -> Result<super::ChainWindow<crate::complex::SimplexId>> {
        self.feature_pairing(feature)?;
        let stage = feature.death.ok_or(invalid("finite death required"))?;
        let source = self
            .simplex_source
            .as_ref()
            .ok_or(invalid("one-simplex source required"))?;
        let lower = self.stage(stage)?.window();
        simplex_window(
            source,
            lower.degree() + 1,
            stage + 1,
            weights,
            lower.source(),
            unit,
            &mut WorkBudget::new(execution)?,
        )
    }
    /// Exact H1 minimum in the birth-stage historical class, with unit edge
    /// costs. This query does not change P or the interval decomposition. Ties
    /// use the smallest bitset; boundary rank above 16 is outside replay scope.
    pub fn optimal_1_cycle(
        &self,
        feature: &Feature,
        execution: &Execution<'_>,
    ) -> Result<MinimumChain<crate::complex::SimplexId>> {
        let mut budget = WorkBudget::new(execution)?;
        self.feature_pairing(feature)?;
        let op = self.stage(feature.birth)?;
        if op.window().degree() != 1 || op.window().weights().iter().any(|&w| w != 1) {
            return Err(invalid("birth-pair cycle requires H1 and unit edge costs"));
        }
        if op.window().boundaries()[1].reduce(&mut budget)?.1.len() > 16 {
            return Err(invalid(
                "exact birth-cycle requires boundary rank at most 16",
            ));
        }
        op.minimum_class_work(&feature.history[0], false, &mut budget)
    }
    /// Exact finite H1/H2 pair volume. Contains the death simplex, uses only
    /// upper simplices strictly after birth and before death, has birth boundary
    /// coefficient one and zero boundary on intervening degree-k simplices.
    /// It need not represent the chosen history class. Costs come from `upper`.
    /// `None` means no finite death (outside the volume domain).
    pub fn optimal_volume(
        &self,
        feature: &Feature,
        upper: &super::ChainWindow<crate::complex::SimplexId>,
        execution: &Execution<'_>,
    ) -> Result<Option<MinimumChain<crate::complex::SimplexId>>> {
        self.pair_volume(feature, upper, None, &mut WorkBudget::new(execution)?)
    }
    /// Exact Definition 6 stable-volume optimization for integer-valued scales
    /// and a positive integer epsilon, both exactly representable as f64.
    /// `filtration_unit` declares scale units independently of upper-chain costs.
    /// Uses the exact cutoff birth+epsilon, includes the death simplex, and
    /// imposes no birth coefficient or history-class constraint. This result
    /// certifies optimization only, not a geometric perturbation guarantee.
    pub fn stable_volume(
        &self,
        feature: &Feature,
        upper: &super::ChainWindow<crate::complex::SimplexId>,
        epsilon: u64,
        filtration_unit: &str,
        execution: &Execution<'_>,
    ) -> Result<Option<MinimumChain<crate::complex::SimplexId>>> {
        let mut budget = WorkBudget::new(execution)?;
        self.check_feature(feature)?;
        if epsilon == 0 || epsilon > (1_u64 << 53) || filtration_unit.trim().is_empty() {
            return Err(invalid(
                "positive exact integer epsilon and filtration unit required",
            ));
        }
        budget.step_by(self.scales.len())?;
        if self
            .scales
            .iter()
            .any(|x| x.fract() != 0. || x.abs() > (1_u64 << 53) as f64)
        {
            return Err(invalid(
                "stable volume requires exact integer filtration scales within 2^53",
            ));
        }
        let birth = self.scales[feature.birth] as i64;
        let cutoff = birth
            .checked_add(epsilon as i64)
            .filter(|x| x.unsigned_abs() <= (1_u64 << 53))
            .ok_or(invalid("exact stable cutoff outside 2^53"))?;
        self.pair_volume(feature, upper, Some(cutoff), &mut budget)
    }
    fn pair_volume<const C: bool>(
        &self,
        feature: &Feature,
        upper: &super::ChainWindow<crate::complex::SimplexId>,
        cutoff: Option<i64>,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Option<MinimumChain<crate::complex::SimplexId>>> {
        self.feature_pairing(feature)?;
        let Some(stage) = feature.death else {
            return Ok(None);
        };
        let window = self.stage(stage)?.window();
        if !matches!(window.degree(), 1 | 2) {
            return Err(invalid("pair volume requires H1 or H2"));
        }
        window.check_upper(upper, budget)?;
        let d = &window.boundaries()[1];
        let after = |id: &crate::complex::SimplexId| {
            cutoff.map_or(id.index() > feature.birth, |c| {
                self.scales[id.index()] >= c as f64
            })
        };
        let allowed: Vec<_> = window.bases()[2]
            .iter()
            .enumerate()
            .filter(|(_, id)| id.index() < stage && after(id))
            .map(|(j, _)| j)
            .collect();
        let mut rows: Vec<_> = window.bases()[1]
            .iter()
            .enumerate()
            .filter(|(_, id)| id.index() < stage && after(id))
            .map(|(i, _)| i)
            .collect();
        let birth = window.bases()[1]
            .iter()
            .position(|id| id.index() == feature.birth)
            .ok_or(invalid("missing birth simplex"))?;
        let death = window.bases()[2]
            .iter()
            .position(|id| id.index() == stage)
            .ok_or(invalid("missing death simplex"))?;
        if cutoff.is_none() {
            rows.insert(0, birth);
        }
        budget.step_by(rows.len())?;
        let mut seed = vec![0; d.ncols()];
        seed[death] = 1;
        let anchor = cutoff.is_some() && rows.iter().all(|&i| !d.bit(i, death));
        let (chain, mass, examined) = if anchor {
            budget.step()?;
            (seed, upper.weights()[death], 1)
        } else {
            let constraint = F2Matrix::from_rows(
                &rows
                    .iter()
                    .map(|&i| allowed.iter().map(|&j| u8::from(d.bit(i, j))).collect())
                    .collect::<Vec<_>>(),
                allowed.len(),
            )?;
            let rhs: Vec<_> = rows
                .iter()
                .map(|&i| u8::from(d.bit(i, death)) ^ u8::from(cutoff.is_none() && i == birth))
                .collect();
            let prepared = constraint.prepare(true, budget)?;
            let free = prepared
                .solve(&rhs, budget)?
                .ok_or(invalid("certified pair has no feasible volume"))?;
            let kernel = prepared.kernel(budget)?;
            if kernel.len() > 16 {
                return Err(invalid(
                    "exact volume requires constraint nullity at most 16",
                ));
            }
            for (&j, &bit) in allowed.iter().zip(&free) {
                seed[j] = bit;
            }
            let basis = kernel
                .iter()
                .map(|v| embed(v, &allowed, d.ncols(), budget))
                .collect::<Result<Vec<_>>>()?;
            minimum_affine(&seed, &basis, upper.weights(), false, budget)?
        };
        let boundary = d.apply_with(&chain, budget)?;
        if chain[death] != 1
            || cutoff.is_none() && boundary[birth] != 1
            || rows
                .iter()
                .any(|&i| boundary[i] != u8::from(cutoff.is_none() && i == birth))
        {
            return Err(invalid("volume failed independent boundary constraints"));
        }
        budget.check()?;
        Ok(Some(MinimumChain {
            window: Arc::new(upper.clone()),
            chain,
            mass,
            boundary,
            witness: None,
            examined,
        }))
    }
    /// Extract and solve borrowed stages of one frozen simplicial source.
    ///
    /// Checks the actual [`crate::filtration::SimplicialStage::inclusion_into`]
    /// relation before extraction, rejecting equal clones from other owners or
    /// mixed certified/bare interpretations. Costs are positive integers in
    /// each stage's current-degree traversal order and inherited on shared IDs.
    /// One budget covers inclusion checks, extraction, validation and solving;
    /// accessors/allocations are not internally preemptible. The owned family
    /// outlives the stages. It certifies only the supplied finite windows, not
    /// completeness of higher dimensions or an infinite filtration.
    pub fn from_stages(
        stages: &[crate::filtration::SimplicialStage<'_>],
        degree: usize,
        weights: Vec<Vec<u64>>,
        source: &str,
        unit: &str,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        let mut budget = WorkBudget::new(execution)?;
        if stages.is_empty() || stages.len() != weights.len() {
            return Err(invalid(
                "one cost vector per nonempty stage sequence is required",
            ));
        }
        for pair in stages.windows(2) {
            budget.step()?;
            let _ = pair[0].inclusion_into(&pair[1])?;
        }
        let mut operators = Vec::new();
        operators
            .try_reserve_exact(stages.len())
            .map_err(|_| crate::Error::AllocationFailed {
                context: "stage operator family",
            })?;
        for (stage, costs) in stages.iter().zip(weights) {
            let window = super::ChainWindow::from_complex_work(
                stage,
                degree,
                None,
                costs,
                source,
                unit,
                &mut budget,
            )?;
            operators.push(HomologyOperator::feasible_work(window, &mut budget)?);
        }
        Self::new_work(
            stages.iter().map(|stage| stage.scale()).collect(),
            operators,
            WeightPolicy::Inherited,
            &mut budget,
        )
    }
}

fn simplex_window<const C: bool>(
    complex: &crate::complex::SimplicialComplex,
    degree: usize,
    count: usize,
    weights: Vec<u64>,
    source: &str,
    unit: &str,
    budget: &mut WorkBudget<'_, C>,
) -> Result<super::ChainWindow<crate::complex::SimplexId>> {
    let next = degree.checked_add(1).ok_or(crate::Error::SizeOverflow {
        operation: "chain degree",
    })?;
    let mut bases: [Vec<_>; 3] = std::array::from_fn(|_| Vec::new());
    for id in complex.cells().take(count) {
        budget.step()?;
        let dim = FilteredComplex::dimension(complex, id);
        if degree > 0 && dim == degree - 1 {
            bases[0].push(id);
        }
        if dim == degree {
            bases[1].push(id);
        }
        if dim == next {
            bases[2].push(id);
        }
    }
    let mut a = F2Matrix::zero(bases[0].len(), bases[1].len())?;
    let mut d = F2Matrix::zero(bases[1].len(), bases[2].len())?;
    for (matrix, rows, columns) in [
        (&mut a, &bases[0], &bases[1]),
        (&mut d, &bases[1], &bases[2]),
    ] {
        let positions: HashMap<_, _> = rows.iter().enumerate().map(|(i, &id)| (id, i)).collect();
        for (j, &id) in columns.iter().enumerate() {
            for (face, coefficient) in FilteredComplex::boundary(complex, id) {
                budget.step()?;
                if coefficient.rem_euclid(2) != 0 {
                    let i = *positions.get(&face).ok_or(invalid("prefix omits a face"))?;
                    matrix.toggle(i, j)?;
                }
            }
        }
    }
    super::ChainWindow::new_work(degree, [a, d], bases, weights, source, unit, budget)
}

// Small rank table; no projection or persistence storage contract is exposed.
struct F2RankTable {
    n: usize,
    values: Vec<usize>,
}
impl F2RankTable {
    fn new(n: usize) -> Result<Self> {
        let size = n.checked_mul(n).ok_or(crate::Error::SizeOverflow {
            operation: "stage rank table",
        })?;
        let mut values = Vec::new();
        values
            .try_reserve_exact(size)
            .map_err(|_| crate::Error::AllocationFailed {
                context: "stage ranks",
            })?;
        values.resize(size, 0);
        Ok(Self { n, values })
    }
    fn set(&mut self, i: usize, j: usize, value: usize) {
        self.values[i * self.n + j] = value;
    }
    fn get(&self, i: usize, j: usize) -> usize {
        self.values[i * self.n + j]
    }
}
fn inclusion<Id: Eq + Hash>(source: &[Id], target: &[Id]) -> Result<Vec<usize>> {
    let positions: HashMap<_, _> = target.iter().enumerate().map(|(i, id)| (id, i)).collect();
    source
        .iter()
        .map(|id| {
            positions
                .get(id)
                .copied()
                .ok_or(invalid("stage bases are not nested by identity"))
        })
        .collect()
}
fn embed<const C: bool>(
    x: &[u8],
    indices: &[usize],
    size: usize,
    budget: &mut WorkBudget<'_, C>,
) -> Result<Vec<u8>> {
    let mut output = zeros(size)?;
    for (&bit, &j) in x.iter().zip(indices) {
        budget.step()?;
        output[j] = bit;
    }
    Ok(output)
}
fn identity<const C: bool>(n: usize, budget: &mut WorkBudget<'_, C>) -> Result<F2Matrix> {
    let mut output = F2Matrix::zero(n, n)?;
    for j in 0..n {
        budget.step()?;
        output.set(j, j)?;
    }
    Ok(output)
}
