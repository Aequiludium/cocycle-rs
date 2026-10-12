//! Finite nested chain windows, coordinate inclusions and rank readouts.
use super::matrix::{Decomposition, invalid, zeros};
use super::{ClassQuery, F2Matrix, HomologyOperator};
use crate::Result;
use crate::execution::{Execution, WorkBudget};
use std::collections::HashMap;
use std::hash::Hash;

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
}

impl OperatorFamily<crate::complex::SimplexId> {
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
