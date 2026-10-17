//! Exact F2 homology projections on ordered, weighted chain windows.
//!
//! This domain implements the feasible construction from `homology-operator`
//! (MIT, revision `a11bf50b8c8320651e9ac015566778d957af9c17`). It has no Python
//! runtime dependency. `A: C[k] -> C[k-1]`, `D: C[k+1] -> C[k]`, and AD=0.
//! The factorized action is P=(I+DU)(I+GA), with generalized inverses G/U;
//! L=I+P is an F2 operator, not a real Hodge Laplacian.
//!
//! Feasibility, exact binary topology and selected integer mass are distinct
//! from global optimization. Opt-in affine searches minimize one class, filling
//! or finite pair volume; they do not minimize projection stretch. Costs are
//! strictly positive u64 values in caller-declared units;
//! checked overflow is an error. Rational and floating cost domains, solver
//! registries and Python snapshot formats are outside this initial integration.
//! Source labels are caller assertions, not authentication of a dataset.
mod matrix;
pub use matrix::F2Matrix;
mod family;
pub use family::{Feature, OperatorFamily, StageInterval, WeightPolicy};

use crate::complex::FilteredComplex;
use crate::execution::{Execution, WorkBudget};
use crate::{Error, Result};
use matrix::{Inverse, check_vector, invalid, pack, unpack};
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::sync::Arc;

/// An owned validated three-degree chain window over F2.
///
/// Bases retain input order and owner-local IDs. `source` names the source
/// snapshot/namespace; `unit` declares positive abstract coordinate costs.
/// A content-equal window does not authenticate that source assertion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChainWindow<Id = usize> {
    degree: usize,
    boundaries: [F2Matrix; 2],
    bases: [Vec<Id>; 3],
    weights: Vec<u64>,
    source: String,
    unit: String,
}

impl<Id: Clone + Eq + Hash> ChainWindow<Id> {
    /// Validate explicit A/D shapes, ordered distinct bases, AD=0 and positive
    /// integer costs. Empty spaces are valid; ordinary H0 has C[-1]=0.
    pub fn new(
        degree: usize,
        boundaries: [F2Matrix; 2],
        bases: [Vec<Id>; 3],
        weights: Vec<u64>,
        source: &str,
        unit: &str,
    ) -> Result<Self> {
        Self::new_work(
            degree,
            boundaries,
            bases,
            weights,
            source,
            unit,
            &mut WorkBudget::unlimited(),
        )
    }
    fn new_work<const C: bool>(
        degree: usize,
        boundaries: [F2Matrix; 2],
        bases: [Vec<Id>; 3],
        weights: Vec<u64>,
        source: &str,
        unit: &str,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Self> {
        budget.check()?;
        for basis in &bases {
            budget.step_by(basis.len())?;
        }
        budget.step_by(weights.len())?;
        let [a, d] = &boundaries;
        if a.n != d.m || [a.m, a.n, d.n] != bases.each_ref().map(|b| b.len()) {
            return Err(invalid("boundary shapes and ordered bases disagree"));
        }
        if degree == 0 && a.m != 0 {
            return Err(invalid("ordinary H0 requires an empty previous space"));
        }
        if bases
            .iter()
            .any(|b| b.iter().collect::<HashSet<_>>().len() != b.len())
        {
            return Err(invalid("basis IDs must be unique within each degree"));
        }
        if weights.len() != a.n || weights.contains(&0) {
            return Err(invalid(
                "one positive integer cost per current coordinate is required",
            ));
        }
        if source.is_empty() || unit.is_empty() {
            return Err(invalid("source snapshot and cost unit must be declared"));
        }
        for j in 0..d.n {
            budget.step()?;
            if a.apply_with(&d.column(j)?, budget)?.iter().any(|&x| x != 0) {
                return Err(invalid("A D must be zero over F2"));
            }
        }
        Ok(Self {
            degree,
            boundaries,
            bases,
            weights,
            source: source.into(),
            unit: unit.into(),
        })
    }

    /// Extract a window from the existing frozen filtered-cell contract.
    ///
    /// `through` selects cells with value <= scale, with ties together; `None`
    /// selects the whole supplied source. Stable IDs and degree-restricted source
    /// traversal order are preserved; signed incidences are reduced modulo two,
    /// including repeated terms. A source skeleton remains the supplied skeleton:
    /// this does not certify missing cofaces, source coverage or integer validity.
    /// Costs are supplied separately from filtration values in current-basis order.
    pub fn from_complex<C: FilteredComplex<CellId = Id>>(
        complex: &C,
        degree: usize,
        through: Option<f64>,
        weights: Vec<u64>,
        source: &str,
        unit: &str,
    ) -> Result<Self>
    where
        Id: Copy,
    {
        Self::from_complex_work(
            complex,
            degree,
            through,
            weights,
            source,
            unit,
            &mut WorkBudget::unlimited(),
        )
    }
    /// As [`Self::from_complex`], with one budget for traversal, incidence
    /// conversion and window validation. Accessors and allocations cannot be
    /// interrupted internally. Failure returns no partial window.
    pub fn from_complex_with<C: FilteredComplex<CellId = Id>>(
        complex: &C,
        degree: usize,
        through: Option<f64>,
        weights: Vec<u64>,
        source: &str,
        unit: &str,
        execution: &Execution<'_>,
    ) -> Result<Self>
    where
        Id: Copy,
    {
        Self::from_complex_work(
            complex,
            degree,
            through,
            weights,
            source,
            unit,
            &mut WorkBudget::new(execution)?,
        )
    }
    fn from_complex_work<C: FilteredComplex<CellId = Id>, const CONTROLLED: bool>(
        complex: &C,
        degree: usize,
        through: Option<f64>,
        weights: Vec<u64>,
        source: &str,
        unit: &str,
        budget: &mut WorkBudget<'_, CONTROLLED>,
    ) -> Result<Self>
    where
        Id: Copy,
    {
        budget.check()?;
        if through.is_some_and(|x| !x.is_finite()) {
            return Err(invalid("stage scale must be finite"));
        }
        let next = degree.checked_add(1).ok_or(Error::SizeOverflow {
            operation: "chain degree",
        })?;
        let mut bases: [Vec<Id>; 3] = std::array::from_fn(|_| Vec::new());
        let mut known = HashMap::new();
        let mut previous_value = f64::NEG_INFINITY;
        for cell in complex.cells() {
            budget.step()?;
            let dim = complex.dimension(cell);
            let value = complex.filtration_value(cell);
            if !value.is_finite() || value < previous_value || known.contains_key(&cell) {
                return Err(invalid(
                    "source traversal must have distinct IDs and finite ordered values",
                ));
            }
            for (face, coefficient) in complex.boundary(cell) {
                budget.step()?;
                if coefficient == 0 {
                    continue;
                }
                let Some(&(face_dim, _)) = known.get(&face) else {
                    return Err(invalid("boundary must precede its coface"));
                };
                if dim == 0 || face_dim != dim - 1 {
                    return Err(invalid("boundary has the wrong degree"));
                }
            }
            known.insert(cell, (dim, value));
            previous_value = value;
            if through.is_some_and(|scale| value > scale) {
                continue;
            }
            if degree > 0 && dim == degree - 1 {
                bases[0].push(cell);
            }
            if dim == degree {
                bases[1].push(cell);
            }
            if dim == next {
                bases[2].push(cell);
            }
        }
        let mut a = F2Matrix::zero(bases[0].len(), bases[1].len())?;
        let mut d = F2Matrix::zero(bases[1].len(), bases[2].len())?;
        for (matrix, rows, columns) in [
            (&mut a, &bases[0], &bases[1]),
            (&mut d, &bases[1], &bases[2]),
        ] {
            budget.step_by(rows.len())?;
            let positions: HashMap<_, _> =
                rows.iter().enumerate().map(|(i, &id)| (id, i)).collect();
            for (j, &cell) in columns.iter().enumerate() {
                budget.step()?;
                for (face, coefficient) in complex.boundary(cell) {
                    budget.step()?;
                    if coefficient.rem_euclid(2) == 0 {
                        continue;
                    }
                    let Some(&i) = positions.get(&face) else {
                        return Err(invalid("selected stage omits a boundary cell"));
                    };
                    matrix.toggle(i, j)?;
                }
            }
        }
        Self::new_work(degree, [a, d], bases, weights, source, unit, budget)
    }
    /// Homology degree of this window.
    pub fn degree(&self) -> usize {
        self.degree
    }
    /// A and D in the declared ordered bases.
    pub fn boundaries(&self) -> &[F2Matrix; 2] {
        &self.boundaries
    }
    /// Previous, current and next degree bases, in that order.
    pub fn bases(&self) -> &[Vec<Id>; 3] {
        &self.bases
    }
    /// Positive current-coordinate costs, in current-basis order.
    pub fn weights(&self) -> &[u64] {
        &self.weights
    }
    /// Caller-declared source namespace/snapshot.
    pub fn source(&self) -> &str {
        &self.source
    }
    /// Caller-declared abstract cost unit.
    pub fn unit(&self) -> &str {
        &self.unit
    }
}

#[derive(Clone, Debug)]
enum Action {
    Factorized { g: Inverse, u: Inverse },
    Explicit(F2Matrix),
}

/// An owned feasible projection result and its original chain/cost context.
///
/// Construction checks projection idempotence, cycle image, boundary killing
/// and same-class preservation. Clones retain the run binding; a fresh solve
/// or import establishes a different binding even when all matrix entries agree.
#[derive(Clone, Debug)]
pub struct HomologyOperator<Id = usize> {
    window: Arc<ChainWindow<Id>>,
    action: Action,
    betti: usize,
    run: Arc<()>,
}

/// An owned representative query bound to its actual operator run and costs.
///
/// The immutable shared context keeps the result readable after dropping the
/// operator. Coordinates refer to the current ordered basis, not vertex labels.
#[derive(Clone, Debug)]
pub struct ClassQuery<Id = usize> {
    window: Arc<ChainWindow<Id>>,
    run: Arc<()>,
    representative: Vec<u8>,
    mass: u64,
}

/// A completed exact integer minimum with its original ordered cost window.
///
/// Construction is private: interrupted searches return an error, never an
/// optimal candidate. Optimality concerns the requested affine set only and
/// does not certify the projection's stretch or any other class. Coordinates
/// and masses remain readable after dropping the originating operator/family.
#[derive(Clone, Debug)]
pub struct MinimumChain<Id = usize> {
    window: Arc<ChainWindow<Id>>,
    chain: Vec<u8>,
    mass: u64,
    boundary: Vec<u8>,
    witness: Option<Vec<u8>>,
    examined: u64,
}
impl<Id> MinimumChain<Id> {
    /// Minimizing chain in the cost window's current-degree basis.
    pub fn chain(&self) -> &[u8] {
        &self.chain
    }
    /// Exact checked integer mass in the declared cost unit.
    pub fn mass(&self) -> u64 {
        self.mass
    }
    /// Boundary of the minimizing chain in the previous-degree basis.
    pub fn boundary(&self) -> &[u8] {
        &self.boundary
    }
    /// For class minimization, next-degree c with D c = input XOR chain.
    /// Filling/volume results have no class witness.
    pub fn boundary_witness(&self) -> Option<&[u8]> {
        self.witness.as_deref()
    }
    /// Original ordered bases, costs, source declaration and unit.
    pub fn window(&self) -> &ChainWindow<Id> {
        &self.window
    }
    /// Number of evaluated candidates. Zero/positive-cost shortcuts may finish
    /// before a full enumeration; a returned result always proves optimality.
    pub fn examined(&self) -> u64 {
        self.examined
    }
}
impl<Id: Clone + Eq + Hash> ClassQuery<Id> {
    /// Selected representative in the originating current basis.
    pub fn representative(&self) -> &[u8] {
        &self.representative
    }
    /// Checked exact integer mass of the selected representative.
    pub fn selected_mass(&self) -> u64 {
        self.mass
    }
    /// Support positions in the originating current ordered basis.
    pub fn support(&self) -> impl Iterator<Item = usize> + '_ {
        self.representative
            .iter()
            .enumerate()
            .filter_map(|(i, &bit)| (bit != 0).then_some(i))
    }
    /// Shared support positions of two records from this same operator run.
    pub fn shared_support(&self, other: &Self) -> Result<Vec<usize>> {
        self.same_class(other)?;
        Ok(self
            .representative
            .iter()
            .zip(&other.representative)
            .enumerate()
            .filter_map(|(i, (&a, &b))| (a != 0 && b != 0).then_some(i))
            .collect())
    }
    /// Union support positions of two records from this same operator run.
    pub fn union_support(&self, other: &Self) -> Result<Vec<usize>> {
        self.same_class(other)?;
        Ok(self
            .representative
            .iter()
            .zip(&other.representative)
            .enumerate()
            .filter_map(|(i, (&a, &b))| (a != 0 || b != 0).then_some(i))
            .collect())
    }
    /// Original window, costs, units and source assertion.
    pub fn window(&self) -> &ChainWindow<Id> {
        &self.window
    }
    /// Compare two records from the same accepted operator run. Foreign runs,
    /// projections, bases or costs are rejected rather than merged implicitly.
    pub fn same_class(&self, other: &Self) -> Result<bool> {
        if !Arc::ptr_eq(&self.run, &other.run) {
            return Err(invalid("query records belong to different operator runs"));
        }
        Ok(self.representative == other.representative)
    }
}

impl<Id: Clone + Eq + Hash> HomologyOperator<Id> {
    /// Construct and independently validate the deterministic feasible action.
    /// This establishes feasibility only; it does not compute an objective.
    pub fn feasible(window: ChainWindow<Id>) -> Result<Self> {
        Self::feasible_with(window, &Execution::default())
    }
    /// As [`Self::feasible`], with one fresh shared cooperative work budget.
    /// Allocation, sorting and source access are not internally preemptible;
    /// work units are arithmetic/coordinate visits, not bytes or time.
    pub fn feasible_with(window: ChainWindow<Id>, execution: &Execution<'_>) -> Result<Self> {
        let mut budget = WorkBudget::new(execution)?;
        Self::feasible_work(window, &mut budget)
    }
    fn feasible_work<const C: bool>(
        window: ChainWindow<Id>,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Self> {
        budget.check()?;
        let g = Inverse::new(&window.boundaries[0], budget)?;
        let u = Inverse::new(&window.boundaries[1], budget)?;
        let betti = window.weights.len() - g.rank() - u.rank();
        let op = Self {
            window: Arc::new(window),
            action: Action::Factorized { g, u },
            betti,
            run: Arc::new(()),
        };
        op.validate(budget)?;
        Ok(op)
    }
    /// Import an explicit P, rechecking feasibility against the actual window.
    /// No solver optimality assertion or external certification is trusted.
    pub fn from_projection(
        window: ChainWindow<Id>,
        projection: F2Matrix,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        let mut budget = WorkBudget::new(execution)?;
        let n = window.weights.len();
        if projection.m != n || projection.n != n {
            return Err(invalid("projection must be current-space square"));
        }
        let rank_a = window.boundaries[0].reduce(&mut budget)?.1.len();
        let rank_d = window.boundaries[1].reduce(&mut budget)?.1.len();
        let op = Self {
            window: Arc::new(window),
            action: Action::Explicit(projection),
            betti: n - rank_a - rank_d,
            run: Arc::new(()),
        };
        op.validate(&mut budget)?;
        Ok(op)
    }
    /// Input context retained by this result, without cloning buffers.
    pub fn window(&self) -> &ChainWindow<Id> {
        &self.window
    }
    /// Exact dimension of ker(A)/im(D), including a valid computed zero.
    pub fn betti(&self) -> usize {
        self.betti
    }
    /// Apply P to any binary chain, including noncycles.
    pub fn project(&self, x: &[u8]) -> Result<Vec<u8>> {
        self.project_work(x, &mut WorkBudget::unlimited())
    }
    /// As [`Self::project`], with fresh cooperative controls for the full action.
    pub fn project_with(&self, x: &[u8], execution: &Execution<'_>) -> Result<Vec<u8>> {
        self.project_work(x, &mut WorkBudget::new(execution)?)
    }
    fn project_work<const C: bool>(
        &self,
        x: &[u8],
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Vec<u8>> {
        check_vector(x, self.window.weights.len())?;
        match &self.action {
            Action::Explicit(p) => p.apply_with(x, budget),
            Action::Factorized { g, u } => {
                let [a, d] = &self.window.boundaries;
                let packed = pack(x, x.len(), budget)?;
                let ga = g.apply_packed(&a.apply_packed(&packed, budget)?, budget)?;
                let r: Vec<_> = packed.iter().zip(ga).map(|(a, b)| a ^ b).collect();
                let du = d.apply_packed(&u.apply_packed(&r, budget)?, budget)?;
                budget.step_by(r.len())?;
                Ok(unpack(
                    &r.into_iter()
                        .zip(du)
                        .map(|(a, b)| a ^ b)
                        .collect::<Vec<_>>(),
                    x.len(),
                ))
            }
        }
    }
    /// Apply L=I+P over F2. This is not a real spectral operator.
    pub fn apply_operator(&self, x: &[u8]) -> Result<Vec<u8>> {
        Ok(x.iter().zip(self.project(x)?).map(|(a, b)| a ^ b).collect())
    }
    /// Materialize P explicitly, paying n-by-n storage only when requested.
    pub fn projection_matrix(&self) -> Result<F2Matrix> {
        let n = self.window.weights.len();
        let mut p = F2Matrix::zero(n, n)?;
        for j in 0..n {
            let mut e = vec![0; n];
            e[j] = 1;
            for (i, bit) in self.project(&e)?.into_iter().enumerate() {
                if bit != 0 {
                    p.set(i, j)?;
                }
            }
        }
        Ok(p)
    }
    /// Whether a binary input is a cycle; malformed coordinates are errors.
    pub fn is_cycle(&self, x: &[u8]) -> Result<bool> {
        Ok(self.window.boundaries[0].apply(x)?.iter().all(|&x| x == 0))
    }
    /// Whether a binary input is a boundary in this actual window.
    pub fn is_boundary(&self, x: &[u8]) -> Result<bool> {
        Ok(self.window.boundaries[1].solve(x)?.is_some())
    }
    /// Find coefficients in the next-degree basis whose boundary is `x`.
    /// `None` means no preimage exists; `Some` may be an empty vector. This is
    /// a witness for class transport, not an optimal filling or certificate.
    pub fn boundary_preimage(&self, x: &[u8]) -> Result<Option<Vec<u8>>> {
        self.boundary_preimage_with(x, &Execution::default())
    }
    /// As [`Self::boundary_preimage`], budgeting reduction and the right-hand
    /// side solve together. The witness uses this result's retained next basis.
    pub fn boundary_preimage_with(
        &self,
        x: &[u8],
        execution: &Execution<'_>,
    ) -> Result<Option<Vec<u8>>> {
        let mut budget = WorkBudget::new(execution)?;
        check_vector(x, self.window.weights.len())?;
        self.window.boundaries[1]
            .prepare(true, &mut budget)?
            .solve(x, &mut budget)
    }
    /// Read a selected representative and its exact integer cost. Noncycles
    /// are rejected. Selected mass is not the shortest mass in the class.
    pub fn class_query(&self, z: &[u8]) -> Result<ClassQuery<Id>> {
        self.class_query_with(z, &Execution::default())
    }
    /// As [`Self::class_query`], budgeting cycle validation, action and mass.
    pub fn class_query_with(&self, z: &[u8], execution: &Execution<'_>) -> Result<ClassQuery<Id>> {
        self.query_work(z, &mut WorkBudget::new(execution)?)
    }
    fn query_work<const C: bool>(
        &self,
        z: &[u8],
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<ClassQuery<Id>> {
        if self.window.boundaries[0]
            .apply_with(z, budget)?
            .iter()
            .any(|&x| x != 0)
        {
            return Err(invalid("class queries require cycles"));
        }
        let representative = self.project_work(z, budget)?;
        budget.step_by(self.window.weights.len())?;
        let mass = self.mass(&representative)?;
        budget.check()?;
        Ok(ClassQuery {
            window: Arc::clone(&self.window),
            run: Arc::clone(&self.run),
            representative,
            mass,
        })
    }
    /// Batch cycle queries, sharing one budget and immutable run binding.
    /// A failed query returns no apparently successful partial batch.
    pub fn class_queries_with(
        &self,
        cycles: &[&[u8]],
        execution: &Execution<'_>,
    ) -> Result<Vec<ClassQuery<Id>>> {
        let mut budget = WorkBudget::new(execution)?;
        let mut output = Vec::new();
        output
            .try_reserve_exact(cycles.len())
            .map_err(|_| Error::AllocationFailed {
                context: "class query batch",
            })?;
        for &z in cycles {
            output.push(self.query_work(z, &mut budget)?);
        }
        budget.check()?;
        Ok(output)
    }
    /// Cost of the difference class in this P, with both operands validated as
    /// cycles. This does not claim a shortest-representative class metric.
    pub fn class_distance(&self, z: &[u8], y: &[u8]) -> Result<u64> {
        let q = self.class_queries_with(&[z, y], &Execution::default())?;
        self.mass(
            &q[0]
                .representative
                .iter()
                .zip(&q[1].representative)
                .map(|(a, b)| a ^ b)
                .collect::<Vec<_>>(),
        )
    }
    /// Minimize mass in z + im(D), independently of this projection's section.
    /// Ties use lexicographic chain coordinates, as in the pinned reference.
    /// Exponential in rank(D); all preparation, enumeration and witness checks
    /// share the supplied cooperative budget. Noncycles and overflow fail.
    pub fn minimum_class_mass(
        &self,
        z: &[u8],
        execution: &Execution<'_>,
    ) -> Result<MinimumChain<Id>> {
        self.minimum_class_work(z, true, &mut WorkBudget::new(execution)?)
    }
    fn minimum_class_work<const C: bool>(
        &self,
        z: &[u8],
        lexicographic: bool,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<MinimumChain<Id>> {
        if self.window.boundaries[0]
            .apply_with(z, budget)?
            .iter()
            .any(|&x| x != 0)
        {
            return Err(invalid("minimum class mass requires a cycle"));
        }
        let d = &self.window.boundaries[1];
        let pivots = d.reduce(budget)?.1;
        let basis = pivots
            .iter()
            .map(|&j| d.column(j))
            .collect::<Result<Vec<_>>>()?;
        let (chain, mass, examined) =
            minimum_affine(z, &basis, &self.window.weights, lexicographic, budget)?;
        budget.step_by(z.len())?;
        let residual: Vec<_> = z.iter().zip(&chain).map(|(a, b)| a ^ b).collect();
        // Re-solve against raw D, independently of the enumeration coefficients.
        let witness = d
            .prepare(true, budget)?
            .solve(&residual, budget)?
            .ok_or(invalid("minimum class witness escaped the boundary space"))?;
        if d.apply_with(&witness, budget)? != residual {
            return Err(invalid("minimum class witness failed"));
        }
        budget.check()?;
        Ok(MinimumChain {
            window: self.window.clone(),
            boundary: self.window.boundaries[0].apply_with(&chain, budget)?,
            chain,
            mass,
            witness: Some(witness),
            examined,
        })
    }
    /// Minimize upper-chain mass subject to D c=z. The upper window must have
    /// degree k+1, A=D and exactly matching ordered previous/current bases and
    /// source. Upper costs/units are explicit and independent of lower costs.
    /// `None` proves infeasibility; empty and zero fillings are valid minima.
    /// Ties use the smallest original-coordinate bitset. As in the reference,
    /// nonzero searches support kernel dimension at most 16 (65536 fillings).
    pub fn minimum_filling(
        &self,
        z: &[u8],
        upper: &ChainWindow<Id>,
        execution: &Execution<'_>,
    ) -> Result<Option<MinimumChain<Id>>> {
        let mut budget = WorkBudget::new(execution)?;
        self.window.check_upper(upper, &mut budget)?;
        if self.window.boundaries[0]
            .apply_with(z, &mut budget)?
            .iter()
            .any(|&x| x != 0)
        {
            return Err(invalid("fixed-boundary filling requires a cycle"));
        }
        let d = &self.window.boundaries[1];
        let Some(seed) = d.prepare(true, &mut budget)?.solve(z, &mut budget)? else {
            budget.check()?;
            return Ok(None);
        };
        let kernel = if z.iter().any(|&x| x != 0) {
            let k = d.kernel_with(&mut budget)?;
            if k.len() > 16 {
                return Err(invalid("exact filling requires nullity at most 16"));
            }
            k
        } else {
            Vec::new()
        };
        let (chain, mass, examined) =
            minimum_affine(&seed, &kernel, &upper.weights, false, &mut budget)?;
        let boundary = upper.boundaries[0].apply_with(&chain, &mut budget)?;
        if boundary != z {
            return Err(invalid("minimum filling failed its boundary equation"));
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
    /// Independently enumerable cycle-basis generators in original coordinates.
    /// This is ker(A), not a claim that these are independent homology classes.
    pub fn cycle_basis(&self) -> Result<Vec<Vec<u8>>> {
        self.window.boundaries[0].kernel_basis()
    }
    pub(super) fn homology_basis_work<const C: bool>(
        &self,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Vec<Vec<u8>>> {
        let columns: Vec<_> = self.window.boundaries[0]
            .kernel_with(budget)?
            .iter()
            .map(|z| self.project_work(z, budget))
            .collect::<Result<_>>()?;
        budget.step_by(self.window.weights.len().saturating_mul(columns.len()))?;
        let matrix = F2Matrix::from_columns(&columns, self.window.weights.len())?;
        let pivots = matrix.reduce(budget)?.1;
        pivots.into_iter().map(|j| matrix.column(j)).collect()
    }
    /// Exhaustively read this P's worst selected/input integer mass ratio.
    /// Returns numerator/denominator exactly, or `None` for no nonzero cycles.
    /// Exponential in dim ker(A); the caller must select a practical work limit.
    /// This is an objective readout for this P, never global optimality.
    pub fn stretch_with(&self, execution: &Execution<'_>) -> Result<Option<(u64, u64)>> {
        let mut budget = WorkBudget::new(execution)?;
        let basis = self.window.boundaries[0].kernel_with(&mut budget)?;
        if basis.is_empty() {
            return Ok(None);
        }
        let projected = basis
            .iter()
            .map(|z| self.project_work(z, &mut budget))
            .collect::<Result<Vec<_>>>()?;
        let mut coefficients = vec![0; basis.len()];
        let mut z = vec![0; self.window.weights.len()];
        let mut p = z.clone();
        let mut best = (0_u64, 1_u64);
        loop {
            budget.step()?;
            let Some(i) = coefficients.iter().position(|&x| x == 0) else {
                break;
            };
            coefficients[..i].fill(0);
            coefficients[i] = 1;
            // Binary increment flips precisely 0..=i. Preserve the existing
            // enumeration/tie order while reusing both linear actions.
            for j in 0..=i {
                budget.step_by(z.len().saturating_mul(2))?;
                for (x, y) in z.iter_mut().zip(&basis[j]) {
                    *x ^= y;
                }
                for (x, y) in p.iter_mut().zip(&projected[j]) {
                    *x ^= y;
                }
            }
            budget.step_by(z.len().saturating_mul(2))?;
            let ratio = (self.mass(&p)?, self.mass(&z)?);
            if u128::from(ratio.0) * u128::from(best.1) > u128::from(best.0) * u128::from(ratio.1) {
                best = ratio;
            }
        }
        Ok(Some(best))
    }
    fn mass(&self, x: &[u8]) -> Result<u64> {
        self.window
            .weights
            .iter()
            .zip(x)
            .filter(|(_, bit)| **bit != 0)
            .try_fold(0_u64, |sum, (&weight, _)| {
                sum.checked_add(weight).ok_or(Error::NumericalFailure {
                    context: "selected integer mass overflow",
                })
            })
    }
    fn validate<const C: bool>(&self, budget: &mut WorkBudget<'_, C>) -> Result<()> {
        let n = self.window.weights.len();
        let [a, d] = &self.window.boundaries;
        // Independently prepare D once for the full residual batch. A solution
        // of D x = z + Pz proves membership; solver inverse factors/caches are
        // not used by this check.
        let membership = d.prepare(true, budget)?;
        for j in 0..n {
            budget.step()?;
            let mut e = vec![0; n];
            e[j] = 1;
            let p = self.project_work(&e, budget)?;
            if self.project_work(&p, budget)? != p
                || a.apply_with(&p, budget)?.iter().any(|&x| x != 0)
            {
                return Err(invalid("projection is not idempotent with cycle image"));
            }
        }
        for j in 0..d.n {
            if self
                .project_work(&d.column(j)?, budget)?
                .iter()
                .any(|&x| x != 0)
            {
                return Err(invalid("projection does not kill boundaries"));
            }
        }
        for z in a.kernel_with(budget)? {
            let p = self.project_work(&z, budget)?;
            let residual: Vec<_> = z.iter().zip(p).map(|(a, b)| a ^ b).collect();
            if membership.solve(&residual, budget)?.is_none() {
                return Err(invalid("projection does not preserve cycle homology"));
            }
        }
        budget.check()
    }
}

impl<Id: Clone + Eq + Hash> ChainWindow<Id> {
    fn check_upper<const C: bool>(
        &self,
        upper: &Self,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<()> {
        budget.step_by(self.bases.iter().map(Vec::len).sum())?;
        if self.degree.checked_add(1) != Some(upper.degree)
            || self.boundaries[1] != upper.boundaries[0]
            || self.bases[1] != upper.bases[0]
            || self.bases[2] != upper.bases[1]
            || self.source != upper.source
        {
            return Err(invalid(
                "upper window must match degree, D, ordered bases and source",
            ));
        }
        Ok(())
    }
}

// Enumerate an independent affine basis without a machine-word dimension cap.
// Every binary increment updates the chain in place; all arithmetic is checked.
fn minimum_affine<const C: bool>(
    seed: &[u8],
    basis: &[Vec<u8>],
    weights: &[u64],
    lexicographic: bool,
    budget: &mut WorkBudget<'_, C>,
) -> Result<(Vec<u8>, u64, u64)> {
    let mut chain = seed.to_vec();
    let mut best = chain.clone();
    let mut best_mass = u64::MAX;
    let mut examined = 0_u64;
    let mut coefficients = vec![0_u8; basis.len()];
    loop {
        budget.step()?;
        budget.step_by(chain.len())?;
        let mass = weights
            .iter()
            .zip(&chain)
            .filter(|(_, b)| **b != 0)
            .try_fold(0_u64, |sum, (&w, _)| {
                sum.checked_add(w).ok_or(Error::NumericalFailure {
                    context: "minimum integer mass overflow",
                })
            })?;
        examined = examined.checked_add(1).ok_or(Error::SizeOverflow {
            operation: "minimum enumeration count",
        })?;
        let earlier = if lexicographic {
            chain < best
        } else {
            chain.iter().rev().cmp(best.iter().rev()).is_lt()
        };
        if mass < best_mass || mass == best_mass && earlier {
            best.clone_from(&chain);
            best_mass = mass;
        }
        // Positive costs prove the unique zero chain globally minimal.
        if best_mass == 0 {
            break;
        }
        let Some(i) = coefficients.iter().position(|&b| b == 0) else {
            break;
        };
        coefficients[..i].fill(0);
        coefficients[i] = 1;
        for vector in &basis[..=i] {
            budget.step_by(chain.len())?;
            for (x, y) in chain.iter_mut().zip(vector) {
                *x ^= y;
            }
        }
    }
    budget.check()?;
    Ok((best, best_mass, examined))
}
