//! Exact F2 homology projections on ordered, weighted chain windows.
//!
//! This domain implements the feasible construction from `homology-operator`
//! (MIT, revision `3daff9bd69033d1ee9c5cba5026bfff94df8d6ac`). It has no Python
//! runtime dependency. `A: C[k] -> C[k-1]`, `D: C[k+1] -> C[k]`, and AD=0.
//! The factorized action is P=(I+DU)(I+GA), with generalized inverses G/U;
//! L=I+P is an F2 operator, not a real Hodge Laplacian.
//!
//! Feasibility, exact binary topology and selected integer mass are distinct
//! from global optimization. No minimum-stretch or shortest-class solver is
//! implied. Costs are strictly positive u64 values in caller-declared units;
//! checked overflow is an error. Rational and floating cost domains, solver
//! registries and Python snapshot formats are outside this initial integration.
//! Source labels are caller assertions, not authentication of a dataset.
mod matrix;
pub use matrix::F2Matrix;
mod family;
pub use family::{OperatorFamily, StageInterval, WeightPolicy};

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
