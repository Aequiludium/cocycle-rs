//! Exact algebraic splittings; these are not real orthogonal projections.
use super::*;

/// Owned algebraic projection and chain homotopy in a PathHomology's Omega bases.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathProjection {
    field: PrimeField,
    mode: PathBoundaryMode,
    annotations: Vec<Vec<PathVector>>,
    projections: Vec<Vec<PathVector>>,
    homotopies: Vec<Vec<PathVector>>,
}
impl PathProjection {
    /// Selected field; no interpretation over another field is established.
    pub fn field(&self) -> PrimeField {
        self.field
    }
    /// Ambient boundary model.
    pub fn mode(&self) -> PathBoundaryMode {
        self.mode
    }
    /// A_n maps Omega_n coordinates to the selected H_n coordinates.
    pub fn annotations(&self, n: usize) -> Option<&[PathVector]> {
        self.annotations.get(n).map(Vec::as_slice)
    }
    /// Columns of P_n=Z_n A_n in Omega_n coordinates.
    pub fn projections(&self, n: usize) -> Option<&[PathVector]> {
        self.projections.get(n).map(Vec::as_slice)
    }
    /// h_n maps Omega_n to Omega_(n+1), including h_d in the upper layer.
    pub fn homotopies(&self, n: usize) -> Option<&[PathVector]> {
        self.homotopies.get(n).map(Vec::as_slice)
    }
}

/// Linear cohomology paired identically with a particular PathHomology basis.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathCohomology {
    field: PrimeField,
    mode: PathBoundaryMode,
    cochains: Vec<Vec<PathVector>>,
}
impl PathCohomology {
    /// Bound coefficient field.
    pub fn field(&self) -> PrimeField {
        self.field
    }
    /// Ambient boundary convention.
    pub fn mode(&self) -> PathBoundaryMode {
        self.mode
    }
    /// Cocycle columns in the dual Omega basis, paired in the original H order.
    pub fn representatives(&self, n: usize) -> Option<&[PathVector]> {
        self.cochains.get(n).map(Vec::as_slice)
    }
}

impl PathHomology {
    /// Build a checked algebraic splitting I-P=B h+h B in degrees 0 through d.
    ///
    /// Matrices use this result's Omega coordinates. This has no inner-product
    /// or harmonic meaning. The positive cap has the exact solver's per-matrix
    /// and retained-plus-temporary entry semantics, not total-memory semantics.
    ///
    /// # Errors
    /// Returns entry/allocation/size errors or a failed exact splitting identity.
    pub fn projection(&self, max_entries: usize) -> Result<PathProjection> {
        self.projection_with(max_entries, &Execution::default())
    }
    /// Build a splitting using one shared budget through all solves and identities.
    ///
    /// # Errors
    /// Returns the errors of [`Self::projection`], cancellation or work exhaustion.
    pub fn projection_with(
        &self,
        max_entries: usize,
        execution: &Execution<'_>,
    ) -> Result<PathProjection> {
        let mut budget = WorkBudget::new(execution)?;
        self.projection_budget(max_entries, &mut budget)
    }
    /// Compute linear cohomology normalized to the stored homology basis.
    ///
    /// Cochains live in Omega dual coordinates. This supplies neither cup
    /// products nor a canonical extension to the whole allowed-path dual space.
    ///
    /// # Errors
    /// Returns the errors of [`Self::projection`].
    pub fn cohomology(&self, max_entries: usize) -> Result<PathCohomology> {
        self.cohomology_with(max_entries, &Execution::default())
    }
    /// Compute dual cocycles with shared projection/transposition work controls.
    ///
    /// # Errors
    /// Returns the errors of [`Self::cohomology`], cancellation or work exhaustion.
    pub fn cohomology_with(
        &self,
        limit: usize,
        execution: &Execution<'_>,
    ) -> Result<PathCohomology> {
        let mut budget = WorkBudget::new(execution)?;
        let projection = self.projection_budget(limit, &mut budget)?;
        let mut cochains = Vec::new();
        for (n, columns) in projection.annotations.iter().enumerate() {
            budget.step()?;
            let mut layer = Vec::new();
            for _ in &self.representatives[n] {
                push(&mut layer, PathVector::new())?;
            }
            for (j, column) in columns.iter().enumerate() {
                for (&i, &value) in column {
                    budget.step()?;
                    layer[i].insert(j, value);
                }
            }
            matrix(&layer, limit)?;
            push(&mut cochains, layer)?;
        }
        Ok(PathCohomology {
            field: self.field,
            mode: self.mode,
            cochains,
        })
    }

    pub(super) fn projection_budget(
        &self,
        limit: usize,
        budget: &mut WorkBudget<'_>,
    ) -> Result<PathProjection> {
        if limit == 0 {
            return Err(Error::InvalidParameter {
                parameter: "max_entries",
                reason: "must be positive",
            });
        }
        count(self.omega_bases.len(), limit)?;
        let field = self.field;
        let (mut images, mut lifts) = (Vec::new(), Vec::new());
        for columns in &self.boundaries {
            budget.step()?;
            let mut solver = Echelon::new(field, limit);
            let (mut image, mut lift) = (Vec::new(), Vec::new());
            for (j, column) in columns.iter().enumerate() {
                if solver.insert(column, j, 0, budget)?.is_none() {
                    push(&mut image, copy_vector(column, budget)?)?;
                    push(&mut lift, BTreeMap::from([(j, 1)]))?;
                }
            }
            matrix(&image, limit)?;
            matrix(&lift, limit)?;
            push(&mut images, image)?;
            push(&mut lifts, lift)?;
        }
        let (mut annotations, mut projections, mut homotopies) =
            (Vec::new(), Vec::new(), Vec::<Vec<PathVector>>::new());
        for n in 0..self.representatives.len() {
            budget.step()?;
            let solver = Echelon::from_columns(&self.omega_bases[n], field, limit, budget)?;
            let mut cycles = Vec::new();
            for v in &self.representatives[n] {
                push(&mut cycles, solver.coordinates(v, budget)?)?;
            }
            matrix(&cycles, limit)?;
            let c = &images[n + 1];
            let mut splitting = Vec::new();
            for v in c.iter().chain(&cycles).chain(&lifts[n]) {
                push(&mut splitting, copy_vector(v, budget)?)?;
                matrix(&splitting, limit)?;
            }
            let solver = Echelon::from_columns(&splitting, field, limit, budget)?;
            let width = self.omega_bases[n].len();
            if splitting.len() != width || solver.rank() != width {
                return Err(invariant("boundary/homology/lift splitting"));
            }
            let (mut alpha, mut h_coeffs) = (Vec::new(), Vec::new());
            for j in 0..width {
                let coordinates = solver.coordinates(&BTreeMap::from([(j, 1)]), budget)?;
                let (mut a, mut h) = (PathVector::new(), PathVector::new());
                for (i, v) in coordinates {
                    budget.step()?;
                    if i < c.len() {
                        h.insert(i, v);
                    } else if i < sum(c.len(), cycles.len())? {
                        a.insert(i - c.len(), v);
                    }
                }
                push(&mut alpha, a)?;
                push(&mut h_coeffs, h)?;
                matrix(&alpha, limit)?;
                matrix(&h_coeffs, limit)?;
            }
            let projection = apply(&cycles, &alpha, field, limit, budget)?;
            let homotopy = apply(&lifts[n + 1], &h_coeffs, field, limit, budget)?;
            if apply(&projection, &projection, field, limit, budget)? != projection
                || apply(&self.boundaries[n], &projection, field, limit, budget)?
                    .iter()
                    .any(|v| !v.is_empty())
                || apply(&alpha, &self.boundaries[n + 1], field, limit, budget)?
                    .iter()
                    .any(|v| !v.is_empty())
            {
                return Err(invariant("projection and cocycle identities"));
            }
            let mut identity = Vec::new();
            for j in 0..cycles.len() {
                budget.step()?;
                push(&mut identity, BTreeMap::from([(j, 1)]))?;
            }
            if apply(&alpha, &cycles, field, limit, budget)? != identity {
                return Err(invariant("cohomology evaluation is identity"));
            }
            let mut rhs = apply(&self.boundaries[n + 1], &homotopy, field, limit, budget)?;
            if n > 0 {
                let lower = apply(
                    &homotopies[n - 1],
                    &self.boundaries[n],
                    field,
                    limit,
                    budget,
                )?;
                for (column, extra) in rhs.iter_mut().zip(lower) {
                    add_scaled(column, &extra, 1, field, 0, limit, budget)?;
                }
            }
            let mut lhs = Vec::new();
            for (j, column) in projection.iter().enumerate() {
                let mut v = BTreeMap::from([(j, 1)]);
                add_scaled(&mut v, column, field.negate(1), field, 0, limit, budget)?;
                push(&mut lhs, v)?;
            }
            matrix(&lhs, limit)?;
            matrix(&rhs, limit)?;
            if lhs != rhs {
                return Err(invariant("I-P = B h + h B"));
            }
            push(&mut annotations, alpha)?;
            push(&mut projections, projection)?;
            push(&mut homotopies, homotopy)?;
        }
        budget.check()?;
        Ok(PathProjection {
            field,
            mode: self.mode,
            annotations,
            projections,
            homotopies,
        })
    }
}
