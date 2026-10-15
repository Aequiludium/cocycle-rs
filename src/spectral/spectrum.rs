//! Nonzero signed pairs from singular vectors of a retained boundary.
use super::*;

/// Approximate nonzero Dirac eigenpairs from one retained boundary map.
///
/// Contains count largest positive singular values of `B[boundary_degree]` in
/// descending order, and their 2*count signed eigenpairs in ascending order.
/// Vectors use the entire frozen truncated Dirac basis. Negative column j is
/// Gamma times positive column 2*count-1-j, with Gamma=(-1)^absolute_degree.
/// Does not compute the kernel or claim a complete Dirac spectrum. Different
/// boundaries must be requested separately; numerical thresholds do not certify rank.
pub struct DiracBoundarySpectrum<I> {
    degrees: RangeInclusive<usize>,
    boundary_degree: usize,
    basis: Vec<I>,
    offsets: Vec<usize>,
    singular_values: Vec<f64>,
    values: Vec<f64>,
    vectors: Vec<f64>,
    residuals: Vec<f64>,
    orthogonality_error: f64,
    gram_residuals: Vec<f64>,
    gram_orthogonality_error: f64,
    iterations: usize,
    gram_actions: usize,
    peak_subspace: usize,
    options: IterativeOptions,
    zero_threshold: f64,
}
impl<I> DiracBoundarySpectrum<I> {
    /// Selected inclusive chain interval; all equations use these same endpoints.
    pub fn degrees(&self) -> RangeInclusive<usize> {
        self.degrees.clone()
    }
    /// Degree j of the retained map `B[j]:C[j]->C[j-1]`.
    pub fn boundary_degree(&self) -> usize {
        self.boundary_degree
    }
    /// Owned graded source-local domain handles.
    pub fn basis(&self) -> &[I] {
        &self.basis
    }
    /// Offsets of all selected degrees, including empty degrees and the final total.
    pub fn degree_offsets(&self) -> &[usize] {
        &self.offsets
    }
    /// Descending positive numerical singular values, each above zero_threshold.
    pub fn singular_values(&self) -> &[f64] {
        &self.singular_values
    }
    /// Ascending eigenvalues [-sigma_max,...,-sigma_min,sigma_min,...,sigma_max].
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    /// N by 2*count row-major vectors in the entire truncated Dirac basis.
    pub fn vectors(&self) -> &[f64] {
        &self.vectors
    }
    /// Actual residuals ||K v-lambda v|| from the full truncated sparse action.
    pub fn residuals(&self) -> &[f64] {
        &self.residuals
    }
    /// Measured ||V^T V-I|| for all returned signed vectors.
    pub fn orthogonality_error(&self) -> f64 {
        self.orthogonality_error
    }
    /// Normal-map residuals, paired with descending singular_values.
    pub fn gram_residuals(&self) -> &[f64] {
        &self.gram_residuals
    }
    /// Measured orthogonality error of the selected right singular Ritz vectors.
    pub fn gram_orthogonality_error(&self) -> f64 {
        self.gram_orthogonality_error
    }
    /// Block updates for the normal-map solve.
    pub fn iterations(&self) -> usize {
        self.iterations
    }
    /// Calls to matrix-free B^T B; final Dirac residual calls are separate.
    pub fn gram_applications(&self) -> usize {
        self.gram_actions
    }
    /// Largest dense Ritz dimension used for the normal map.
    pub fn peak_subspace_dimension(&self) -> usize {
        self.peak_subspace
    }
    /// Retained block-solver policy.
    pub fn options(&self) -> IterativeOptions {
        self.options
    }
    /// Explicit numerical exclusion threshold in singular-value units.
    pub fn zero_threshold(&self) -> f64 {
        self.zero_threshold
    }
}

impl<I: Copy + Eq + Hash> DiracOperator<I> {
    /// Approximate count largest positive singular modes of one retained `B[j]`.
    ///
    /// Uses matrix-free B^T B and the same block LOBPCG/Ritz backend as Hodge.
    /// From Bv=sigma*u and B^T u=sigma*v, form w+=(u,v)/sqrt(2) and
    /// w-=Gamma*w+. The chain condition makes these eigenpairs of the entire
    /// same truncated K, rather than just a two-degree surrogate. Both signed
    /// residuals are verified with that actual K before any result is returned.
    ///
    /// No full K, Hodge block or SVD matrix is materialized. Conditioning is
    /// squared by the normal map; tiny modes may fail residual/threshold checks.
    /// IterativeOptions tolerances have sigma² units on B^T B. The signed
    /// tolerance is absolute_tolerance/|sigma| + relative_tolerance*|sigma|,
    /// in Dirac eigenvalue units; this reflects normal-map conditioning.
    /// This does not compute all nonzero modes, certify rank or infer the kernel.
    /// A random block cannot certify extremal index coverage or missing multiplicity.
    /// One Execution budget spans starts, Gram actions, Ritz solve and verification.
    /// # Errors
    /// Rejects j outside the retained maps (a<j<=b), invalid count/options or
    /// threshold, workspace/output caps, cancellation/work exhaustion, allocation,
    /// arithmetic/nonconvergence, requested modes at/below zero_threshold, and
    /// signed residual or orthogonality failures. Zero modes are not padded.
    pub fn boundary_spectrum(
        &self,
        degree: usize,
        count: usize,
        zero_threshold: f64,
        options: &IterativeOptions,
        execution: &Execution<'_>,
    ) -> Result<DiracBoundarySpectrum<I>> {
        let mut budget = WorkBudget::new(execution)?;
        if degree <= *self.degrees.start() || degree > *self.degrees.end() {
            return Err(Error::InvalidParameter {
                parameter: "Dirac boundary degree",
                reason: "require a < degree <= b in the same retained interval",
            });
        }
        if !zero_threshold.is_finite() || zero_threshold < 0. {
            return Err(Error::InvalidParameter {
                parameter: "Dirac zero threshold",
                reason: "must be finite and nonnegative",
            });
        }
        let index = degree - self.degrees.start();
        let lower = self.offsets[index - 1]..self.offsets[index];
        let upper = self.offsets[index]..self.offsets[index + 1];
        let n = self.basis.len();
        let columns = count.checked_mul(2).ok_or(Error::SizeOverflow {
            operation: "Dirac eigenvector columns",
        })?;
        let entries = n.checked_mul(columns).ok_or(Error::SizeOverflow {
            operation: "Dirac eigenvector entries",
        })?;
        if entries > options.max_basis_entries {
            return Err(Error::InvalidParameter {
                parameter: "Dirac eigenvector entries",
                reason: "exceeds the explicit workspace entry limit",
            });
        }
        let gram = crate::eigensolver::lobpcg(
            upper.len(),
            count,
            true,
            options,
            &mut budget,
            |x, budget| {
                budget.step_by(lower.len())?;
                let mut bx = buffer(lower.len(), 0.)?;
                for (j, terms) in self.boundaries[upper.clone()].iter().enumerate() {
                    for &(row, c) in terms {
                        budget.step()?;
                        bx[row - lower.start] = accumulate(bx[row - lower.start], c as f64, x[j])?;
                    }
                }
                budget.step_by(upper.len())?;
                let mut result = buffer(upper.len(), 0.)?;
                for (j, terms) in self.boundaries[upper.clone()].iter().enumerate() {
                    for &(row, c) in terms {
                        budget.step()?;
                        result[j] = accumulate(result[j], c as f64, bx[row - lower.start])?;
                    }
                }
                Ok(result)
            },
        )?;
        let mut left = Vec::new();
        let mut sigmas = Vec::new();
        for j in 0..count {
            budget.step_by(lower.len())?;
            let mut u = buffer(lower.len(), 0.)?;
            for (i, terms) in self.boundaries[upper.clone()].iter().enumerate() {
                for &(row, c) in terms {
                    budget.step()?;
                    u[row - lower.start] =
                        accumulate(u[row - lower.start], c as f64, gram.vectors[i * count + j])?;
                }
            }
            budget.step_by(u.len())?;
            let sigma = u.iter().fold(0_f64, |a, &b| a.hypot(b));
            if !sigma.is_finite() || sigma <= zero_threshold || gram.values[j] <= 0. {
                return Err(Error::NumericalFailure {
                    context: "Dirac requested nonzero singular mode",
                });
            }
            for x in &mut u {
                *x /= sigma;
            }
            left.try_reserve(1).map_err(|_| allocation())?;
            left.push(u);
            sigmas.try_reserve(1).map_err(|_| allocation())?;
            sigmas.push(sigma);
        }
        let mut order: Vec<usize> = Vec::new();
        order.try_reserve_exact(count).map_err(|_| allocation())?;
        budget.step_by(count)?;
        order.extend(0..count);
        order.sort_unstable_by(|&i, &j| sigmas[j].total_cmp(&sigmas[i]));
        budget.step_by(entries)?;
        let mut vectors = buffer(entries, 0.)?;
        let mut values = buffer(columns, 0.)?;
        let mut singular_values = buffer(count, 0.)?;
        let mut gram_residuals = buffer(count, 0.)?;
        for (j, &original) in order.iter().enumerate() {
            let positive = columns - 1 - j;
            singular_values[j] = sigmas[original];
            gram_residuals[j] = gram.residuals[original];
            values[j] = -sigmas[original];
            values[positive] = sigmas[original];
            budget.step_by(lower.len() + upper.len())?;
            for row in lower.clone() {
                vectors[row * columns + positive] = left[original][row - lower.start];
            }
            for row in upper.clone() {
                vectors[row * columns + positive] =
                    gram.vectors[(row - upper.start) * count + original];
            }
            let length = (0..n).fold(0_f64, |a, row| a.hypot(vectors[row * columns + positive]));
            if !length.is_finite() || length == 0. {
                return Err(Error::NumericalFailure {
                    context: "Dirac pair normalization",
                });
            }
            budget.step_by(n)?;
            for row in 0..n {
                vectors[row * columns + positive] /= length;
            }
            for (slot, range) in self.offsets.windows(2).enumerate() {
                let sign = if (self.degrees.start() + slot).is_multiple_of(2) {
                    1.
                } else {
                    -1.
                };
                budget.step_by(range[1] - range[0])?;
                for row in range[0]..range[1] {
                    vectors[row * columns + j] = sign * vectors[row * columns + positive];
                }
            }
        }
        let mut residuals = buffer(columns, 0.)?;
        for j in 0..columns {
            budget.step_by(n)?;
            let mut v = buffer(n, 0.)?;
            for row in 0..n {
                v[row] = vectors[row * columns + j];
            }
            let mut image = self.action_budget(&v, &mut budget)?;
            budget.step_by(n)?;
            for row in 0..n {
                image[row] -= values[j] * v[row];
            }
            let residual = image.iter().fold(0_f64, |a, &b| a.hypot(b));
            let target = options.absolute_tolerance / values[j].abs()
                + options.relative_tolerance * values[j].abs();
            if !residual.is_finite() || !target.is_finite() || residual > target {
                return Err(Error::NumericalFailure {
                    context: "Dirac signed eigenpair residual",
                });
            }
            residuals[j] = residual;
        }
        let mut orthogonality_error = 0_f64;
        for i in 0..columns {
            for j in 0..columns {
                budget.step_by(n)?;
                let dot: f64 = (0..n)
                    .map(|row| vectors[row * columns + i] * vectors[row * columns + j])
                    .sum();
                orthogonality_error = orthogonality_error.hypot(dot - if i == j { 1. } else { 0. });
            }
        }
        if !orthogonality_error.is_finite() || orthogonality_error > 1e-8 {
            return Err(Error::NumericalFailure {
                context: "Dirac paired orthogonality",
            });
        }
        let mut basis = Vec::new();
        basis.try_reserve_exact(n).map_err(|_| allocation())?;
        budget.step_by(n)?;
        basis.extend_from_slice(&self.basis);
        budget.step_by(self.offsets.len())?;
        let mut offsets = buffer(self.offsets.len(), 0)?;
        offsets.copy_from_slice(&self.offsets);
        budget.check()?;
        Ok(DiracBoundarySpectrum {
            degrees: self.degrees.clone(),
            boundary_degree: degree,
            basis,
            offsets,
            singular_values,
            values,
            vectors,
            residuals,
            orthogonality_error,
            gram_residuals,
            gram_orthogonality_error: gram.orthogonality_error,
            iterations: gram.iterations,
            gram_actions: gram.actions,
            peak_subspace: gram.peak_subspace,
            options: *options,
            zero_threshold,
        })
    }
}
