//! Bounded corrected weak column reduction, independently transcribed from
//! Botnan-Dong Algorithms 1-3. Keep explicit dense transforms for this baseline;
//! compact representations and automatic dispatch are separate work.
use super::{SpectralOptions, allocation};
use crate::execution::WorkBudget;
use crate::{Error, Result};

use super::schur::dense;

/// Disjoint signed supports of ker(B2), with no dense R/E/V/Z buffers.
#[derive(Clone, Debug)]
pub(super) struct KernelProjection {
    groups: Vec<usize>,
    signs: Vec<f64>,
    sizes: Vec<usize>,
}
impl KernelProjection {
    pub(super) fn new(
        rows: &[Vec<(usize, f64)>],
        columns: usize,
        max_entries: usize,
        max_incidence: usize,
        budget: &mut WorkBudget<'_>,
    ) -> Result<Self> {
        let mut adjacency: Vec<Vec<(usize, f64)>> = buffer(columns, Vec::new(), budget)?;
        let mut forced = buffer(columns, false, budget)?;
        let mut entries = 0usize;
        for row in rows {
            budget.step_by(row.len())?;
            if row.len() > 2 || row.iter().any(|&(_, v)| v != 1. && v != -1.) {
                return Err(ineligible(
                    "outside rows must have at most two signed unit entries",
                ));
            }
            if row.len() == 1 {
                forced[row[0].0] = true;
            } else if row.len() == 2 {
                let (a, b) = (row[0], row[1]);
                entries = entries.checked_add(2).ok_or(Error::SizeOverflow {
                    operation: "compact projection adjacency",
                })?;
                if entries > max_incidence {
                    return Err(ineligible("compact adjacency exceeds incidence cap"));
                }
                for (u, v) in [(a.0, b.0), (b.0, a.0)] {
                    adjacency[u].try_reserve(1).map_err(|_| allocation())?;
                    adjacency[u].push((v, -a.1 * b.1));
                }
            }
        }
        let mut groups = buffer(columns, usize::MAX, budget)?;
        let mut signs = buffer(columns, 0., budget)?;
        let mut sizes = Vec::new();
        let mut pending = Vec::new();
        let mut support = Vec::new();
        pending
            .try_reserve_exact(columns)
            .map_err(|_| allocation())?;
        support
            .try_reserve_exact(columns)
            .map_err(|_| allocation())?;
        sizes.try_reserve_exact(columns).map_err(|_| allocation())?;
        for root in 0..columns {
            budget.step()?;
            if signs[root] != 0. {
                continue;
            }
            support.clear();
            pending.push(root);
            signs[root] = 1.;
            let mut blocked = false;
            while let Some(u) = pending.pop() {
                budget.step()?;
                support.push(u);
                blocked |= forced[u];
                for &(v, relation) in &adjacency[u] {
                    budget.step()?;
                    let wanted = signs[u] * relation;
                    if signs[v] == 0. {
                        signs[v] = wanted;
                        pending.push(v);
                    } else {
                        blocked |= signs[v] != wanted;
                    }
                }
            }
            if !blocked {
                for &i in &support {
                    budget.step()?;
                    groups[i] = sizes.len();
                }
                sizes.push(support.len());
            }
        }
        let stored = columns
            .checked_add(sizes.len())
            .ok_or(Error::SizeOverflow {
                operation: "compact projection entries",
            })?;
        if stored > max_entries {
            return Err(ineligible("compact projection exceeds entry cap"));
        }
        budget.check()?;
        Ok(Self {
            groups,
            signs,
            sizes,
        })
    }
    pub(super) fn entries(&self) -> usize {
        self.signs.len() + self.sizes.len()
    }
    pub(super) fn rank(&self) -> usize {
        self.signs.len() - self.sizes.len()
    }
    pub(super) fn copy(&self, budget: &mut WorkBudget<'_>) -> Result<Self> {
        let mut groups = buffer(self.groups.len(), 0, budget)?;
        let mut signs = buffer(self.signs.len(), 0., budget)?;
        let mut sizes = buffer(self.sizes.len(), 0, budget)?;
        groups.copy_from_slice(&self.groups);
        signs.copy_from_slice(&self.signs);
        sizes.copy_from_slice(&self.sizes);
        Ok(Self {
            groups,
            signs,
            sizes,
        })
    }
    pub(super) fn apply(&self, y: &mut [f64], budget: &mut WorkBudget<'_>) -> Result<()> {
        let mut sums = buffer(self.sizes.len(), 0., budget)?;
        for (i, &value) in y.iter().enumerate() {
            budget.step()?;
            let group = self.groups[i];
            if group != usize::MAX {
                sums[group] += self.signs[i] * value;
            }
        }
        for (sum, &size) in sums.iter_mut().zip(&self.sizes) {
            budget.step()?;
            *sum /= size as f64;
            if !sum.is_finite() {
                return Err(super::numerical("compact kernel projection"));
            }
        }
        for (i, value) in y.iter_mut().enumerate() {
            budget.step()?;
            *value = if self.groups[i] == usize::MAX {
                0.
            } else {
                self.signs[i] * sums[self.groups[i]]
            };
        }
        Ok(())
    }
    // Equivalent signed constraints, used when extending an existing kernel.
    pub(super) fn constraints(
        &self,
        budget: &mut WorkBudget<'_>,
    ) -> Result<Vec<Vec<(usize, f64)>>> {
        let mut roots = buffer(self.sizes.len(), usize::MAX, budget)?;
        let mut rows = Vec::new();
        rows.try_reserve_exact(self.groups.len())
            .map_err(|_| allocation())?;
        for (i, &group) in self.groups.iter().enumerate() {
            budget.step()?;
            let mut row = Vec::new();
            row.try_reserve_exact(2).map_err(|_| allocation())?;
            if group == usize::MAX {
                row.push((i, 1.));
            } else if roots[group] == usize::MAX {
                roots[group] = i;
            } else {
                let root = roots[group];
                row.push((root, -self.signs[root] * self.signs[i]));
                row.push((i, 1.));
            }
            if !row.is_empty() {
                rows.push(row);
            }
        }
        Ok(rows)
    }
}
fn buffer<T: Clone>(n: usize, value: T, budget: &mut WorkBudget<'_>) -> Result<Vec<T>> {
    budget.step_by(n)?;
    let mut out = Vec::new();
    out.try_reserve_exact(n).map_err(|_| allocation())?;
    out.resize(n, value);
    Ok(out)
}

// The nonzero-column index map is separate from the original zero columns.
// R, diagonal E and V are all consumed when selecting kernel directions.
pub(super) struct WeakReduction {
    pub(super) nonzero: Vec<usize>,
    pub(super) original_zero: Vec<usize>,
    pub(super) zero_columns: Vec<usize>,
    pub(super) e: Vec<f64>,
    pub(super) v: Vec<f64>,
    pub(super) r: Vec<f64>,
}

pub(super) fn weak_reduction(
    d: &[f64],
    rows: usize,
    columns: usize,
    options: &SpectralOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<WeakReduction> {
    check_shape(d, rows, columns)?;
    let mut adjacency = Vec::new();
    adjacency
        .try_reserve_exact(columns)
        .map_err(|_| allocation())?;
    let mut active = Vec::new();
    active
        .try_reserve_exact(columns)
        .map_err(|_| allocation())?;
    active.resize(columns, false);
    let mut forced = Vec::new();
    forced
        .try_reserve_exact(columns)
        .map_err(|_| allocation())?;
    forced.resize(columns, false);
    for _ in 0..columns {
        budget.step()?;
        adjacency.push(Vec::<(usize, f64)>::new());
    }
    for i in 0..rows {
        budget.step()?;
        let mut entries = [(0, 0.); 2];
        let mut count = 0;
        for j in 0..columns {
            budget.step()?;
            let value = d[i * columns + j];
            if value == 0. {
                continue;
            }
            if value != 1. && value != -1. {
                return Err(ineligible("require signed unit outside incidence"));
            }
            if count == 2 {
                return Err(ineligible("outside rows must have at most two nonzeros"));
            }
            entries[count] = (j, value);
            count += 1;
            active[j] = true;
        }
        if count == 1 {
            forced[entries[0].0] = true;
        } else if count == 2 {
            let (u, a) = entries[0];
            let (v, b) = entries[1];
            adjacency[u].try_reserve(1).map_err(|_| allocation())?;
            adjacency[v].try_reserve(1).map_err(|_| allocation())?;
            adjacency[u].push((v, -a * b));
            adjacency[v].push((u, -a * b));
        }
    }
    let mut nonzero = Vec::new();
    let mut original_zero = Vec::new();
    let mut local = Vec::new();
    local.try_reserve_exact(columns).map_err(|_| allocation())?;
    for (j, &is_active) in active.iter().enumerate() {
        budget.step()?;
        local.push(nonzero.len());
        let indices = if is_active {
            &mut nonzero
        } else {
            &mut original_zero
        };
        indices.try_reserve(1).map_err(|_| allocation())?;
        indices.push(j);
    }
    let c = nonzero.len();
    let mut e = dense(c, 1, options, budget)?;
    e.fill(1.);
    let mut v = dense(c, c, options, budget)?;
    for i in 0..c {
        budget.step()?;
        v[i * c + i] = 1.;
    }
    let mut signs = dense(columns, 1, options, budget)?;
    let mut pending = Vec::new();
    let mut support = Vec::new();
    let mut zero_columns = Vec::new();
    for root in 0..columns {
        budget.step()?;
        if !active[root] || signs[root] != 0. {
            continue;
        }
        pending.clear();
        support.clear();
        pending.try_reserve(1).map_err(|_| allocation())?;
        pending.push(root);
        signs[root] = 1.;
        let mut consistent = true;
        let mut singular = false;
        while let Some(u) = pending.pop() {
            budget.step()?;
            support.try_reserve(1).map_err(|_| allocation())?;
            support.push(u);
            singular |= forced[u];
            for &(other, relation) in &adjacency[u] {
                budget.step()?;
                let wanted = signs[u] * relation;
                if signs[other] == 0. {
                    signs[other] = wanted;
                    pending.try_reserve(1).map_err(|_| allocation())?;
                    pending.push(other);
                } else {
                    consistent &= signs[other] == wanted;
                }
            }
        }
        if consistent && !singular {
            budget.check()?;
            support.sort_unstable();
            budget.check()?;
            // Every traversal starts at an active vertex, so support is nonempty.
            let last = local[*support.last().ok_or(Error::InternalInvariant {
                reason: "empty non-branching component",
            })?];
            for &j in &support {
                budget.step()?;
                e[local[j]] = signs[j];
                v[local[j] * c + last] = 1.;
            }
            zero_columns.try_reserve(1).map_err(|_| allocation())?;
            zero_columns.push(last);
        }
    }
    let mut r = dense(rows, c, options, budget)?;
    for row in 0..rows {
        for col in 0..c {
            for i in 0..c {
                budget.step()?;
                r[row * c + col] += d[row * columns + nonzero[i]] * e[i] * v[i * c + col];
            }
        }
    }
    budget.check()?;
    Ok(WeakReduction {
        nonzero,
        original_zero,
        zero_columns,
        e,
        v,
        r,
    })
}

impl WeakReduction {
    pub(super) fn kernel(
        &self,
        columns: usize,
        options: &SpectralOptions,
        budget: &mut WorkBudget<'_>,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        let c = self.nonzero.len();
        let rows = self.r.len() / c.max(1);
        // Unlike the printed Algorithm 2, R includes all reoriented columns.
        // With signed unit rows and a 0/1 transform, these zero tests are exact.
        for col in 0..c {
            budget.step_by(self.zero_columns.len())?;
            let selected = self.zero_columns.contains(&col);
            let mut zero = true;
            for row in 0..rows {
                budget.step()?;
                zero &= self.r[row * c + col] == 0.;
            }
            if selected != zero {
                return Err(Error::InternalInvariant {
                    reason: "weak reduction zero-column mismatch",
                });
            }
        }
        let count = self
            .original_zero
            .len()
            .checked_add(self.zero_columns.len())
            .ok_or(Error::SizeOverflow {
                operation: "non-branching kernel width",
            })?;
        let mut z = dense(columns, count, options, budget)?;
        let mut metric = dense(count, 1, options, budget)?;
        for (j, &original) in self.original_zero.iter().enumerate() {
            budget.step()?;
            z[original * count + j] = 1.;
            metric[j] = 1.;
        }
        for (j, &col) in self.zero_columns.iter().enumerate() {
            let k = self.original_zero.len() + j;
            for (i, &original) in self.nonzero.iter().enumerate() {
                budget.step()?;
                let value = self.e[i] * self.v[i * c + col];
                z[original * count + k] = value;
                metric[k] += value * value;
            }
        }
        budget.check()?;
        Ok((z, metric))
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn upper(
    b1: &[f64],
    d: &[f64],
    n: usize,
    m: usize,
    p: usize,
    options: &SpectralOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    check_shape(b1, n, p)?;
    for &value in b1 {
        budget.step()?;
        if value != 0. && value != 1. && value != -1. {
            return Err(ineligible("require signed unit retained incidence"));
        }
    }
    let reduction = weak_reduction(d, m, p, options, budget)?;
    let (z, metric) = reduction.kernel(p, options, budget)?;
    let k = metric.len();
    for row in 0..m {
        for col in 0..k {
            let mut dot = 0.;
            for j in 0..p {
                budget.step()?;
                dot += d[row * p + j] * z[j * k + col];
            }
            if dot != 0. {
                return Err(Error::InternalInvariant {
                    reason: "non-branching kernel violates outside constraint",
                });
            }
        }
    }
    let mut factor = dense(n, k, options, budget)?;
    for row in 0..n {
        for col in 0..k {
            let mut sum = 0.;
            for j in 0..p {
                budget.step()?;
                sum += b1[row * p + j] * z[j * k + col];
            }
            factor[row * k + col] = sum / metric[col].sqrt();
        }
    }
    let mut matrix = dense(n, n, options, budget)?;
    for row in 0..n {
        for col in 0..=row {
            let mut sum = 0.;
            for j in 0..k {
                budget.step()?;
                sum += factor[row * k + j] * factor[col * k + j];
            }
            matrix[row * n + col] = sum;
            matrix[col * n + row] = sum;
        }
    }
    budget.check()?;
    Ok(matrix)
}

fn check_shape(values: &[f64], rows: usize, columns: usize) -> Result<()> {
    let expected = rows.checked_mul(columns).ok_or(Error::SizeOverflow {
        operation: "non-branching boundary shape",
    })?;
    if values.len() != expected {
        return Err(Error::ShapeMismatch {
            input: "non-branching boundary",
            expected,
            actual: values.len(),
        });
    }
    Ok(())
}

fn ineligible(reason: &'static str) -> Error {
    Error::InvalidParameter {
        parameter: "non-branching boundary",
        reason,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::execution::Execution;
    use std::sync::atomic::{AtomicBool, Ordering};

    #[test]
    fn interruption_during_reduction_and_assembly_returns_no_result() {
        let options = SpectralOptions::default();
        let b = [-1., 0., 0., 1.];
        let d = [1., -1.];
        let mut measured = WorkBudget::new(&Execution::default().max_work(u64::MAX)).unwrap();
        upper(&b, &d, 2, 1, 2, &options, &mut measured).unwrap();
        let work = measured.used();
        for at in [1, work / 4, work / 2, work - 1] {
            let flag = AtomicBool::new(false);
            let execution = Execution::default().max_work(u64::MAX).cancellation(&flag);
            let mut budget = WorkBudget::new(&execution).unwrap();
            budget.cancel_at_work(at);
            assert_eq!(
                upper(&b, &d, 2, 1, 2, &options, &mut budget),
                Err(Error::Cancelled)
            );
            assert!(flag.load(Ordering::Relaxed));
            let mut budget = WorkBudget::new(&Execution::default().max_work(at)).unwrap();
            assert_eq!(
                upper(&b, &d, 2, 1, 2, &options, &mut budget),
                Err(Error::WorkLimitExceeded { limit: at })
            );
        }
    }

    #[test]
    fn algebraic_sign_contradiction_and_invalid_entries_are_distinct() {
        let options = SpectralOptions::default();
        let mut budget = WorkBudget::new(&Execution::default()).unwrap();
        assert_eq!(
            upper(
                &[1., 0., 0., 1.],
                &[1., 1., 1., -1.],
                2,
                2,
                2,
                &options,
                &mut budget
            )
            .unwrap(),
            [0.; 4]
        );
        for d in [&[1., 1., 1.][..], &[1., 2., 0.], &[f64::NAN, 0., 0.]] {
            let mut budget = WorkBudget::new(&Execution::default()).unwrap();
            assert!(matches!(
                upper(&[1., 0., 0.], d, 1, 1, 3, &options, &mut budget),
                Err(Error::InvalidParameter {
                    parameter: "non-branching boundary",
                    ..
                })
            ));
        }
        let mut budget = WorkBudget::new(&Execution::default()).unwrap();
        assert!(matches!(
            weak_reduction(&[], usize::MAX, 2, &options, &mut budget),
            Err(Error::SizeOverflow { .. })
        ));
    }
}
