//! Bounded real event reduction and new-chain evaluation, Algorithms 1-2 of
//! Kirchner--Vittorietti, arXiv:2611.09610v1, with checked rational pivots.
use super::{
    LaplacianComponent, PersistentLaplacian, SpectralOptions, allocation, basis_ids, finish,
    numerical, schur, validate_pair,
};
use crate::execution::{Execution, WorkBudget};
use crate::filtration::SimplicialStage;
use crate::{Error, Result};

/// One requested full operator and its evaluated upper-term representative.
///
/// Scheduling is exact over Q (and therefore R for integer incidence), while
/// orthogonalization, matrices and spectra use f64. Integer overflow is an error;
/// no floating-point pivot threshold or finite-field rank prunes a pair.
#[derive(Clone, Debug)]
pub struct RelevantPairLaplacian {
    operator: PersistentLaplacian,
    representative: Option<(usize, usize)>,
    representative_order: usize,
}
impl RelevantPairLaplacian {
    /// Owned full operator in the requested K basis, including its own lower term.
    pub fn operator(&self) -> &PersistentLaplacian {
        &self.operator
    }

    /// Evaluated upper representative as indices into the supplied stage grid.
    ///
    /// The order is maximal endpoint, then maximal start, among the relevant
    /// pairs below this request. `None` means the upper term is zero. Its matrix
    /// is embedded in the requested basis by appending zero rows and columns;
    /// the representative's full operator must never be reused as a full result.
    pub fn representative(&self) -> Option<(usize, usize)> {
        self.representative
    }

    /// Requested coordinates occupied by the representative's upper matrix.
    /// Remaining coordinates are zero extension; `None` representatives use
    /// an empty range. This embedding survives release of the source/grid.
    pub fn upper_embedding(&self) -> std::ops::Range<usize> {
        0..self.representative_order
    }
}

/// Compute requested full persistent Laplacians using real new-chain events.
///
/// `stages` is a nondecreasing inclusive grid from one immutable source, and each
/// `(s,t)` in `requests` indexes that grid with `s <= t`. Output preserves request
/// order, duplicates, original scales, bases, coverage and skeleton certificates.
/// Equal births enter together; repeated grid scales share their earliest event
/// index. The grid may omit births. The empty-stage sentinel remains internal.
///
/// The final supplied stage bounds the signed boundary. Checked i128 rational
/// reduction retains R=B V with invertible upper-triangular V, including zero
/// columns. Algorithm 2 combines two event buckets; it is not a barcode-endpoint
/// filter. Only selected predecessors needed by requests are evaluated. At each,
/// all eligible V columns are orthonormalized together in final-chain coordinates
/// before applying B. No separately computed nonorthogonal Grams are summed.
/// Each request receives zero extension and its own lower term.
///
/// This explicit bounded alternative does not change default dispatch or cache a
/// filtration globally. Preparation and the schedule can be dense/quadratic.
/// `max_dense_entries` caps each rational boundary/transformation buffer, each
/// floating buffer, the triangular grid size and cumulative returned matrix
/// entries. A rational entry occupies 32 bytes; this is not a total memory quota.
/// A single execution budget covers validation, reduction, scheduling, numerical
/// assembly and returned metadata. Spectrum calls start separate budgets.
///
/// # Errors
/// Rejects invalid indices/order, mixed sources/certificates, missing skeletons,
/// invalid options, storage/allocation limits, checked rational overflow,
/// cancellation, work exhaustion or failed orthogonality/boundary residuals.
/// Overflow never triggers an unchecked or over-budget exhaustive fallback.
pub fn relevant_pair_laplacians(
    stages: &[SimplicialStage<'_>],
    requests: &[(usize, usize)],
    degree: usize,
    options: SpectralOptions,
    execution: &Execution<'_>,
) -> Result<Vec<RelevantPairLaplacian>> {
    let mut budget = WorkBudget::new(execution)?;
    options.validate()?;
    degree.checked_add(1).ok_or(Error::SizeOverflow {
        operation: "adjacent spectral degree",
    })?;
    for &(s, t) in requests {
        budget.step()?;
        if s > t || t >= stages.len() {
            return Err(Error::InvalidParameter {
                parameter: "spectral pair indices",
                reason: "require start <= end < stage count",
            });
        }
    }
    let mut result = Vec::new();
    result
        .try_reserve_exact(requests.len())
        .map_err(|_| allocation())?;
    let Some(target) = stages.last() else {
        return Ok(result);
    };
    for (i, stage) in stages.iter().enumerate() {
        validate_pair(
            stage,
            target,
            degree,
            LaplacianComponent::Full,
            &options,
            &mut budget,
        )?;
        if i > 0 && stages[i - 1].scale() > stage.scale() {
            return Err(Error::InvalidParameter {
                parameter: "spectral stage grid",
                reason: "scales must be nondecreasing",
            });
        }
    }
    let grid = stages.len();
    let pair_count = grid
        .checked_add(1)
        .and_then(|n| grid.checked_mul(n))
        .map(|n| n / 2)
        .ok_or_else(overflow)?;
    schur::check_entries(pair_count, &options)?;
    let ids = basis_ids(target, degree, &mut budget)?;
    let upper_ids = basis_ids(target, degree + 1, &mut budget)?;
    let n = ids.len();
    let p = upper_ids.len();
    let mut sizes = Vec::new();
    sizes.try_reserve_exact(grid).map_err(|_| allocation())?;
    for stage in stages {
        budget.step_by(n)?;
        sizes.push(ids.partition_point(|id| id.index() < stage.len()));
    }
    let total = requests.iter().try_fold(0usize, |total, &(s, _)| {
        sizes[s]
            .checked_mul(sizes[s])
            .and_then(|n| total.checked_add(n))
            .ok_or_else(overflow)
    })?;
    schur::check_entries(total, &options)?;
    if requests.is_empty() {
        return Ok(result);
    }
    let entries = n.checked_mul(p).ok_or_else(overflow)?;
    let square = p.checked_mul(p).ok_or_else(overflow)?;
    schur::check_entries(entries, &options)?;
    schur::check_entries(square, &options)?;
    let mut boundary = Vec::new();
    boundary
        .try_reserve_exact(entries)
        .map_err(|_| allocation())?;
    boundary.resize(entries, Rational::ZERO);
    for (j, &id) in upper_ids.iter().enumerate() {
        for term in target.boundary(id).ok_or(Error::InternalInvariant {
            reason: "stage boundary missing",
        })? {
            budget.step()?;
            let i = ids
                .binary_search(&term.face)
                .map_err(|_| Error::InternalInvariant {
                    reason: "q-face absent",
                })?;
            boundary[j * n + i] = Rational {
                num: i128::from(term.coefficient),
                den: 1,
            };
        }
    }
    let (reduced, transform) = reduce(&boundary, n, p, &mut budget)?;
    let mut events = Vec::new();
    events.try_reserve_exact(p).map_err(|_| allocation())?;
    for j in 0..p {
        budget.step_by(n.checked_add(grid).ok_or_else(overflow)?)?;
        let low = (0..n).rev().find(|&i| reduced[j * n + i].num != 0);
        // Stage 0 is the internal empty stage; real grid indices are shifted.
        let start = low.map_or(0, |i| {
            stages.partition_point(|stage| stage.len() <= ids[i].index()) + 1
        });
        let end = stages.partition_point(|stage| stage.len() <= upper_ids[j].index()) + 1;
        events.push((start, end));
    }
    drop(reduced);
    let selected = select_pairs(&events, grid, pair_count, &mut budget)?;
    let mut solved: Vec<((usize, usize), Vec<f64>, f64)> = Vec::new();
    solved
        .try_reserve_exact(requests.len().min(selected.len()))
        .map_err(|_| allocation())?;
    for &(s, t) in requests {
        budget.step_by(selected.len())?;
        let representative = selected
            .iter()
            .rev()
            .copied()
            .find(|&(a, b)| a <= s && b <= t);
        let mut matrix = schur::dense(sizes[s], sizes[s], &options, &mut budget)?;
        let mut residual = 0.0;
        if let Some(pair) = representative {
            budget.step_by(solved.len())?;
            let index = if let Some(index) = solved.iter().position(|(p, _, _)| *p == pair) {
                index
            } else {
                let (up, error) = new_chain_upper(
                    &boundary,
                    &transform,
                    &events,
                    n,
                    p,
                    pair.0 + 1,
                    pair.1 + 1,
                    sizes[pair.0],
                    &options,
                    &mut budget,
                )?;
                solved.push((pair, up, error));
                solved.len() - 1
            };
            let (_, old, error) = &solved[index];
            residual = *error;
            let count = sizes[pair.0];
            for i in 0..count {
                budget.step_by(count)?;
                matrix[i * sizes[s]..i * sizes[s] + count]
                    .copy_from_slice(&old[i * count..(i + 1) * count]);
            }
        }
        // Match the original pair's Gram scale, not its predecessor's scale.
        let mut scale = 0.0_f64;
        for i in 0..sizes[t] {
            let mut diagonal = 0.0_f64;
            for (j, &(_, end)) in events.iter().enumerate() {
                budget.step()?;
                if end <= t + 1 && boundary[j * n + i].num != 0 {
                    diagonal += 1.0;
                }
            }
            scale = scale.max(diagonal);
        }
        let operator = finish(
            &stages[s],
            &stages[t],
            degree,
            &ids[..sizes[s]],
            &ids[..sizes[t]],
            matrix,
            residual,
            options,
            LaplacianComponent::Full,
            scale,
            &mut budget,
        )?;
        result.push(RelevantPairLaplacian {
            operator,
            representative,
            representative_order: representative.map_or(0, |(s, _)| sizes[s]),
        });
    }
    budget.check()?;
    Ok(result)
}

fn select_pairs(
    events: &[(usize, usize)],
    grid: usize,
    capacity: usize,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<(usize, usize)>> {
    let mut selected = Vec::new();
    selected
        .try_reserve_exact(capacity)
        .map_err(|_| allocation())?;
    for t in 1..=grid {
        for s in 1..=t {
            budget.step_by(events.len())?;
            if events.iter().any(|&(a, b)| a <= s && b == t)
                && events.iter().any(|&(a, b)| a == s && b <= t)
            {
                selected.push((s - 1, t - 1));
            }
        }
    }
    Ok(selected)
}

#[allow(clippy::too_many_arguments)]
fn new_chain_upper(
    boundary: &[Rational],
    transform: &[Rational],
    events: &[(usize, usize)],
    n: usize,
    p: usize,
    s: usize,
    t: usize,
    count: usize,
    options: &SpectralOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<(Vec<f64>, f64)> {
    let mut q = schur::dense(p, p, options, budget)?;
    let mut rank = 0;
    let mut residual = 0.0_f64;
    for (j, &(a, b)) in events.iter().enumerate() {
        budget.step()?;
        if a > s || b > t {
            continue;
        }
        for i in 0..p {
            budget.step()?;
            q[rank * p + i] = transform[j * p + i].as_f64();
        }
        budget.step_by(p)?;
        let original_norm = (0..p).map(|i| q[rank * p + i].powi(2)).sum::<f64>().sqrt();
        // Reorthogonalize the whole collected space, including cycle columns.
        for _ in 0..2 {
            for k in 0..rank {
                budget.step_by(p)?;
                let dot: f64 = (0..p).map(|i| q[rank * p + i] * q[k * p + i]).sum();
                for i in 0..p {
                    budget.step()?;
                    q[rank * p + i] -= dot * q[k * p + i];
                }
            }
        }
        budget.step_by(p)?;
        let norm = (0..p).map(|i| q[rank * p + i].powi(2)).sum::<f64>().sqrt();
        if !norm.is_finite() || norm <= options.relative_tolerance * original_norm {
            return Err(numerical("new-chain orthogonalization"));
        }
        for i in 0..p {
            q[rank * p + i] /= norm;
        }
        for k in 0..=rank {
            budget.step_by(p)?;
            let dot: f64 = (0..p).map(|i| q[rank * p + i] * q[k * p + i]).sum();
            residual = residual.max((dot - f64::from(k == rank)).abs());
        }
        rank += 1;
    }
    let mut factors = schur::dense(n, rank, options, budget)?;
    let mut inside = 0.0_f64;
    let mut outside = 0.0_f64;
    for i in 0..n {
        for k in 0..rank {
            budget.step_by(p)?;
            let value: f64 = (0..p)
                .map(|j| boundary[j * n + i].as_f64() * q[k * p + j])
                .sum();
            factors[i * rank + k] = value;
            if i < count {
                inside = inside.max(value.abs());
            } else {
                outside = outside.max(value.abs());
            }
        }
    }
    residual = residual.max(outside / inside.max(1.0));
    if !residual.is_finite() || residual > options.residual_tolerance {
        return Err(numerical("new-chain basis residual"));
    }
    let mut matrix = schur::dense(count, count, options, budget)?;
    for i in 0..count {
        for j in 0..=i {
            budget.step_by(rank)?;
            let dot: f64 = (0..rank)
                .map(|k| factors[i * rank + k] * factors[j * rank + k])
                .sum();
            matrix[i * count + j] = dot;
            matrix[j * count + i] = dot;
        }
    }
    Ok((matrix, residual))
}

fn reduce(
    boundary: &[Rational],
    rows: usize,
    columns: usize,
    budget: &mut WorkBudget<'_>,
) -> Result<(Vec<Rational>, Vec<Rational>)> {
    let mut reduced = Vec::new();
    budget.check()?;
    reduced
        .try_reserve_exact(boundary.len())
        .map_err(|_| allocation())?;
    reduced.extend_from_slice(boundary);
    let mut transform = Vec::new();
    let entries = columns.checked_mul(columns).ok_or_else(overflow)?;
    transform
        .try_reserve_exact(entries)
        .map_err(|_| allocation())?;
    transform.resize(entries, Rational::ZERO);
    let mut owners = Vec::new();
    owners.try_reserve_exact(rows).map_err(|_| allocation())?;
    owners.resize(rows, None);
    budget.check()?;
    for j in 0..columns {
        transform[j * columns + j] = Rational { num: 1, den: 1 };
        loop {
            budget.step_by(rows)?;
            let Some(low) = (0..rows).rev().find(|&i| reduced[j * rows + i].num != 0) else {
                break;
            };
            let Some(i) = owners[low] else {
                owners[low] = Some(j);
                break;
            };
            let factor = reduced[j * rows + low].divide(reduced[i * rows + low])?;
            for k in 0..rows {
                budget.step()?;
                reduced[j * rows + k] =
                    reduced[j * rows + k].subtract(factor.multiply(reduced[i * rows + k])?)?;
            }
            for k in 0..columns {
                budget.step()?;
                transform[j * columns + k] = transform[j * columns + k]
                    .subtract(factor.multiply(transform[i * columns + k])?)?;
            }
        }
    }
    Ok((reduced, transform))
}

fn overflow() -> Error {
    Error::SizeOverflow {
        operation: "exact real event reduction",
    }
}

/// Normalized checked fractions; deliberately private and bounded, not a field
/// extension of the finite-field persistence reducer.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Rational {
    num: i128,
    den: i128,
}
impl Rational {
    const ZERO: Self = Self { num: 0, den: 1 };
    fn new(num: i128, den: i128) -> Result<Self> {
        if den == 0 {
            return Err(overflow());
        }
        let abs = num.checked_abs().ok_or_else(overflow)?;
        let divisor = gcd(abs, den.checked_abs().ok_or_else(overflow)?);
        let sign = if den < 0 { -1 } else { 1 };
        Ok(Self {
            num: (num / divisor).checked_mul(sign).ok_or_else(overflow)?,
            den: den / divisor * sign,
        })
    }
    fn multiply(self, other: Self) -> Result<Self> {
        let a = gcd(self.num.checked_abs().ok_or_else(overflow)?, other.den);
        let b = gcd(other.num.checked_abs().ok_or_else(overflow)?, self.den);
        Self::new(
            (self.num / a)
                .checked_mul(other.num / b)
                .ok_or_else(overflow)?,
            (self.den / b)
                .checked_mul(other.den / a)
                .ok_or_else(overflow)?,
        )
    }
    fn divide(self, other: Self) -> Result<Self> {
        self.multiply(Self::new(other.den, other.num)?)
    }
    fn subtract(self, other: Self) -> Result<Self> {
        let divisor = gcd(self.den, other.den);
        let left = self
            .num
            .checked_mul(other.den / divisor)
            .ok_or_else(overflow)?;
        let right = other
            .num
            .checked_mul(self.den / divisor)
            .ok_or_else(overflow)?;
        Self::new(
            left.checked_sub(right).ok_or_else(overflow)?,
            self.den
                .checked_mul(other.den / divisor)
                .ok_or_else(overflow)?,
        )
    }
    fn as_f64(self) -> f64 {
        self.num as f64 / self.den as f64
    }
}
fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        let remainder = a % b;
        a = b;
        b = remainder;
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checked_rationals_do_not_turn_overflow_into_a_zero_pivot() {
        assert_eq!(Rational::new(2, -4).unwrap(), Rational { num: -1, den: 2 });
        assert!(Rational::new(i128::MIN, 1).is_err());
        assert!(
            Rational::new(i128::MAX, 1)
                .unwrap()
                .multiply(Rational::new(2, 1).unwrap())
                .is_err()
        );
        assert!(Rational::ZERO.divide(Rational::ZERO).is_err());
        let a = Rational::new(i128::MAX, 2).unwrap();
        let b = Rational::new(2, i128::MAX).unwrap();
        assert_eq!(a.multiply(b).unwrap(), Rational::new(1, 1).unwrap());
    }

    #[test]
    fn fractional_transform_retains_zero_reduced_column() {
        let b = [1, 0, 0, 2, 1, 1].map(|num| Rational { num, den: 1 });
        let (r, v) = reduce(
            &b,
            2,
            3,
            &mut WorkBudget::new(&Execution::default()).unwrap(),
        )
        .unwrap();
        assert_eq!(&r[4..], &[Rational::ZERO; 2]);
        assert_eq!(
            &v[6..],
            &[
                Rational::new(-1, 1).unwrap(),
                Rational::new(-1, 2).unwrap(),
                Rational::new(1, 1).unwrap()
            ]
        );
    }
}
