//! Bounded ordinary q=1 HeCS reference solve, independently transcribed from
//! Savostianov--Tudisco--Guglielmi, arXiv:2401.15492v1, Algorithms 5.1/6.1.
use super::*;

/// Policy for the ordinary triangular upper-Hodge HeCS/CGLS reference solve.
///
/// Uses dense f64 factors/active matrices, without a runtime dependency. Each
/// dense matrix buffer is bounded separately; this is not an aggregate memory quota.
/// No sparse scaling, persistent adaptation or end-to-end speedup is promised.
#[derive(Clone, Copy, Debug)]
pub struct HecsOptions {
    /// Maximum entries per dense matrix buffer; default [`DEFAULT_MAX_DENSE_ENTRIES`].
    /// The m by m original matrix is checked before any dense allocation.
    pub max_dense_entries: usize,
    /// Maximum CGLS updates; zero permits an already satisfactory zero iterate.
    pub max_iterations: usize,
    /// Positive finite target for `||Ux-b||2 / max(1,||b||2)`; default 1e-9.
    pub residual_tolerance: f64,
    /// Peel free triangles before heavy selection, as in paper Section 6.2.
    /// Default false selects directly by Algorithm 6.1.
    pub peel_core: bool,
}
impl Default for HecsOptions {
    fn default() -> Self {
        Self {
            max_dense_entries: DEFAULT_MAX_DENSE_ENTRIES,
            max_iterations: 4000,
            residual_tolerance: 1e-9,
            peel_core: false,
        }
    }
}

/// Owned result of a successful ordinary q=1 upper-Hodge HeCS solve.
///
/// Coordinates, weights and original basis handles survive source release;
/// retain the source description to interpret owner-local handles geometrically
/// or reconstruct the full equation. The solution uses the Moore-Penrose lift
/// and lies numerically in the original boundary image. Factor/active matrices
/// are bounded dense reference outputs, not ambient inverses or eigenpairs.
pub struct HecsSolution<Id> {
    basis: Vec<Id>,
    triangles: Vec<Id>,
    weights: Vec<f64>,
    selected: Vec<usize>,
    pairs: Vec<[usize; 2]>,
    edge_order: Vec<usize>,
    triangle_order: Vec<usize>,
    factor: Vec<f64>,
    active: Vec<f64>,
    values: Vec<f64>,
    residuals: Vec<f64>,
    options: HecsOptions,
}
impl<Id> HecsSolution<Id> {
    /// Original edge coordinates, preserving the complete lower skeleton.
    pub fn basis(&self) -> &[Id] {
        &self.basis
    }
    /// Original triangle handles in weight/input order.
    pub fn triangle_basis(&self) -> &[Id] {
        &self.triangles
    }
    /// Explicit positive face weights: `U=B diag(weights) B^T`.
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
    /// Original triangle indexes in heavy acceptance order (initial peel first).
    pub fn selected_triangles(&self) -> &[usize] {
        &self.selected
    }
    /// `[original edge index, original triangle index]` in collapse order.
    pub fn collapse_pairs(&self) -> &[[usize; 2]] {
        &self.pairs
    }
    /// Collapse edges followed by unused original edge indexes in input order.
    pub fn edge_permutation(&self) -> &[usize] {
        &self.edge_order
    }
    /// Collapse triangles followed by unused original triangle indexes.
    pub fn triangle_permutation(&self) -> &[usize] {
        &self.triangle_order
    }
    /// Exact real boundary rank certified by unit-triangular image coverage.
    pub fn rank(&self) -> usize {
        self.pairs.len()
    }
    /// Row-major m by r weighted factor, in `edge_permutation()` row order.
    /// Columns follow `collapse_pairs()`. Its leading r by r minor is triangular.
    pub fn factor(&self) -> &[f64] {
        &self.factor
    }
    /// Row-major r by r `C+ U_permuted C+^T`, checked numerically positive.
    pub fn active_matrix(&self) -> &[f64] {
        &self.active
    }
    /// Original edge-coordinate solution; a fresh original residual is checked.
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    /// Fresh normalized original residuals, including iteration zero and final.
    pub fn residual_history(&self) -> &[f64] {
        &self.residuals
    }
    /// Number of completed CGLS updates.
    pub fn iterations(&self) -> usize {
        self.residuals.len() - 1
    }
    /// Numeric/resource policy used for this solve.
    pub fn options(&self) -> HecsOptions {
        self.options
    }
}

impl<Id: Copy + Eq + Hash> HodgeLaplacian<Id> {
    /// Solve an explicitly weighted ordinary q=1 upper-Hodge equation by HeCS.
    ///
    /// Requires this snapshot to have degree 1 and selection Upper. Each edge
    /// must have two distinct endpoints with coefficients -1/+1; each triangle
    /// has three signed unit edges on three vertices. Parallel edges and repeated
    /// faces are rejected, including arbitrary cellular chains that merely have
    /// three nonzeros per column. Orientations and source traversal are retained.
    /// Higher source degrees do not contribute to this q=1 equation.
    ///
    /// `weights` follows next_basis(), is finite/strictly positive, and defines
    /// `U=B diag(weights) B^T`; filtration values are never used as weights.
    /// `rhs` follows basis(). Heavy ties use the original triangle index; each
    /// greedy collapse picks the smallest free edge. Selected-image coverage is
    /// proved exactly over checked integers before forming a reduced problem.
    ///
    /// One fresh execution budget covers topology checks, selection, factor,
    /// solve, residual verification and output copies. Previously completed
    /// snapshot preparation is separate. Every success meets the fresh original
    /// residual target. The compatible-RHS check uses that same tolerance;
    /// this is numerical evidence, not an exact RHS membership certificate.
    ///
    /// # Errors
    /// Rejects unsupported windows/topology, shape/weight/options errors, dense
    /// limits before allocation, integer overflow, selected rank loss, numerically
    /// ambiguous pivots, incompatible RHS, breakdown/nonconvergence, arithmetic
    /// failure, allocation failure, cancellation and work exhaustion. No ridge,
    /// ambient inverse, fallback or apparently successful partial result is used.
    pub fn solve_hecs(
        &self,
        weights: &[f64],
        rhs: &[f64],
        options: &HecsOptions,
        execution: &Execution<'_>,
    ) -> Result<HecsSolution<Id>> {
        let mut budget = WorkBudget::new(execution)?;
        self.hecs_budget(weights, rhs, options, &mut budget)
    }

    fn hecs_budget(
        &self,
        weights: &[f64],
        rhs: &[f64],
        options: &HecsOptions,
        budget: &mut WorkBudget<'_>,
    ) -> Result<HecsSolution<Id>> {
        let m = self.basis.len();
        let t = self.upper.len();
        if self.degree != 1 || self.selection != HodgeSelection::Upper {
            return Err(parameter(
                "HeCS window",
                "requires ordinary degree 1, Upper selection",
            ));
        }
        dense_limit(m, m, options.max_dense_entries)?;
        if !options.residual_tolerance.is_finite() || options.residual_tolerance <= 0. {
            return Err(parameter(
                "HeCS residual tolerance",
                "must be finite and positive",
            ));
        }
        shape(weights.len(), t, "HeCS weights")?;
        shape(rhs.len(), m, "HeCS right-hand side")?;
        for (i, &w) in weights.iter().enumerate() {
            budget.step()?;
            if !w.is_finite() {
                return Err(Error::NonFiniteValue {
                    field: "HeCS weights",
                    index: Some(i),
                });
            }
            if w <= 0. {
                return Err(parameter("HeCS weights", "must be strictly positive"));
            }
        }
        for (i, &v) in rhs.iter().enumerate() {
            budget.step()?;
            if !v.is_finite() {
                return Err(Error::NonFiniteValue {
                    field: "HeCS right-hand side",
                    index: Some(i),
                });
            }
        }
        let incidence = simplicial(self, budget)?;
        let (selected, pairs) = heavy(&incidence, m, weights, options.peel_core, budget)?;
        cover(&self.upper, m, &pairs, budget)?;
        let r = pairs.len();
        let mut edge_order = Vec::new();
        for pair in &pairs {
            budget.step()?;
            push(&mut edge_order, pair[0])?;
        }
        for e in 0..m {
            budget.step_by(pairs.len())?;
            if !pairs.iter().any(|p| p[0] == e) {
                push(&mut edge_order, e)?;
            }
        }
        let mut triangle_order = Vec::new();
        for pair in &pairs {
            budget.step()?;
            push(&mut triangle_order, pair[1])?;
        }
        for j in 0..t {
            budget.step_by(pairs.len())?;
            if !pairs.iter().any(|p| p[1] == j) {
                push(&mut triangle_order, j)?;
            }
        }
        let mut f = dense(m, r, options, budget)?;
        for (j, pair) in pairs.iter().enumerate() {
            for &(e, sign) in &self.upper[pair[1]] {
                budget.step()?;
                f[e * r + j] = sign as f64 * weights[pair[1]].sqrt();
            }
        }
        let mut u = dense(m, m, options, budget)?;
        for (terms, &weight) in self.upper.iter().zip(weights) {
            for &(a, left) in terms {
                for &(b, right) in terms {
                    budget.step()?;
                    u[a * m + b] += left as f64 * right as f64 * weight;
                    finite(u[a * m + b])?;
                }
            }
        }
        let ft = transpose(&f, m, r, budget)?;
        let gram = product(&ft, &f, [r, m, r], options, budget)?;
        let l = cholesky(&gram, r, options, budget)?;
        let mut lift = dense(m, r, options, budget)?;
        for e in 0..m {
            let row = solve_cholesky(&l, &f[e * r..(e + 1) * r], budget)?;
            budget.step_by(r)?;
            lift[e * r..(e + 1) * r].copy_from_slice(&row);
        }
        let cp = transpose(&lift, m, r, budget)?;
        let active = product(
            &product(&cp, &u, [r, m, m], options, budget)?,
            &lift,
            [r, m, r],
            options,
            budget,
        )?;
        cholesky(&active, r, options, budget)?;
        let h = action(&cp, r, rhs, budget)?;
        let projected = action(&f, m, &h, budget)?;
        if relative_residual(&projected, rhs, budget)? > options.residual_tolerance {
            return Err(parameter(
                "HeCS right-hand side",
                "outside the real boundary image at the requested tolerance",
            ));
        }
        let (values, residuals) = cgls(&active, &h, &u, rhs, &lift, options, budget)?;
        let mut factor = dense(m, r, options, budget)?;
        for (i, &e) in edge_order.iter().enumerate() {
            budget.step_by(r)?;
            factor[i * r..(i + 1) * r].copy_from_slice(&f[e * r..(e + 1) * r]);
        }
        let result = HecsSolution {
            basis: copy(&self.basis, budget)?,
            triangles: copy(&self.next_basis, budget)?,
            weights: copy(weights, budget)?,
            selected,
            pairs,
            edge_order,
            triangle_order,
            factor,
            active,
            values,
            residuals,
            options: *options,
        };
        budget.check()?;
        Ok(result)
    }
}

fn parameter(parameter: &'static str, reason: &'static str) -> Error {
    Error::InvalidParameter { parameter, reason }
}
fn numerical(context: &'static str) -> Error {
    Error::NumericalFailure { context }
}
fn shape(actual: usize, expected: usize, input: &'static str) -> Result<()> {
    if actual != expected {
        Err(Error::ShapeMismatch {
            input,
            actual,
            expected,
        })
    } else {
        Ok(())
    }
}
fn copy<T: Copy>(input: &[T], budget: &mut WorkBudget<'_>) -> Result<Vec<T>> {
    budget.step_by(input.len())?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(input.len())
        .map_err(|_| allocation())?;
    result.extend_from_slice(input);
    budget.check()?;
    Ok(result)
}
fn dense_limit(rows: usize, columns: usize, limit: usize) -> Result<usize> {
    let entries = rows.checked_mul(columns).ok_or(size_error())?;
    if entries > limit {
        return Err(parameter(
            "HeCS dense entries",
            "exceeds the per-buffer limit",
        ));
    }
    Ok(entries)
}
fn dense(
    rows: usize,
    columns: usize,
    options: &HecsOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    zeros(
        dense_limit(rows, columns, options.max_dense_entries)?,
        budget,
    )
}

fn simplicial<Id>(op: &HodgeLaplacian<Id>, budget: &mut WorkBudget<'_>) -> Result<Vec<[usize; 3]>> {
    let mut edges = Vec::new();
    for terms in &op.lower {
        budget.step()?;
        if terms.len() != 2 || terms[0].1.abs() != 1 || terms[1].1 != -terms[0].1 {
            return Err(parameter(
                "HeCS topology",
                "requires signed simplicial edges",
            ));
        }
        push(&mut edges, [terms[0].0, terms[1].0])?;
    }
    let mut sorted = copy(&edges, budget)?;
    sorted.sort_unstable();
    budget.check()?;
    for pair in sorted.windows(2) {
        budget.step()?;
        if pair[0] == pair[1] {
            return Err(parameter("HeCS topology", "parallel edges are unsupported"));
        }
    }
    let mut incidence = Vec::new();
    for terms in &op.upper {
        budget.step()?;
        if terms.len() != 3 || terms.iter().any(|t| t.1.abs() != 1) {
            return Err(parameter(
                "HeCS topology",
                "requires signed simplicial triangles",
            ));
        }
        let mut vertices = Vec::new();
        for &(e, _) in terms {
            for &v in &edges[e] {
                budget.step()?;
                if !vertices.contains(&v) {
                    push(&mut vertices, v)?;
                }
            }
        }
        if vertices.len() != 3 {
            return Err(parameter(
                "HeCS topology",
                "triangle support must have three vertices",
            ));
        }
        push(&mut incidence, [terms[0].0, terms[1].0, terms[2].0])?;
    }
    let mut sorted = copy(&incidence, budget)?;
    sorted.sort_unstable();
    budget.check()?;
    for pair in sorted.windows(2) {
        budget.step()?;
        if pair[0] == pair[1] {
            return Err(parameter(
                "HeCS topology",
                "repeated triangles are unsupported",
            ));
        }
    }
    Ok(incidence)
}

fn collapse(
    incidence: &[[usize; 3]],
    m: usize,
    selected: &[usize],
    budget: &mut WorkBudget<'_>,
) -> Result<(Vec<[usize; 2]>, Vec<usize>)> {
    let mut remaining = copy(selected, budget)?;
    remaining.sort_unstable();
    budget.check()?;
    let mut pairs = Vec::new();
    while !remaining.is_empty() {
        let mut free = None;
        for e in 0..m {
            let mut adjacent = None;
            let mut count = 0;
            for (i, &t) in remaining.iter().enumerate() {
                budget.step()?;
                if incidence[t].contains(&e) {
                    count += 1;
                    adjacent = Some((i, t));
                }
                if count == 2 {
                    break;
                }
            }
            if count == 1 {
                free = adjacent.map(|(i, t)| (e, i, t));
                break;
            }
        }
        let Some((e, i, t)) = free else { break };
        budget.step_by(remaining.len())?;
        remaining.remove(i);
        push(&mut pairs, [e, t])?;
    }
    Ok((pairs, remaining))
}
fn heavy(
    incidence: &[[usize; 3]],
    m: usize,
    weights: &[f64],
    peel: bool,
    budget: &mut WorkBudget<'_>,
) -> Result<(Vec<usize>, Vec<[usize; 2]>)> {
    let mut candidates = Vec::new();
    for t in 0..weights.len() {
        budget.step()?;
        push(&mut candidates, t)?;
    }
    let mut selected = Vec::new();
    if peel {
        let (initial, core) = collapse(incidence, m, &candidates, budget)?;
        for p in initial {
            budget.step()?;
            push(&mut selected, p[1])?;
        }
        candidates = core;
    }
    budget.check()?;
    candidates.sort_unstable_by(|&a, &b| weights[b].total_cmp(&weights[a]).then(a.cmp(&b)));
    budget.check()?;
    for t in candidates {
        budget.step()?;
        push(&mut selected, t)?;
        if !collapse(incidence, m, &selected, budget)?.1.is_empty() {
            selected.pop();
        }
    }
    let (pairs, _) = collapse(incidence, m, &selected, budget)?;
    Ok((selected, pairs))
}
fn coefficient(column: &[(usize, i64)], row: usize) -> i128 {
    column
        .iter()
        .find(|&&(r, _)| r == row)
        .map_or(0, |&(_, v)| i128::from(v))
}
fn cover(
    columns: &[Vec<(usize, i64)>],
    m: usize,
    pairs: &[[usize; 2]],
    budget: &mut WorkBudget<'_>,
) -> Result<()> {
    let r = pairs.len();
    let mut x = Vec::new();
    x.try_reserve_exact(r).map_err(|_| allocation())?;
    x.resize(r, 0_i128);
    for column in columns {
        for (i, &[e, t]) in pairs.iter().enumerate() {
            let mut value = coefficient(column, e);
            for j in 0..i {
                budget.step()?;
                value = value
                    .checked_sub(
                        coefficient(&columns[pairs[j][1]], e)
                            .checked_mul(x[j])
                            .ok_or(size_error())?,
                    )
                    .ok_or(size_error())?;
            }
            budget.step()?;
            x[i] = value
                .checked_mul(coefficient(&columns[t], e))
                .ok_or(size_error())?;
        }
        for e in 0..m {
            let mut value = 0_i128;
            for (p, &v) in pairs.iter().zip(&x) {
                budget.step()?;
                value = value
                    .checked_add(
                        coefficient(&columns[p[1]], e)
                            .checked_mul(v)
                            .ok_or(size_error())?,
                    )
                    .ok_or(size_error())?;
            }
            budget.step()?;
            if value != coefficient(column, e) {
                return Err(parameter(
                    "HeCS rank",
                    "selected triangles do not span the original real boundary image",
                ));
            }
        }
    }
    budget.check()
}

fn dot(a: &[f64], b: &[f64], budget: &mut WorkBudget<'_>) -> Result<f64> {
    let mut value = 0.;
    for (&a, &b) in a.iter().zip(b) {
        budget.step()?;
        value += a * b;
    }
    finite(value)?;
    Ok(value)
}
fn norm(a: &[f64], budget: &mut WorkBudget<'_>) -> Result<f64> {
    let mut value: f64 = 0.;
    for &a in a {
        budget.step()?;
        value = value.hypot(a);
    }
    finite(value)?;
    Ok(value)
}
fn action(a: &[f64], rows: usize, x: &[f64], budget: &mut WorkBudget<'_>) -> Result<Vec<f64>> {
    let mut y = zeros(rows, budget)?;
    for (i, v) in y.iter_mut().enumerate() {
        *v = dot(&a[i * x.len()..(i + 1) * x.len()], x, budget)?;
    }
    Ok(y)
}
fn transpose(a: &[f64], rows: usize, cols: usize, budget: &mut WorkBudget<'_>) -> Result<Vec<f64>> {
    let mut result = zeros(a.len(), budget)?;
    for i in 0..rows {
        for j in 0..cols {
            budget.step()?;
            result[j * rows + i] = a[i * cols + j];
        }
    }
    Ok(result)
}
fn product(
    a: &[f64],
    b: &[f64],
    shape: [usize; 3],
    options: &HecsOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    let [rows, inner, cols] = shape;
    let mut result = dense(rows, cols, options, budget)?;
    for i in 0..rows {
        for j in 0..cols {
            let mut value = 0.;
            for k in 0..inner {
                budget.step()?;
                value += a[i * inner + k] * b[k * cols + j];
            }
            finite(value)?;
            result[i * cols + j] = value;
        }
    }
    Ok(result)
}
fn cholesky(
    a: &[f64],
    n: usize,
    options: &HecsOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    let mut l = dense(n, n, options, budget)?;
    let mut scale: f64 = 0.;
    for i in 0..n {
        budget.step()?;
        scale = scale.max(a[i * n + i]);
    }
    for i in 0..n {
        for j in 0..=i {
            let value = a[i * n + j] - dot(&l[i * n..i * n + j], &l[j * n..j * n + j], budget)?;
            budget.step()?;
            finite(value)?;
            l[i * n + j] = if i == j {
                if value <= 64. * n as f64 * f64::EPSILON * scale {
                    return Err(numerical("HeCS ambiguous active pivot"));
                }
                value.sqrt()
            } else {
                value / l[j * n + j]
            };
            finite(l[i * n + j])?;
        }
    }
    Ok(l)
}
fn solve_cholesky(l: &[f64], b: &[f64], budget: &mut WorkBudget<'_>) -> Result<Vec<f64>> {
    let n = b.len();
    let mut x = copy(b, budget)?;
    for i in 0..n {
        budget.step()?;
        x[i] = (x[i] - dot(&l[i * n..i * n + i], &x[..i], budget)?) / l[i * n + i];
        finite(x[i])?;
    }
    for i in (0..n).rev() {
        let mut sum = 0.;
        for j in i + 1..n {
            budget.step()?;
            sum += l[j * n + i] * x[j];
        }
        budget.step()?;
        x[i] = (x[i] - sum) / l[i * n + i];
        finite(x[i])?;
    }
    Ok(x)
}
fn relative_residual(actual: &[f64], rhs: &[f64], budget: &mut WorkBudget<'_>) -> Result<f64> {
    let mut difference: f64 = 0.;
    for (&a, &b) in actual.iter().zip(rhs) {
        budget.step()?;
        difference = difference.hypot(a - b);
    }
    let residual = difference / norm(rhs, budget)?.max(1.);
    finite(residual)?;
    Ok(residual)
}
fn cgls(
    a: &[f64],
    h: &[f64],
    u: &[f64],
    b: &[f64],
    lift: &[f64],
    options: &HecsOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let n = h.len();
    let at = transpose(a, n, n, budget)?;
    let mut x = zeros(n, budget)?;
    let mut residual = copy(h, budget)?;
    let mut s = action(&at, n, &residual, budget)?;
    let mut p = copy(&s, budget)?;
    let mut gamma = dot(&s, &s, budget)?;
    let mut history = Vec::new();
    for iteration in 0..=options.max_iterations {
        let values = action(lift, b.len(), &x, budget)?;
        let fresh = relative_residual(&action(u, b.len(), &values, budget)?, b, budget)?;
        budget.step()?;
        push(&mut history, fresh)?;
        if fresh <= options.residual_tolerance {
            return Ok((values, history));
        }
        if iteration == options.max_iterations {
            return Err(numerical("HeCS CGLS iteration limit"));
        }
        let q = action(a, n, &p, budget)?;
        let denominator = dot(&q, &q, budget)?;
        if gamma == 0. || denominator == 0. {
            return Err(numerical("HeCS CGLS breakdown"));
        }
        let alpha = gamma / denominator;
        finite(alpha)?;
        for (x, &p) in x.iter_mut().zip(&p) {
            budget.step()?;
            *x += alpha * p;
            finite(*x)?;
        }
        for (r, &q) in residual.iter_mut().zip(&q) {
            budget.step()?;
            *r -= alpha * q;
            finite(*r)?;
        }
        s = action(&at, n, &residual, budget)?;
        let next_gamma = dot(&s, &s, budget)?;
        let beta = next_gamma / gamma;
        finite(beta)?;
        for (p, &s) in p.iter_mut().zip(&s) {
            budget.step()?;
            *p = s + beta * *p;
            finite(*p)?;
        }
        gamma = next_gamma;
    }
    Err(numerical("HeCS CGLS iteration limit"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicBool, Ordering};

    fn interrupted<T>(mut operation: impl FnMut(&mut WorkBudget<'_>) -> Result<T>) {
        let flag = AtomicBool::new(false);
        let execution = Execution::default().max_work(u64::MAX).cancellation(&flag);
        let mut complete = WorkBudget::new(&execution).unwrap();
        assert!(operation(&mut complete).is_ok());
        let mut stopped = WorkBudget::new(&execution).unwrap();
        stopped.cancel_at_work(complete.used() / 2);
        assert!(matches!(operation(&mut stopped), Err(Error::Cancelled)));
        assert!(flag.load(Ordering::Relaxed));
    }

    #[test]
    fn cancellation_inside_selection_factor_triangular_and_iterative_work() {
        let options = HecsOptions::default();
        interrupted(|budget| heavy(&[[0, 1, 2], [0, 3, 4]], 5, &[1., 2.], false, budget));
        let matrix = [4., 1., 0., 1., 4., 1., 0., 1., 4.];
        interrupted(|budget| cholesky(&matrix, 3, &options, budget));
        let l = cholesky(
            &matrix,
            3,
            &options,
            &mut WorkBudget::new(&Execution::default()).unwrap(),
        )
        .unwrap();
        interrupted(|budget| solve_cholesky(&l, &[1., 2., 3.], budget));
        interrupted(|budget| {
            cgls(
                &[1., 0., 0., 2.],
                &[1., 4.],
                &[1., 0., 0., 2.],
                &[1., 4.],
                &[1., 0., 0., 1.],
                &options,
                budget,
            )
        });
    }

    #[test]
    fn ambiguous_pivots_and_integer_certificate_overflow_are_errors() {
        let options = HecsOptions::default();
        let mut budget = WorkBudget::new(&Execution::default()).unwrap();
        assert!(cholesky(&[1., 0., 0., 1e-30], 2, &options, &mut budget).is_err());
        // Test the checked arithmetic leaf directly, outside public simplicial
        // input: multiplying this reconstructed coefficient must overflow i128.
        let columns = vec![
            vec![(0, 1), (1, i64::MAX)],
            vec![(1, 1), (2, i64::MAX)],
            vec![(2, 1), (3, i64::MAX)],
            vec![(0, 1)],
        ];
        assert!(matches!(
            cover(&columns, 4, &[[0, 0], [1, 1], [2, 2]], &mut budget),
            Err(Error::SizeOverflow { .. })
        ));
    }
}
