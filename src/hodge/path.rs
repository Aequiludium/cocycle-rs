//! Numerical path Omega spaces recomputed over the reals, never lifted from Fp.
use crate::complex::{DirectedGraph, Path, PathBoundaryMode, PathComplex, path_boundary_budget};
use crate::execution::{Execution, WorkBudget};
use crate::{Error, Result};
use std::collections::BTreeMap;

/// Bounded numerical real chain with orthonormal Omega bases through d+1.
///
/// Rank depends on the recorded tolerance. Ambiguous pivots and failed residuals
/// return errors; this is not an exact real-rank or arbitrary-condition certificate.
#[derive(Clone, Debug)]
pub struct PathRealChain {
    mode: PathBoundaryMode,
    tolerance: f64,
    limit: usize,
    boundary_residual: f64,
    paths: Vec<Vec<Path>>,
    bases: Vec<Vec<Vec<f64>>>,
    boundaries: Vec<Vec<Vec<f64>>>,
}
impl PathRealChain {
    /// Boundary convention; allowed words remain regular in both models.
    pub fn mode(&self) -> PathBoundaryMode {
        self.mode
    }
    /// Relative pivot/residual tolerance, strictly between zero and 1e-3.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }
    /// Largest measured forbidden-face, coordinate-boundary or boundary-square residual.
    pub fn boundary_residual(&self) -> f64 {
        self.boundary_residual
    }
    /// Owned allowed-path identities in a retained degree.
    pub fn paths(&self, n: usize) -> Option<&[Path]> {
        self.paths.get(n).map(Vec::as_slice)
    }
    /// Orthonormal real Omega columns in allowed-path coordinates.
    pub fn bases(&self, n: usize) -> Option<&[Vec<f64>]> {
        self.bases.get(n).map(Vec::as_slice)
    }
    /// Coordinate boundaries between orthonormal Omega bases; B0 is zero.
    pub fn boundaries(&self, n: usize) -> Option<&[Vec<f64>]> {
        self.boundaries.get(n).map(Vec::as_slice)
    }
    /// Own the two adjacent factors of the standard real path Laplacian.
    ///
    /// # Errors
    /// Rejects a degree outside 0..=d or failed reservations. No square matrix
    /// or spectrum is materialized. This terminal starts a fresh execution budget.
    pub fn laplacian(&self, n: usize) -> Result<PathRealLaplacian> {
        self.laplacian_with(n, &Execution::default())
    }
    /// Snapshot with cooperative controls across all copied path/factor coordinates.
    ///
    /// # Errors
    /// Returns the errors of [`Self::laplacian`], cancellation or work exhaustion.
    pub fn laplacian_with(&self, n: usize, execution: &Execution<'_>) -> Result<PathRealLaplacian> {
        let mut budget = WorkBudget::new(execution)?;
        if n.checked_add(1)
            .is_none_or(|upper| upper >= self.bases.len())
        {
            return Err(Error::InvalidParameter {
                parameter: "path Laplacian degree",
                reason: "requires both adjacent boundary layers",
            });
        }
        Ok(PathRealLaplacian {
            mode: self.mode,
            dimension: n,
            source_time: 0.0,
            target_time: 0.0,
            paths: copy_paths(&self.paths[n], &mut budget)?,
            basis: copy_matrix(&self.bases[n], self.limit, &mut budget)?,
            down: copy_matrix(&self.boundaries[n], self.limit, &mut budget)?,
            upper: copy_matrix(&self.boundaries[n + 1], self.limit, &mut budget)?,
            limit: self.limit,
            tolerance: self.tolerance,
            boundary_residual: self.boundary_residual,
        })
    }
}

/// Owned factorized real operator D^T D+U U^T in an orthonormal source Omega basis.
#[derive(Clone, Debug)]
pub struct PathRealLaplacian {
    mode: PathBoundaryMode,
    dimension: usize,
    source_time: f64,
    target_time: f64,
    paths: Vec<Path>,
    basis: Vec<Vec<f64>>,
    down: Vec<Vec<f64>>,
    upper: Vec<Vec<f64>>,
    limit: usize,
    tolerance: f64,
    boundary_residual: f64,
}
impl PathRealLaplacian {
    /// Selected ambient model.
    pub fn mode(&self) -> PathBoundaryMode {
        self.mode
    }
    /// Computed chain degree.
    pub fn dimension(&self) -> usize {
        self.dimension
    }
    /// Actual source/target events, or (0,0) for a static chain snapshot.
    pub fn times(&self) -> (f64, f64) {
        (self.source_time, self.target_time)
    }
    /// Relative tolerance used for real Omega and persistent-domain construction.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }
    /// Largest measured chain/domain boundary residual; not an exact error bound.
    pub fn boundary_residual(&self) -> f64 {
        self.boundary_residual
    }
    /// Ambient target-path identities used to interpret the source basis columns.
    pub fn paths(&self) -> &[Path] {
        &self.paths
    }
    /// Orthonormal source Omega columns in these path coordinates.
    pub fn basis(&self) -> &[Vec<f64>] {
        &self.basis
    }
    /// Down-boundary columns; the lower term is D^T D.
    pub fn down(&self) -> &[Vec<f64>] {
        &self.down
    }
    /// Induced upper-boundary columns; the upper term is U U^T.
    pub fn upper(&self) -> &[Vec<f64>] {
        &self.upper
    }
    /// Apply the two rectangular factors without constructing a square matrix.
    ///
    /// # Errors
    /// Rejects shape mismatch, nonfinite input, failed reservations or arithmetic
    /// overflow. Work starts a fresh terminal; no eigenpair or exact rank is claimed.
    pub fn apply(&self, v: &[f64]) -> Result<Vec<f64>> {
        self.apply_with(v, &Execution::default())
    }
    /// Apply with shared coordinate/arithmetic work and cancellation controls.
    ///
    /// # Errors
    /// Returns the errors of [`Self::apply`], cancellation or work exhaustion.
    pub fn apply_with(&self, v: &[f64], execution: &Execution<'_>) -> Result<Vec<f64>> {
        let mut budget = WorkBudget::new(execution)?;
        if v.len() != self.basis.len() {
            return Err(Error::ShapeMismatch {
                input: "path action",
                expected: self.basis.len(),
                actual: v.len(),
            });
        }
        for (i, x) in v.iter().enumerate() {
            budget.step()?;
            if !x.is_finite() {
                return Err(Error::NonFiniteValue {
                    field: "path action",
                    index: Some(i),
                });
            }
        }
        let rows = self.down.first().map_or(0, Vec::len);
        let image = product(&self.down, v, rows, self.limit, &mut budget)?;
        let mut result = Vec::new();
        for c in &self.down {
            push(&mut result, dot(c, &image, &mut budget)?)?;
        }
        for c in &self.upper {
            let scalar = dot(c, v, &mut budget)?;
            for (x, y) in result.iter_mut().zip(c) {
                budget.step()?;
                *x += scalar * y;
                finite(*x)?;
            }
        }
        Ok(result)
    }
}

fn numerical() -> Error {
    Error::NumericalFailure {
        context: "real path Omega rank, residual or arithmetic",
    }
}
fn finite(v: f64) -> Result<f64> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(numerical())
    }
}
fn count(n: usize, limit: usize) -> Result<()> {
    if n > limit {
        Err(Error::InvalidParameter {
            parameter: "max_dense_entries",
            reason: "real path matrix/vector or degree limit exceeded",
        })
    } else {
        Ok(())
    }
}
fn cells(rows: usize, cols: usize, limit: usize) -> Result<()> {
    count(rows, limit)?;
    count(cols, limit)?;
    count(
        rows.checked_mul(cols).ok_or(Error::SizeOverflow {
            operation: "real path matrix",
        })?,
        limit,
    )
}
fn push<T>(v: &mut Vec<T>, value: T) -> Result<()> {
    v.try_reserve(1).map_err(|_| Error::AllocationFailed {
        context: "real path computation",
    })?;
    v.push(value);
    Ok(())
}
fn zeros(n: usize, limit: usize, budget: &mut WorkBudget<'_>) -> Result<Vec<f64>> {
    count(n, limit)?;
    budget.step_by(n)?;
    let mut v = Vec::new();
    v.try_reserve_exact(n)
        .map_err(|_| Error::AllocationFailed {
            context: "real path vector",
        })?;
    v.resize(n, 0.0);
    Ok(v)
}
fn copy_matrix(
    matrix: &[Vec<f64>],
    limit: usize,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<Vec<f64>>> {
    cells(matrix.first().map_or(0, Vec::len), matrix.len(), limit)?;
    let mut output = Vec::new();
    for column in matrix {
        budget.step()?;
        let mut c = zeros(column.len(), limit, budget)?;
        c.copy_from_slice(column);
        push(&mut output, c)?;
    }
    Ok(output)
}
fn copy_paths(paths: &[Path], budget: &mut WorkBudget<'_>) -> Result<Vec<Path>> {
    let mut output = Vec::new();
    for path in paths {
        budget.step_by(path.len())?;
        let mut c = Vec::new();
        c.try_reserve_exact(path.len())
            .map_err(|_| Error::AllocationFailed {
                context: "real path identities",
            })?;
        c.extend_from_slice(path);
        push(&mut output, c)?;
    }
    Ok(output)
}
fn dot(a: &[f64], b: &[f64], budget: &mut WorkBudget<'_>) -> Result<f64> {
    if a.len() != b.len() {
        return Err(numerical());
    }
    let mut sum = 0.0;
    for (x, y) in a.iter().zip(b) {
        budget.step()?;
        sum += x * y;
        finite(sum)?;
    }
    Ok(sum)
}
fn norm(v: &[f64], budget: &mut WorkBudget<'_>) -> Result<f64> {
    Ok(dot(v, v, budget)?.sqrt())
}
fn close(a: &[f64], b: &[f64], tolerance: f64, budget: &mut WorkBudget<'_>) -> Result<bool> {
    if a.len() != b.len() {
        return Ok(false);
    }
    let threshold = finite(100.0 * tolerance * norm(a, budget)?.max(norm(b, budget)?).max(1.0))?;
    for (x, y) in a.iter().zip(b) {
        budget.step()?;
        if finite(x - y)?.abs() > threshold {
            return Ok(false);
        }
    }
    Ok(true)
}
fn coordinate_residual(a: &[f64], b: &[f64], budget: &mut WorkBudget<'_>) -> Result<f64> {
    if a.len() != b.len() {
        return Err(numerical());
    }
    let mut error = 0.0f64;
    for (x, y) in a.iter().zip(b) {
        budget.step()?;
        error = error.max(finite(x - y)?.abs());
    }
    Ok(error)
}

fn product(
    columns: &[Vec<f64>],
    v: &[f64],
    rows: usize,
    limit: usize,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    if columns.len() != v.len() {
        return Err(numerical());
    }
    let mut out = zeros(rows, limit, budget)?;
    for (c, &s) in columns.iter().zip(v) {
        budget.step()?;
        if c.len() != rows {
            return Err(numerical());
        }
        for (x, y) in out.iter_mut().zip(c) {
            budget.step()?;
            *x += s * y;
            finite(*x)?;
        }
    }
    Ok(out)
}

fn nullspace(
    columns: &[Vec<f64>],
    rows: usize,
    tolerance: f64,
    limit: usize,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<Vec<f64>>> {
    let width = columns.len();
    cells(rows, width, limit)?;
    cells(width, width, limit)?;
    let mut a = Vec::new();
    for _ in 0..rows {
        budget.step()?;
        push(&mut a, zeros(width, limit, budget)?)?;
    }
    let mut scale = 1.0f64;
    for (j, c) in columns.iter().enumerate() {
        budget.step()?;
        if c.len() != rows {
            return Err(numerical());
        }
        for (i, &v) in c.iter().enumerate() {
            budget.step()?;
            finite(v)?;
            scale = scale.max(v.abs());
            a[i][j] = v;
        }
    }
    let threshold = finite(tolerance * scale)?;
    let mut pivots = Vec::new();
    for j in 0..width {
        budget.step()?;
        let start = pivots.len();
        let mut best = None;
        for (i, row) in a.iter().enumerate().skip(start) {
            budget.step()?;
            if best.is_none_or(|k: usize| row[j].abs() > a[k][j].abs()) {
                best = Some(i);
            }
        }
        let Some(row) = best else { break };
        let value = a[row][j].abs();
        if value > threshold / 100.0 && value < 100.0 * threshold {
            return Err(numerical());
        }
        if value <= threshold {
            continue;
        }
        a.swap(row, start);
        let value = a[start][j];
        for x in &mut a[start] {
            budget.step()?;
            *x = finite(*x / value)?;
        }
        let mut pivot = zeros(width, limit, budget)?;
        pivot.copy_from_slice(&a[start]);
        for (i, row) in a.iter_mut().enumerate() {
            budget.step()?;
            if i != start {
                let s = row[j];
                for (x, y) in row.iter_mut().zip(&pivot) {
                    budget.step()?;
                    *x = finite(*x - s * y)?;
                }
            }
        }
        push(&mut pivots, j)?;
    }
    let mut basis: Vec<Vec<f64>> = Vec::new();
    for j in 0..width {
        budget.step()?;
        if pivots.contains(&j) {
            continue;
        }
        let mut v = zeros(width, limit, budget)?;
        v[j] = 1.0;
        for (i, &pivot) in pivots.iter().enumerate() {
            budget.step()?;
            v[pivot] = -a[i][j];
        }
        for _ in 0..2 {
            for q in &basis {
                let s = dot(q, &v, budget)?;
                for (x, y) in v.iter_mut().zip(q) {
                    budget.step()?;
                    *x = finite(*x - s * y)?;
                }
            }
        }
        let length = norm(&v, budget)?;
        if length < 100.0 * threshold || length == 0.0 {
            return Err(numerical());
        }
        for x in &mut v {
            budget.step()?;
            *x = finite(*x / length)?;
        }
        let residual = product(columns, &v, rows, limit, budget)?;
        let zero = zeros(rows, limit, budget)?;
        if !close(&residual, &zero, finite(tolerance * scale)?, budget)? {
            return Err(numerical());
        }
        for q in &basis {
            if dot(q, &v, budget)?.abs() > 100.0 * tolerance {
                return Err(numerical());
            }
        }
        push(&mut basis, v)?;
    }
    Ok(basis)
}

/// Recompute real Omega and orthonormal coordinate boundaries through degree d+1.
///
/// This uses integer ambient boundaries, numerical row elimination and twice
/// reorthogonalized kernels, never a modular basis. Positive max_dense_entries
/// bounds each matrix/vector, layer count and kernel width squared, including
/// zero entries; it is not total memory. Tolerance must lie strictly in (0,1e-3).
///
/// # Errors
/// Rejects insufficient graph truncation, invalid limits/tolerance, size or
/// allocation failure, ambiguous numerical rank or failed chain residuals.
pub fn compute_real_path_chain(
    source: &PathComplex,
    d: usize,
    mode: PathBoundaryMode,
    max_dense_entries: usize,
    tolerance: f64,
) -> Result<PathRealChain> {
    compute_real_path_chain_with(
        source,
        d,
        mode,
        max_dense_entries,
        tolerance,
        &Execution::default(),
    )
}
/// Construct with one execution budget across path boundaries, dense algebra and residuals.
///
/// # Errors
/// Returns the errors of [`compute_real_path_chain`], cancellation or work exhaustion.
pub fn compute_real_path_chain_with(
    source: &PathComplex,
    d: usize,
    mode: PathBoundaryMode,
    limit: usize,
    tolerance: f64,
    execution: &Execution<'_>,
) -> Result<PathRealChain> {
    let mut budget = WorkBudget::new(execution)?;
    real_chain(source, d, mode, limit, tolerance, &mut budget)
}
fn real_chain(
    source: &PathComplex,
    d: usize,
    mode: PathBoundaryMode,
    limit: usize,
    tolerance: f64,
    budget: &mut WorkBudget<'_>,
) -> Result<PathRealChain> {
    if limit == 0 || !tolerance.is_finite() || tolerance <= 0.0 || tolerance >= 1e-3 {
        return Err(Error::InvalidParameter {
            parameter: "real path limits/tolerance",
            reason: "positive dense cap and finite tolerance strictly in (0,1e-3) required",
        });
    }
    let layers = d.checked_add(2).ok_or(Error::SizeOverflow {
        operation: "real path chain degrees",
    })?;
    count(layers, limit)?;
    if source.truncated_at().is_some_and(|n| n < layers - 1) {
        return Err(Error::InvalidParameter {
            parameter: "path dimension",
            reason: "real H_d requires Omega_(d+1)",
        });
    }
    let mut paths = Vec::new();
    for n in 0..layers {
        budget.step()?;
        cells(source.paths(n).len(), source.paths(n).len(), limit)?;
        push(&mut paths, copy_paths(source.paths(n), budget)?)?;
    }
    let (mut bases, mut boundaries) = (Vec::<Vec<Vec<f64>>>::new(), Vec::<Vec<Vec<f64>>>::new());
    let mut boundary_residual = 0.0f64;
    for (n, layer) in paths.iter().enumerate() {
        budget.step()?;
        count(layer.len(), limit)?;
        let lower = if n == 0 {
            &[][..]
        } else {
            paths[n - 1].as_slice()
        };
        let mut index = BTreeMap::new();
        for (i, w) in lower.iter().enumerate() {
            budget.step_by(w.len())?;
            index.insert(w.as_slice(), i);
        }
        let mut ambient = Vec::new();
        for w in layer {
            cells(w.len(), w.len(), limit)?;
            push(&mut ambient, path_boundary_budget(w, mode, budget)?)?;
        }
        let mut bad_index = BTreeMap::new();
        for column in &ambient {
            for (face, _) in column {
                budget.step_by(face.len())?;
                if !index.contains_key(face.as_slice()) {
                    let i = bad_index.len();
                    bad_index.entry(face.as_slice()).or_insert(i);
                }
            }
        }
        cells(lower.len(), layer.len(), limit)?;
        cells(bad_index.len(), layer.len(), limit)?;
        let (mut allowed, mut forbidden) = (Vec::new(), Vec::new());
        for column in &ambient {
            let mut a = zeros(lower.len(), limit, budget)?;
            let mut f = zeros(bad_index.len(), limit, budget)?;
            for (face, value) in column {
                budget.step_by(face.len())?;
                if let Some(&i) = index.get(face.as_slice()) {
                    a[i] = *value as f64;
                } else {
                    f[bad_index[face.as_slice()]] = *value as f64;
                }
            }
            push(&mut allowed, a)?;
            push(&mut forbidden, f)?;
        }
        let basis = nullspace(&forbidden, bad_index.len(), tolerance, limit, budget)?;
        let mut boundary = Vec::new();
        for v in &basis {
            let bad = product(&forbidden, v, bad_index.len(), limit, budget)?;
            let bad_zero = zeros(bad_index.len(), limit, budget)?;
            boundary_residual =
                boundary_residual.max(coordinate_residual(&bad, &bad_zero, budget)?);
            let image = product(&allowed, v, lower.len(), limit, budget)?;
            let mut coordinates = Vec::new();
            if n > 0 {
                for q in &bases[n - 1] {
                    push(&mut coordinates, dot(q, &image, budget)?)?;
                }
                let reconstructed =
                    product(&bases[n - 1], &coordinates, lower.len(), limit, budget)?;
                boundary_residual =
                    boundary_residual.max(coordinate_residual(&reconstructed, &image, budget)?);
                if !close(&reconstructed, &image, tolerance, budget)? {
                    return Err(numerical());
                }
                let rows = if n == 1 { 0 } else { bases[n - 2].len() };
                let zero = zeros(rows, limit, budget)?;
                let square = product(&boundaries[n - 1], &coordinates, rows, limit, budget)?;
                boundary_residual =
                    boundary_residual.max(coordinate_residual(&square, &zero, budget)?);
                if !close(&square, &zero, tolerance, budget)? {
                    return Err(numerical());
                }
            }
            push(&mut boundary, coordinates)?;
        }
        cells(
            basis.len(),
            if n == 0 { 0 } else { bases[n - 1].len() },
            limit,
        )?;
        push(&mut bases, basis)?;
        push(&mut boundaries, boundary)?;
    }
    budget.check()?;
    Ok(PathRealChain {
        mode,
        tolerance,
        limit,
        boundary_residual,
        paths,
        bases,
        boundaries,
    })
}

/// Numerical real persistent path Laplacian for one exact source/target event pair.
///
/// Events are sorted zero plus finite nonnegative edge births. `pair` contains
/// indices s<=t. Budgets are (max_paths,max_dense_entries). The upper domain
/// consists of target Omega_(d+1) vectors whose full boundary lies in source
/// Omega_d, with its induced real orthonormal metric. No eigenproblem is solved.
///
/// # Errors
/// Rejects invalid input/pairs/budgets, insufficient numerical rank separation,
/// failed residuals, overflow or allocation failure, without a partial operator.
pub fn compute_persistent_path_laplacian(
    n: usize,
    edges: &[(usize, usize, f64)],
    d: usize,
    pair: (usize, usize),
    mode: PathBoundaryMode,
    budgets: (usize, usize),
    tolerance: f64,
) -> Result<PathRealLaplacian> {
    compute_persistent_path_laplacian_with(
        n,
        edges,
        d,
        pair,
        mode,
        budgets,
        tolerance,
        &Execution::default(),
    )
}
/// Share work through validation, both real chains, embeddings and the induced domain.
///
/// # Errors
/// Returns the errors of [`compute_persistent_path_laplacian`], cancellation or work exhaustion.
#[allow(clippy::too_many_arguments)]
pub fn compute_persistent_path_laplacian_with(
    n: usize,
    edges: &[(usize, usize, f64)],
    d: usize,
    pair: (usize, usize),
    mode: PathBoundaryMode,
    budgets: (usize, usize),
    tolerance: f64,
    execution: &Execution<'_>,
) -> Result<PathRealLaplacian> {
    let (path_limit, limit) = budgets;
    let mut budget = WorkBudget::new(execution)?;
    if path_limit == 0 || limit == 0 {
        return Err(Error::InvalidParameter {
            parameter: "real path limits",
            reason: "must be positive",
        });
    }
    let upper = d.checked_add(1).ok_or(Error::SizeOverflow {
        operation: "real path chain degrees",
    })?;
    let path_count = n.checked_add(edges.len()).ok_or(Error::SizeOverflow {
        operation: "directed path count",
    })?;
    crate::complex::path_count_cap(path_count, path_limit)?;
    count(edges.len(), limit)?;
    let (mut pairs, mut times) = (Vec::new(), Vec::new());
    push(&mut times, 0.0)?;
    for (i, &(a, b, t)) in edges.iter().enumerate() {
        budget.step()?;
        if !t.is_finite() || t < 0.0 {
            return Err(Error::InvalidGraph {
                edge: i,
                reason: "birth must be finite and nonnegative",
            });
        }
        push(&mut pairs, (a, b))?;
        push(&mut times, crate::canonical_zero(t))?;
    }
    let validated = DirectedGraph::new_budget(n, pairs, &mut budget)?;
    budget.check()?;
    times.sort_by(f64::total_cmp);
    budget.check()?;
    times.dedup();
    count(times.len(), limit)?;
    if pair.0 > pair.1 || pair.1 >= times.len() {
        return Err(Error::InvalidParameter {
            parameter: "path event pair",
            reason: "expected source<=target in the complete event list",
        });
    }
    let final_complex = validated.allowed_paths_budget(upper, path_limit, &mut budget)?;
    let mut chains = Vec::new();
    for event in [pair.0, pair.1] {
        let mut pairs = Vec::new();
        for &(a, b, t) in edges {
            budget.step()?;
            if t <= times[event] {
                push(&mut pairs, (a, b))?;
            }
        }
        let source = DirectedGraph::new_budget(n, pairs, &mut budget)?.allowed_paths_budget(
            upper,
            path_limit,
            &mut budget,
        )?;
        push(
            &mut chains,
            real_chain(&source, d, mode, limit, tolerance, &mut budget)?,
        )?;
    }
    let source = &chains[0];
    let target = &chains[1];
    let paths = final_complex.paths(d);
    let mut index = BTreeMap::new();
    for (i, w) in paths.iter().enumerate() {
        budget.step_by(w.len())?;
        index.insert(w.as_slice(), i);
    }
    cells(paths.len(), source.bases[d].len(), limit)?;
    let mut basis = Vec::new();
    for column in &source.bases[d] {
        let mut mapped = zeros(paths.len(), limit, &mut budget)?;
        for (w, &x) in source.paths[d].iter().zip(column) {
            budget.step_by(w.len())?;
            let &i = index.get(w.as_slice()).ok_or(numerical())?;
            mapped[i] = x;
        }
        push(&mut basis, mapped)?;
    }
    cells(paths.len(), target.bases[d + 1].len(), limit)?;
    let (mut residuals, mut coordinates, mut images) = (Vec::new(), Vec::new(), Vec::new());
    for column in &target.boundaries[d + 1] {
        let local_image = product(
            &target.bases[d],
            column,
            target.paths[d].len(),
            limit,
            &mut budget,
        )?;
        let mut image = zeros(paths.len(), limit, &mut budget)?;
        for (path, value) in target.paths[d].iter().zip(local_image) {
            budget.step_by(path.len())?;
            image[index[path.as_slice()]] = value;
        }
        let mut c = Vec::new();
        for q in &basis {
            push(&mut c, dot(q, &image, &mut budget)?)?;
        }
        let projected = product(&basis, &c, paths.len(), limit, &mut budget)?;
        let mut residual = zeros(paths.len(), limit, &mut budget)?;
        for ((r, &x), y) in residual.iter_mut().zip(&image).zip(projected) {
            budget.step()?;
            *r = finite(x - y)?;
        }
        push(&mut residuals, residual)?;
        push(&mut coordinates, c)?;
        push(&mut images, image)?;
    }
    let domain = nullspace(&residuals, paths.len(), tolerance, limit, &mut budget)?;
    let mut up = Vec::new();
    let mut boundary_residual = source.boundary_residual.max(target.boundary_residual);
    for v in &domain {
        let column = product(&coordinates, v, basis.len(), limit, &mut budget)?;
        let image = product(&images, v, paths.len(), limit, &mut budget)?;
        let reconstructed = product(&basis, &column, paths.len(), limit, &mut budget)?;
        boundary_residual =
            boundary_residual.max(coordinate_residual(&image, &reconstructed, &mut budget)?);
        if !close(&image, &reconstructed, tolerance, &mut budget)? {
            return Err(numerical());
        }
        push(&mut up, column)?;
    }
    let operator = PathRealLaplacian {
        mode,
        dimension: d,
        source_time: times[pair.0],
        target_time: times[pair.1],
        paths: copy_paths(paths, &mut budget)?,
        basis,
        down: copy_matrix(&source.boundaries[d], limit, &mut budget)?,
        upper: up,
        limit,
        tolerance,
        boundary_residual,
    };
    budget.check()?;
    Ok(operator)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn ambiguous_rank_and_overflow_are_errors() {
        let mut budget = WorkBudget::new(&Execution::default()).unwrap();
        assert!(
            nullspace(
                &[vec![1.0, 0.0], vec![0.0, 1e-10]],
                2,
                1e-10,
                100,
                &mut budget
            )
            .is_err()
        );
        assert!(dot(&[f64::MAX], &[2.0], &mut budget).is_err());
        assert!(cells(usize::MAX, 2, usize::MAX).is_err());
    }
}
