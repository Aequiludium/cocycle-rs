//! Adaptive sparse/packed F2 arithmetic with stable coordinate columns.
use crate::execution::WorkBudget;
use crate::{Error, Result};

/// An owned binary matrix with explicit shape, including empty dimensions.
///
/// Entries use column-vector arithmetic over F2. Internal sparse/packed storage
/// is private; callers retain their declared coordinate order.
#[derive(Clone, Debug)]
pub struct F2Matrix {
    pub(super) m: usize,
    pub(super) n: usize,
    words: usize,
    data: Vec<u64>,
    sparse: Option<Vec<Vec<usize>>>,
    nonzeros: usize,
}
impl PartialEq for F2Matrix {
    fn eq(&self, other: &Self) -> bool {
        self.m == other.m
            && self.n == other.n
            && match (&self.sparse, &other.sparse) {
                (Some(a), Some(b)) => a == b,
                (None, None) => self.data == other.data,
                _ => (0..self.m).all(|i| self.indices(i).eq(other.indices(i))),
            }
    }
}
impl Eq for F2Matrix {}

pub(super) fn invalid(reason: &'static str) -> Error {
    Error::InvalidParameter {
        parameter: "homology operator",
        reason,
    }
}
pub(super) fn zeros<T: Clone + Default>(n: usize) -> Result<Vec<T>> {
    let mut v = Vec::new();
    v.try_reserve_exact(n)
        .map_err(|_| Error::AllocationFailed {
            context: "F2 matrix",
        })?;
    v.resize(n, T::default());
    Ok(v)
}
impl F2Matrix {
    /// Create a zero matrix. Shape multiplication is checked before allocation.
    pub fn zero(nrows: usize, ncols: usize) -> Result<Self> {
        let words = ncols.div_ceil(64);
        let size = nrows.checked_mul(words).ok_or(Error::SizeOverflow {
            operation: "F2 matrix storage",
        })?;
        Ok(Self {
            m: nrows,
            n: ncols,
            words,
            data: if words >= 8 { Vec::new() } else { zeros(size)? },
            sparse: if words >= 8 {
                Some(zeros(nrows)?)
            } else {
                None
            },
            nonzeros: 0,
        })
    }
    /// Read equally sized binary rows; ncols preserves the shape of no rows.
    /// Nonbinary entries and inconsistent row lengths are errors.
    pub fn from_rows(rows: &[Vec<u8>], ncols: usize) -> Result<Self> {
        let mut matrix = Self::zero(rows.len(), ncols)?;
        for (i, row) in rows.iter().enumerate() {
            check_vector(row, ncols)?;
            for (j, &bit) in row.iter().enumerate() {
                if bit != 0 {
                    matrix.set(i, j)?;
                }
            }
        }
        Ok(matrix)
    }
    /// Read binary column vectors, preserving nrows for no columns.
    pub fn from_columns(columns: &[Vec<u8>], nrows: usize) -> Result<Self> {
        let mut matrix = Self::zero(nrows, columns.len())?;
        for (j, column) in columns.iter().enumerate() {
            check_vector(column, nrows)?;
            for (i, &bit) in column.iter().enumerate() {
                if bit != 0 {
                    matrix.set(i, j)?;
                }
            }
        }
        Ok(matrix)
    }
    /// Read nonzero row positions per column without dense coordinate buffers.
    /// Empty columns are retained. Duplicate or out-of-range indices are errors;
    /// each column may list its distinct indices in any order.
    pub fn from_sparse_columns(columns: &[Vec<usize>], nrows: usize) -> Result<Self> {
        let mut matrix = Self::zero(nrows, columns.len())?;
        for (j, column) in columns.iter().enumerate() {
            for &i in column {
                if i >= nrows || matrix.bit(i, j) {
                    return Err(invalid(
                        "sparse column indices must be distinct and in range",
                    ));
                }
                matrix.set(i, j)?;
            }
        }
        Ok(matrix)
    }
    /// Number of rows in the declared basis.
    pub fn nrows(&self) -> usize {
        self.m
    }
    /// Number of columns in the declared basis.
    pub fn ncols(&self) -> usize {
        self.n
    }
    /// Read one logical row; an out-of-range position is an error.
    pub fn row(&self, row: usize) -> Result<Vec<u8>> {
        if row >= self.m {
            return Err(invalid("matrix row is out of range"));
        }
        let mut output = zeros(self.n)?;
        for j in self.indices(row) {
            output[j] = 1;
        }
        Ok(output)
    }
    /// Read one logical column; an out-of-range position is an error.
    pub fn column(&self, column: usize) -> Result<Vec<u8>> {
        if column >= self.n {
            return Err(invalid("matrix column is out of range"));
        }
        Ok((0..self.m).map(|i| u8::from(self.bit(i, column))).collect())
    }
    /// Apply this matrix to a binary column vector.
    pub fn apply(&self, vector: &[u8]) -> Result<Vec<u8>> {
        self.apply_with(vector, &mut WorkBudget::unlimited())
    }
    pub(super) fn apply_with<const C: bool>(
        &self,
        vector: &[u8],
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Vec<u8>> {
        let packed = pack(vector, self.n, budget)?;
        Ok(unpack(&self.apply_packed(&packed, budget)?, self.m))
    }
    pub(super) fn apply_packed<const C: bool>(
        &self,
        vector: &[u64],
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Vec<u64>> {
        let mut output = zeros::<u64>(self.m.div_ceil(64))?;
        for i in 0..self.m {
            let bit = if let Some(rows) = &self.sparse {
                budget.step_by(rows[i].len().saturating_add(1))?;
                rows[i]
                    .iter()
                    .fold(0, |p, &j| p ^ ((vector[j / 64] >> (j % 64)) & 1))
            } else {
                budget.step_by(self.words.saturating_add(1))?;
                self.packed_row(i)
                    .iter()
                    .zip(vector)
                    .fold(0, |p, (a, b)| p ^ u64::from((a & b).count_ones() & 1))
            };
            output[i / 64] |= bit << (i % 64);
        }
        Ok(output)
    }
    /// Exact binary rank, without reordering columns.
    pub fn rank(&self) -> Result<usize> {
        Ok(self.reduce(&mut WorkBudget::unlimited())?.1.len())
    }
    /// Kernel generators in increasing free-coordinate order.
    pub fn kernel_basis(&self) -> Result<Vec<Vec<u8>>> {
        self.kernel_with(&mut WorkBudget::unlimited())
    }
    pub(super) fn kernel_with<const C: bool>(
        &self,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Vec<Vec<u8>>> {
        self.prepare(false, budget)?.kernel(budget)
    }
    /// Solve M x = b with free coordinates zero. None means inconsistent;
    /// a valid solution of dimension zero is Some(vec![]).
    pub fn solve(&self, b: &[u8]) -> Result<Option<Vec<u8>>> {
        check_vector(b, self.m)?;
        let mut budget = WorkBudget::unlimited();
        self.prepare(true, &mut budget)?.solve(b, &mut budget)
    }
    pub(super) fn bit(&self, i: usize, j: usize) -> bool {
        self.sparse.as_ref().map_or_else(
            || self.data[i * self.words + j / 64] >> (j % 64) & 1 != 0,
            |rows| rows[i].binary_search(&j).is_ok(),
        )
    }
    pub(super) fn set(&mut self, i: usize, j: usize) -> Result<()> {
        if !self.bit(i, j) {
            self.toggle(i, j)?;
        }
        Ok(())
    }
    pub(super) fn toggle(&mut self, i: usize, j: usize) -> Result<()> {
        if let Some(rows) = &mut self.sparse {
            match rows[i].binary_search(&j) {
                Ok(k) => {
                    rows[i].remove(k);
                    self.nonzeros -= 1;
                }
                Err(k) => {
                    rows[i]
                        .try_reserve(1)
                        .map_err(|_| Error::AllocationFailed {
                            context: "sparse F2 row",
                        })?;
                    rows[i].insert(k, j);
                    self.nonzeros += 1;
                }
            }
            self.densify()?;
        } else {
            self.data[i * self.words + j / 64] ^= 1 << (j % 64);
        }
        Ok(())
    }
    fn densify(&mut self) -> Result<()> {
        if self.sparse.is_some() && self.nonzeros > self.m * self.words / 2 {
            let mut data = zeros::<u64>(self.m * self.words)?;
            for (i, row) in self.sparse.as_ref().into_iter().flatten().enumerate() {
                for &j in row {
                    data[i * self.words + j / 64] |= 1 << (j % 64);
                }
            }
            self.data = data;
            self.sparse = None;
        }
        Ok(())
    }
    fn packed_row(&self, i: usize) -> &[u64] {
        &self.data[i * self.words..(i + 1) * self.words]
    }
    fn indices(&self, i: usize) -> impl Iterator<Item = usize> + '_ {
        let sparse = self.sparse.as_ref().map(|rows| rows[i].iter().copied());
        let dense = (self.sparse.is_none()).then(|| {
            self.packed_row(i)
                .iter()
                .copied()
                .enumerate()
                .flat_map(|(w, mut bits)| {
                    std::iter::from_fn(move || {
                        if bits == 0 {
                            return None;
                        }
                        let j = w * 64 + bits.trailing_zeros() as usize;
                        bits &= bits - 1;
                        Some(j)
                    })
                })
        });
        sparse
            .into_iter()
            .flatten()
            .chain(dense.into_iter().flatten())
    }
    fn swap_rows(&mut self, i: usize, j: usize) {
        if let Some(rows) = &mut self.sparse {
            rows.swap(i, j);
        } else {
            for w in 0..self.words {
                self.data.swap(i * self.words + w, j * self.words + w);
            }
        }
    }
    fn xor_row<const C: bool>(
        &mut self,
        i: usize,
        j: usize,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<()> {
        if let Some(rows) = &mut self.sparse {
            let (a, b) = (&rows[i], &rows[j]);
            budget.step_by(a.len().saturating_add(b.len()))?;
            let mut merged = Vec::new();
            merged
                .try_reserve(a.len().saturating_add(b.len()))
                .map_err(|_| Error::AllocationFailed {
                    context: "sparse F2 elimination",
                })?;
            let (mut x, mut y) = (0, 0);
            while x < a.len() || y < b.len() {
                if y == b.len() || (x < a.len() && a[x] < b[y]) {
                    merged.push(a[x]);
                    x += 1;
                } else if x == a.len() || b[y] < a[x] {
                    merged.push(b[y]);
                    y += 1;
                } else {
                    x += 1;
                    y += 1;
                }
            }
            self.nonzeros = self.nonzeros - rows[i].len() + merged.len();
            rows[i] = merged;
            self.densify()?;
        } else {
            budget.step_by(self.words)?;
            for w in 0..self.words {
                self.data[i * self.words + w] ^= self.data[j * self.words + w];
            }
        }
        Ok(())
    }
    pub(super) fn product<const C: bool>(
        &self,
        other: &Self,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Self> {
        if self.n != other.m {
            return Err(invalid("F2 product dimensions disagree"));
        }
        let mut output = Self::zero(self.m, other.n)?;
        let mut row = zeros::<u64>(other.words)?;
        for i in 0..self.m {
            budget.step_by(other.words.saturating_add(1))?;
            row.fill(0);
            for j in self.indices(i) {
                if let Some(rows) = &other.sparse {
                    budget.step_by(rows[j].len().saturating_add(1))?;
                    for &k in &rows[j] {
                        row[k / 64] ^= 1 << (k % 64);
                    }
                } else {
                    budget.step_by(other.words.saturating_add(1))?;
                    for (x, y) in row.iter_mut().zip(other.packed_row(j)) {
                        *x ^= y;
                    }
                }
            }
            if output.sparse.is_none() {
                output.data[i * other.words..(i + 1) * other.words].copy_from_slice(&row);
            } else {
                for (w, &bits) in row.iter().enumerate() {
                    let mut bits = bits;
                    while bits != 0 {
                        output.set(i, w * 64 + bits.trailing_zeros() as usize)?;
                        bits &= bits - 1;
                    }
                }
            }
        }
        Ok(output)
    }
    pub(super) fn reduce<const C: bool>(
        &self,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<(Self, Vec<usize>)> {
        let prepared = self.prepare(false, budget)?;
        Ok((prepared.reduced, prepared.pivots))
    }
    pub(super) fn prepare<const C: bool>(
        &self,
        full: bool,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Decomposition> {
        budget.check()?;
        let mut reduced = self.clone();
        let mut transform = if full {
            Some(Self::zero(self.m, self.m)?)
        } else {
            None
        };
        if let Some(e) = &mut transform {
            for i in 0..self.m {
                budget.step()?;
                e.set(i, i)?;
            }
        }
        let mut pivots = Vec::new();
        for column in 0..self.n {
            budget.step()?;
            let row = pivots.len();
            let mut selected = None;
            for i in row..self.m {
                budget.step()?;
                if reduced.bit(i, column) {
                    selected = Some(i);
                    break;
                }
            }
            let Some(selected) = selected else {
                continue;
            };
            budget.step_by(reduced.words)?;
            reduced.swap_rows(row, selected);
            if let Some(e) = &mut transform {
                e.swap_rows(row, selected);
            }
            for i in 0..self.m {
                budget.step()?;
                if i != row && reduced.bit(i, column) {
                    reduced.xor_row(i, row, budget)?;
                    if let Some(e) = &mut transform {
                        e.xor_row(i, row, budget)?;
                    }
                }
            }
            pivots.push(column);
            if pivots.len() == self.m {
                break;
            }
        }
        Ok(Decomposition {
            reduced,
            transform,
            pivots,
        })
    }
}
pub(super) fn check_vector(x: &[u8], n: usize) -> Result<()> {
    if x.len() != n {
        return Err(Error::ShapeMismatch {
            input: "F2 vector",
            expected: n,
            actual: x.len(),
        });
    }
    if x.iter().any(|&bit| bit > 1) {
        return Err(invalid("F2 coordinates must be 0 or 1"));
    }
    Ok(())
}
pub(super) fn pack<const C: bool>(
    x: &[u8],
    n: usize,
    budget: &mut WorkBudget<'_, C>,
) -> Result<Vec<u64>> {
    check_vector(x, n)?;
    let mut output = zeros::<u64>(n.div_ceil(64))?;
    for (j, &bit) in x.iter().enumerate() {
        budget.step()?;
        output[j / 64] |= u64::from(bit) << (j % 64);
    }
    Ok(output)
}
pub(super) fn unpack(x: &[u64], n: usize) -> Vec<u8> {
    (0..n)
        .map(|j| ((x[j / 64] >> (j % 64)) & 1) as u8)
        .collect()
}
/// A request-local stable RREF; rank/kernel-only work has no row transform.
pub(super) struct Decomposition {
    reduced: F2Matrix,
    transform: Option<F2Matrix>,
    pivots: Vec<usize>,
}
impl Decomposition {
    pub(super) fn kernel<const C: bool>(
        &self,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Vec<Vec<u8>>> {
        let mut vectors = Vec::new();
        for free in 0..self.reduced.n {
            budget.step()?;
            if self.pivots.binary_search(&free).is_ok() {
                continue;
            }
            let mut x = zeros::<u8>(self.reduced.n)?;
            x[free] = 1;
            for (i, &pivot) in self.pivots.iter().enumerate() {
                budget.step()?;
                x[pivot] = u8::from(self.reduced.bit(i, free));
            }
            vectors.push(x);
        }
        Ok(vectors)
    }
    pub(super) fn solve<const C: bool>(
        &self,
        b: &[u8],
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Option<Vec<u8>>> {
        let e = self
            .transform
            .as_ref()
            .ok_or(invalid("coordinate solve needs a full decomposition"))?;
        let y = e.apply_with(b, budget)?;
        if y[self.pivots.len()..].iter().any(|&bit| bit != 0) {
            return Ok(None);
        }
        let mut x = zeros::<u8>(self.reduced.n)?;
        for (i, &pivot) in self.pivots.iter().enumerate() {
            budget.step()?;
            x[pivot] = y[i];
        }
        Ok(Some(x))
    }
    pub(super) fn inverse<const C: bool>(&self, budget: &mut WorkBudget<'_, C>) -> Result<Inverse> {
        let mut e = self
            .transform
            .as_ref()
            .ok_or(invalid("inverse needs a full decomposition"))?
            .clone();
        // Continue through the right block of [M|I]. Null coefficient rows
        // change upper inverse rows and therefore the action on noncycles.
        let mut row = self.pivots.len();
        for column in 0..e.n {
            budget.step()?;
            let mut selected = None;
            for i in row..e.m {
                budget.step()?;
                if e.bit(i, column) {
                    selected = Some(i);
                    break;
                }
            }
            let Some(selected) = selected else {
                continue;
            };
            e.swap_rows(row, selected);
            for i in 0..e.m {
                budget.step()?;
                if i != row && e.bit(i, column) {
                    e.xor_row(i, row, budget)?;
                }
            }
            row += 1;
            if row == e.m {
                break;
            }
        }
        let mut rows = F2Matrix::zero(self.pivots.len(), e.n)?;
        for i in 0..self.pivots.len() {
            for j in e.indices(i) {
                budget.step()?;
                rows.set(i, j)?;
            }
        }
        Ok(Inverse {
            rows,
            pivots: self.pivots.clone(),
            size: self.reduced.n,
        })
    }
}
/// Private generalized inverse, storing only nonzero pivot rows.
#[derive(Clone, Debug)]
pub(super) struct Inverse {
    rows: F2Matrix,
    pivots: Vec<usize>,
    size: usize,
}
impl Inverse {
    pub(super) fn new<const C: bool>(
        matrix: &F2Matrix,
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Self> {
        matrix.prepare(true, budget)?.inverse(budget)
    }
    pub(super) fn apply_packed<const C: bool>(
        &self,
        x: &[u64],
        budget: &mut WorkBudget<'_, C>,
    ) -> Result<Vec<u64>> {
        let coefficients = self.rows.apply_packed(x, budget)?;
        let mut output = zeros::<u64>(self.size.div_ceil(64))?;
        for (i, &pivot) in self.pivots.iter().enumerate() {
            budget.step()?;
            output[pivot / 64] |= ((coefficients[i / 64] >> (i % 64)) & 1) << (pivot % 64);
        }
        Ok(output)
    }
    pub(super) fn rank(&self) -> usize {
        self.pivots.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sparse_fill_switch_preserves_coordinates_and_equality() {
        let mut columns = vec![vec![]; 513];
        columns[0] = vec![0];
        columns[64] = vec![1];
        columns[512] = vec![2];
        let matrix = F2Matrix::from_sparse_columns(&columns, 3).unwrap();
        assert!(matrix.sparse.is_some());
        let mut dense = matrix.clone();
        for j in 0..20 {
            dense.toggle(0, j).unwrap();
        }
        for j in 0..20 {
            dense.toggle(0, j).unwrap();
        }
        assert!(dense.sparse.is_none());
        assert_eq!(matrix, dense);
        assert_eq!(matrix.apply(&vec![1; 513]).unwrap(), vec![1, 1, 1]);
        assert_eq!(
            matrix.reduce(&mut WorkBudget::unlimited()).unwrap().0,
            dense.reduce(&mut WorkBudget::unlimited()).unwrap().0
        );
        assert!(F2Matrix::from_sparse_columns(&[vec![1, 1]], 3).is_err());
        assert!(F2Matrix::from_sparse_columns(&[vec![3]], 3).is_err());
    }
    #[test]
    fn prepared_solves_match_independent_small_enumeration() {
        for code in 0..64 {
            let rows: Vec<_> = (0..2)
                .map(|i| (0..3).map(|j| ((code >> (3 * i + j)) & 1) as u8).collect())
                .collect();
            let matrix = F2Matrix::from_rows(&rows, 3).unwrap();
            let prepared = matrix.prepare(true, &mut WorkBudget::unlimited()).unwrap();
            assert!(
                matrix
                    .prepare(false, &mut WorkBudget::unlimited())
                    .unwrap()
                    .transform
                    .is_none()
            );
            for b in 0..4 {
                let rhs = vec![(b & 1) as u8, ((b >> 1) & 1) as u8];
                let possible = (0..8).any(|x| {
                    (0..2).all(|i| {
                        (0..3).fold(0, |p, j| p ^ (rows[i][j] & ((x >> j) & 1) as u8)) == rhs[i]
                    })
                });
                let solution = prepared.solve(&rhs, &mut WorkBudget::unlimited()).unwrap();
                assert_eq!(solution.is_some(), possible);
                if let Some(x) = solution {
                    assert_eq!(matrix.apply(&x).unwrap(), rhs);
                }
            }
        }
    }
}
