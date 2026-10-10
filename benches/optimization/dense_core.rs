//! Benchmark-only dense critical-set baseline; never linked into the crate.
//! Extracted without formula changes from local research commit
//! 167026e7ee8dbaafe8abaf5e325149af6e93f660 (big_steps.rs).

use std::collections::{BTreeMap, BTreeSet};

type Matrix = Vec<Vec<u8>>; // Columns, then rows; all entries in F2.

#[derive(Clone)]
struct Cell {
    vertices: Vec<usize>,
    value: f64,
}

fn subset(a: &[usize], b: &[usize]) -> bool {
    a.iter().all(|v| b.binary_search(v).is_ok())
}

fn identity(n: usize) -> Matrix {
    (0..n)
        .map(|j| (0..n).map(|i| u8::from(i == j)).collect())
        .collect()
}

fn low(column: &[u8]) -> Option<usize> {
    column.iter().rposition(|&x| x != 0)
}

fn add_column(matrix: &mut Matrix, source: usize, destination: usize) {
    for row in 0..matrix[destination].len() {
        matrix[destination][row] ^= matrix[source][row];
    }
}

struct Reduced {
    r: Matrix,
    v: Matrix,
    u: Matrix,
}

// Algorithm 1, with alpha=1. In U the elementary operation acts on ROWS.
fn reduce(d: &Matrix, rows: usize) -> Reduced {
    let mut r = d.clone();
    let mut v = identity(d.len());
    let mut u = identity(d.len());
    let mut owner = vec![None; rows];
    for j in 0..d.len() {
        while let Some(pivot) = low(&r[j]) {
            if let Some(i) = owner[pivot] {
                add_column(&mut r, i, j);
                add_column(&mut v, i, j);
                for column in &mut u {
                    column[i] ^= column[j];
                }
            } else {
                owner[pivot] = Some(j);
                break;
            }
        }
    }
    Reduced { r, v, u }
}

fn anti_transpose(d: &Matrix, rows: usize) -> Matrix {
    (0..rows)
        .map(|j| {
            (0..d.len())
                .map(|i| d[d.len() - 1 - i][rows - 1 - j])
                .collect()
        })
        .collect()
}

struct Block {
    rows: Vec<usize>,
    columns: Vec<usize>,
    #[cfg(test)]
    d: Matrix,
    primal: Reduced,
    dual: Reduced,
}

struct Analysis {
    blocks: Vec<Block>,
    pairs: Vec<(usize, usize)>,
    essential: Vec<usize>,
}

fn validate(cells: &[Cell]) -> Result<(), String> {
    let mut seen = BTreeSet::new();
    for c in cells {
        if c.vertices.is_empty()
            || !c.value.is_finite()
            || c.vertices.windows(2).any(|w| w[0] >= w[1])
            || !seen.insert(c.vertices.clone())
        {
            return Err("invalid simplex, duplicate or nonfinite value".into());
        }
    }
    for c in cells.iter().filter(|c| c.vertices.len() > 1) {
        for removed in 0..c.vertices.len() {
            let mut face = c.vertices.clone();
            face.remove(removed);
            let f = cells
                .iter()
                .find(|f| f.vertices == face)
                .ok_or("missing face")?;
            if f.value > c.value {
                return Err("face value exceeds coface value".into());
            }
        }
    }
    Ok(())
}

fn analyze(cells: &[Cell]) -> Analysis {
    let mut order: Vec<_> = (0..cells.len()).collect();
    order.sort_by(|&i, &j| {
        cells[i]
            .value
            .total_cmp(&cells[j].value)
            .then(cells[i].vertices.len().cmp(&cells[j].vertices.len()))
            .then(cells[i].vertices.cmp(&cells[j].vertices))
    });
    let top = cells.iter().map(|c| c.vertices.len()).max().unwrap_or(0);
    let mut blocks = Vec::new();
    let mut pairs = Vec::new();
    let mut positive = BTreeSet::new();
    for p in 0..=top {
        let rows: Vec<_> = order
            .iter()
            .copied()
            .filter(|&i| cells[i].vertices.len() == p)
            .collect();
        let columns: Vec<_> = order
            .iter()
            .copied()
            .filter(|&i| cells[i].vertices.len() == p + 1)
            .collect();
        let d: Matrix = columns
            .iter()
            .map(|&j| {
                rows.iter()
                    .map(|&i| u8::from(subset(&cells[i].vertices, &cells[j].vertices)))
                    .collect()
            })
            .collect();
        let primal = reduce(&d, rows.len());
        let dual = reduce(&anti_transpose(&d, rows.len()), columns.len());
        for (j, c) in primal.r.iter().enumerate() {
            if let Some(i) = low(c) {
                pairs.push((rows[i], columns[j]));
            } else {
                positive.insert(columns[j]);
            }
        }
        blocks.push(Block {
            rows,
            columns,
            #[cfg(test)]
            d,
            primal,
            dual,
        });
    }
    for &(birth, _) in &pairs {
        positive.remove(&birth);
    }
    pairs.sort_unstable();
    Analysis {
        blocks,
        pairs,
        essential: positive.into_iter().collect(),
    }
}

fn support(matrix: &Reduced, focus: usize, values: &[f64], target: f64) -> Vec<usize> {
    let start = values[focus];
    if target == start {
        return vec![focus];
    }
    (0..values.len())
        .filter(|&j| {
            if target > start {
                start <= values[j] && values[j] <= target && matrix.u[j][focus] != 0
            } else {
                target <= values[j] && values[j] <= start && matrix.v[focus][j] != 0
            }
        })
        .collect()
}

fn critical(cells: &[Cell], a: &Analysis, id: usize, target: f64, birth: bool) -> Vec<usize> {
    let p = cells[id].vertices.len() - 1;
    let mut result: Vec<usize> = if birth {
        let block = &a.blocks[p + 1];
        let ids: Vec<_> = block.rows.iter().rev().copied().collect();
        let focus = ids.iter().position(|&i| i == id).unwrap();
        let values: Vec<_> = ids.iter().map(|&i| -cells[i].value).collect();
        support(&block.dual, focus, &values, -target)
            .into_iter()
            .map(|j| ids[j])
            .collect()
    } else {
        let block = &a.blocks[p];
        let focus = block.columns.iter().position(|&i| i == id).unwrap();
        let values: Vec<_> = block.columns.iter().map(|&i| cells[i].value).collect();
        support(&block.primal, focus, &values, target)
            .into_iter()
            .map(|j| block.columns[j])
            .collect()
    };
    result.sort_unstable();
    result
}
