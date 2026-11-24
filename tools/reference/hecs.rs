//! Unoptimized ordinary q=1 HeCS parity worker for arXiv:2401.15492v1.
//! Standalone experimental code: no library API, persistent adaptation or solver
//! wrapper is changed. Input/limits are documented in tools/README.md.
#![forbid(unsafe_code)]

use std::collections::BTreeSet;
use std::io::{self, Read};
use std::time::Instant;

type Matrix = Vec<Vec<f64>>;
type Pairs = Vec<[usize; 2]>;

fn transpose(a: &Matrix, columns: usize) -> Matrix {
    (0..columns)
        .map(|j| a.iter().map(|row| row[j]).collect())
        .collect()
}

fn product(a: &Matrix, b: &Matrix, columns: usize) -> Matrix {
    a.iter()
        .map(|row| {
            (0..columns)
                .map(|j| row.iter().zip(b).map(|(x, r)| x * r[j]).sum())
                .collect()
        })
        .collect()
}

fn action(a: &Matrix, x: &[f64]) -> Vec<f64> {
    a.iter()
        .map(|row| row.iter().zip(x).map(|(a, b)| a * b).sum())
        .collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn norm(x: &[f64]) -> f64 {
    x.iter().fold(0., |value, &entry| value.hypot(entry))
}

fn eye(n: usize) -> Matrix {
    (0..n)
        .map(|i| (0..n).map(|j| f64::from(i == j)).collect())
        .collect()
}

fn finite(a: &Matrix) -> Result<(), &'static str> {
    if a.iter().flatten().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err("nonfinite arithmetic")
    }
}

fn collapse_scan(incidence: &[[usize; 3]], m: usize, selected: &[usize]) -> (Pairs, Vec<usize>) {
    let mut remaining: BTreeSet<_> = selected.iter().copied().collect();
    let mut pairs = Vec::new();
    while !remaining.is_empty() {
        let free = (0..m).find_map(|e| {
            let mut adjacent = remaining.iter().filter(|&&t| incidence[t].contains(&e));
            let &t = adjacent.next()?;
            adjacent.next().is_none().then_some([e, t])
        });
        let Some(pair) = free else { break };
        remaining.remove(&pair[1]);
        pairs.push(pair);
    }
    (pairs, remaining.into_iter().collect())
}

// Round-3 candidate changes only free-edge storage, preserving the pop policy.
fn collapse_adjacency(
    incidence: &[[usize; 3]],
    m: usize,
    selected: &[usize],
) -> (Pairs, Vec<usize>) {
    let mut adjacent = vec![BTreeSet::new(); m];
    let mut remaining: BTreeSet<_> = selected.iter().copied().collect();
    for &t in selected {
        for &e in &incidence[t] {
            adjacent[e].insert(t);
        }
    }
    let mut free: BTreeSet<_> = (0..m).filter(|&e| adjacent[e].len() == 1).collect();
    let mut pairs = Vec::new();
    while let Some(e) = free.pop_first() {
        let Some(&t) = adjacent[e].first() else {
            continue;
        };
        remaining.remove(&t);
        pairs.push([e, t]);
        for &edge in &incidence[t] {
            adjacent[edge].remove(&t);
            if adjacent[edge].len() == 1 {
                free.insert(edge);
            } else {
                free.remove(&edge);
            }
        }
    }
    (pairs, remaining.into_iter().collect())
}

fn collapse(
    incidence: &[[usize; 3]],
    m: usize,
    selected: &[usize],
    adjacency: bool,
) -> (Pairs, Vec<usize>) {
    if adjacency {
        collapse_adjacency(incidence, m, selected)
    } else {
        collapse_scan(incidence, m, selected)
    }
}

fn heavy(
    incidence: &[[usize; 3]],
    m: usize,
    weights: &[f64],
    peel: bool,
    adjacency: bool,
) -> (Vec<usize>, Pairs) {
    let all: Vec<_> = (0..weights.len()).collect();
    let (initial, core) = collapse(incidence, m, &all, adjacency);
    let mut selected: Vec<_> = if peel {
        initial.iter().map(|p| p[1]).collect()
    } else {
        Vec::new()
    };
    let mut candidates = if peel { core } else { all };
    candidates.sort_by(|&a, &b| weights[b].total_cmp(&weights[a]).then(a.cmp(&b)));
    for t in candidates {
        selected.push(t);
        if !collapse(incidence, m, &selected, adjacency).1.is_empty() {
            selected.pop();
        }
    }
    let (pairs, remaining) = collapse(incidence, m, &selected, adjacency);
    debug_assert!(remaining.is_empty());
    (selected, pairs)
}

/// Selected unit-incidence columns have a unit triangular minor. Solve on it
/// over checked integers, then test ALL ambient rows of ALL original columns.
/// This proves real image coverage exactly, without using an F2 rank or a
/// numerical threshold. A maximal collapsible selection need not cover it.
fn covers(b: &[Vec<i128>], pairs: &Pairs, t: usize) -> Result<bool, &'static str> {
    let r = pairs.len();
    for column in 0..t {
        let mut coefficients = vec![0_i128; r];
        for i in 0..r {
            let [edge, triangle] = pairs[i];
            let mut value = b[edge][column];
            for j in 0..i {
                value = value
                    .checked_sub(
                        b[edge][pairs[j][1]]
                            .checked_mul(coefficients[j])
                            .ok_or("coverage certificate overflow")?,
                    )
                    .ok_or("coverage certificate overflow")?;
            }
            let pivot = b[edge][triangle];
            if pivot != 1 && pivot != -1 {
                return Err("invalid collapse pivot");
            }
            coefficients[i] = value
                .checked_mul(pivot)
                .ok_or("coverage certificate overflow")?;
        }
        for row in b {
            let value = pairs
                .iter()
                .zip(&coefficients)
                .try_fold(0_i128, |v, (p, c)| {
                    v.checked_add(row[p[1]].checked_mul(*c)?)
                });
            if value.ok_or("coverage certificate overflow")? != row[column] {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn cholesky(a: &Matrix) -> Result<Matrix, &'static str> {
    let n = a.len();
    let mut l = vec![vec![0.; n]; n];
    let scale = (0..n).map(|i| a[i][i]).fold(0., f64::max);
    for i in 0..n {
        for j in 0..=i {
            let value = a[i][j] - dot(&l[i][..j], &l[j][..j]);
            if i == j {
                if !value.is_finite() || value <= 64. * n as f64 * f64::EPSILON * scale {
                    return Err("numerically ambiguous active pivot");
                }
                l[i][j] = value.sqrt();
            } else {
                l[i][j] = value / l[j][j];
            }
        }
    }
    finite(&l)?;
    Ok(l)
}

fn solve_cholesky(l: &Matrix, b: &[f64]) -> Vec<f64> {
    let mut x = b.to_vec();
    for i in 0..x.len() {
        x[i] = (x[i] - dot(&l[i][..i], &x[..i])) / l[i][i];
    }
    for i in (0..x.len()).rev() {
        let value: f64 = (i + 1..x.len()).map(|j| l[j][i] * x[j]).sum();
        x[i] = (x[i] - value) / l[i][i];
    }
    x
}

struct Solve {
    status: &'static str,
    x: Vec<f64>,
    trace: Vec<String>,
    solve_ns: u128,
    action_ns: u128,
    actions: usize,
}

fn cgls(a: &Matrix, h: &[f64], u: &Matrix, b: &[f64], lift: &Matrix, n: usize) -> Solve {
    let start = Instant::now();
    let at = transpose(a, n);
    let mut x = vec![0.; n];
    let mut residual = h.to_vec();
    let mut s = action(&at, &residual);
    let mut p = s.clone();
    let mut gamma = dot(&s, &s);
    let mut trace = Vec::new();
    let mut action_ns = 0;
    let mut actions = 0;
    let status = (0..=4000)
        .find_map(|iteration| {
            let original_x = action(lift, &x);
            let original: Vec<_> = action(u, &original_x)
                .iter()
                .zip(b)
                .map(|(a, b)| a - b)
                .collect();
            let fresh = norm(&original) / norm(b).max(1.);
            let mut entry = format!(
                "\"iteration\":{iteration},\"x\":{x:?},\"residual\":{},\"normal_residual\":{},\"original_residual\":{fresh}",
                norm(&residual), gamma.sqrt()
            );
            let stop = if !fresh.is_finite() || !gamma.is_finite() {
                Some("nonfinite")
            } else if fresh <= 1e-9 {
                Some("converged")
            } else if iteration == 4000 {
                Some("iteration_limit")
            } else {
                None
            };
            if stop.is_some() {
                trace.push(format!("{{{entry}}}"));
                return stop;
            }
            let action_start = Instant::now();
            let q = action(a, &p);
            action_ns += action_start.elapsed().as_nanos();
            actions += 1;
            let denominator = dot(&q, &q);
            if !denominator.is_finite() {
                return Some("nonfinite");
            }
            if gamma == 0. || denominator == 0. {
                trace.push(format!("{{{entry}}}"));
                return Some("breakdown");
            }
            let alpha = gamma / denominator;
            for (x, p) in x.iter_mut().zip(&p) {
                *x += alpha * p;
            }
            for (r, q) in residual.iter_mut().zip(q) {
                *r -= alpha * q;
            }
            let action_start = Instant::now();
            s = action(&at, &residual);
            action_ns += action_start.elapsed().as_nanos();
            actions += 1;
            let next_gamma = dot(&s, &s);
            let beta = next_gamma / gamma;
            entry.push_str(&format!(",\"alpha\":{alpha},\"beta\":{beta},\"gamma\":{gamma}"));
            trace.push(format!("{{{entry}}}"));
            for (p, s) in p.iter_mut().zip(&s) {
                *p = s + beta * *p;
            }
            gamma = next_gamma;
            None
        })
        .unwrap_or("iteration_limit");
    Solve {
        status,
        x: action(lift, &x),
        trace,
        solve_ns: start.elapsed().as_nanos(),
        action_ns,
        actions,
    }
}

fn solve_json(solve: Solve, setup_ns: u128) -> Result<String, &'static str> {
    if solve.status == "nonfinite" || !solve.x.iter().all(|v| v.is_finite()) {
        return Err("nonfinite CGLS arithmetic; no successful partial result");
    }
    Ok(format!(
        "{{\"status\":\"{}\",\"original_x\":{:?},\"trace\":[{}],\"setup_ns\":{setup_ns},\"solve_ns\":{},\"iteration_action_ns\":{},\"iteration_actions\":{}}}",
        solve.status,
        solve.x,
        solve.trace.join(","),
        solve.solve_ns,
        solve.action_ns,
        solve.actions
    ))
}

fn run(input: &str, adjacency: bool) -> Result<String, Box<dyn std::error::Error>> {
    run_route(input, adjacency, false)
}

fn weighted_gram(b: &[Vec<i128>], weights: &[f64]) -> Matrix {
    b.iter()
        .map(|left| {
            b.iter()
                .map(|right| {
                    left.iter()
                        .zip(right)
                        .zip(weights)
                        .map(|((&a, &b), &w)| a as f64 * b as f64 * w)
                        .sum()
                })
                .collect()
        })
        .collect()
}

fn run_route(
    input: &str,
    adjacency: bool,
    identity_only: bool,
) -> Result<String, Box<dyn std::error::Error>> {
    let total = Instant::now();
    let mut tokens = input.split_whitespace();
    let mut take = || tokens.next().ok_or("missing input");
    let n: usize = take()?.parse()?;
    let m: usize = take()?.parse()?;
    let t: usize = take()?.parse()?;
    let peel: usize = take()?.parse()?;
    if n > 256 || m > 512 || t > 1024 || peel > 1 {
        return Err("experiment limits: n<=256, m<=512, t<=1024, peel=0/1".into());
    }
    let mut edges = Vec::new();
    for _ in 0..m {
        let edge = [take()?.parse::<usize>()?, take()?.parse::<usize>()?];
        if edge[0] >= edge[1] || edge[1] >= n || edges.contains(&edge) {
            return Err("edges must be unique, increasing vertex pairs".into());
        }
        edges.push(edge);
    }
    let mut triangles = Vec::new();
    let mut weights = Vec::new();
    let mut incidence = Vec::new();
    let mut b = vec![vec![0_i128; t]; m];
    for j in 0..t {
        let face = [
            take()?.parse::<usize>()?,
            take()?.parse::<usize>()?,
            take()?.parse::<usize>()?,
        ];
        let weight: f64 = take()?.parse()?;
        if face[0] >= face[1] || face[1] >= face[2] || face[2] >= n || triangles.contains(&face) {
            return Err("triangles must be unique, increasing vertex triples".into());
        }
        if !weight.is_finite() || weight <= 0. {
            return Err("triangle weights must be finite and strictly positive".into());
        }
        let mut adjacent = [0; 3];
        for removed in 0..3 {
            let edge: Vec<_> = face
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != removed)
                .map(|(_, &v)| v)
                .collect();
            let e = edges
                .iter()
                .position(|x| x.as_slice() == edge)
                .ok_or("missing lower skeleton edge")?;
            b[e][j] = if removed % 2 == 0 { 1 } else { -1 };
            adjacent[removed] = e;
        }
        triangles.push(face);
        incidence.push(adjacent);
        weights.push(weight);
    }
    let rhs = (0..m)
        .map(|_| Ok(take()?.parse::<f64>()?))
        .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?;
    if tokens.next().is_some() || !rhs.iter().all(|x| x.is_finite()) {
        return Err("trailing input or nonfinite RHS".into());
    }
    if identity_only {
        let u = weighted_gram(&b, &weights);
        finite(&u)?;
        let solve = cgls(&u, &rhs, &u, &rhs, &eye(m), m);
        if solve.status != "converged" {
            return Err(format!("unpreconditioned CGLS {}", solve.status).into());
        }
        return Ok(format!(
            "{{\"values\":{:?},\"iterations\":{}}}",
            solve.x,
            solve.trace.len().saturating_sub(1)
        ));
    }
    let preparation_ns = total.elapsed().as_nanos();
    let selection_start = Instant::now();
    let (selected, pairs) = heavy(&incidence, m, &weights, peel == 1, adjacency);
    let selection_ns = selection_start.elapsed().as_nanos();
    let factor_start = Instant::now();
    let rank = pairs.len();
    let edge_permutation: Vec<_> = pairs
        .iter()
        .map(|p| p[0])
        .chain((0..m).filter(|e| !pairs.iter().any(|p| p[0] == *e)))
        .collect();
    let triangle_permutation: Vec<_> = pairs
        .iter()
        .map(|p| p[1])
        .chain((0..t).filter(|t| !pairs.iter().any(|p| p[1] == *t)))
        .collect();
    let f: Matrix = b
        .iter()
        .map(|row| {
            pairs
                .iter()
                .map(|p| row[p[1]] as f64 * weights[p[1]].sqrt())
                .collect()
        })
        .collect();
    let factor: Matrix = edge_permutation.iter().map(|&e| f[e].clone()).collect();
    let u = weighted_gram(&b, &weights);
    finite(&u)?;
    let coverage = covers(&b, &pairs, t)?;
    let mut active = Vec::new();
    let mut hecs = "{\"status\":\"rank_loss\"}".to_string();
    let mut factor_ns = factor_start.elapsed().as_nanos();
    let mut factor_apply_ns = 0;
    let factor_nnz = factor.iter().flatten().filter(|&&x| x != 0.).count();
    if coverage {
        let ft = transpose(&f, rank);
        let gram = product(&ft, &f, rank);
        let l = cholesky(&gram)?;
        let cp = transpose(&f.iter().map(|row| solve_cholesky(&l, row)).collect(), rank);
        let lift = transpose(&cp, m);
        active = product(&product(&cp, &u, m), &lift, rank);
        finite(&active)?;
        cholesky(&active)?; // reject ambiguous active-space positivity
        let apply_start = Instant::now();
        let h = action(&cp, &rhs);
        let probe = action(&lift, &h);
        factor_apply_ns = apply_start.elapsed().as_nanos();
        if !probe.iter().all(|v| v.is_finite()) {
            return Err("nonfinite factor application".into());
        }
        factor_ns = factor_start.elapsed().as_nanos();
        let setup_ns = selection_ns + factor_start.elapsed().as_nanos();
        hecs = solve_json(cgls(&active, &h, &u, &rhs, &lift, rank), setup_ns)?;
    }
    let mut solves = Vec::new();
    for name in ["identity", "diagonal"] {
        let start = Instant::now();
        let scale: Vec<_> = (0..m)
            .map(|i| {
                if name == "diagonal" && u[i][i] > 0. {
                    1. / u[i][i].sqrt()
                } else {
                    1.
                }
            })
            .collect();
        let a: Matrix = (0..m)
            .map(|i| (0..m).map(|j| scale[i] * u[i][j] * scale[j]).collect())
            .collect();
        let mut lift = eye(m);
        for i in 0..m {
            lift[i][i] = scale[i];
        }
        let h: Vec<_> = scale.iter().zip(&rhs).map(|(s, b)| s * b).collect();
        let setup_ns = start.elapsed().as_nanos();
        solves.push(format!(
            "\"{name}\":{}",
            solve_json(cgls(&a, &h, &u, &rhs, &lift, m), setup_ns)?
        ));
    }
    solves.push(format!("\"hecs\":{hecs}"));
    Ok(format!(
        "{{\"selected\":{selected:?},\"pairs\":{pairs:?},\"edge_permutation\":{edge_permutation:?},\"triangle_permutation\":{triangle_permutation:?},\"factor\":{factor:?},\"factor_nnz\":{factor_nnz},\"active\":{active:?},\"preparation_ns\":{preparation_ns},\"selection_ns\":{selection_ns},\"factor_ns\":{factor_ns},\"factor_apply_ns\":{factor_apply_ns},\"total_ns\":{},\"solves\":{{{}}}}}",
        total.elapsed().as_nanos(),
        solves.join(",")
    ))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if !args.is_empty() && args != ["--adjacency"] && args != ["--identity-only"] {
        return Err("optional argument is --adjacency or --identity-only".into());
    }
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    if args == ["--identity-only"] {
        println!("{}", run_route(&input, false, true)?);
    } else {
        println!("{}", run(&input, !args.is_empty())?);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_invalid_and_unbounded_input() {
        for input in [
            "257 0 0 0",
            "2 1 0 0 1 0 0",
            "3 3 1 0 0 1 0 2 1 2 0 1 2 -1 0 0 0",
            "2 1 0 0 0 1 NaN",
        ] {
            assert!(run(input, false).is_err(), "{input}");
        }
    }

    #[test]
    fn empty_zero_rhs_and_inconsistent_rhs() {
        assert!(
            run("0 0 0 0", false)
                .unwrap()
                .contains("\"status\":\"converged\"")
        );
        let zero = cgls(
            &vec![vec![1., 0.], vec![0., 0.]],
            &[0., 0.],
            &eye(2),
            &[0., 0.],
            &eye(2),
            2,
        );
        assert_eq!(zero.status, "converged");
        assert_eq!(zero.trace.len(), 1);
        let inconsistent = cgls(
            &vec![vec![1., 0.], vec![0., 0.]],
            &[0., 1.],
            &eye(2),
            &[0., 1.],
            &eye(2),
            2,
        );
        assert_eq!(inconsistent.status, "breakdown");
    }

    #[test]
    fn coverage_checks_all_ambient_rows() {
        assert!(covers(&[vec![1, 2], vec![-1, -2]], &vec![[0, 0]], 2).unwrap());
        assert!(!covers(&[vec![1, 2], vec![-1, 0]], &vec![[0, 0]], 2).unwrap());
    }

    #[test]
    fn adjacency_preserves_all_shell_collapse_certificates() {
        let incidence = [[0, 1, 3], [0, 2, 4], [1, 2, 5], [3, 4, 5]];
        for mask in 0..16 {
            let selected: Vec<_> = (0..4).filter(|i| mask & (1 << i) != 0).collect();
            let scan = collapse_scan(&incidence, 7, &selected);
            let adjacency = collapse_adjacency(&incidence, 7, &selected);
            assert_eq!(scan, adjacency);
            assert_eq!(scan.1.is_empty(), mask != 15);
        }
    }

    #[test]
    fn rejects_ambiguous_numeric_pivot_without_regularization() {
        assert!(cholesky(&vec![vec![1., 0.], vec![0., 1e-30]]).is_err());
        assert!(cholesky(&vec![vec![f64::INFINITY]]).is_err());
    }
}
