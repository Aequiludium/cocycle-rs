//! Matched unaccelerated/alternative spectral workflows; benchmark-only adapter.
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::persistent_laplacian::{
    HarmonicRestriction, IterativeOptions, OperatorOptions, PersistentLaplacian,
    PersistentLaplacianOperator, PersistentSpectrumTarget, SpectralOptions,
    fixed_endpoint_eigenpairs, gram_eigenpairs, gram_spectrum, relevant_pair_eigenpairs,
    relevant_pair_laplacians,
};
pub use cocycle::{Error, Result};
use std::io::{self, Read};
use std::time::Instant;

#[allow(dead_code)]
#[path = "../../src/eigensolver.rs"]
mod eigensolver;
#[allow(dead_code)]
#[path = "../../src/execution/mod.rs"]
mod execution;

type Trial = std::result::Result<(Vec<f64>, u128, u128), Box<dyn std::error::Error>>;

fn canonical_zero(value: f64) -> f64 {
    if value == 0. { 0. } else { value }
}

fn options() -> SpectralOptions {
    SpectralOptions {
        max_dense_entries: 4_000_000,
        ..Default::default()
    }
}

fn action(route: &str, n: usize, actions: usize) -> Trial {
    let start = Instant::now();
    let mut cells: Vec<_> = (0..n)
        .map(|v| Simplex::new(vec![v], 0.))
        .collect::<Result<_>>()?;
    cells.push(Simplex::new(vec![n], 1.)?);
    for v in 0..n {
        cells.push(Simplex::new(vec![v, n], 2.)?);
    }
    let source = SimplicialComplex::new(cells)?;
    let k = source.stage(0.)?;
    let l = source.stage(2.)?;
    let execution = Execution::default();
    let x: Vec<_> = (0..n).map(|i| (i % 17) as f64 / 17.).collect();
    let mean = x.iter().sum::<f64>() / n as f64;
    let dense = if route == "dense" {
        Some(PersistentLaplacian::compute(
            &k,
            &l,
            0,
            options(),
            &execution,
        )?)
    } else {
        None
    };
    let sparse = if route == "factorized" {
        Some(PersistentLaplacianOperator::new(
            &k,
            &l,
            0,
            OperatorOptions {
                spectral: options(),
                ..Default::default()
            },
            &execution,
        )?)
    } else {
        None
    };
    if dense.is_none() && sparse.is_none() {
        return Err("action route must be dense or factorized".into());
    }
    let prepare_ns = start.elapsed().as_nanos();
    let start = Instant::now();
    let mut values = Vec::new();
    for _ in 0..actions {
        values = if let Some(operator) = &dense {
            operator.action(&x, &execution)?
        } else {
            sparse.as_ref().unwrap().action(&x, &execution)?
        };
        if values
            .iter()
            .zip(&x)
            .any(|(a, b)| (a - (b - mean)).abs() > 1e-9)
        {
            return Err("star action failed independent I-11^T/n control".into());
        }
        std::hint::black_box(&values);
    }
    Ok((values, prepare_ns, start.elapsed().as_nanos()))
}

fn full_modes(h: &[f64], n: usize, rank: usize) -> Result<(Vec<f64>, Vec<f64>)> {
    let mut budget = execution::WorkBudget::new(&execution::Execution::default())?;
    let (values, vectors, _) = eigensolver::jacobi(h, n, 0., 4. * f64::EPSILON, 256, &mut budget)?;
    let mut order: Vec<_> = (0..n).collect();
    order.sort_by(|&i, &j| values[i].total_cmp(&values[j]));
    let mut positive = Vec::new();
    let mut basis = vec![0.; n * rank];
    for (j, &index) in order[n - rank..].iter().enumerate() {
        positive.push(values[index]);
        for i in 0..n {
            basis[i * rank + j] = vectors[i * n + index];
        }
    }
    Ok((positive, basis))
}

fn spectral(route: &str, kind: &str, n: usize) -> Trial {
    if !n.is_power_of_two() || n < 8 {
        return Err("spectral control requires a power of two >= 8".into());
    }
    // Independent Hadamard columns: exact eigenvalues 1, 4, 9, all others zero.
    let start = Instant::now();
    let rank = 3;
    let entry = |i: usize, j: usize| {
        let sign = if (i & j).count_ones() % 2 == 0 {
            1.
        } else {
            -1.
        };
        sign / (n as f64).sqrt()
    };
    let mut a = vec![0.; rank * n];
    for j in 0..rank {
        for i in 0..n {
            a[j * n + i] = (j + 1) as f64 * entry(i, j);
        }
    }
    let kernel: Vec<_> = if kind == "harmonic" {
        (0..n)
            .flat_map(|i| (rank..n).map(move |j| entry(i, j)))
            .collect()
    } else {
        vec![]
    };
    let make_h = || {
        (0..n)
            .flat_map(|i| {
                let a = &a;
                (0..n).map(move |j| (0..rank).map(|k| a[k * n + i] * a[k * n + j]).sum())
            })
            .collect::<Vec<f64>>()
    };
    let supplied_h = if kind == "harmonic" { make_h() } else { vec![] };
    let prepare_ns = start.elapsed().as_nanos();
    let start = Instant::now();
    let (mut values, vectors) = if route == "full" {
        if kind == "harmonic" {
            full_modes(&supplied_h, n, rank)?
        } else {
            full_modes(&make_h(), n, rank)?
        }
    } else if kind == "gram" && route == "both-grams" {
        let result = gram_spectrum(&a, rank, n, options(), &Execution::default())?;
        (result.positive, result.lifted_vectors)
    } else if kind == "gram" && route == "single-gram" {
        let result = gram_eigenpairs(&a, rank, n, options(), &Execution::default())?;
        (result.eigenvalues[n - rank..].to_vec(), result.vectors)
    } else if kind == "harmonic" && route == "restriction" {
        let result = HarmonicRestriction::compute(
            &supplied_h,
            n,
            &kernel,
            n - rank,
            options(),
            &Execution::default(),
        )?;
        (result.eigenvalues().to_vec(), result.vectors().to_vec())
    } else {
        return Err("invalid spectral route".into());
    };
    if values.len() != rank
        || values
            .iter()
            .zip([1., 4., 9.])
            .any(|(a, b)| (a - b).abs() > 1e-9)
    {
        return Err("spectral control eigenvalues disagree".into());
    }
    for i in 0..n {
        for j in 0..n {
            let projector: f64 = (0..rank)
                .map(|k| vectors[i * rank + k] * vectors[j * rank + k])
                .sum();
            let expected: f64 = (0..rank).map(|k| entry(i, k) * entry(j, k)).sum();
            if (projector - expected).abs() > 1e-9 {
                return Err("spectral control projector disagrees".into());
            }
            values.push(projector);
        }
    }
    Ok((values, prepare_ns, start.elapsed().as_nanos()))
}

fn next<T: std::str::FromStr>(
    words: &mut std::str::SplitWhitespace<'_>,
) -> std::result::Result<T, Box<dyn std::error::Error>> {
    words
        .next()
        .ok_or("missing fixture field")?
        .parse()
        .map_err(|_| "invalid fixture field".into())
}

fn pairs(input: &str, kind: &str, route: &str, spectra: bool, single: bool) -> Trial {
    let mut words = input.split_whitespace();
    let mut values = Vec::new();
    let mut prepare_ns = 0;
    let mut solve_ns = 0;
    while words.clone().next().is_some() {
        let start = Instant::now();
        let q: usize = next(&mut words)?;
        let count: usize = next(&mut words)?;
        let cells_count: usize = next(&mut words)?;
        let end: f64 = next(&mut words)?;
        let scales = (0..count)
            .map(|_| next(&mut words))
            .collect::<std::result::Result<Vec<f64>, _>>()?;
        let mut cells = Vec::new();
        for _ in 0..cells_count {
            let scale = next(&mut words)?;
            let size = next(&mut words)?;
            let vertices = (0..size)
                .map(|_| next(&mut words))
                .collect::<std::result::Result<Vec<usize>, _>>()?;
            cells.push(Simplex::new(vertices, scale)?);
        }
        let source = SimplicialComplex::new(cells)?;
        let stages = scales
            .iter()
            .map(|&s| source.stage(s))
            .collect::<Result<Vec<_>>>()?;
        let target = source.stage(end)?;
        let requests: Vec<_> = if single {
            if stages.is_empty() {
                vec![]
            } else {
                vec![(stages.len() - 1, stages.len() - 1)]
            }
        } else {
            (0..stages.len())
                .flat_map(|t| (0..=t).map(move |s| (s, t)))
                .collect()
        };
        let execution = Execution::default();
        prepare_ns += start.elapsed().as_nanos();
        let start = Instant::now();
        if kind == "partial-pairs" {
            let iteration = IterativeOptions {
                absolute_tolerance: 1e-10,
                relative_tolerance: 1e-10,
                ..Default::default()
            };
            let operator_options = OperatorOptions {
                spectral: options(),
                max_projection_entries: 4_000_000,
                ..Default::default()
            };
            let nonempty = |stage: &cocycle::filtration::SimplicialStage<'_>| {
                stage.simplices().iter().any(|s| s.dimension() == q)
            };
            if route == "relevant" || route == "exhaustive" {
                let requests: Vec<_> = requests
                    .into_iter()
                    .filter(|&(s, _)| nonempty(&stages[s]))
                    .collect();
                if route == "relevant" {
                    for result in relevant_pair_eigenpairs(
                        &stages,
                        &requests,
                        q,
                        1,
                        PersistentSpectrumTarget::Largest,
                        operator_options,
                        &iteration,
                        true,
                        &execution,
                    )? {
                        values.extend(result.values);
                    }
                } else {
                    for (s, t) in requests {
                        values.push(
                            *PersistentLaplacian::compute(
                                &stages[s],
                                &stages[t],
                                q,
                                options(),
                                &execution,
                            )?
                            .spectrum(&execution)?
                            .eigenvalues
                            .last()
                            .ok_or("empty requested spectrum")?,
                        );
                    }
                }
            } else {
                let starts: Vec<_> = stages.iter().copied().filter(nonempty).collect();
                if route == "direct-fixed" {
                    for stage in &starts {
                        values.push(
                            *PersistentLaplacian::compute(
                                stage,
                                &target,
                                q,
                                options(),
                                &execution,
                            )?
                            .spectrum(&execution)?
                            .eigenvalues
                            .last()
                            .ok_or("empty requested spectrum")?,
                        );
                    }
                } else if route == "fixed-cold" || route == "fixed-warm" {
                    for result in fixed_endpoint_eigenpairs(
                        &starts,
                        &target,
                        q,
                        1,
                        PersistentSpectrumTarget::Largest,
                        operator_options,
                        &iteration,
                        route == "fixed-warm",
                        &execution,
                    )? {
                        values.extend(result.values);
                    }
                } else {
                    return Err("invalid partial batch route".into());
                }
            }
            solve_ns += start.elapsed().as_nanos();
            continue;
        }
        if kind == "pairs" && route == "batch" {
            let operators = relevant_pair_laplacians(&stages, &requests, q, options(), &execution)?;
            for operator in &operators {
                append_operator(operator.operator(), &mut values, spectra, &execution)?;
            }
            solve_ns += start.elapsed().as_nanos();
            continue;
        }
        let operators = if kind == "pairs" && route == "exhaustive" {
            requests
                .iter()
                .map(|&(s, t)| {
                    PersistentLaplacian::compute(&stages[s], &stages[t], q, options(), &execution)
                })
                .collect::<Result<Vec<_>>>()?
        } else if kind == "nonbranching" {
            stages
                .iter()
                .map(|stage| match route {
                    "schur" => {
                        PersistentLaplacian::compute(stage, &target, q, options(), &execution)
                    }
                    "specialized" => PersistentLaplacian::compute_non_branching(
                        stage,
                        &target,
                        q,
                        options(),
                        &execution,
                    ),
                    _ => Err(Error::InvalidParameter {
                        parameter: "benchmark route",
                        reason: "unknown non-branching route",
                    }),
                })
                .collect::<Result<Vec<_>>>()?
        } else {
            return Err("invalid pair route".into());
        };
        for operator in operators {
            append_operator(&operator, &mut values, spectra, &execution)?;
        }
        solve_ns += start.elapsed().as_nanos();
    }
    Ok((values, prepare_ns, solve_ns))
}

fn partial(route: &str, kind: &str, n: usize) -> Trial {
    let start = Instant::now();
    let mut cells = Vec::new();
    for v in 0..n {
        cells.push(Simplex::new(
            vec![v],
            if kind == "compact" { v as f64 } else { 0. },
        )?);
    }
    let (s, t, positive) = match kind {
        "partial" => {
            for j in 0..3 {
                cells.push(Simplex::new(vec![2 * j, 2 * j + 1], 1.)?);
            }
            (1., 1., true)
        }
        "compact" => {
            for v in 0..n - 1 {
                cells.push(Simplex::new(vec![v, v + 1], n as f64)?);
            }
            (1., n as f64, true)
        }
        "partial-star" | "partial-star-positive" => {
            cells.push(Simplex::new(vec![n], 1.)?);
            for v in 0..n {
                cells.push(Simplex::new(vec![v, n], 2.)?);
            }
            (0., 2., kind == "partial-star-positive")
        }
        _ => return Err("unknown partial fixture".into()),
    };
    let source = SimplicialComplex::new(cells)?;
    let (k, l) = (source.stage(s)?, source.stage(t)?);
    let execution = Execution::default();
    let prepare_ns = start.elapsed().as_nanos();
    let start = Instant::now();
    let count = if kind == "compact" { 1 } else { 3 };
    let values = if route == "dense" {
        let full = PersistentLaplacian::compute(&k, &l, 0, options(), &execution)?
            .spectrum(&execution)?
            .eigenvalues;
        if positive {
            full.into_iter()
                .filter(|v| *v > 1e-8)
                .take(count)
                .collect::<Vec<_>>()
        } else {
            full.into_iter().rev().take(count).collect::<Vec<_>>()
        }
    } else {
        let policy = OperatorOptions {
            spectral: options(),
            max_projection_entries: 4_000_000,
            ..Default::default()
        };
        let op = if route == "compact" {
            PersistentLaplacianOperator::new_non_branching(&k, &l, 0, policy, &execution)?
        } else if route == "factorized" {
            PersistentLaplacianOperator::new(&k, &l, 0, policy, &execution)?
        } else {
            return Err("invalid partial route".into());
        };
        op.eigenpairs(
            count,
            if positive {
                PersistentSpectrumTarget::SmallestPositive
            } else {
                PersistentSpectrumTarget::Largest
            },
            None,
            &IterativeOptions {
                absolute_tolerance: 1e-10,
                relative_tolerance: 1e-10,
                ..Default::default()
            },
            &execution,
        )?
        .values
    };
    let expected = if kind.starts_with("partial-star") {
        1.
    } else {
        2.
    };
    if values.len() != count || values.iter().any(|v| (v - expected).abs() > 1e-8) {
        return Err("independent partial-spectrum expectation failed".into());
    }
    Ok((values, prepare_ns, start.elapsed().as_nanos()))
}

fn append_operator(
    operator: &PersistentLaplacian,
    values: &mut Vec<f64>,
    spectra: bool,
    execution: &Execution<'_>,
) -> Result<()> {
    values.push(operator.basis().len() as f64);
    values.extend_from_slice(operator.matrix());
    if spectra {
        values.extend(operator.spectrum(execution)?.eigenvalues);
    }
    Ok(())
}

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() < 2 {
        return Err("expected kind route [size repetitions | matrices|spectra all|single]".into());
    }
    let (kind, route) = (args[0].as_str(), args[1].as_str());
    let (values, prepare_ns, solve_ns) = match kind {
        "partial" | "partial-star" | "partial-star-positive" | "compact" if args.len() == 3 => {
            let n = args[2].parse()?;
            if !(8..=512).contains(&n) {
                return Err("partial dimension limit".into());
            }
            partial(route, kind, n)?
        }
        "action" if args.len() == 4 => {
            let n = args[2].parse()?;
            let repeats = args[3].parse()?;
            if n == 0 || n > 1024 || repeats == 0 || repeats > 1000 {
                return Err("action limits".into());
            }
            action(route, n, repeats)?
        }
        "gram" | "harmonic" if args.len() == 3 => {
            let n = args[2].parse()?;
            if n > 128 {
                return Err("spectral control dimension limit".into());
            }
            spectral(route, kind, n)?
        }
        "pairs" | "nonbranching" | "partial-pairs" if args.len() == 4 => {
            if !["matrices", "spectra"].contains(&args[2].as_str())
                || !["all", "single"].contains(&args[3].as_str())
            {
                return Err("pair output/request mode".into());
            }
            let mut input = String::new();
            io::stdin().read_to_string(&mut input)?;
            pairs(
                &input,
                kind,
                route,
                args[2] == "spectra",
                args[3] == "single",
            )?
        }
        _ => return Err("invalid benchmark mode".into()),
    };
    println!("{{\"values\":{values:?},\"prepare_ns\":{prepare_ns},\"solve_ns\":{solve_ns}}}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn star_routes_match_analytic_action() {
        let dense = action("dense", 5, 2).unwrap().0;
        let sparse = action("factorized", 5, 2).unwrap().0;
        for (a, b) in dense.iter().zip(sparse) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn spectral_routes_preserve_analytic_projector() {
        for (kind, route) in [("gram", "both-grams"), ("harmonic", "restriction")] {
            let direct = spectral("full", kind, 8).unwrap().0;
            let selected = spectral(route, kind, 8).unwrap().0;
            assert_eq!(direct.len(), selected.len());
            for (a, b) in direct.iter().zip(selected) {
                assert!((a - b).abs() < 1e-9);
            }
        }
    }
}
