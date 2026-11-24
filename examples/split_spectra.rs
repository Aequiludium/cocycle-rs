//! Bounded real Gram and harmonic-restriction example; --oracle is a test protocol.
use cocycle::execution::Execution;
use cocycle::persistent_laplacian::{
    HarmonicRestriction, SpectralOptions, gram_spectrum, split_spectrum,
};
use std::io::{self, Read};
use std::str::{FromStr, SplitWhitespace};
use std::sync::atomic::AtomicBool;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() {
        let options = SpectralOptions::default();
        let execution = Execution::default();
        let triangle = gram_spectrum(
            &[-1., -1., 0., 1., 0., -1., 0., 1., 1.],
            3,
            3,
            options,
            &execution,
        )?;
        println!("triangle edge spectrum: {:?}", triangle.reconstructed);
        let restriction = HarmonicRestriction::compute(
            &[1., -1., -1., 1.],
            2,
            &[1., 1.],
            1,
            options,
            &execution,
        )?;
        let modes = restriction.smallest_eigenpairs(&[1., -1.], 1, 1, 120, 1e-9, &execution)?;
        println!("connected edge positive spectrum: {:?}", modes.values);
        println!("original-coordinate residuals: {:?}", modes.residuals);
        println!(
            "index coverage certified: {}",
            modes.is_index_coverage_certified()
        );
        return Ok(());
    }
    if args[0] != "--oracle" {
        return Err(
            "usage: split_spectra [--oracle [--max-work N] [--max-entries N] [--cancel]]".into(),
        );
    }
    let mut options = SpectralOptions::default();
    let mut execution = Execution::default();
    let cancel = AtomicBool::new(args.iter().any(|a| a == "--cancel"));
    execution = execution.cancellation(&cancel);
    let mut index = 1;
    while index < args.len() {
        match args[index].as_str() {
            "--cancel" => index += 1,
            "--max-work" => {
                execution =
                    execution.max_work(args.get(index + 1).ok_or("missing work limit")?.parse()?);
                index += 2;
            }
            "--max-entries" => {
                options.max_dense_entries = args
                    .get(index + 1)
                    .ok_or("missing storage limit")?
                    .parse()?;
                index += 2;
            }
            _ => return Err("unknown option".into()),
        }
    }
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let mut words = input.split_whitespace();
    let started = Instant::now();
    let operation = words.next().ok_or("missing operation")?;
    match operation {
        "GRAM" => {
            let m = next(&mut words)?;
            let n = next(&mut words)?;
            let a = matrix(&mut words, m, n)?;
            let g = gram_spectrum(&a, m, n, options, &execution)?;
            println!("PRIMAL {:?}\nDUAL {:?}", g.primal, g.dual);
            println!(
                "PRIMAL_SPECTRUM {:?}\nDUAL_SPECTRUM {:?}",
                g.primal_spectrum, g.dual_spectrum
            );
            println!("DUAL_IS_SMALLER {}", g.dual_is_smaller);
            println!(
                "POSITIVE {:?}\nRECONSTRUCTED {:?}\nLIFTED {:?}",
                g.positive, g.reconstructed, g.lifted_vectors
            );
            println!(
                "RANK {}\nTHRESHOLD {:?}\nAMBIGUOUS {}\nRESIDUAL {:?}\nORTHOGONALITY {:?}",
                g.numerical_rank,
                g.threshold,
                g.rank_ambiguous,
                g.lifted_residual,
                g.orthogonality_error
            );
        }
        "SPLIT" => {
            let m = next(&mut words)?;
            let n = next(&mut words)?;
            let p = next(&mut words)?;
            let a = matrix(&mut words, m, n)?;
            let f = matrix(&mut words, n, p)?;
            let s = split_spectrum(&a, m, n, &f, p, options, &execution)?;
            println!(
                "MATRIX {:?}\nVALUES {:?}\nAMBIGUOUS {}\nCHAIN_RESIDUAL {:?}",
                s.matrix, s.eigenvalues, s.rank_ambiguous, s.chain_residual
            );
        }
        "RESTRICT" | "SOLVE" => {
            let n = next(&mut words)?;
            let k = next(&mut words)?;
            let (block, count, iterations, tolerance) = if operation == "SOLVE" {
                (
                    next(&mut words)?,
                    next(&mut words)?,
                    next(&mut words)?,
                    next(&mut words)?,
                )
            } else {
                (0, 0, 0, 1e-9)
            };
            let h = matrix(&mut words, n, n)?;
            let kernel = matrix(&mut words, n, k)?;
            let r = HarmonicRestriction::compute(&h, n, &kernel, k, options, &execution)?;
            println!("SETUP_SECONDS {:?}", started.elapsed().as_secs_f64());
            println!(
                "KERNEL {:?}\nCOMPLEMENT {:?}\nMATRIX {:?}\nVALUES {:?}\nVECTORS {:?}",
                r.kernel(),
                r.complement(),
                r.matrix(),
                r.eigenvalues(),
                r.vectors()
            );
            println!(
                "THRESHOLD {:?}\nAMBIGUOUS {}\nKERNEL_RESIDUAL {:?}\nRESIDUAL {:?}\nORTHOGONALITY {:?}",
                r.zero_threshold(),
                r.rank_ambiguous(),
                r.kernel_residual(),
                r.residual(),
                r.orthogonality_error()
            );
            if operation == "SOLVE" {
                let initial = matrix(&mut words, n, block)?;
                let solve_started = Instant::now();
                let pairs = r.smallest_eigenpairs(
                    &initial, block, count, iterations, tolerance, &execution,
                )?;
                println!("SOLVE_SECONDS {:?}", solve_started.elapsed().as_secs_f64());
                println!(
                    "PAIR_VALUES {:?}\nPAIR_VECTORS {:?}\nPAIR_RESIDUALS {:?}",
                    pairs.values, pairs.vectors, pairs.residuals
                );
                println!(
                    "LEAKAGE {:?}\nPAIR_ORTHOGONALITY {:?}\nITERATIONS {}\nCERTIFIED {}",
                    pairs.range_leakage,
                    pairs.orthogonality_error,
                    pairs.iterations,
                    pairs.is_index_coverage_certified()
                );
            }
        }
        _ => return Err("unknown operation".into()),
    }
    if words.next().is_some() {
        return Err("trailing protocol data".into());
    }
    println!("TOTAL_SECONDS {:?}", started.elapsed().as_secs_f64());
    Ok(())
}

fn next<T: FromStr>(words: &mut SplitWhitespace<'_>) -> Result<T, Box<dyn std::error::Error>>
where
    T::Err: std::error::Error + 'static,
{
    Ok(words.next().ok_or("incomplete protocol input")?.parse()?)
}
fn matrix(
    words: &mut SplitWhitespace<'_>,
    m: usize,
    n: usize,
) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    let entries = m.checked_mul(n).ok_or("protocol shape overflow")?;
    (0..entries).map(|_| next(words)).collect()
}
