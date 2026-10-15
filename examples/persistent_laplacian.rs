//! A runnable path example and a whitespace protocol for independent verification.
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::persistent_laplacian::{
    PersistentLaplacian, SpectralOptions, fixed_endpoint_laplacians, fixed_endpoint_up_laplacians,
};
use std::io::{self, Read};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        let source = SimplicialComplex::new(vec![
            Simplex::new(vec![0], 0.)?,
            Simplex::new(vec![2], 0.)?,
            Simplex::new(vec![1], 1.)?,
            Simplex::new(vec![0, 1], 2.)?,
            Simplex::new(vec![1, 2], 2.)?,
        ])?;
        let k = source.stage(0.)?;
        let l = source.stage(2.)?;
        let result = PersistentLaplacian::compute(
            &k,
            &l,
            0,
            SpectralOptions::default(),
            &Execution::default(),
        )?;
        println!(
            "basis: {:?}",
            result
                .basis()
                .iter()
                .map(Simplex::vertices)
                .collect::<Vec<_>>()
        );
        println!("matrix: {:?}", result.matrix());
        println!(
            "spectrum: {:?}",
            result.spectrum(&Execution::default())?.eigenvalues
        );
        println!(
            "action on [1,-1]: {:?}",
            result.action(&[1., -1.], &Execution::default())?
        );
        return Ok(());
    }
    let (benchmark, components) = match args.as_slice() {
        [a] if a == "--oracle" => (None, false),
        [a] if a == "--components" => (None, true),
        [a, b] if a == "--benchmark" && (b == "direct" || b == "stepwise") => {
            (Some(b.as_str()), false)
        }
        _ => {
            return Err(
                "usage: persistent_laplacian [--oracle | --components | --benchmark direct|stepwise]".into(),
            );
        }
    };
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let mut words = input.split_whitespace();
    while let Some(first) = words.next() {
        let degree: usize = first.parse()?;
        let start_count: usize = next(&mut words)?;
        let cell_count: usize = next(&mut words)?;
        let end: f64 = next(&mut words)?;
        let scales: Vec<f64> = (0..start_count)
            .map(|_| next(&mut words))
            .collect::<Result<_, _>>()?;
        let mut cells = Vec::new();
        for _ in 0..cell_count {
            let value: f64 = next(&mut words)?;
            let length: usize = next(&mut words)?;
            let vertices: Vec<usize> = (0..length)
                .map(|_| next(&mut words))
                .collect::<Result<_, _>>()?;
            cells.push(Simplex::new(vertices, value)?);
        }
        let source = SimplicialComplex::new(cells)?;
        let target = source.stage(end)?;
        let starts = scales
            .iter()
            .map(|s| source.stage(*s))
            .collect::<cocycle::Result<Vec<_>>>()?;
        let options = SpectralOptions {
            max_dense_entries: if benchmark.is_some() {
                65_536
            } else {
                4_000_000
            },
            ..Default::default()
        };
        let execution = Execution::default();
        if let Some(route) = benchmark {
            let begin = Instant::now();
            let results = if route == "direct" {
                starts
                    .iter()
                    .map(|k| PersistentLaplacian::compute(k, &target, degree, options, &execution))
                    .collect::<cocycle::Result<Vec<_>>>()?
            } else {
                fixed_endpoint_laplacians(&starts, &target, degree, options, &execution)?
            };
            let elapsed = begin.elapsed().as_secs_f64();
            let checksum: f64 = results.iter().flat_map(|r| r.matrix()).sum();
            std::hint::black_box(&results);
            let fallbacks = results.iter().filter(|r| r.used_direct_fallback()).count();
            println!("TIME {elapsed:.17e} CHECKSUM {checksum:.17e} FALLBACK {fallbacks}");
            continue;
        }
        let steps = fixed_endpoint_laplacians(&starts, &target, degree, options, &execution)?;
        let upper_steps = if components {
            Some(fixed_endpoint_up_laplacians(
                &starts, &target, degree, options, &execution,
            )?)
        } else {
            None
        };
        for (index, (k, step)) in starts.iter().zip(steps).enumerate() {
            let r = PersistentLaplacian::compute(k, &target, degree, options, &execution)?;
            let spectrum = r.spectrum(&execution)?;
            println!(
                "PAIR {} {} {} {}",
                degree,
                k.scale(),
                target.scale(),
                r.basis().len()
            );
            print!("BASIS");
            for s in r.basis() {
                print!(" {}", s.vertices().len());
                for v in s.vertices() {
                    print!(" {v}");
                }
            }
            println!();
            print_values("MATRIX", r.matrix());
            print_values("STEPWISE", step.matrix());
            print_values("SPECTRUM", &spectrum.eigenvalues);
            let x: Vec<_> = (0..r.basis().len()).map(|i| (i as f64 + 1.) / 7.).collect();
            print_values("ACTION", &r.action(&x, &execution)?);
            println!(
                "RESIDUAL {} {} {}",
                r.solver_residual(),
                spectrum.residual,
                step.solver_residual()
            );
            if let Some(upper_steps) = &upper_steps {
                let up = PersistentLaplacian::compute_up(k, &target, degree, options, &execution)?;
                let down = PersistentLaplacian::compute_down(k, degree, options, &execution)?;
                let upper_spectrum = up.spectrum(&execution)?;
                let lower_spectrum = down.spectrum(&execution)?;
                print_values("UP", up.matrix());
                print_values("DOWN", down.matrix());
                print_values("UP_STEPWISE", upper_steps[index].matrix());
                print_values("UP_SPECTRUM", &upper_spectrum.eigenvalues);
                print_values("DOWN_SPECTRUM", &lower_spectrum.eigenvalues);
                print_values("UP_ACTION", &up.action(&x, &execution)?);
                print_values("DOWN_ACTION", &down.action(&x, &execution)?);
                print_values(
                    "COMPONENT_RESIDUAL",
                    &[
                        up.solver_residual(),
                        upper_spectrum.residual,
                        lower_spectrum.residual,
                        upper_steps[index].solver_residual(),
                    ],
                );
                println!(
                    "FALLBACK {} {}",
                    u8::from(step.used_direct_fallback()),
                    u8::from(upper_steps[index].used_direct_fallback())
                );
            }
        }
    }
    Ok(())
}
fn next<T: std::str::FromStr>(
    words: &mut std::str::SplitWhitespace<'_>,
) -> Result<T, Box<dyn std::error::Error>>
where
    T::Err: std::error::Error + 'static,
{
    Ok(words.next().ok_or("truncated fixture")?.parse()?)
}
fn print_values(name: &str, values: &[f64]) {
    print!("{name}");
    for v in values {
        print!(" {v:.17e}");
    }
    println!();
}
