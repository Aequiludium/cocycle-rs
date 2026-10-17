//! Standalone exact F2 operator worker; fixture I/O and export are outside timing.
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::homology_operators::{ChainWindow, F2Matrix, HomologyOperator, OperatorFamily};
use std::hint::black_box;
use std::io::{self, Read};
use std::time::Instant;

fn geometry() -> Result<(), Box<dyn std::error::Error>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let mut tokens = input.split_whitespace();
    let mut take = || tokens.next().ok_or("truncated geometry fixture");
    let degree: usize = take()?.parse()?;
    let count: usize = take()?.parse()?;
    let mut simplices = Vec::new();
    for _ in 0..count {
        let scale: f64 = take()?.parse()?;
        let size: usize = take()?.parse()?;
        let vertices = (0..size)
            .map(|_| Ok(take()?.parse()?))
            .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?;
        simplices.push(Simplex::new(vertices, scale)?);
    }
    if tokens.next().is_some() {
        return Err("extra geometry fixture tokens".into());
    }
    let source = SimplicialComplex::new(simplices)?;
    let controls = Execution::default().max_work(100_000_000);
    let weights = vec![
        1;
        source
            .simplices()
            .iter()
            .filter(|s| s.dimension() == degree)
            .count()
    ];
    let family = OperatorFamily::from_simplices(
        &source,
        degree,
        weights,
        "geometry-fixture",
        "unit",
        &controls,
    )?;
    let mut records = Vec::new();
    for f in family.features(&controls)? {
        let histories = (f.birth()..f.death().unwrap_or(count))
            .map(|i| {
                Ok(family
                    .resolve_feature(&f, i, &controls)?
                    .ok_or("inactive history")?
                    .representative()
                    .to_vec())
            })
            .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?;
        let minimum = family
            .stage(f.birth())?
            .minimum_class_mass(&histories[0], &controls)?;
        let cycle = if degree == 1 {
            format!("{:?}", family.optimal_1_cycle(&f, &controls)?.chain())
        } else {
            "null".into()
        };
        let mut volumes = Vec::new();
        if let Some(death) = f.death() {
            for weighted in [false, true] {
                let window = family.stage(death)?.window();
                let costs = window.bases()[2]
                    .iter()
                    .map(|id| {
                        if weighted {
                            (id.index() % 5 + 1) as u64
                        } else {
                            1
                        }
                    })
                    .collect();
                let upper = family.filling_window(&f, costs, "upper cost", &controls)?;
                let volume = family
                    .optimal_volume(&f, &upper, &controls)?
                    .ok_or("missing volume")?;
                let stable = family
                    .stable_volume(&f, &upper, 1, "filtration time", &controls)?
                    .ok_or("missing stable volume")?;
                let mut z = vec![0; window.weights().len()];
                for (bit, id) in histories[0]
                    .iter()
                    .zip(&family.stage(f.birth())?.window().bases()[1])
                {
                    z[window.bases()[1]
                        .iter()
                        .position(|x| x == id)
                        .ok_or("missing inclusion")?] = *bit;
                }
                let filling = family
                    .stage(death)?
                    .minimum_filling(&z, &upper, &controls)?
                    .ok_or("missing filling")?;
                volumes.push(format!("{{\"volume\":{:?},\"mass\":{},\"boundary\":{:?},\"stable\":{:?},\"stable_mass\":{},\"filling\":{:?},\"filling_mass\":{}}}",volume.chain(),volume.mass(),volume.boundary(),stable.chain(),stable.mass(),filling.chain(),filling.mass()));
            }
        }
        let death = f.death().map_or("null".into(), |d| d.to_string());
        let locations = family.locate_chain(f.birth(), &histories[0])?;
        records.push(format!("{{\"birth\":{},\"death\":{death},\"occurrence\":{},\"histories\":{histories:?},\"minimum\":{:?},\"minimum_mass\":{},\"cycle\":{cycle},\"locations\":{locations:?},\"volumes\":[{}]}}",f.birth(),f.occurrence(),minimum.chain(),minimum.mass(),volumes.join(",")));
    }
    println!(
        "{{\"state\":\"Computed\",\"features\":[{}]}}",
        records.join(",")
    );
    Ok(())
}

fn rss(field: &str) -> Option<u64> {
    let text = std::fs::read_to_string("/proc/self/status").ok()?;
    text.lines().find_map(|line| {
        line.strip_prefix(field)?
            .split_whitespace()
            .next()?
            .parse()
            .ok()
    })
}
fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let mode = args.get(1).map(String::as_str).unwrap_or("solve");
    if mode == "geometry" {
        return geometry();
    }
    if mode == "family" {
        let a = F2Matrix::from_rows(&[vec![1, 1, 0], vec![1, 0, 1], vec![0, 1, 1]], 3)?;
        let ring = ChainWindow::new(
            1,
            [a.clone(), F2Matrix::zero(3, 0)?],
            [vec![0, 1, 2], vec![3, 4, 5], vec![]],
            vec![2, 3, 5],
            "triangle",
            "cost",
        )?;
        let fill = ChainWindow::new(
            1,
            [a, F2Matrix::from_rows(&[vec![1], vec![1], vec![1]], 1)?],
            [vec![0, 1, 2], vec![3, 4, 5], vec![6]],
            vec![2, 3, 5],
            "triangle",
            "cost",
        )?;
        let r = HomologyOperator::feasible(ring)?;
        let family = OperatorFamily::new(
            vec![0., 1., 2.],
            vec![r.clone(), r, HomologyOperator::feasible(fill)?],
        )?;
        let ranks: Vec<_> = (0..3)
            .map(|i| {
                (i..3)
                    .map(|j| family.transport_rank(i, j))
                    .collect::<cocycle::Result<Vec<_>>>()
            })
            .collect::<cocycle::Result<_>>()?;
        let bars: Vec<_> = family
            .barcode()?
            .iter()
            .map(|b| {
                format!(
                    "[{},{},{}]",
                    b.birth,
                    b.death.map_or("null".into(), |d| d.to_string()),
                    b.multiplicity
                )
            })
            .collect();
        println!("{{\"ranks\":{ranks:?},\"bars\":[{}]}}", bars.join(","));
        return Ok(());
    }
    if mode == "family-bench" {
        let n: usize = args.get(2).ok_or("missing family size")?.parse()?;
        let stages: usize = args.get(3).ok_or("missing stage count")?.parse()?;
        let iterations: usize = args.get(4).map(String::as_str).unwrap_or("1").parse()?;
        if stages == 0 || stages > n || iterations == 0 {
            return Err("invalid family dimensions".into());
        }
        let operators = (0..stages)
            .map(|j| {
                let columns: Vec<_> = (0..j)
                    .map(|k| {
                        let mut x = vec![0; n];
                        x[k] = 1;
                        x[k + 1] = 1;
                        x
                    })
                    .collect();
                HomologyOperator::feasible(ChainWindow::new(
                    0,
                    [F2Matrix::zero(0, n)?, F2Matrix::from_columns(&columns, n)?],
                    [vec![], (0..n).collect(), (0..j).map(|k| n + k).collect()],
                    vec![1; n],
                    "shared-merge-tree",
                    "cost",
                )?)
            })
            .collect::<cocycle::Result<Vec<_>>>()?;
        let family = OperatorFamily::new((0..stages).map(|j| (j / 2) as f64).collect(), operators)?;
        let baseline = rss("VmRSS:");
        let started = Instant::now();
        let mut bars = Vec::new();
        for _ in 0..iterations {
            bars = black_box(family.barcode()?);
        }
        let elapsed = started.elapsed().as_nanos();
        let peak = rss("VmHWM:");
        let bars: Vec<_> = bars
            .iter()
            .map(|b| {
                format!(
                    "[{},{},{}]",
                    b.birth,
                    b.death.map_or("null".into(), |d| d.to_string()),
                    b.multiplicity
                )
            })
            .collect();
        let nullable = |v: Option<u64>| v.map_or("null".into(), |x| x.to_string());
        println!(
            "{{\"state\":\"Computed\",\"bars\":[{}],\"elapsed_ns\":{elapsed},\"iterations\":{iterations},\"baseline_rss_kib\":{},\"peak_rss_kib\":{}}}",
            bars.join(","),
            nullable(baseline),
            nullable(peak)
        );
        return Ok(());
    }
    if mode != "solve" && mode != "query" {
        return Err("unknown timing mode".into());
    }
    let iterations: usize = args.get(2).map(String::as_str).unwrap_or("1").parse()?;
    if iterations == 0 {
        return Err("positive iterations required".into());
    }
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let values: Vec<u64> = input
        .split_whitespace()
        .map(str::parse)
        .collect::<Result<_, _>>()?;
    let mut pos = 0;
    let mut take = || {
        let value = values.get(pos).copied().ok_or("truncated fixture");
        pos += 1;
        value
    };
    let m = usize::try_from(take()?)?;
    let n = usize::try_from(take()?)?;
    let p = usize::try_from(take()?)?;
    let q = usize::try_from(take()?)?;
    let mut rows = |r: usize, c: usize| -> Result<Vec<Vec<u8>>, Box<dyn std::error::Error>> {
        (0..r)
            .map(|_| (0..c).map(|_| Ok(u8::try_from(take()?)?)).collect())
            .collect()
    };
    let a = F2Matrix::from_rows(&rows(m, n)?, n)?;
    let d = F2Matrix::from_rows(&rows(n, p)?, p)?;
    let weights = (0..n).map(|_| take()).collect::<Result<Vec<_>, _>>()?;
    let vectors = (0..q)
        .map(|_| {
            (0..n)
                .map(|_| Ok(u8::try_from(take()?)?))
                .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    if pos != values.len() {
        return Err("unexpected fixture tokens".into());
    }
    let window = ChainWindow::new(
        1,
        [a, d],
        [(0..m).collect(), (0..n).collect(), (0..p).collect()],
        weights,
        "shared-fixture",
        "cost",
    )?;
    let baseline_rss = rss("VmRSS:");
    let baseline_hwm = rss("VmHWM:");
    let op = HomologyOperator::feasible(window.clone())?;
    let cycle_vectors: Vec<&[u8]> = vectors
        .iter()
        .filter_map(|z| match op.is_cycle(z) {
            Ok(true) => Some(Ok(z.as_slice())),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        })
        .collect::<cocycle::Result<_>>()?;
    let started = Instant::now();
    for _ in 0..iterations {
        if mode == "solve" {
            black_box(HomologyOperator::feasible(black_box(window.clone()))?);
        } else {
            black_box(op.class_queries_with(black_box(&cycle_vectors), &Execution::default())?);
        }
    }
    let elapsed = started.elapsed().as_nanos();
    // Capture resource data before large logical exports can raise high water.
    let peak = rss("VmHWM:");
    let projection = op.projection_matrix()?;
    let projection: Vec<_> = (0..n)
        .map(|i| projection.row(i))
        .collect::<cocycle::Result<_>>()?;
    let actions: Vec<_> = vectors
        .iter()
        .map(|x| op.project(x))
        .collect::<cocycle::Result<_>>()?;
    let masses: Vec<_> = vectors
        .iter()
        .map(|x| {
            if op.is_cycle(x)? {
                Ok(Some(op.class_query(x)?.selected_mass()))
            } else {
                Ok(None)
            }
        })
        .collect::<cocycle::Result<_>>()?;
    let masses: Vec<_> = masses
        .iter()
        .map(|v| v.map_or("null".into(), |m| m.to_string()))
        .collect();
    let nullable = |v: Option<u64>| v.map_or("null".into(), |x| x.to_string());
    println!(
        "{{\"state\":\"Computed\",\"betti\":{},\"projection\":{projection:?},\"actions\":{actions:?},\"masses\":[{}],\"elapsed_ns\":{elapsed},\"iterations\":{iterations},\"baseline_rss_kib\":{},\"baseline_hwm_kib\":{},\"peak_rss_kib\":{}}}",
        op.betti(),
        masses.join(","),
        nullable(baseline_rss),
        nullable(baseline_hwm),
        nullable(peak)
    );
    Ok(())
}
fn main() {
    if let Err(error) = run() {
        eprintln!("homology worker: {error}");
        std::process::exit(1);
    }
}
