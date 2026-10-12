//! Standalone exact F2 operator worker; fixture I/O and export are outside timing.
use cocycle::execution::Execution;
use cocycle::homology_operators::{ChainWindow, F2Matrix, HomologyOperator, OperatorFamily};
use std::hint::black_box;
use std::io::{self, Read};
use std::time::Instant;

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
