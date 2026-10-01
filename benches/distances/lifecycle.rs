//! Persistent-process experiment. Compiled against generated copies, never src/.
#![forbid(unsafe_code)]
#![allow(dead_code)]
pub use cocycle::{Error, Result, diagram, filtration};
#[path = "src/diagram_distances/mod.rs"]
mod diagram_distances;
#[path = "src/execution/mod.rs"]
mod execution;
#[path = "pool.rs"]
mod lifecycle;

use diagram::{Coverage, IntervalEnd, PersistenceDiagram, PersistenceInterval};
use diagram_distances::{Kind, bottleneck as bn, wasserstein as w};
use execution::WorkBudget;
use std::{
    env, fs,
    io::{self, Write},
    time::Instant,
};
type WPrepared = (Vec<w::LifecyclePoint>, f64, bool);
type Op = (usize, usize, usize);

fn make(points: &[[f64; 2]]) -> Result<PersistenceDiagram> {
    PersistenceDiagram::new(
        0,
        Coverage::Complete,
        points
            .iter()
            .map(|p| {
                PersistenceInterval::new(
                    0,
                    p[0],
                    if p[1] == f64::INFINITY {
                        IntervalEnd::Essential
                    } else {
                        IntervalEnd::Finite(p[1])
                    },
                )
            })
            .collect::<Result<Vec<_>>>()?,
    )
}
fn public(
    a: &PersistenceDiagram,
    b: &PersistenceDiagram,
    dimension: usize,
    metric: &str,
) -> Result<f64> {
    match metric {
        "bottleneck" => cocycle::diagram_distances::bottleneck_distance(a, b, dimension),
        "w1" => cocycle::diagram_distances::wasserstein_1_infinity(a, b, dimension),
        "w2" => cocycle::diagram_distances::wasserstein_2_euclidean(a, b, dimension),
        _ => panic!("invalid metric"),
    }
}
fn kind(metric: &str) -> Kind {
    match metric {
        "bottleneck" => Kind::Bottleneck,
        "w1" => Kind::W1,
        "w2" => Kind::W2,
        _ => panic!("invalid metric"),
    }
}
fn coordinates(view: &diagram::DiagramDimension<'_>, count: usize) -> Vec<[f64; 2]> {
    let _clock = lifecycle::Clock::new(1);
    let mut points = Vec::with_capacity(count);
    for i in view.iter() {
        if let IntervalEnd::Finite(d) = i.end() {
            points.push([i.birth(), d]);
        }
    }
    points
}
fn candidate<const C: bool>(
    a: &PersistenceDiagram,
    b: &PersistenceDiagram,
    dimension: usize,
    metric: &str,
    pa: Option<&bn::Prepared<'_>>,
    pb: Option<&bn::Prepared<'_>>,
    wa: Option<&WPrepared>,
    wb: Option<&WPrepared>,
    budget: &mut WorkBudget<'_, C>,
) -> Result<f64> {
    diagram_distances::lifecycle_distance(
        a,
        b,
        dimension,
        kind(metric),
        budget,
        |a, b, counts, budget| {
            if metric == "bottleneck" {
                // Charge the extraction scan on uncached controlled calls as R0 does.
                let ac;
                let bc;
                let ap;
                let bp;
                let pa = if let Some(p) = pa {
                    p
                } else {
                    budget.step_by(a.len())?;
                    ac = coordinates(a, counts.0);
                    ap = bn::lifecycle_prepare_with(&ac, budget)?;
                    &ap
                };
                let pb = if let Some(p) = pb {
                    p
                } else {
                    budget.step_by(b.len())?;
                    bc = coordinates(b, counts.1);
                    bp = bn::lifecycle_prepare_with(&bc, budget)?;
                    &bp
                };
                bn::lifecycle_solve(
                    pa,
                    pb,
                    bn::Options::default(),
                    &mut bn::Diagnostics::default(),
                    budget,
                )
            } else {
                let ap;
                let bp;
                let wa = if let Some(p) = wa {
                    p
                } else {
                    ap = w::lifecycle_prepare(a, counts.0, budget)?;
                    &ap
                };
                let wb = if let Some(p) = wb {
                    p
                } else {
                    bp = w::lifecycle_prepare(b, counts.1, budget)?;
                    &bp
                };
                w::lifecycle_solve(
                    wa,
                    wb,
                    if metric == "w1" {
                        w::Metric::W1
                    } else {
                        w::Metric::W2
                    },
                    &mut w::Stats::default(),
                    budget,
                )
            }
        },
    )
}
fn rss() -> String {
    fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines().find_map(|l| {
                l.strip_prefix("VmRSS:")?
                    .split_whitespace()
                    .next()?
                    .parse::<usize>()
                    .ok()
            })
        })
        .map_or("null".into(), |x| (x * 1024).to_string())
}
fn emit(line: &str, handshake: bool) {
    println!("{line}");
    io::stdout().flush().unwrap();
    if handshake {
        let mut ack = String::new();
        io::stdin().read_line(&mut ack).unwrap();
        assert_eq!(ack.trim(), "ack");
    }
}
fn same(a: f64, b: f64) -> bool {
    a == b || (a.is_finite() && b.is_finite() && (a - b).abs() <= 1e-10 * b.abs().max(1.0))
}

fn run(
    diagrams: &[PersistenceDiagram],
    ops: &[Op],
    fixed: usize,
    metric: &str,
    variant: &str,
    detailed: bool,
    trace: bool,
    handshake: bool,
    reference: Option<Vec<f64>>,
) {
    assert!(!ops.is_empty());
    assert!(["public", "r0", "r1", "r2", "r3", "r4", "r5", "bounded"].contains(&variant));
    // The ordinary public R0 is an equivalence control, not an independent oracle.
    let mut expectations = std::collections::HashMap::new();
    let expected: Vec<_> = reference.unwrap_or_else(|| {
        ops.iter()
            .map(|&(a, b, _)| {
                *expectations
                    .entry((a, b))
                    .or_insert_with(|| public(&diagrams[a], &diagrams[b], 0, metric).unwrap())
            })
            .collect()
    });
    assert_eq!(expected.len(), ops.len());
    lifecycle::configure(
        ["r3", "r4", "r5", "bounded"].contains(&variant),
        detailed,
        if variant == "bounded" {
            4096
        } else {
            usize::MAX
        },
    );
    // Warm only the first two small operations to fault in generated code. All
    // preparation and retained buffers are dropped before the measured sequence.
    for &(a, b, _) in ops.iter().take(2) {
        lifecycle::begin();
        if variant == "public" {
            std::hint::black_box(public(&diagrams[a], &diagrams[b], 0, metric).unwrap());
        } else {
            std::hint::black_box(
                candidate(
                    &diagrams[a],
                    &diagrams[b],
                    0,
                    metric,
                    None,
                    None,
                    None,
                    None,
                    &mut WorkBudget::unlimited(),
                )
                .unwrap(),
            );
        }
    }
    lifecycle::clear();
    lifecycle::reset_stats();
    if handshake {
        emit("{\"type\":\"ready\"}", true);
    }
    let compute_started_unix_ns = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let start = Instant::now();
    let mut frequency = vec![0usize; diagrams.len()];
    for &(a, b, _) in ops {
        frequency[a] += 1;
        frequency[b] += 1;
    }
    let mask: Vec<_> = frequency
        .iter()
        .enumerate()
        .map(|(i, &f)| match variant {
            "public" | "r0" => false,
            "r1" => i == fixed,
            _ => f > 1,
        })
        .collect();
    let prep_start = Instant::now();
    let coords: Vec<_> = diagrams
        .iter()
        .zip(&mask)
        .map(|(d, &yes)| {
            if yes && metric == "bottleneck" {
                let v = d.dimension(0).unwrap();
                Some(coordinates(&v, v.len()))
            } else {
                None
            }
        })
        .collect();
    let bcache: Vec<_> = coords
        .iter()
        .map(|p| p.as_ref().map(|p| bn::lifecycle_prepare(p).unwrap()))
        .collect();
    let wcache: Vec<_> = diagrams
        .iter()
        .zip(&mask)
        .map(|(d, &yes)| {
            if yes && metric != "bottleneck" {
                let v = d.dimension(0).unwrap();
                Some(w::lifecycle_prepare(&v, v.len(), &mut WorkBudget::unlimited()).unwrap())
            } else {
                None
            }
        })
        .collect();
    let prep_ns = prep_start.elapsed().as_nanos();
    let prepared_bytes: usize = coords
        .iter()
        .flatten()
        .map(|p| p.capacity() * 16)
        .sum::<usize>()
        + bcache.iter().flatten().map(|p| p.bytes()).sum::<usize>()
        + wcache
            .iter()
            .flatten()
            .map(|p| p.0.capacity() * std::mem::size_of::<w::LifecyclePoint>())
            .sum::<usize>();
    let rss_pre = if trace { rss() } else { "null".into() };
    let mut output = Vec::new();
    let mut previous_batch = None;
    let mut checked = 0;
    let mut checksum = 0.0;
    let mut peak_retained = 0;
    let mut loop_ns = 0u128;
    // Timed loops exclude trace I/O and RSS reads. Formal non-trace total uses
    // one whole-sequence clock; individual clocks belong only to poison traces.
    let loop_start = Instant::now();
    for (index, &(a, b, batch)) in ops.iter().enumerate() {
        let step = trace.then(Instant::now);
        if previous_batch != Some(batch) {
            if variant == "r3" {
                lifecycle::clear();
            }
            if variant != "r5" {
                output = Vec::new();
            } else {
                output.clear();
            }
            previous_batch = Some(batch);
        }
        lifecycle::begin();
        let value = if variant == "public" {
            public(&diagrams[a], &diagrams[b], 0, metric)
        } else {
            candidate(
                &diagrams[a],
                &diagrams[b],
                0,
                metric,
                bcache[a].as_ref(),
                bcache[b].as_ref(),
                wcache[a].as_ref(),
                wcache[b].as_ref(),
                &mut WorkBudget::unlimited(),
            )
        }
        .unwrap();
        output.push(std::hint::black_box(value));
        assert!(same(value, expected[index]), "R0 mismatch at {index}");
        checked += 1;
        checksum += value;
        let elapsed = step.map_or(0, |s| s.elapsed().as_nanos());
        if trace {
            loop_ns += elapsed;
            let (_, hits, growths, capacity) = lifecycle::stats();
            peak_retained = peak_retained.max(capacity);
            emit(
                &format!(
                    "{{\"type\":\"trace\",\"index\":{index},\"pair\":[{a},{b}],\"batch\":{batch},\"elapsed_ns\":{elapsed},\"retained_bytes\":{capacity},\"hits\":{hits},\"growths\":{growths},\"rss_bytes\":{}}}",
                    rss()
                ),
                handshake,
            );
        }
    }
    if !trace {
        loop_ns = loop_start.elapsed().as_nanos();
    }
    let total_ns = if trace {
        prep_ns + loop_ns
    } else {
        start.elapsed().as_nanos()
    };
    let compute_finished_unix_ns = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let (phases, hits, growths, retained) = lifecycle::stats();
    peak_retained = peak_retained.max(retained);
    emit(
        &format!(
            "{{\"type\":\"result\",\"protocol\":\"cocycle-distance-lifecycle-v1\",\"metric\":\"{metric}\",\"variant\":\"{variant}\",\"operations\":{checked},\"total_ns\":{total_ns},\"compute_started_unix_ns\":{compute_started_unix_ns},\"compute_finished_unix_ns\":{compute_finished_unix_ns},\"prep_ns\":{prep_ns},\"loop_ns\":{loop_ns},\"prepared_bytes\":{prepared_bytes},\"retained_bytes\":{retained},\"peak_retained_bytes\":{peak_retained},\"output_bytes\":{},\"hits\":{hits},\"growths\":{growths},\"phases_ns\":{phases:?},\"checksum\":{checksum},\"rss_pre_bytes\":{rss_pre},\"rss_bytes\":{}}}",
            output.capacity() * 8,
            rss()
        ),
        handshake,
    );
}
fn selftest() {
    lifecycle::configure(true, false, 4096);
    lifecycle::begin();
    lifecycle::put("test-slot", Vec::<usize>::with_capacity(128));
    assert_eq!(lifecycle::empty::<usize>("test-slot").capacity(), 0);
    lifecycle::begin();
    assert!(lifecycle::empty::<usize>("test-slot").capacity() >= 128);
    lifecycle::put("test-slot", Vec::<usize>::with_capacity(1024));
    lifecycle::begin();
    assert_eq!(lifecycle::empty::<usize>("test-slot").capacity(), 0);
    lifecycle::clear();
    let a = make(&[[0., 2.], [0., 2.]]).unwrap();
    let b = make(&[[0.25, 2.25]]).unwrap();
    let empty = make(&[]).unwrap();
    let essential = make(&[[0., f64::INFINITY]]).unwrap();
    let truncated = PersistenceDiagram::new(0, Coverage::Through(3.), vec![]).unwrap();
    let noncontiguous = PersistenceDiagram::with_dimensions(
        diagram::ComputedDimensions::new(vec![0, 2]).unwrap(),
        Coverage::Complete,
        vec![],
    )
    .unwrap();
    let semantic = [
        (&a, &empty, 0),
        (&empty, &empty, 0),
        (&a, &essential, 0),
        (&a, &a, 1),
        (&a, &truncated, 0),
        (&noncontiguous, &noncontiguous, 1),
        (&noncontiguous, &noncontiguous, 2),
    ];
    for retain in [false, true] {
        lifecycle::configure(retain, false, usize::MAX);
        for metric in ["bottleneck", "w1", "w2"] {
            for &(x, y, d) in &semantic {
                lifecycle::begin();
                let actual = candidate(
                    x,
                    y,
                    d,
                    metric,
                    None,
                    None,
                    None,
                    None,
                    &mut WorkBudget::unlimited(),
                );
                let expected = public(x, y, d, metric);
                match (actual, expected) {
                    (Ok(x), Ok(y)) => assert!(same(x, y)),
                    (Err(x), Err(y)) => assert_eq!(format!("{x:?}"), format!("{y:?}")),
                    _ => panic!("semantic mismatch"),
                }
            }
            for limit in [0, 1, 8, 32, 64] {
                lifecycle::begin();
                let mut budget =
                    WorkBudget::new(&execution::Execution::default().max_work(limit)).unwrap();
                assert!(
                    candidate(&a, &b, 0, metric, None, None, None, None, &mut budget).is_err()
                        || limit == 64
                );
                lifecycle::begin();
                assert!(same(
                    candidate(
                        &a,
                        &b,
                        0,
                        metric,
                        None,
                        None,
                        None,
                        None,
                        &mut WorkBudget::unlimited()
                    )
                    .unwrap(),
                    public(&a, &b, 0, metric).unwrap()
                ));
            }
            let flag = std::sync::atomic::AtomicBool::new(true);
            assert!(matches!(
                WorkBudget::new(&execution::Execution::default().cancellation(&flag)),
                Err(Error::Cancelled)
            ));
        }
    }
    // Checked wide scaling and an actual numerical failure, followed by retry.
    let wide = make(&[[2f64.powi(500), 2f64.powi(501)]]).unwrap();
    let tiny = make(&[[0., f64::from_bits(1)]]).unwrap();
    let bad = make(&[[0., f64::from_bits(1)], [2f64.powi(500), 2f64.powi(501)]]).unwrap();
    for metric in ["w1", "bottleneck", "w2", "w1"] {
        for (x, y) in [(&wide, &wide), (&tiny, &empty), (&a, &b)] {
            lifecycle::begin();
            let v = candidate(
                x,
                y,
                0,
                metric,
                None,
                None,
                None,
                None,
                &mut WorkBudget::unlimited(),
            );
            let p = public(x, y, 0, metric);
            match (v, p) {
                (Ok(x), Ok(y)) => assert!(same(x, y)),
                (Err(_), Err(_)) => {}
                _ => panic!("scale/failure mismatch"),
            }
        }
    }
    // Hand-derived independent values: diagonal cost, multiplicity, W2 norm.
    for (metric, value) in [("bottleneck", 1.), ("w1", 2.), ("w2", 2.)] {
        assert!(same(
            candidate(
                &a,
                &empty,
                0,
                metric,
                None,
                None,
                None,
                None,
                &mut WorkBudget::unlimited()
            )
            .unwrap(),
            value
        ));
    }
    let ac = coordinates(&a.dimension(0).unwrap(), a.len());
    let bc = coordinates(&b.dimension(0).unwrap(), b.len());
    let ap = bn::lifecycle_prepare(&ac).unwrap();
    let bp = bn::lifecycle_prepare(&bc).unwrap();
    let wa = w::lifecycle_prepare(
        &a.dimension(0).unwrap(),
        a.len(),
        &mut WorkBudget::unlimited(),
    )
    .unwrap();
    let wb = w::lifecycle_prepare(
        &b.dimension(0).unwrap(),
        b.len(),
        &mut WorkBudget::unlimited(),
    )
    .unwrap();
    let ww = w::lifecycle_prepare(
        &wide.dimension(0).unwrap(),
        wide.len(),
        &mut WorkBudget::unlimited(),
    )
    .unwrap();
    let wt = w::lifecycle_prepare(
        &bad.dimension(0).unwrap(),
        bad.len(),
        &mut WorkBudget::unlimited(),
    )
    .unwrap();
    for metric in ["w2", "bottleneck", "w1", "bottleneck", "w2"] {
        lifecycle::begin();
        let actual = candidate(
            &a,
            &b,
            0,
            metric,
            Some(&ap),
            Some(&bp),
            Some(&wa),
            Some(&wb),
            &mut WorkBudget::unlimited(),
        )
        .unwrap();
        assert!(same(actual, public(&a, &b, 0, metric).unwrap()));
        if metric != "bottleneck" {
            lifecycle::begin();
            assert!(same(
                candidate(
                    &wide,
                    &wide,
                    0,
                    metric,
                    None,
                    None,
                    Some(&ww),
                    Some(&ww),
                    &mut WorkBudget::unlimited()
                )
                .unwrap(),
                0.
            ));
            lifecycle::begin();
            assert!(
                candidate(
                    &bad,
                    &empty,
                    0,
                    metric,
                    None,
                    None,
                    Some(&wt),
                    None,
                    &mut WorkBudget::unlimited()
                )
                .is_err()
            );
            lifecycle::begin();
            assert!(same(
                candidate(
                    &b,
                    &a,
                    0,
                    metric,
                    None,
                    None,
                    Some(&wb),
                    Some(&wa),
                    &mut WorkBudget::unlimited()
                )
                .unwrap(),
                actual
            ));
        }
    }
    // Use library-produced contexts; external callers cannot forge these types.
    use cocycle::persistence::PersistenceExt;
    let source = cocycle::filtration::RipsBuilder::from_points(
        cocycle::geometry::PointCloudView::new(&[0., 1., 3.], 3, 1).unwrap(),
    );
    let field2 = source.persistence().compute().unwrap().into_data();
    let field3 = source
        .persistence()
        .field(cocycle::algebra::PrimeField::new(3).unwrap())
        .compute()
        .unwrap()
        .into_data();
    assert!(diagram_distances::check_context(&field2, &field3).is_err());
    assert!(diagram_distances::check_context(&field2, &field2).is_ok());
    assert!(cocycle::diagram_distances::bottleneck_distance_results(&field2, &field3, 0).is_err());
    let supplied = cocycle::complex::SimplicialComplex::new(vec![
        cocycle::complex::Simplex::new(vec![0], -2.).unwrap(),
    ])
    .unwrap()
    .persistence()
    .compute()
    .unwrap()
    .into_data();
    assert!(diagram_distances::check_context(&supplied, &supplied).is_err());
    let dense_a = make(
        &(0..16)
            .map(|i| {
                [
                    (i * 37 % 97) as f64 / 256.,
                    32. + (i * 37 % 97) as f64 / 256.,
                ]
            })
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let dense_b = make(
        &(0..16)
            .map(|i| {
                [
                    (i * 43 % 89) as f64 / 256.,
                    32.125 + (i * 43 % 89) as f64 / 256.,
                ]
            })
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let da = coordinates(&dense_a.dimension(0).unwrap(), 16);
    let db = coordinates(&dense_b.dimension(0).unwrap(), 16);
    let dpa = bn::lifecycle_prepare(&da).unwrap();
    let dpb = bn::lifecycle_prepare(&db).unwrap();
    lifecycle::clear();
    lifecycle::configure(true, false, usize::MAX);
    let mut failed_with_scratch = false;
    for limit in [256, 512, 1024, 2048, 4096, 8192] {
        lifecycle::begin();
        let result = candidate(
            &dense_a,
            &dense_b,
            0,
            "bottleneck",
            Some(&dpa),
            Some(&dpb),
            None,
            None,
            &mut WorkBudget::new(&execution::Execution::default().max_work(limit)).unwrap(),
        );
        failed_with_scratch |= result.is_err() && lifecycle::stats().3 > 0;
        lifecycle::begin();
        assert!(same(
            candidate(
                &dense_a,
                &dense_b,
                0,
                "bottleneck",
                Some(&dpa),
                Some(&dpb),
                None,
                None,
                &mut WorkBudget::unlimited()
            )
            .unwrap(),
            public(&dense_a, &dense_b, 0, "bottleneck").unwrap()
        ));
    }
    assert!(
        failed_with_scratch,
        "must exercise failure after scratch allocation"
    );
    lifecycle::begin();
    assert!(matches!(
        lifecycle::filled::<usize>("bn-left", usize::MAX, 0),
        Err(Error::AllocationFailed { .. })
    ));
    let sparse_a = make(
        &(0..16)
            .map(|i| [i as f64 * 8., i as f64 * 8. + 1.])
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let sparse_b = make(
        &(0..16)
            .map(|i| [i as f64 * 8. + 0.125, i as f64 * 8. + 1.125])
            .collect::<Vec<_>>(),
    )
    .unwrap();
    let pairs = [
        (&dense_a, &dense_b),
        (&sparse_a, &sparse_b),
        (&a, &b),
        (&wide, &wide),
    ];
    for reverse in [false, true, false] {
        for metric in ["w1", "w2", "bottleneck", "w2", "w1"] {
            for index in 0..pairs.len() {
                let (x, y) = pairs[if reverse {
                    pairs.len() - 1 - index
                } else {
                    index
                }];
                lifecycle::begin();
                assert!(same(
                    candidate(
                        x,
                        y,
                        0,
                        metric,
                        None,
                        None,
                        None,
                        None,
                        &mut WorkBudget::unlimited()
                    )
                    .unwrap(),
                    public(x, y, 0, metric).unwrap()
                ));
            }
        }
    }
    emit("{\"type\":\"selftest\",\"status\":\"passed\"}", false);
}
fn read(path: &str) -> (Vec<PersistenceDiagram>, Vec<Op>, usize) {
    let data = fs::read(path).unwrap();
    assert_eq!(&data[..8], b"COCLIF1\0");
    let mut offset = 8;
    let mut next = || {
        let x = u64::from_le_bytes(data[offset..offset + 8].try_into().unwrap());
        offset += 8;
        x
    };
    let nd = next() as usize;
    let no = next() as usize;
    let fixed = next() as usize;
    assert!(fixed < nd);
    let diagrams = (0..nd)
        .map(|_| {
            let n = next() as usize;
            let points = (0..n)
                .map(|_| [f64::from_bits(next()), f64::from_bits(next())])
                .collect::<Vec<_>>();
            make(&points).unwrap()
        })
        .collect();
    let ops = (0..no)
        .map(|_| {
            let a = next() as usize;
            let b = next() as usize;
            let batch = next() as usize;
            assert!(a < nd && b < nd);
            (a, b, batch)
        })
        .collect();
    assert_eq!(offset, data.len());
    (diagrams, ops, fixed)
}
fn main() {
    let args: Vec<_> = env::args().collect();
    if args.len() == 2 && args[1] == "--selftest" {
        selftest();
        return;
    }
    if args.len() == 4 && args[1] == "--reference" {
        let (diagrams, ops, _) = read(&args[2]);
        let values: Vec<_> = ops
            .iter()
            .map(|&(a, b, _)| public(&diagrams[a], &diagrams[b], 0, &args[3]).unwrap())
            .collect();
        println!("{{\"values\":{values:?}}}");
        return;
    }
    assert_eq!(
        args.len(),
        7,
        "input metric variant detailed trace handshake"
    );
    let (diagrams, ops, fixed) = read(&args[1]);
    let reference = if args[5] == "1" {
        let bytes = fs::read(format!("{}.expected.bin", args[1])).unwrap();
        assert_eq!(bytes.len(), ops.len() * 8);
        Some(
            bytes
                .chunks_exact(8)
                .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
                .collect(),
        )
    } else {
        None
    };
    run(
        &diagrams,
        &ops,
        fixed,
        &args[2],
        &args[3],
        args[4] == "1",
        args[5] == "1",
        args[6] == "1",
        reference,
    );
}
