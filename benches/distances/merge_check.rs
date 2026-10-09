//! Before/after scalar and opt-in prepared calls. Parsing precedes timing;
//! preparation, K queries and preparation destruction are inside the timer.
use cocycle::diagram::{Coverage, IntervalEnd, PersistenceDiagram, PersistenceInterval};
use cocycle::diagram_distances::{
    bottleneck_distance, wasserstein_1_infinity, wasserstein_2_euclidean,
};
#[cfg(prepared_candidate)]
use cocycle::{diagram_distances::PreparedDiagram, execution::Execution};
use std::{env, fs, hint::black_box, time::Instant};

fn main() -> cocycle::Result<()> {
    let args: Vec<_> = env::args().collect();
    assert_eq!(args.len(), 5);
    let bytes = fs::read(&args[1]).unwrap();
    assert_eq!(&bytes[..8], b"COCDST1\0");
    let n = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
    let m = u64::from_le_bytes(bytes[16..24].try_into().unwrap()) as usize;
    assert_eq!(bytes.len(), 24 + (n + m) * 16);
    let bars: Vec<_> = bytes[24..]
        .chunks_exact(16)
        .map(|p| {
            PersistenceInterval::new(
                0,
                f64::from_le_bytes(p[..8].try_into().unwrap()),
                IntervalEnd::Finite(f64::from_le_bytes(p[8..].try_into().unwrap())),
            )
            .unwrap()
        })
        .collect();
    let a = PersistenceDiagram::new(0, Coverage::Complete, bars[..n].to_vec())?;
    let b = PersistenceDiagram::new(0, Coverage::Complete, bars[n..].to_vec())?;
    drop(bars);
    drop(bytes);
    let distance = match args[2].as_str() {
        "bottleneck" => bottleneck_distance,
        "w1" => wasserstein_1_infinity,
        "w2" => wasserstein_2_euclidean,
        _ => panic!("unknown metric"),
    };
    let k: usize = args[3].parse().unwrap();
    assert!(k > 0);
    let expected = distance(&a, &b, 0)?;
    for _ in 0..2 {
        assert_eq!(distance(&a, &b, 0)?, expected);
    }
    let start = Instant::now();
    let mut value = 0.;
    if args[4] == "prepared" {
        #[cfg(prepared_candidate)]
        {
            let controls = Execution::default();
            let prepare = match args[2].as_str() {
                "bottleneck" => PreparedDiagram::bottleneck_with,
                "w1" => PreparedDiagram::wasserstein_1_infinity_with,
                "w2" => PreparedDiagram::wasserstein_2_euclidean_with,
                _ => unreachable!(),
            };
            let pa = prepare(&a, 0, &controls)?;
            let pb = prepare(&b, 0, &controls)?;
            for _ in 0..k {
                value = black_box(pa.distance_with(&pb, &controls)?);
                assert_eq!(value, expected);
            }
            drop(pa);
            drop(pb);
        }
        #[cfg(not(prepared_candidate))]
        panic!("prepared mode unavailable in baseline");
    } else {
        assert_eq!(args[4], "scalar");
        for _ in 0..k {
            value = black_box(distance(black_box(&a), black_box(&b), 0)?);
            assert_eq!(value, expected);
        }
    }
    let elapsed_ms = start.elapsed().as_secs_f64() * 1000.;
    println!("{{\"value\":{value},\"elapsed_ms\":{elapsed_ms}}}");
    Ok(())
}
