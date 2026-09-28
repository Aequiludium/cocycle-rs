//! Warmed public calls for the distance preparation reproducer.
#![forbid(unsafe_code)]

use cocycle::diagram::{Coverage, IntervalEnd, PersistenceDiagram, PersistenceInterval};
use cocycle::diagram_distances::{
    bottleneck_distance, wasserstein_1_infinity, wasserstein_2_euclidean,
};
use std::{env, fs, hint::black_box, time::Instant};

fn memory(field: &str) -> usize {
    fs::read_to_string("/proc/self/status")
        .unwrap()
        .lines()
        .find_map(|line| {
            line.strip_prefix(field)?
                .split_whitespace()
                .next()?
                .parse()
                .ok()
        })
        .unwrap()
}

fn main() -> cocycle::Result<()> {
    let args: Vec<_> = env::args().collect();
    assert_eq!(args.len(), 5, "fixture metric iterations variant");
    let bytes = fs::read(&args[1]).unwrap();
    assert_eq!(&bytes[..8], b"COCDST1\0");
    let n = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
    let m = u64::from_le_bytes(bytes[16..24].try_into().unwrap()) as usize;
    assert_eq!(bytes.len(), 24 + (n + m) * 16);
    let intervals: Vec<_> = bytes[24..]
        .chunks_exact(16)
        .map(|p| {
            let birth = f64::from_le_bytes(p[..8].try_into().unwrap());
            let death = f64::from_le_bytes(p[8..].try_into().unwrap());
            PersistenceInterval::new(0, birth, IntervalEnd::Finite(death)).unwrap()
        })
        .collect();
    let a = PersistenceDiagram::new(0, Coverage::Complete, intervals[..n].to_vec())?;
    let b = PersistenceDiagram::new(0, Coverage::Complete, intervals[n..].to_vec())?;
    drop(intervals);
    drop(bytes);
    let distance = match args[2].as_str() {
        "bottleneck" => bottleneck_distance,
        "w1" => wasserstein_1_infinity,
        "w2" => wasserstein_2_euclidean,
        _ => panic!("unknown metric"),
    };
    let iterations: u64 = args[3].parse().unwrap();
    assert!(iterations > 0);
    #[cfg(cocycle_distance_bench)]
    let variant = args[4].parse().unwrap();
    #[cfg(cocycle_distance_bench)]
    cocycle::diagram_distances::select_preparation_experiment(variant, true);
    let first_start = Instant::now();
    let expected = black_box(distance(black_box(&a), black_box(&b), 0)?);
    let first_ns = first_start.elapsed().as_nanos();
    for _ in 0..2 {
        assert_eq!(
            black_box(distance(black_box(&a), black_box(&b), 0)?),
            expected
        );
    }
    let rss_before = memory("VmRSS:");
    let hwm_before = memory("VmHWM:");
    #[cfg(cocycle_distance_bench)]
    cocycle::diagram_distances::select_preparation_experiment(variant, false);
    let started = Instant::now();
    let mut value = 0.0;
    for _ in 0..iterations {
        value = black_box(distance(black_box(&a), black_box(&b), 0)?);
    }
    let elapsed_ns = started.elapsed().as_nanos();
    let peak_rss = memory("VmHWM:");
    assert_eq!(value, expected);
    #[cfg(cocycle_distance_bench)]
    let counters = cocycle::diagram_distances::preparation_experiment_phases();
    #[cfg(not(cocycle_distance_bench))]
    let counters = [0_u64; 10];
    println!(
        "{{\"elapsed_ns\":{elapsed_ns},\"iterations\":{iterations},\"value\":{value},\"first_ns\":{first_ns},\"rss_before_kib\":{rss_before},\"hwm_before_kib\":{hwm_before},\"peak_rss_kib\":{peak_rss},\"counters\":{counters:?}}}"
    );
    Ok(())
}
