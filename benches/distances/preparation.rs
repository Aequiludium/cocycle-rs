//! Warmed public calls for the distance preparation reproducer.
#![forbid(unsafe_code)]

use cocycle::diagram::{Coverage, IntervalEnd, PersistenceDiagram, PersistenceInterval};
use cocycle::diagram_distances::{
    bottleneck_distance, wasserstein_1_infinity, wasserstein_2_euclidean,
};
use std::{env, fs, hint::black_box, time::Instant};

fn parse_fixture(bytes: &[u8]) -> cocycle::Result<(PersistenceDiagram, PersistenceDiagram)> {
    if bytes.len() < 24 || &bytes[..8] != b"COCDST1\0" {
        return Err(cocycle::Error::InvalidParameter {
            parameter: "fixture",
            reason: "expected a complete COCDST1 header",
        });
    }
    // The complete header above establishes both fixed-width count slices.
    let count = |range: std::ops::Range<usize>| {
        usize::try_from(u64::from_le_bytes(bytes[range].try_into().unwrap())).map_err(|_| {
            cocycle::Error::SizeOverflow {
                operation: "fixture count",
            }
        })
    };
    let n = count(8..16)?;
    let m = count(16..24)?;
    let expected = n
        .checked_add(m)
        .and_then(|count| count.checked_mul(16))
        .and_then(|size| size.checked_add(24))
        .ok_or(cocycle::Error::SizeOverflow {
            operation: "fixture length",
        })?;
    if bytes.len() != expected {
        return Err(cocycle::Error::ShapeMismatch {
            input: "fixture bytes",
            expected,
            actual: bytes.len(),
        });
    }
    let mut intervals = Vec::new();
    intervals
        .try_reserve_exact(n + m)
        .map_err(|_| cocycle::Error::AllocationFailed {
            context: "fixture intervals",
        })?;
    for point in bytes[24..].chunks_exact(16) {
        let birth = f64::from_le_bytes(point[..8].try_into().unwrap());
        let death = f64::from_le_bytes(point[8..].try_into().unwrap());
        intervals.push(PersistenceInterval::new(
            0,
            birth,
            IntervalEnd::Finite(death),
        )?);
    }
    Ok((
        PersistenceDiagram::new(0, Coverage::Complete, intervals[..n].to_vec())?,
        PersistenceDiagram::new(0, Coverage::Complete, intervals[n..].to_vec())?,
    ))
}

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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args().collect();
    assert_eq!(args.len(), 5, "fixture metric iterations variant");
    let bytes = fs::read(&args[1])?;
    let (a, b) = parse_fixture(&bytes)?;
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

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(n: u64, m: u64, points: &[[f64; 2]]) -> Vec<u8> {
        let mut bytes = b"COCDST1\0".to_vec();
        bytes.extend(n.to_le_bytes());
        bytes.extend(m.to_le_bytes());
        for [birth, death] in points {
            bytes.extend(birth.to_le_bytes());
            bytes.extend(death.to_le_bytes());
        }
        bytes
    }

    #[test]
    fn rejects_every_truncated_header_and_wrong_magic() {
        let mut bytes = fixture(0, 0, &[]);
        for end in 0..24 {
            assert!(matches!(
                parse_fixture(&bytes[..end]),
                Err(cocycle::Error::InvalidParameter { .. })
            ));
        }
        bytes[0] = b'X';
        assert!(matches!(
            parse_fixture(&bytes),
            Err(cocycle::Error::InvalidParameter { .. })
        ));
    }

    #[test]
    fn rejects_count_overflow_before_allocation() {
        for (n, m) in [(u64::MAX, 1), (u64::MAX / 16, 0), (u64::MAX / 16 + 1, 0)] {
            assert!(matches!(
                parse_fixture(&fixture(n, m, &[])),
                Err(cocycle::Error::SizeOverflow { .. })
            ));
        }
    }

    #[test]
    fn rejects_truncated_payload_and_trailing_bytes() {
        let mut bytes = fixture(1, 0, &[[0., 2.]]);
        for end in 24..bytes.len() {
            assert!(matches!(
                parse_fixture(&bytes[..end]),
                Err(cocycle::Error::ShapeMismatch { .. })
            ));
        }
        bytes.push(0);
        assert!(matches!(
            parse_fixture(&bytes),
            Err(cocycle::Error::ShapeMismatch { .. })
        ));
    }

    #[test]
    fn rejects_invalid_endpoints() {
        for point in [[2., 1.], [1., 1.], [f64::NAN, 2.], [0., f64::INFINITY]] {
            assert!(parse_fixture(&fixture(1, 0, &[point])).is_err());
        }
    }

    #[test]
    fn keeps_empty_sides_and_interval_multiplicity() {
        for (n, m) in [(0, 0), (2, 0), (0, 2), (1, 1)] {
            let points = vec![[0., 2.]; (n + m) as usize];
            let (a, b) = parse_fixture(&fixture(n, m, &points)).unwrap();
            assert_eq!(a.len(), n as usize);
            assert_eq!(b.len(), m as usize);
            let expected = if n == m { 0. } else { 1. };
            assert_eq!(bottleneck_distance(&a, &b, 0).unwrap(), expected);
            assert_eq!(wasserstein_1_infinity(&a, &b, 0).unwrap(), 2. * expected);
            assert_eq!(wasserstein_2_euclidean(&a, &b, 0).unwrap(), 2. * expected);
        }
    }
}
