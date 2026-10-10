// Inserted only into an artifact source copy by reproduce_distance_preparation.py.
// The ordinary library never compiles these experimental entry points.
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};

static PREPARATION: AtomicUsize = AtomicUsize::new(0);
static VERIFY: AtomicBool = AtomicBool::new(false);
static PHASES: [AtomicU64; 10] = [const { AtomicU64::new(0) }; 10];

pub(crate) fn select_experiment(variant: usize, verify: bool) {
    assert!(variant < 7);
    PREPARATION.store(variant, Ordering::Relaxed);
    VERIFY.store(verify, Ordering::Relaxed);
    for value in &PHASES {
        value.store(0, Ordering::Relaxed);
    }
}

pub(crate) fn experiment_phases() -> [u64; 10] {
    std::array::from_fn(|i| PHASES[i].load(Ordering::Relaxed))
}

fn reserved(count: usize) -> Result<Vec<[f64; 2]>> {
    let mut points = Vec::new();
    points.try_reserve_exact(count).map_err(|_| allocation())?;
    Ok(points)
}

fn coordinates<const CONTROLLED: bool>(
    view: &DiagramDimension<'_>,
    count: usize,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<Vec<[f64; 2]>> {
    let mut points = reserved(count)?;
    for interval in view.iter() {
        budget.step()?;
        if let crate::diagram::IntervalEnd::Finite(death) = interval.end() {
            points.push([interval.birth(), death]);
        }
    }
    Ok(points)
}

pub(crate) fn from_dimensions<const CONTROLLED: bool>(
    first: &DiagramDimension<'_>,
    second: &DiagramDimension<'_>,
    counts: (usize, usize),
    metric: Metric,
    stats: &mut Stats,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<f64> {
    let started = std::time::Instant::now();
    let variant = PREPARATION.load(Ordering::Relaxed);
    let mut copies = if variant == 0 || variant == 6 {
        None
    } else if variant >= 4 {
        Some((reserved(counts.0)?, reserved(counts.1)?))
    } else {
        Some((
            coordinates(first, counts.0, budget)?,
            coordinates(second, counts.1, budget)?,
        ))
    };
    let (a, b, scale) = if variant == 1 || variant == 2 {
        let (left, right) = copies.as_ref().unwrap();
        numeric::prepare_coordinates(left, right, counts, budget)?
    } else {
        prepare_dimensions(first, second, counts, budget)?
    };
    if variant == 6 {
        copies = Some((reserved(counts.0)?, reserved(counts.1)?));
    }
    if VERIFY.swap(false, Ordering::Relaxed) {
        let (expected_a, expected_b, expected_scale) =
            prepare_dimensions(first, second, counts, budget)?;
        assert_eq!(scale.to_bits(), expected_scale.to_bits());
        for (actual, expected) in [(&a, &expected_a), (&b, &expected_b)] {
            assert_eq!(actual.len(), expected.len());
            for (p, q) in actual.iter().zip(expected) {
                let bits = |p: &Point| {
                    [
                        p.coordinates[0].to_bits(),
                        p.coordinates[1].to_bits(),
                        p.midpoint.to_bits(),
                        p.half.to_bits(),
                    ]
                };
                assert_eq!(bits(p), bits(q));
            }
        }
    }
    if variant == 2 || variant == 5 {
        drop(copies.take());
    }
    PHASES[0].fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
    let started = std::time::Instant::now();
    // Every runtime-selected path reaches this same non-inlined solver body.
    let result = solve_prepared(a, b, scale, metric, Options::default(), stats, budget);
    PHASES[1].fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
    let started = std::time::Instant::now();
    std::hint::black_box(&copies);
    drop(copies);
    PHASES[2].fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
    PHASES[3].fetch_add(1, Ordering::Relaxed);
    result
}

// Growth diagnostics: inserted only for the separately timed growth mode.
fn trace_edge_push(values: &mut Vec<graph::Edge>, value: graph::Edge) -> Result<()> {
    let capacity = values.capacity();
    let pointer = values.as_ptr();
    let old_bytes = values.len() * std::mem::size_of::<graph::Edge>();
    push(values, value)?;
    if capacity != values.capacity() {
        PHASES[7].fetch_add(1, Ordering::Relaxed);
        if capacity != 0 && pointer != values.as_ptr() {
            PHASES[8].fetch_add(1, Ordering::Relaxed);
            PHASES[9].fetch_add(old_bytes as u64, Ordering::Relaxed);
        }
    }
    Ok(())
}
