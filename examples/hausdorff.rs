//! Direct point-set comparison, identified sample matrices and oracle fixtures.
use cocycle::geometry::{PointCloudView, hausdorff_distance, hausdorff_distance_matrix};
use cocycle::{Error, Result};

fn main() -> Result<()> {
    match std::env::args().nth(1).as_deref() {
        Some("--fixtures") => fixtures(),
        Some("--measure") => measure(),
        None => workflow(),
        _ => Err(Error::InvalidParameter {
            parameter: "argument",
            reason: "expected no argument, --fixtures or --measure",
        }),
    }
}

fn workflow() -> Result<()> {
    let a = PointCloudView::new(&[0., 2.], 2, 1)?;
    let b = PointCloudView::new(&[1., 3.], 2, 1)?;
    let c = PointCloudView::new(&[4., 6.], 2, 1)?;
    println!("sample-a to sample-b: {}", hausdorff_distance(a, b)?);
    let matrix = hausdorff_distance_matrix(&[("sample-a", a), ("sample-b", b), ("sample-c", c)])?;
    println!("both axes: {:?}", matrix.sample_ids());
    println!("lower triangle: {:?}", matrix.values());
    let nearest = (1..matrix.len())
        .map(|j| (j, matrix.get(0, j).unwrap()))
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .unwrap();
    assert_eq!(nearest, (1, 1.));
    println!(
        "nearest to {}: {} at {}",
        matrix.sample_ids()[0],
        matrix.sample_ids()[nearest.0],
        nearest.1
    );
    Ok(())
}

fn emit(name: &str, first: &[f64], second: &[f64], dimension: usize) -> Result<()> {
    let first_view = PointCloudView::new(first, first.len() / dimension, dimension)?;
    let second_view = PointCloudView::new(second, second.len() / dimension, dimension)?;
    let distance = match hausdorff_distance(first_view, second_view) {
        Ok(value) => format!("{value:?}"),
        Err(Error::NumericalFailure { .. }) => String::from("\"numerical_failure\""),
        Err(error) => return Err(error),
    };
    println!(
        "{{\"name\":\"{name}\",\"dimension\":{dimension},\"first\":{first:?},\"second\":{second:?},\"distance\":{distance}}}"
    );
    Ok(())
}

fn fixtures() -> Result<()> {
    emit("single_345", &[0., 0.], &[3., 4.], 2)?;
    emit("asymmetric", &[0., 2.], &[0.], 1)?;
    emit("duplicates", &[0., 0., 2.], &[2., 0.], 1)?;
    emit("translated", &[0., 2.], &[1., 3.], 1)?;
    emit("outlier", &[0., 1.], &[0., 1., 100.], 1)?;
    let mut state = 110_u64;
    let mut next = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        state >> 32
    };
    for dimension in [1, 2, 3, 8] {
        for case in 0..20 {
            let m = (next() % 7 + 1) as usize;
            let n = (next() % 7 + 1) as usize;
            let a: Vec<_> = (0..m * dimension)
                .map(|_| (next() % 17) as f64 - 8.)
                .collect();
            let b: Vec<_> = (0..n * dimension)
                .map(|_| (next() % 17) as f64 - 8.)
                .collect();
            emit(&format!("random_{dimension}_{case}"), &a, &b, dimension)?;
        }
    }
    emit("stress_large", &[0.], &[1e200], 1)?;
    emit("stress_tiny", &[0.], &[1e-200], 1)?;
    emit("stress_subnormal", &[0.], &[f64::from_bits(1)], 1)?;
    emit(
        "stress_irrelevant_overflow",
        &[-f64::MAX, f64::MAX],
        &[f64::MAX, -f64::MAX],
        1,
    )?;
    emit("stress_final_overflow", &[-f64::MAX], &[f64::MAX], 1)?;
    let a = PointCloudView::new(&[0.], 1, 1)?;
    let b = PointCloudView::new(&[1.], 1, 1)?;
    let c = PointCloudView::new(&[3.], 1, 1)?;
    let matrix = hausdorff_distance_matrix(&[("c", c), ("a", a), ("b", b), ("same", a)])?;
    println!(
        "{{\"matrix_ids\":{:?},\"matrix_values\":{:?}}}",
        matrix.sample_ids(),
        matrix.values()
    );
    Ok(())
}

fn measure() -> Result<()> {
    // Resource snapshot only. Generation/validation precede the timers; no
    // reference-engine timing or before/after speedup is implied.
    for (n, dimension, repeats) in [(32, 2, 20), (128, 3, 5), (512, 16, 2)] {
        let a: Vec<_> = (0..n * dimension).map(|i| (i % 97) as f64).collect();
        let b: Vec<_> = a.iter().map(|x| x + 0.25).collect();
        let first = PointCloudView::new(&a, n, dimension)?;
        let second = PointCloudView::new(&b, n, dimension)?;
        let started = std::time::Instant::now();
        let mut value = 0.;
        for _ in 0..repeats {
            value = std::hint::black_box(hausdorff_distance(
                std::hint::black_box(first),
                std::hint::black_box(second),
            )?);
        }
        println!(
            "{{\"operation\":\"pair\",\"points_per_cloud\":{n},\"dimension\":{dimension},\"repeats\":{repeats},\"seconds_per_call\":{},\"distance\":{value}}}",
            started.elapsed().as_secs_f64() / f64::from(repeats)
        );
    }
    let coordinates: Vec<Vec<f64>> = (0..12)
        .map(|s| {
            (0..96)
                .map(|i| f64::from(i) + f64::from(s) * 0.25)
                .collect()
        })
        .collect();
    let labels: Vec<_> = (0..12).map(|s| format!("sample-{s}")).collect();
    let samples: Vec<_> = labels
        .iter()
        .zip(&coordinates)
        .map(|(id, data)| Ok((id.as_str(), PointCloudView::new(data, 32, 3)?)))
        .collect::<Result<_>>()?;
    let started = std::time::Instant::now();
    let matrix = std::hint::black_box(hausdorff_distance_matrix(&samples)?);
    println!(
        "{{\"operation\":\"matrix\",\"samples\":12,\"points_per_cloud\":32,\"dimension\":3,\"seconds_per_call\":{},\"stored_values\":{},\"numeric_bytes\":{},\"copied_id_bytes\":{}}}",
        started.elapsed().as_secs_f64(),
        matrix.values().len(),
        std::mem::size_of_val(matrix.values()),
        matrix.sample_ids().iter().map(String::len).sum::<usize>()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn direct_comparison_to_nearest_sample() {
        super::workflow().unwrap();
    }
}
