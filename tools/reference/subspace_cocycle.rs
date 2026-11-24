//! Unoptimized controlled sequence-subspace rewrite; outside the product API.
//! Frozen NumPy inputs, not new random inputs, drive all parity comparisons.
pub use cocycle::{Error, Result};
#[allow(dead_code)]
#[path = "../../src/eigensolver.rs"]
mod eigensolver;
#[allow(dead_code)]
#[path = "../../src/execution/mod.rs"]
mod execution;

fn canonical_zero(value: f64) -> f64 {
    if value == 0. { 0. } else { value }
}

fn read<T: std::str::FromStr>(words: &mut std::str::SplitWhitespace<'_>) -> T
where
    T::Err: std::fmt::Debug,
{
    words
        .next()
        .expect("missing fixture field")
        .parse()
        .expect("fixture field")
}

type ControlResult<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn norm(values: &[f64]) -> f64 {
    values.iter().fold(0., |a, &v| a.hypot(v))
}

fn filter(
    a: &[f64],
    y: &[f64],
    n: usize,
    b: usize,
    degrees: &[usize],
    parameters: [f64; 3],
) -> ControlResult<Vec<f64>> {
    let [lower, upper, anchor] = parameters;
    if !anchor.is_finite()
        || !lower.is_finite()
        || !upper.is_finite()
        || !(anchor < lower && lower < upper)
        || degrees.iter().any(|&m| m > 64)
    {
        return Err("filter bounds/degree".into());
    }
    let center = (lower + upper) / 2.;
    let radius = (upper - lower) / 2.;
    let sigma1 = radius / (anchor - center);
    let mut output = y.to_vec();
    for column in 0..b {
        let mut previous: Vec<f64> = (0..n).map(|i| y[i * b + column]).collect();
        let mut current = previous.clone();
        let mut sigma = sigma1;
        for degree in 0..degrees[column] {
            let mut next = vec![0.; n];
            let next_sigma = if degree == 0 {
                sigma1
            } else {
                1. / (2. / sigma1 - sigma)
            };
            for row in 0..n {
                let image: f64 = (0..n).map(|j| a[row * n + j] * current[j]).sum();
                let image = (image - center * current[row]) / radius;
                next[row] = if degree == 0 {
                    sigma1 * image
                } else {
                    2. * next_sigma * image - next_sigma * sigma * previous[row]
                };
            }
            previous = current;
            current = next;
            sigma = next_sigma;
        }
        for row in 0..n {
            output[row * b + column] = current[row];
        }
    }
    if output.iter().any(|v| !v.is_finite()) {
        return Err("nonfinite filter output".into());
    }
    Ok(output)
}

fn order_fields(
    fields: &[f64],
    samples: usize,
    p: usize,
    truncation: usize,
) -> ControlResult<(Vec<usize>, Vec<f64>, Vec<f64>)> {
    if p == 0
        || p > 32
        || truncation == 0
        || truncation > p
        || samples > 64
        || fields.len() != samples * p * p
        || fields.iter().any(|x| !x.is_finite())
    {
        return Err("FFT field shape/truncation/finiteness".into());
    }
    let width = truncation * truncation;
    let mut real = vec![0.; samples * width];
    let mut imag = real.clone();
    // Direct DFT is the simple reference representation of the same unshifted,
    // unnormalized FFT2 descriptor. It makes no FFT-complexity/speedup claim.
    for sample in 0..samples {
        for u in 0..truncation {
            for v in 0..truncation {
                let index = sample * width + u * truncation + v;
                for x in 0..p {
                    for y in 0..p {
                        let angle = -std::f64::consts::TAU * (u * x + v * y) as f64 / p as f64;
                        let value = fields[sample * p * p + x * p + y];
                        real[index] += value * angle.cos();
                        imag[index] += value * angle.sin();
                    }
                }
            }
        }
    }
    if real.iter().chain(&imag).any(|v| !v.is_finite()) {
        return Err("nonfinite FFT descriptor".into());
    }
    let mut order = if samples == 0 { vec![] } else { vec![0] };
    while order.len() < samples {
        let previous = order[order.len() - 1];
        let mut best = None;
        for j in 0..samples {
            if order.contains(&j) {
                continue;
            }
            let mut distance = 0_f64;
            for k in 0..width {
                distance = distance.hypot(real[previous * width + k] - real[j * width + k]);
                distance = distance.hypot(imag[previous * width + k] - imag[j * width + k]);
            }
            if best.is_none_or(|(_, d)| distance < d) {
                best = Some((j, distance));
            }
        }
        order.push(best.ok_or("ordering progress")?.0);
    }
    Ok((order, real, imag))
}

fn minimum_singular_value(input: &[f64], n: usize, b: usize) -> ControlResult<f64> {
    let scale = input.iter().fold(0_f64, |s, x| s.max(x.abs()));
    if scale == 0. {
        return Ok(0.);
    }
    let mut a: Vec<f64> = input.iter().map(|v| v / scale).collect();
    // One-sided Jacobi avoids squaring the 1e-13 singular-value rank threshold.
    for _ in 0..128 {
        let mut changed = false;
        for p in 0..b {
            for q in p + 1..b {
                let alpha: f64 = (0..n).map(|i| a[i * b + p].powi(2)).sum();
                let beta: f64 = (0..n).map(|i| a[i * b + q].powi(2)).sum();
                let gamma: f64 = (0..n).map(|i| a[i * b + p] * a[i * b + q]).sum();
                if gamma.abs() <= 4. * f64::EPSILON * alpha.sqrt() * beta.sqrt() {
                    continue;
                }
                let angle = 0.5 * (2. * gamma).atan2(beta - alpha);
                let (s, c) = angle.sin_cos();
                for i in 0..n {
                    let (x, y) = (a[i * b + p], a[i * b + q]);
                    a[i * b + p] = c * x - s * y;
                    a[i * b + q] = s * x + c * y;
                }
                changed = true;
            }
        }
        if !changed {
            return Ok((0..b)
                .map(|j| (0..n).fold(0_f64, |s, i| s.hypot(a[i * b + j])))
                .fold(f64::INFINITY, f64::min)
                * scale);
        }
    }
    Err("rank SVD convergence".into())
}

fn thin_qr(input: &[f64], n: usize, b: usize) -> ControlResult<Vec<f64>> {
    let mut a = input.to_vec();
    let mut reflectors = Vec::new();
    for j in 0..b {
        let mut v: Vec<f64> = (j..n).map(|i| a[i * b + j]).collect();
        let length = norm(&v);
        if length == 0. || !length.is_finite() {
            return Err("QR rank".into());
        }
        v[0] += length.copysign(v[0]);
        let length = norm(&v);
        for x in &mut v {
            *x /= length;
        }
        for column in j..b {
            let dot: f64 = (j..n).map(|i| v[i - j] * a[i * b + column]).sum();
            for i in j..n {
                a[i * b + column] -= 2. * v[i - j] * dot;
            }
        }
        reflectors.push(v);
    }
    let mut q = vec![0.; n * b];
    for j in 0..b {
        q[j * b + j] = 1.;
    }
    for j in (0..b).rev() {
        let v = &reflectors[j];
        for column in 0..b {
            let dot: f64 = (j..n).map(|i| v[i - j] * q[i * b + column]).sum();
            for i in j..n {
                q[i * b + column] -= 2. * v[i - j] * dot;
            }
        }
    }
    Ok(q)
}

fn orthogonal_ritz(
    a: &[f64],
    filtered: &[f64],
    locked: &[f64],
    n: usize,
    width: usize,
    l: usize,
    budget: &mut execution::WorkBudget<'_>,
) -> ControlResult<(Vec<f64>, Vec<f64>, Vec<f64>, String)> {
    let b = width + l;
    let mut normalized = filtered.to_vec();
    for j in 0..width {
        let length = (0..n).fold(0_f64, |s, i| s.hypot(filtered[i * width + j]));
        if length == 0. || !length.is_finite() {
            return Err("filtered rank".into());
        }
        for i in 0..n {
            normalized[i * width + j] /= length;
        }
    }
    let mut joined = vec![0.; n * b];
    for i in 0..n {
        joined[i * b..i * b + l].copy_from_slice(&locked[i * l..(i + 1) * l]);
        joined[i * b + l..(i + 1) * b].copy_from_slice(&normalized[i * width..(i + 1) * width]);
    }
    if minimum_singular_value(&joined, n, b)? < 1e-13 {
        return Err("filtered rank".into());
    }
    let full_q = thin_qr(&joined, n, b)?;
    let q: Vec<f64> = (0..n)
        .flat_map(|i| full_q[i * b + l..(i + 1) * b].to_vec())
        .collect();
    let image: Vec<f64> = (0..n)
        .flat_map(|i| {
            (0..width)
                .map(|j| (0..n).map(|k| a[i * n + k] * q[k * width + j]).sum())
                .collect::<Vec<_>>()
        })
        .collect();
    let mut quotient = vec![0.; width * width];
    for i in 0..width {
        for j in i..width {
            let left: f64 = (0..n)
                .map(|k| q[k * width + i] * image[k * width + j])
                .sum();
            let right: f64 = (0..n)
                .map(|k| q[k * width + j] * image[k * width + i])
                .sum();
            quotient[i * width + j] = (left + right) / 2.;
            quotient[j * width + i] = quotient[i * width + j];
        }
    }
    let (values, rotations, _) =
        eigensolver::jacobi(&quotient, width, 0., 4. * f64::EPSILON, 64, budget)?;
    let mut order: Vec<usize> = (0..width).collect();
    order.sort_by(|&i, &j| values[i].total_cmp(&values[j]).then(i.cmp(&j)));
    let values: Vec<f64> = order.iter().map(|&i| values[i]).collect();
    let mut ritz = vec![0.; n * width];
    for i in 0..n {
        for (j, &index) in order.iter().enumerate() {
            ritz[i * width + j] = (0..width)
                .map(|k| q[i * width + k] * rotations[k * width + index])
                .sum();
        }
    }
    let residuals = fresh_residuals(a, &ritz, &values, n);
    let trace = format!(
        "{{\"locked_shape\":[{n},{l}],\"active_shape\":[{n},{width}],\"qr_shape\":[{n},{width}],\"ritz_shape\":[{width},{width}],\"q\":{q:?},\"quotient\":{quotient:?},\"ritz\":{ritz:?},\"values\":{values:?},\"fresh_residuals\":{residuals:?}}}"
    );
    Ok((values, ritz, residuals, trace))
}

fn fresh_residuals(a: &[f64], vectors: &[f64], values: &[f64], n: usize) -> Vec<f64> {
    let b = values.len();
    (0..b)
        .map(|j| {
            (0..n).fold(0_f64, |s, i| {
                let image: f64 = (0..n).map(|k| a[i * n + k] * vectors[k * b + j]).sum();
                s.hypot(image - values[j] * vectors[i * b + j])
            })
        })
        .collect()
}

fn transfer(
    block: &[f64],
    old: &[(Vec<u64>, i32)],
    new: &[(Vec<u64>, i32)],
    fresh: &[f64],
    limits: [usize; 3],
    projector: Option<&[f64]>,
) -> ControlResult<(Vec<f64>, usize)> {
    let (old_n, n) = (old.len(), new.len());
    let [k, width, max_entries] = limits;
    if block.len() != old_n * k
        || width > n
        || k > 16
        || width > 16
        || n * width.max(k) > max_entries
        || fresh.len() != n * width
        || block.iter().chain(fresh).any(|x| !x.is_finite())
    {
        return Err("transfer shape/finiteness/storage".into());
    }
    if width == 0 {
        return Ok((Vec::new(), 0));
    }
    if let Some(p) = projector {
        if p.len() != n * n
            || n * n > max_entries
            || p.iter().any(|v| !v.is_finite())
            || (0..n).any(|i| {
                (0..n).any(|j| {
                    (p[i * n + j] - p[j * n + i]).abs() > 1e-10
                        || ((0..n).map(|r| p[i * n + r] * p[r * n + j]).sum::<f64>() - p[i * n + j])
                            .abs()
                            > 1e-10
                })
            })
        {
            return Err("new range projector".into());
        }
    }
    let project = |v: Vec<f64>| -> Vec<f64> {
        if let Some(p) = projector {
            (0..n)
                .map(|i| (0..n).map(|j| p[i * n + j] * v[j]).sum())
                .collect()
        } else {
            v
        }
    };
    let mut mapped = vec![vec![0.; n]; k];
    for (j, (simplex, sign)) in new.iter().enumerate() {
        if let Some(i) = old.iter().position(|(s, _)| s == simplex) {
            for column in 0..k {
                mapped[column][j] = (*sign * old[i].1) as f64 * block[i * k + column];
            }
        }
    }
    for column in &mut mapped {
        *column = project(column.clone());
    }
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let mut retained_rank = 0;
    for candidate in 0..k + width + n {
        let v = if candidate < k {
            mapped[candidate].clone()
        } else if candidate < k + width {
            fresh[(candidate - k) * n..(candidate - k + 1) * n].to_vec()
        } else {
            let mut v = vec![0.; n];
            v[candidate - k - width] = 1.;
            v
        };
        let mut v = project(v);
        let length = norm(&v);
        if length != 0. {
            for x in &mut v {
                *x /= length;
            }
        }
        for _ in 0..2 {
            for q in &basis {
                let dot: f64 = q.iter().zip(&v).map(|(a, b)| a * b).sum();
                for i in 0..n {
                    v[i] -= dot * q[i];
                }
            }
        }
        let length = norm(&v);
        if length > 1e-13 {
            for x in &mut v {
                *x /= length;
            }
            basis.push(v);
            if candidate < k {
                retained_rank += 1;
            }
        }
        if basis.len() == width {
            break;
        }
    }
    if basis.len() != width {
        return Err("projected rank insufficient; caller must choose cold/full route".into());
    }
    let output = (0..n)
        .flat_map(|i| basis.iter().map(|v| v[i]).collect::<Vec<_>>())
        .collect();
    Ok((output, retained_rank))
}

struct SolvedSubspace {
    values: Vec<f64>,
    vectors: Vec<f64>,
    retained: Vec<f64>,
    residuals: Vec<f64>,
    actions: usize,
    trace: Vec<String>,
}

impl SolvedSubspace {
    fn json(&self) -> String {
        format!(
            "{{\"values\":{:?},\"vectors\":{:?},\"retained\":{:?},\"residuals\":{:?},\"actions\":{},\"trace\":[{}]}}",
            self.values,
            self.vectors,
            self.retained,
            self.residuals,
            self.actions,
            self.trace.join(",")
        )
    }
}

fn chfsi(
    a: &[f64],
    initial: &[f64],
    n: usize,
    b: usize,
    controls: [usize; 3],
    parameters: [f64; 3],
    tolerance: f64,
) -> ControlResult<SolvedSubspace> {
    let [count, degree, max_iterations] = controls;
    if n == 0
        || n > 64
        || b == 0
        || b > n
        || b > 16
        || count == 0
        || count > b
        || a.len() != n * n
        || initial.len() != n * b
        || degree == 0
        || degree > 64
        || max_iterations == 0
        || max_iterations > 200
        || !tolerance.is_finite()
        || tolerance <= 0.
        || a.iter().chain(initial).any(|v| !v.is_finite())
        || (0..n).any(|i| (0..n).any(|j| (a[i * n + j] - a[j * n + i]).abs() > 1e-14))
    {
        return Err("ChFSI input".into());
    }
    if minimum_singular_value(initial, n, b)? <= 1e-13 {
        return Err("initial rank".into());
    }
    let mut active = thin_qr(initial, n, b)?;
    let mut locked = Vec::new();
    let mut locked_values: Vec<f64> = Vec::new();
    let mut trace = Vec::new();
    let mut actions = 0;
    let mut budget = execution::WorkBudget::new(&execution::Execution::default())?;
    for iteration in 0..max_iterations {
        let l = locked_values.len();
        let width = b - l;
        let filtered = filter(a, &active, n, width, &vec![degree; width], parameters)?;
        let (values, ritz, errors, step) =
            orthogonal_ritz(a, &filtered, &locked, n, width, l, &mut budget)?;
        actions += (degree + 2) * width;
        let mut take = 0;
        while take < width.min(count - l) && errors[take] <= tolerance * values[take].abs().max(1.)
        {
            take += 1;
        }
        trace.push(format!("{{\"iteration\":{iteration},\"newly_locked\":{take},\"filter_input\":{active:?},\"filter_output\":{filtered:?},\"step\":{step}}}"));
        let next_l = l + take;
        let mut next_locked = vec![0.; n * next_l];
        let mut remaining = vec![0.; n * (width - take)];
        for i in 0..n {
            next_locked[i * next_l..i * next_l + l].copy_from_slice(&locked[i * l..(i + 1) * l]);
            next_locked[i * next_l + l..(i + 1) * next_l]
                .copy_from_slice(&ritz[i * width..i * width + take]);
            remaining[i * (width - take)..(i + 1) * (width - take)]
                .copy_from_slice(&ritz[i * width + take..(i + 1) * width]);
        }
        locked = next_locked;
        locked_values.extend_from_slice(&values[..take]);
        if next_l == count {
            let residuals = fresh_residuals(a, &locked, &locked_values, n);
            actions += count;
            if residuals
                .iter()
                .zip(&locked_values)
                .any(|(&r, v)| !r.is_finite() || r > tolerance * v.abs().max(1.))
            {
                return Err("fresh final residual".into());
            }
            let mut retained = vec![0.; n * b];
            for i in 0..n {
                retained[i * b..i * b + next_l]
                    .copy_from_slice(&locked[i * next_l..(i + 1) * next_l]);
                retained[i * b + next_l..(i + 1) * b]
                    .copy_from_slice(&remaining[i * (b - next_l)..(i + 1) * (b - next_l)]);
            }
            return Ok(SolvedSubspace {
                values: locked_values,
                vectors: locked,
                retained,
                residuals,
                actions,
                trace,
            });
        }
        active = remaining;
    }
    Err("ChFSI convergence".into())
}

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let text = std::fs::read_to_string(std::env::args().nth(1).ok_or("fixture path")?)?;
    let mut words = text.split_whitespace();
    let mode = words.next().ok_or("mode")?;
    if mode == "transfer" {
        let old_source = words.next().ok_or("source")?;
        let new_source = words.next().ok_or("source")?;
        let old_degree: usize = read(&mut words);
        let new_degree: usize = read(&mut words);
        let old_metric = words.next().ok_or("metric")?;
        let new_metric = words.next().ok_or("metric")?;
        if old_source != new_source
            || old_degree != new_degree
            || old_metric != "euclidean"
            || new_metric != "euclidean"
        {
            return Err("incompatible source/degree/inner product".into());
        }
        let old_n: usize = read(&mut words);
        let n: usize = read(&mut words);
        let k: usize = read(&mut words);
        let width: usize = read(&mut words);
        let max_entries: usize = read(&mut words);
        let has_projection: usize = read(&mut words);
        if old_degree > 16 || n > 256 || old_n > 256 || k > 16 || width > 16 || has_projection > 1 {
            return Err("transfer control limit".into());
        }
        let mut bases = Vec::new();
        for length in [old_n, n] {
            let mut basis: Vec<(Vec<u64>, i32)> = Vec::new();
            for _ in 0..length {
                let sign: i32 = read(&mut words);
                let simplex: Vec<u64> = (0..=old_degree).map(|_| read(&mut words)).collect();
                if ![-1, 1].contains(&sign)
                    || simplex.windows(2).any(|w| w[0] >= w[1])
                    || basis.iter().any(|(s, _)| s == &simplex)
                {
                    return Err("invalid oriented simplex basis".into());
                }
                basis.push((simplex, sign));
            }
            bases.push(basis);
        }
        let block: Vec<f64> = (0..old_n * k).map(|_| read(&mut words)).collect();
        let fresh: Vec<f64> = (0..n * width).map(|_| read(&mut words)).collect();
        let projector: Vec<f64> = if has_projection == 1 {
            (0..n * n).map(|_| read(&mut words)).collect()
        } else {
            Vec::new()
        };
        if words.next().is_some() {
            return Err("transfer trailing fields".into());
        }
        let (mapped, rank) = transfer(
            &block,
            &bases[0],
            &bases[1],
            &fresh,
            [k, width, max_entries],
            (has_projection == 1).then_some(projector.as_slice()),
        )?;
        println!("{{\"mapped\":{mapped:?},\"rank\":{rank}}}");
        return Ok(());
    }
    if mode == "order" {
        let samples: usize = read(&mut words);
        let p: usize = read(&mut words);
        let truncation: usize = read(&mut words);
        if p > 32 || samples > 64 {
            return Err("ordering size limit".into());
        }
        let fields: Vec<f64> = words
            .map(str::parse)
            .collect::<std::result::Result<_, _>>()?;
        let (order, real, imag) = order_fields(&fields, samples, p, truncation)?;
        println!(
            "{{\"order\":{order:?},\"descriptors_real\":{real:?},\"descriptors_imag\":{imag:?}}}"
        );
        return Ok(());
    }
    let n: usize = read(&mut words);
    let b: usize = read(&mut words);
    if n == 0 || n > 256 || b == 0 || b > n {
        return Err("control shape: require 0 < b <= n <= 256".into());
    }
    if mode == "chfsi" || mode == "step" || mode == "sequence" {
        let count: usize = read(&mut words);
        let degree: usize = read(&mut words);
        let max_iterations: usize = read(&mut words);
        let tolerance: f64 = read(&mut words);
        let parameters: [f64; 3] = [read(&mut words), read(&mut words), read(&mut words)];
        if n > 64 || b > 16 {
            return Err("ChFSI control limit".into());
        }
        let mut budget = execution::WorkBudget::new(&execution::Execution::default())?;
        if mode == "step" {
            if count >= b {
                return Err("step locked shape".into());
            }
            let a: Vec<f64> = (0..n * n).map(|_| read(&mut words)).collect();
            let active: Vec<f64> = (0..n * (b - count)).map(|_| read(&mut words)).collect();
            let locked: Vec<f64> = (0..n * count).map(|_| read(&mut words)).collect();
            if words.next().is_some()
                || a.iter()
                    .chain(&active)
                    .chain(&locked)
                    .any(|v| !v.is_finite())
            {
                return Err("step input".into());
            }
            let filtered = filter(
                &a,
                &active,
                n,
                b - count,
                &vec![degree; b - count],
                parameters,
            )?;
            let (_, _, _, step) =
                orthogonal_ritz(&a, &filtered, &locked, n, b - count, count, &mut budget)?;
            println!("{step}");
        } else if mode == "chfsi" {
            let a: Vec<f64> = (0..n * n).map(|_| read(&mut words)).collect();
            let initial: Vec<f64> = (0..n * b).map(|_| read(&mut words)).collect();
            if words.next().is_some() {
                return Err("ChFSI trailing fields".into());
            }
            println!(
                "{}",
                chfsi(
                    &a,
                    &initial,
                    n,
                    b,
                    [count, degree, max_iterations],
                    parameters,
                    tolerance
                )?
                .json()
            );
        } else {
            let samples: usize = read(&mut words);
            let p: usize = read(&mut words);
            let truncation: usize = read(&mut words);
            let route: usize = read(&mut words);
            if samples == 0 || samples > 64 || p > 32 || route > 2 {
                return Err("sequence input".into());
            }
            let fields: Vec<f64> = (0..samples * p * p).map(|_| read(&mut words)).collect();
            let matrices: Vec<f64> = (0..samples * n * n).map(|_| read(&mut words)).collect();
            let initial: Vec<f64> = (0..n * b).map(|_| read(&mut words)).collect();
            if words.next().is_some() {
                return Err("sequence trailing fields".into());
            }
            let (sorted, _, _) = order_fields(&fields, samples, p, truncation)?;
            let order: Vec<usize> = if route == 2 {
                sorted
            } else {
                (0..samples).collect()
            };
            let mut previous = initial.clone();
            let mut results: Vec<String> = (0..samples).map(|_| String::new()).collect();
            for &index in &order {
                let a = &matrices[index * n * n..(index + 1) * n * n];
                let start = if route == 0 { &initial } else { &previous };
                let result = chfsi(
                    a,
                    start,
                    n,
                    b,
                    [count, degree, max_iterations],
                    parameters,
                    tolerance,
                )?;
                results[index] = result.json();
                previous = result.retained;
            }
            println!("{{\"order\":{order:?},\"output\":[{}]}}", results.join(","));
        }
    } else if mode == "filter" {
        let lower: f64 = read(&mut words);
        let upper: f64 = read(&mut words);
        let anchor: f64 = read(&mut words);
        let degrees: Vec<usize> = (0..b).map(|_| read(&mut words)).collect();
        if !anchor.is_finite()
            || !lower.is_finite()
            || !upper.is_finite()
            || !(anchor < lower && lower < upper)
            || degrees.iter().any(|&m| m > 64)
        {
            return Err("filter bounds/degree".into());
        }
        let a: Vec<f64> = (0..n * n).map(|_| read(&mut words)).collect();
        let y: Vec<f64> = (0..n * b).map(|_| read(&mut words)).collect();
        if a.iter().chain(&y).any(|v| !v.is_finite()) || words.next().is_some() {
            return Err("filter values/trailing fields".into());
        }
        let output = filter(&a, &y, n, b, &degrees, [lower, upper, anchor])?;
        println!("{{\"filtered\":{output:?}}}");
    } else if mode == "lobpcg" {
        let count: usize = read(&mut words);
        let warm: usize = read(&mut words);
        let a: Vec<f64> = (0..n * n).map(|_| read(&mut words)).collect();
        let initial: Vec<f64> = if warm == 1 {
            (0..n * b).map(|_| read(&mut words)).collect()
        } else {
            Vec::new()
        };
        if warm > 1 || words.next().is_some() || a.iter().any(|v| !v.is_finite()) {
            return Err("LOBPCG fixture".into());
        }
        let options = eigensolver::IterativeOptions {
            block_size: b,
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-11,
            ..Default::default()
        };
        let result = eigensolver::lobpcg_with_initial(
            n,
            count,
            false,
            &options,
            (warm == 1).then_some(initial.as_slice()),
            &mut execution::WorkBudget::new(&execution::Execution::default())?,
            |x, budget| {
                budget.step_by(n * n)?;
                Ok((0..n)
                    .map(|i| (0..n).map(|j| a[i * n + j] * x[j]).sum())
                    .collect())
            },
        )?;
        println!(
            "{{\"values\":{:?},\"vectors\":{:?},\"residuals\":{:?},\"actions\":{},\"iterations\":{}}}",
            result.values, result.vectors, result.residuals, result.actions, result.iterations
        );
    } else {
        return Err("unknown control mode".into());
    }
    Ok(())
}

#[cfg(test)]
mod control_tests {
    use super::*;

    #[test]
    fn singular_rank_uses_unsquared_threshold_and_qr_is_orthonormal() {
        for small in [1e-14, 1e-12] {
            let a = [1., 0., 0., small, 0., 0.];
            let actual = minimum_singular_value(&a, 3, 2).unwrap();
            assert!((actual / small - 1.).abs() < 1e-12);
        }
        let a = [1., 1., 1., -1., 0., 1.];
        let q = thin_qr(&a, 3, 2).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                let dot: f64 = (0..3).map(|k| q[k * 2 + i] * q[k * 2 + j]).sum();
                assert!((dot - f64::from(i == j)).abs() < 1e-13);
            }
        }
    }

    #[test]
    fn known_cluster_solves_and_retains_fresh_original_residuals() {
        let n = 8;
        let mut a = vec![0.; n * n];
        for (i, value) in [0., 0., 1., 3., 4., 5., 6., 7.].iter().enumerate() {
            a[i * n + i] = *value;
        }
        let mut initial = vec![0.; n * 4];
        for j in 0..4 {
            initial[j * 4 + j] = (j + 1) as f64;
        }
        let result = chfsi(&a, &initial, n, 4, [3, 8, 200], [2., 7.1, -0.1], 1e-11).unwrap();
        for (value, expected) in result.values.iter().zip([0., 0., 1.]) {
            assert!((value - expected).abs() < 1e-12);
        }
        assert!(result.residuals.iter().all(|v| *v < 1e-12));
        assert_eq!(result.trace.len(), 1);
        assert_eq!(result.actions, 8 * 4 + 2 * 4 + 3);
        for i in 0..n {
            for j in 0..n {
                let projector: f64 = (0..3)
                    .map(|k| result.vectors[i * 3 + k] * result.vectors[j * 3 + k])
                    .sum();
                assert!((projector - f64::from(i == j && i < 3)).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn dft_descriptor_keeps_start_ties_and_large_distances() {
        let fields: Vec<f64> = [0., 10000., 10000., 30000.]
            .into_iter()
            .flat_map(|x| vec![x; 16])
            .collect();
        let (order, real, _) = order_fields(&fields, 4, 4, 2).unwrap();
        assert_eq!(order, [0, 1, 2, 3]);
        assert_eq!(real[4], 160000.);
        assert!(order_fields(&[], 0, 4, 2).unwrap().0.is_empty());
    }

    #[test]
    fn transfer_preserves_large_integer_identity_and_projection_rank_failure() {
        let id = (1_u64 << 53) + 1;
        let old = vec![(vec![id], 1)];
        let new = vec![(vec![id + 1], 1), (vec![id], -1)];
        let (mapped, rank) = transfer(&[1.], &old, &new, &[0., 1.], [1, 1, 20], None).unwrap();
        assert_eq!(mapped, [0., -1.]);
        assert_eq!(rank, 1);
        let p = [1., 0., 0., 0.];
        let result = transfer(
            &[1., 0., 0., 1.],
            &new,
            &new,
            &[1., 0., 0., 1.],
            [2, 2, 20],
            Some(&p),
        );
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("rank insufficient")
        );
    }
}
