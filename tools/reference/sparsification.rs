//! Optional bounded rewrite of reproduced sampling and Section 4.3 moments.
//! This worker is outside the library; no weighted-complex API is introduced.
//! Input: n p degree probes draws, row-major integer B, Rademacher Z (p by
//! probes), p probabilities, draw indexes. Frozen draws avoid conflating RNGs.
//! Output: numerical leverage scores, odd moments, weights and sampled Gram.

use cocycle::complex::FilteredComplex;
use cocycle::execution::Execution;
use cocycle::hodge::{HodgeLaplacian, HodgeSelection, SpectrumOptions};
use std::io::{self, Read};
use std::str::{FromStr, SplitWhitespace};

// Embed B as a two-degree integer chain window to reuse the existing checked
// Hodge eigensolver. These formal coordinates are not geometric vertex labels.
struct Boundary {
    rows: usize,
    columns: Vec<Vec<(usize, i32)>>,
}

impl FilteredComplex for Boundary {
    type CellId = usize;
    fn cells(&self) -> impl Iterator<Item = usize> + '_ {
        0..self.rows + self.columns.len()
    }
    fn dimension(&self, cell: usize) -> usize {
        usize::from(cell >= self.rows)
    }
    fn filtration_value(&self, _cell: usize) -> f64 {
        0.
    }
    fn boundary(&self, cell: usize) -> impl Iterator<Item = (usize, i32)> + '_ {
        let entries: &[(usize, i32)] = if cell < self.rows {
            &[]
        } else {
            &self.columns[cell - self.rows]
        };
        entries.iter().copied()
    }
}

fn number<T: FromStr>(tokens: &mut SplitWhitespace<'_>) -> Result<T, &'static str> {
    tokens
        .next()
        .ok_or("missing number")?
        .parse()
        .map_err(|_| "invalid number")
}

fn multiply(h: &[f64], x: &[f64], n: usize, probes: usize) -> Vec<f64> {
    let mut out = vec![0.; x.len()];
    for i in 0..n {
        for j in 0..n {
            for z in 0..probes {
                out[i * probes + z] += h[i * n + j] * x[j * probes + z];
            }
        }
    }
    out
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut input = String::new();
    io::stdin().take(2_000_001).read_to_string(&mut input)?;
    if input.len() > 2_000_000 {
        return Err("input cap exceeded".into());
    }
    let mut tokens = input.split_whitespace();
    let n: usize = number(&mut tokens)?;
    let p: usize = number(&mut tokens)?;
    let degree: usize = number(&mut tokens)?;
    let probes: usize = number(&mut tokens)?;
    let count: usize = number(&mut tokens)?;
    if n > 128
        || p > 128
        || !(1..=32).contains(&degree)
        || !(1..=512).contains(&probes)
        || !(1..=2048).contains(&count)
    {
        return Err("worker dimension/degree/probe/draw cap exceeded".into());
    }
    let mut b = vec![0.; n * p];
    let mut source = Boundary {
        rows: n,
        columns: vec![Vec::new(); p],
    };
    for i in 0..n {
        for j in 0..p {
            let value: i32 = number(&mut tokens)?;
            if value.abs_diff(0) > 16 {
                return Err("integer coefficient cap exceeded".into());
            }
            b[i * p + j] = f64::from(value);
            if value != 0 {
                source.columns[j].push((i, value));
            }
        }
    }
    let mut z = vec![0.; p * probes];
    for value in &mut z {
        *value = number::<f64>(&mut tokens)?;
        if *value != 1. && *value != -1. {
            return Err("expected Rademacher probes".into());
        }
    }
    let mut probabilities = vec![0.; p];
    for probability in &mut probabilities {
        *probability = number::<f64>(&mut tokens)?;
        if !probability.is_finite() || *probability < 0. {
            return Err("invalid probability".into());
        }
    }
    let draws: Vec<usize> = (0..count)
        .map(|_| number(&mut tokens))
        .collect::<Result<_, _>>()?;
    if tokens.next().is_some() {
        return Err("trailing input".into());
    }
    let hodge = HodgeLaplacian::new(&source, 0, HodgeSelection::Upper)?;
    let spectrum = hodge.spectrum(0..n, &SpectrumOptions::default(), &Execution::default())?;
    let scale = spectrum.values().last().copied().unwrap_or(0.);
    let threshold = 1e-10 * scale.max(1.);
    let mut scores = vec![0.; p];
    let mut rank = 0;
    for (mode, &value) in spectrum.values().iter().enumerate() {
        if value > threshold {
            rank += 1;
            for j in 0..p {
                let projection: f64 = (0..n)
                    .map(|i| spectrum.vectors()[i * n + mode] * b[i * p + j])
                    .sum();
                scores[j] += projection * projection / value;
            }
        }
    }
    let mut h = vec![0.; p * p];
    if rank != 0 {
        for i in 0..p {
            for j in 0..p {
                h[i * p + j] = (0..n)
                    .map(|row| b[row * p + i] * b[row * p + j])
                    .sum::<f64>()
                    / scale;
            }
        }
    }
    let mut moments = vec![0.; (degree + 1) * p];
    let mut previous = z.clone();
    let mut current = multiply(&h, &z, p, probes);
    for m in 1..=degree {
        if m % 2 == 1 {
            for i in 0..p {
                moments[m * p + i] = 2.
                    * (0..probes)
                        .map(|j| z[i * probes + j] * current[i * probes + j])
                        .sum::<f64>()
                    / probes as f64;
            }
        }
        let mut next = multiply(&h, &current, p, probes);
        for (value, old) in next.iter_mut().zip(&previous) {
            *value = 2. * *value - old;
        }
        previous = current;
        current = next;
    }
    let mut weights = vec![0.; p];
    if rank == 0 {
        if probabilities.iter().any(|&x| x != 0.) {
            return Err("zero rank requires zero probabilities".into());
        }
        // Fixed-size zero placeholders for draws are ignored in this zero case.
    } else {
        if (probabilities.iter().sum::<f64>() - 1.).abs() > 1e-12
            || source
                .columns
                .iter()
                .zip(&probabilities)
                .any(|(column, &probability)| !column.is_empty() && probability == 0.)
        {
            return Err("probability sum/support mismatch".into());
        }
        for &draw in &draws {
            if draw >= p || probabilities[draw] == 0. {
                return Err("invalid draw index/support".into());
            }
            weights[draw] += 1. / (count as f64 * probabilities[draw]);
        }
    }
    let mut sampled = vec![0.; n * n];
    for i in 0..n {
        for j in 0..n {
            sampled[i * n + j] = (0..p)
                .map(|k| b[i * p + k] * weights[k] * b[j * p + k])
                .sum();
        }
    }
    if scores
        .iter()
        .chain(&moments)
        .chain(&weights)
        .chain(&sampled)
        .any(|x| !x.is_finite())
    {
        return Err("nonfinite numeric result".into());
    }
    println!(
        "{{\"rank\":{rank},\"scores\":{scores:?},\"scale\":{scale},\"moments\":{moments:?},\"weights\":{weights:?},\"sampled\":{sampled:?}}}"
    );
    Ok(())
}
