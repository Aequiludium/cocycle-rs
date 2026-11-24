//! Optional bounded rewrite of reproduced sampling and Section 4.3 moments.
//! This worker is outside the library; no weighted-complex API is introduced.
//! Input: n p degree probes draws, row-major integer B, Rademacher Z (p by
//! probes), p probabilities, draw indexes. Frozen draws avoid conflating RNGs.
//! Output: numerical leverage scores, odd moments, weights and sampled Gram.
//! Native input: `native n p degree probes draws seed leverage|uniform`, B, Z.
//! Native draws use SplitMix64 high-53-bit uniforms and an ordered inverse CDF.
//! `--batch` accepts one bounded request per line, without caching solver state.

use cocycle::complex::FilteredComplex;
use cocycle::execution::Execution;
use cocycle::hodge::{HodgeLaplacian, HodgeSelection, SpectrumOptions};
use std::io::{self, BufRead, Read};
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

fn random_word(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e3779b97f4a7c15);
    let mut word = *state;
    word = (word ^ (word >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    word = (word ^ (word >> 27)).wrapping_mul(0x94d049bb133111eb);
    word ^ (word >> 31)
}

fn native_draws(probabilities: &[f64], count: usize, seed: u64) -> Vec<usize> {
    let mut state = seed;
    let mass: f64 = probabilities.iter().sum();
    // evaluate validates finite positive mass and support before calling here.
    let last_positive = probabilities.iter().rposition(|&x| x > 0.).unwrap();
    (0..count)
        .map(|_| {
            let uniform = (random_word(&mut state) >> 11) as f64 / 9_007_199_254_740_992.;
            let target = uniform * mass;
            let mut cumulative = 0.;
            for (index, &probability) in probabilities.iter().enumerate() {
                cumulative += probability;
                if target < cumulative {
                    return index;
                }
            }
            // The last bin absorbs at most floating-point CDF endpoint roundoff.
            last_positive
        })
        .collect()
}

fn evaluate(input: &str) -> Result<String, Box<dyn std::error::Error>> {
    evaluate_route(input, None)
}

fn evaluate_route(
    input: &str,
    benchmark: Option<&str>,
) -> Result<String, Box<dyn std::error::Error>> {
    let mut tokens = input.split_whitespace();
    let first = tokens.next().ok_or("missing request")?;
    let native = first == "native";
    let n: usize = if native {
        number(&mut tokens)?
    } else {
        first.parse()?
    };
    let p: usize = number(&mut tokens)?;
    let degree: usize = number(&mut tokens)?;
    let probes: usize = number(&mut tokens)?;
    let count: usize = number(&mut tokens)?;
    let seed = if native {
        Some(number::<u64>(&mut tokens)?)
    } else {
        None
    };
    let route = if native {
        tokens.next().ok_or("missing sampling route")?
    } else {
        "replay"
    };
    if !matches!(route, "leverage" | "uniform" | "replay") || (native && route == "replay") {
        return Err("unknown native sampling route".into());
    }
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
    if !native {
        for probability in &mut probabilities {
            *probability = number::<f64>(&mut tokens)?;
            if !probability.is_finite() || *probability < 0. {
                return Err("invalid probability".into());
            }
        }
    }
    let mut draws: Vec<usize> = (0..if native { 0 } else { count })
        .map(|_| number(&mut tokens))
        .collect::<Result<_, _>>()?;
    if tokens.next().is_some() {
        return Err("trailing input".into());
    }
    if benchmark == Some("exact") {
        let mut values = vec![0.; n * n];
        for i in 0..n {
            for j in 0..n {
                values[i * n + j] = (0..p).map(|k| b[i * p + k] * b[j * p + k]).sum();
            }
        }
        return Ok(format!("{{\"values\":{values:?}}}"));
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
    let score_sum: f64 = scores.iter().sum();
    if (score_sum - rank as f64).abs() > 1e-8 * (rank as f64).max(1.) {
        return Err("leverage trace/rank inconsistency".into());
    }
    if rank == 0 && b.iter().any(|&x| x != 0.) {
        return Err("nonzero boundary classified as zero rank".into());
    }
    if native && rank != 0 {
        for (index, probability) in probabilities.iter_mut().enumerate() {
            *probability = if route == "uniform" {
                1. / p as f64
            } else {
                scores[index] / score_sum
            };
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
        draws.clear();
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
        if let Some(seed) = seed {
            // Positive rank and checked support establish a positive probability mass.
            draws = native_draws(&probabilities, count, seed);
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
    let seed_json = seed.map_or_else(|| "null".to_owned(), |value| value.to_string());
    if benchmark == Some("sampled") {
        return Ok(format!("{{\"values\":{sampled:?}}}"));
    }
    Ok(format!(
        "{{\"rank\":{rank},\"scores\":{scores:?},\"score_sum\":{score_sum},\"scale\":{scale},\"moments\":{moments:?},\"probabilities\":{probabilities:?},\"draws\":{draws:?},\"weights\":{weights:?},\"sampled\":{sampled:?},\"seed\":{seed_json},\"route\":\"{route}\"}}"
    ))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if let [flag, route] = args.as_slice()
        && flag == "--benchmark"
        && (route == "exact" || route == "sampled")
    {
        let mut input = String::new();
        io::stdin().take(2_000_001).read_to_string(&mut input)?;
        if input.len() > 2_000_000 {
            return Err("input cap exceeded".into());
        }
        println!("{}", evaluate_route(&input, Some(route))?);
    } else if args == ["--batch"] {
        let mut reader = io::stdin().lock();
        let mut requests = 0;
        loop {
            let mut line = String::new();
            let length = reader.by_ref().take(2_000_001).read_line(&mut line)?;
            if length == 0 {
                break;
            }
            if length > 2_000_000 {
                return Err("input cap exceeded".into());
            }
            requests += 1;
            let result = evaluate(&line).map_err(|error| format!("request {requests}: {error}"))?;
            println!("{result}");
        }
        if requests == 0 {
            return Err("missing request".into());
        }
    } else if args.is_empty() {
        let mut input = String::new();
        io::stdin().take(2_000_001).read_to_string(&mut input)?;
        if input.len() > 2_000_000 {
            return Err("input cap exceeded".into());
        }
        println!("{}", evaluate(&input)?);
    } else {
        return Err("expected no arguments or --batch".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn values(output: &str, field: &str) -> Vec<f64> {
        output
            .split(&format!("\"{field}\":["))
            .nth(1)
            .unwrap()
            .split(']')
            .next()
            .unwrap()
            .split(',')
            .filter(|x| !x.is_empty())
            .map(|x| x.trim().parse().unwrap())
            .collect()
    }

    #[test]
    fn splitmix64_golden_words() {
        let mut state = 0;
        assert_eq!(random_word(&mut state), 0xe220a8397b1dcdaf);
        assert_eq!(random_word(&mut state), 0x6e789e6aa1b965f4);
    }

    #[test]
    fn inverse_cdf_never_selects_zero_mass_columns() {
        assert_eq!(native_draws(&[0., 1., 0.], 20, 0), vec![1; 20]);
        assert_eq!(native_draws(&[0.5, 0.5], 2, 0), vec![1, 0]);
    }

    #[test]
    fn native_nonuniform_one_row_has_constant_unbiased_gram() {
        // B=[1,2], G=5, tau=[1/5,4/5]. Every reweighted draw contributes 5/k.
        let out = evaluate("native 1 2 3 1 30 0 leverage 1 2 1 1").unwrap();
        let scores = values(&out, "scores");
        assert!((scores[0] - 0.2).abs() < 1e-14);
        assert!((scores[1] - 0.8).abs() < 1e-14);
        assert!((values(&out, "sampled")[0] - 5.).abs() < 1e-13);
        assert_eq!(values(&out, "draws").len(), 30);
        assert_eq!(
            out,
            evaluate("native 1 2 3 1 30 0 leverage 1 2 1 1").unwrap()
        );
    }

    #[test]
    fn zero_rank_has_no_draws_or_probability_normalization() {
        for request in [
            "native 3 0 3 1 4 0 leverage",
            "native 0 0 3 1 4 0 uniform",
            "native 1 2 3 1 4 0 leverage 0 0 1 -1",
        ] {
            let out = evaluate(request).unwrap();
            assert!(out.contains("\"rank\":0"));
            assert!(values(&out, "draws").is_empty());
            assert!(values(&out, "probabilities").iter().all(|&x| x == 0.));
        }
    }

    #[test]
    fn scalar_odd_moments_follow_direct_polynomials() {
        let out = evaluate("native 1 1 5 1 1 0 leverage 1 1").unwrap();
        assert_eq!(values(&out, "moments"), vec![0., 2., 0., 2., 0., 2.]);
    }

    #[test]
    fn replay_rejects_missing_support_and_invalid_native_headers() {
        for request in [
            "1 2 1 1 1 1 2 1 1 1 0 0",
            "native 1 1 1 1 1 0 unknown 1 1",
            "native 1 1 1 1 0 0 leverage 1 1",
            "native 1 1 1 1 1 0 leverage 1 0",
            "native 1 1 1 1 1 0 leverage 1 1 trailing",
            "native 129 1 1 1 1 0 leverage",
            "native 1 1 1 1 1 0 leverage 17 1",
        ] {
            assert!(evaluate(request).is_err(), "accepted {request}");
        }
    }
}
