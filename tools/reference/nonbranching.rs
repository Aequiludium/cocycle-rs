//! Standalone identity-metric #147 experiment; no production backend dispatch.
//!
//! Input: n m p, then row-major B1 (n*p) and D (m*p), signed unit entries.
//! Args: baseline|compact, repetitions. Every dense buffer is capped at 65,536
//! entries. Baseline retains corrected paper R/E/V/Z; compact omits these dense
//! diagnostics and uses the same traversal and induced metric. Both return the
//! same factor action and dense upper matrix. The Python controller owns the
//! literal-paper replay and independent rational oracle.
#![forbid(unsafe_code)]

use std::hint::black_box;
use std::io::{self, Read};
use std::time::Instant;

const CAP: usize = 65_536;

fn size(rows: usize, columns: usize) -> Result<usize, String> {
    rows.checked_mul(columns)
        .filter(|&n| n <= CAP)
        .ok_or_else(|| "dense entry limit exceeded".into())
}

#[derive(Debug)]
struct Component {
    support: Vec<usize>,
    signs: Vec<i64>,
    kind: &'static str,
}

struct Reduction {
    components: Vec<Component>,
    original_zero: Vec<usize>,
    directions: Vec<Vec<(usize, i64)>>,
    r: Vec<i64>,
    e: Vec<i64>,
    v: Vec<i64>,
    z: Vec<i64>,
}

fn reduce(d: &[i64], m: usize, p: usize, compact: bool) -> Result<Reduction, String> {
    let mut adjacency = vec![Vec::new(); p];
    let mut active = vec![false; p];
    let mut forced = vec![false; p];
    for row in 0..m {
        let indices: Vec<_> = (0..p).filter(|&j| d[row * p + j] != 0).collect();
        for &j in &indices {
            active[j] = true;
        }
        match indices.as_slice() {
            [] => (),
            &[u] => forced[u] = true,
            &[u, v] => {
                let sign = -d[row * p + u] * d[row * p + v];
                adjacency[u].push((v, sign));
                adjacency[v].push((u, sign));
            }
            _ => return Err("branching outside row".into()),
        }
    }
    let original_zero: Vec<_> = (0..p).filter(|&j| !active[j]).collect();
    let nonzero: Vec<_> = (0..p).filter(|&j| active[j]).collect();
    let c = nonzero.len();
    let mut components = Vec::new();
    let mut flags = vec![0i64; p];
    let mut directions: Vec<_> = original_zero.iter().map(|&j| vec![(j, 1)]).collect();
    for root in 0..p {
        if !active[root] || flags[root] != 0 {
            continue;
        }
        flags[root] = 1;
        let mut pending = vec![root];
        let mut support = Vec::new();
        let mut consistent = true;
        let mut singular = false;
        while let Some(u) = pending.pop() {
            support.push(u);
            singular |= forced[u];
            for &(other, relation) in &adjacency[u] {
                let wanted = flags[u] * relation;
                if flags[other] == 0 {
                    flags[other] = wanted;
                    pending.push(other);
                } else {
                    consistent &= flags[other] == wanted;
                }
            }
        }
        support.sort_unstable();
        let signs: Vec<_> = support.iter().map(|&j| flags[j]).collect();
        let kind = if singular {
            "row-singular"
        } else if consistent {
            "regulable"
        } else {
            "irregular"
        };
        if kind == "regulable" {
            directions.push(support.iter().copied().zip(signs.iter().copied()).collect());
        }
        components.push(Component {
            support,
            signs,
            kind,
        });
    }
    let (mut r, mut e, mut v, mut z) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    if !compact {
        r = vec![0; size(m, c)?];
        e = vec![0; size(c, c)?];
        v = vec![0; size(c, c)?];
        z = vec![0; size(p, directions.len())?];
        let mut local = vec![0; p];
        for (i, &j) in nonzero.iter().enumerate() {
            local[j] = i;
            e[i * c + i] = 1;
            v[i * c + i] = 1;
        }
        for component in &components {
            if component.kind == "regulable" {
                let last = local[*component.support.last().unwrap()];
                for (&j, &flag) in component.support.iter().zip(&component.signs) {
                    let i = local[j];
                    e[i * c + i] = flag;
                    v[i * c + last] = 1;
                }
            }
        }
        // Explicit corrected paper product, rather than setting a column to
        // zero and leaving the other reoriented columns unchanged.
        for row in 0..m {
            for col in 0..c {
                r[row * c + col] = (0..c)
                    .map(|i| d[row * p + nonzero[i]] * e[i * c + i] * v[i * c + col])
                    .sum();
            }
        }
        for (col, direction) in directions.iter().enumerate() {
            for &(row, sign) in direction {
                z[row * directions.len() + col] = sign;
            }
        }
    }
    Ok(Reduction {
        components,
        original_zero,
        directions,
        r,
        e,
        v,
        z,
    })
}

fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if !(args.len() == 2 || (args.len() == 3 && args[2] == "--summary"))
        || !["baseline", "compact"].contains(&args[0].as_str())
    {
        return Err("usage: nonbranching baseline|compact repetitions [--summary]".into());
    }
    let compact = args[0] == "compact";
    let repetitions: u32 = args[1].parse().map_err(|_| "invalid repetitions")?;
    if repetitions == 0 || repetitions > 1_000_000 {
        return Err("repetitions must be in 1..=1000000".into());
    }
    let prepare = Instant::now();
    let mut input = String::new();
    io::stdin()
        .read_to_string(&mut input)
        .map_err(|e| e.to_string())?;
    let mut tokens = input.split_whitespace();
    let mut dimension = || -> Result<usize, String> {
        tokens
            .next()
            .ok_or("missing dimension")?
            .parse()
            .map_err(|_| "invalid dimension".into())
    };
    let (n, m, p) = (dimension()?, dimension()?, dimension()?);
    if n > CAP || m > CAP || p > CAP {
        return Err("dimension limit exceeded".into());
    }
    let bp = size(n, p)?;
    let dp = size(m, p)?;
    let nn = size(n, n)?;
    let mut read_matrix = |length| -> Result<Vec<i64>, String> {
        (0..length)
            .map(|_| {
                let value: i64 = tokens
                    .next()
                    .ok_or("missing matrix entry")?
                    .parse()
                    .map_err(|_| "invalid matrix entry")?;
                if !(-1..=1).contains(&value) {
                    return Err("only signed unit incidence is supported".into());
                }
                Ok(value)
            })
            .collect()
    };
    let b = read_matrix(bp)?;
    let d = read_matrix(dp)?;
    if tokens.next().is_some() {
        return Err("trailing input".into());
    }
    let prepare_ns = prepare.elapsed().as_nanos();
    let mut elapsed = [0u128; 4];
    let total = Instant::now();
    let mut last = None;
    for _ in 0..repetitions {
        let start = Instant::now();
        let reduction = reduce(black_box(&d), m, p, compact)?;
        elapsed[0] += start.elapsed().as_nanos();
        let k = reduction.directions.len();
        let start = Instant::now();
        let mut factor = vec![0.; size(n, k)?];
        let metric: Vec<_> = reduction.directions.iter().map(Vec::len).collect();
        for row in 0..n {
            for col in 0..k {
                let numerator: i64 = if compact {
                    reduction.directions[col]
                        .iter()
                        .map(|&(j, s)| b[row * p + j] * s)
                        .sum()
                } else {
                    (0..p)
                        .map(|j| b[row * p + j] * reduction.z[j * k + col])
                        .sum()
                };
                factor[row * k + col] = numerator as f64 / (metric[col] as f64).sqrt();
            }
        }
        elapsed[1] += start.elapsed().as_nanos();
        let start = Instant::now();
        let mut action = vec![0.; n];
        for col in 0..k {
            let dot: f64 = (0..n).map(|row| factor[row * k + col]).sum();
            for row in 0..n {
                action[row] += factor[row * k + col] * dot;
            }
        }
        black_box(&action);
        elapsed[2] += start.elapsed().as_nanos();
        let start = Instant::now();
        let mut upper = vec![0.; nn];
        for i in 0..n {
            for j in 0..n {
                upper[i * n + j] = (0..k)
                    .map(|col| factor[i * k + col] * factor[j * k + col])
                    .sum();
            }
        }
        black_box(&upper);
        elapsed[3] += start.elapsed().as_nanos();
        last = Some((reduction, metric, factor, action, upper));
    }
    let total_ns = total.elapsed().as_nanos();
    let (reduction, metric, factor, action, upper) = last.unwrap();
    let components: Vec<_> = reduction
        .components
        .iter()
        .map(|c| {
            format!(
                "{{\"support\":{:?},\"signs\":{:?},\"kind\":\"{}\"}}",
                c.support, c.signs, c.kind
            )
        })
        .collect();
    let diagnostic_entries =
        reduction.r.len() + reduction.e.len() + reduction.v.len() + reduction.z.len();
    let incidence_nnz = b.iter().chain(&d).filter(|&&x| x != 0).count();
    let kernel_support_entries: usize = reduction.directions.iter().map(Vec::len).sum();
    let factor_nnz = factor.iter().filter(|&&x| x != 0.).count();
    let summary = args.len() == 3;
    let (r, e, v, z) = if summary {
        (&[][..], &[][..], &[][..], &[][..])
    } else {
        (
            reduction.r.as_slice(),
            reduction.e.as_slice(),
            reduction.v.as_slice(),
            reduction.z.as_slice(),
        )
    };
    let export = Instant::now();
    let encoded = format!(
        "{{\"R\":{:?},\"E\":{:?},\"V\":{:?},\"Z\":{:?},\"original_zero\":{:?},\"components\":[{}],\"metric\":{:?},\"factor\":{:?},\"action\":{:?},\"upper\":{:?},\"phase_ns\":{:?},\"total_ns\":{},\"diagnostic_entries\":{},\"input_entries\":{},\"factor_entries\":{},\"returned_entries\":{},\"prepare_ns\":{},\"incidence_nnz\":{},\"kernel_support_entries\":{},\"factor_nnz\":{}",
        r,
        e,
        v,
        z,
        reduction.original_zero,
        components.join(","),
        metric,
        factor,
        action,
        upper,
        elapsed,
        total_ns,
        diagnostic_entries,
        b.len() + d.len(),
        factor.len(),
        upper.len(),
        prepare_ns,
        incidence_nnz,
        kernel_support_entries,
        factor_nnz
    );
    let serialization_ns = export.elapsed().as_nanos();
    println!("{encoded},\"serialization_ns\":{serialization_ns}}}");
    Ok(())
}
