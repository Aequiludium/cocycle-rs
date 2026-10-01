// Appended to a generated copy of cocycle.rs by the resource controller.
// The ordinary private options and raw-input boundary stay in that owner.
use std::{
    io::{self, Write},
    sync::Barrier,
    thread,
};

fn read_resource_input(
    path: &str,
) -> std::result::Result<(Vec<[f64; 2]>, usize), Box<dyn std::error::Error>> {
    let bytes = fs::read(path)?;
    if bytes.len() < 24 || &bytes[..8] != b"COCDST1\0" {
        return Err("invalid resource fixture".into());
    }
    let n = usize::try_from(u64::from_le_bytes(bytes[8..16].try_into()?))?;
    let m = usize::try_from(u64::from_le_bytes(bytes[16..24].try_into()?))?;
    if n.checked_add(m)
        .and_then(|x| x.checked_mul(16))
        .and_then(|x| x.checked_add(24))
        != Some(bytes.len())
    {
        return Err("invalid resource fixture length".into());
    }
    let points = bytes[24..]
        .chunks_exact(16)
        .map(|b| {
            [
                f64::from_le_bytes(b[..8].try_into().unwrap()),
                f64::from_le_bytes(b[8..].try_into().unwrap()),
            ]
        })
        .collect();
    Ok((points, n))
}

fn resource_run() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args().collect();
    if args.len() == 4 && args[1] == "--reference" {
        let (points, n) = read_resource_input(&args[2])?;
        let (value, _, _) = compute(&points[..n], &points[n..], &args[3], "public")?;
        println!("{value:.17e}");
        return Ok(());
    }
    if args.len() != 9 {
        return Err("fixture metric variant threads jobs expected tolerance warm|cold".into());
    }
    let (points, n) = read_resource_input(&args[1])?;
    let count: usize = args[4].parse()?;
    let jobs: usize = args[5].parse()?;
    let expected: f64 = args[6].parse()?;
    let tolerance: f64 = args[7].parse()?;
    let warm = match args[8].as_str() {
        "warm" => true,
        "cold" => false,
        _ => return Err("invalid warmup mode".into()),
    };
    if count == 0 || jobs == 0 || !expected.is_finite() || !tolerance.is_finite() || tolerance < 0.0
    {
        return Err("invalid resource parameters".into());
    }
    let ready = Barrier::new(count + 1);
    let go = Barrier::new(count + 1);
    let results = thread::scope(
        |scope| -> std::result::Result<_, Box<dyn std::error::Error>> {
            let mut handles = Vec::new();
            for _ in 0..count {
                handles.push(scope.spawn(|| -> Result<_> {
                    let mut durations = Vec::with_capacity(jobs);
                    if warm {
                        let (value, _, _) =
                            compute(&points[..n], &points[n..], &args[2], &args[3])?;
                        assert!(
                            (value - expected).abs() <= tolerance,
                            "warmup scalar mismatch"
                        );
                    }
                    ready.wait();
                    go.wait();
                    let mut last = (
                        bottleneck::Diagnostics::default(),
                        wasserstein::Stats::default(),
                    );
                    let mut checksum = 0.0;
                    for _ in 0..jobs {
                        let start = Instant::now();
                        let (value, bs, ws) =
                            compute(&points[..n], &points[n..], &args[2], &args[3])?;
                        let elapsed = start.elapsed().as_nanos();
                        assert!(
                            value.is_finite() && (value - expected).abs() <= tolerance,
                            "scalar mismatch"
                        );
                        durations.push(elapsed);
                        checksum += std::hint::black_box(value);
                        last = (bs, ws);
                    }
                    Ok((durations, checksum, last))
                }));
            }
            ready.wait();
            println!("{{\"type\":\"ready\",\"pid\":{}}}", std::process::id());
            io::stdout().flush()?;
            let mut command = String::new();
            io::stdin().read_line(&mut command)?;
            if command.trim() != "go" {
                return Err("missing go handshake".into());
            }
            let start = Instant::now();
            go.wait();
            let mut results = Vec::new();
            for handle in handles {
                results.push(handle.join().map_err(|_| "resource thread panicked")??);
            }
            Ok((results, start.elapsed().as_nanos()))
        },
    )?;
    let (results, elapsed) = results;
    let durations: Vec<_> = results.iter().map(|r| &r.0).collect();
    let checksum: f64 = results.iter().map(|r| r.1).sum();
    let (bs, ws) = &results[0].2;
    println!(
        concat!(
            "{{\"type\":\"resource-result\",\"protocol\":\"cocycle-distance-resources-v2\",",
            "\"metric\":\"{}\",\"variant\":\"{}\",\"threads\":{},\"jobs_per_thread\":{},",
            "\"native_group_ns\":{},\"durations_ns\":{:?},\"checksum\":{},",
            "\"peak_rss_kib\":{},\"stats\":{{\"bottleneck_route\":\"{:?}\",",
            "\"dense_solves\":{},\"sparse_solves\":{},\"candidate_pairs\":{},\"positive_edges\":{},",
            "\"components\":{},\"duplicate_groups\":{},\"graph_bytes\":{},\"residual_bytes\":{},\"scratch_bytes\":{},",
            "\"direct_cost_fallbacks\":{},\"greedy_certificates\":{},\"tiny_components\":{},\"augmentations\":{},\"preparation_bytes\":{}}}}}"
        ),
        args[2],
        args[3],
        count,
        jobs,
        elapsed,
        durations,
        checksum,
        nullable(memory("VmHWM:")),
        bs.route,
        ws.dense_solves,
        ws.sparse_solves,
        ws.candidate_pairs,
        ws.positive_edges,
        ws.components,
        ws.duplicate_groups,
        ws.peak_graph_storage_bytes,
        ws.peak_residual_storage_bytes,
        ws.peak_sparse_scratch_bytes,
        ws.direct_cost_fallbacks,
        ws.greedy_certificates,
        ws.tiny_components,
        ws.augmentations,
        ws.preparation_bytes
    );
    io::stdout().flush()?;
    let mut ack = String::new();
    io::stdin().read_line(&mut ack)?;
    if ack.trim() != "ack" {
        return Err("missing final observer handshake".into());
    }
    Ok(())
}

fn main() {
    if let Err(error) = resource_run() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}
