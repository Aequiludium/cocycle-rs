//! Native critical-set and persistence context worker; generated oracle modules
//! are copied from pinned repository revisions by the benchmark controller.
use cocycle::complex::{Simplex, SimplexId, SimplicialComplex, WeightedEdge, WeightedGraph};
use cocycle::diagram::{IntervalEnd, PersistenceDiagram};
use cocycle::filtration::FlagFiltration;
use cocycle::optimization::CriticalSetWorkspace;
use cocycle::persistence::{
    PersistenceBuilder, PersistenceExt, RepresentativeRequest, RepresentativeSelection,
};
use cocycle::{Error, Result};
use std::collections::BTreeSet;
use std::hint::black_box;
use std::time::Instant;

mod dense;
mod reference;

fn memory(field: &str) -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .unwrap()
        .lines()
        .find_map(|line| {
            line.strip_prefix(field)
                .and_then(|value| value.split_whitespace().next()?.parse().ok())
        })
        .unwrap()
}

fn intervals(diagram: &PersistenceDiagram) -> Vec<(usize, f64, Option<f64>)> {
    diagram
        .intervals()
        .map(|bar| {
            (
                bar.dimension(),
                bar.birth(),
                match bar.end() {
                    IntervalEnd::Finite(d) => Some(d),
                    IntervalEnd::Essential => None,
                    IntervalEnd::RightCensored { .. } => {
                        panic!("complete supplied source required")
                    }
                },
            )
        })
        .collect()
}

fn normalize(
    source: &SimplicialComplex,
    pairs: &[(usize, usize)],
    essential: &[usize],
    q: usize,
) -> Vec<(usize, f64, Option<f64>)> {
    let cells = source.simplices();
    let mut result: Vec<_> = pairs
        .iter()
        .filter_map(|&(b, d)| {
            (cells[b].dimension() <= q && cells[b].value() < cells[d].value()).then_some((
                cells[b].dimension(),
                cells[b].value(),
                Some(cells[d].value()),
            ))
        })
        .collect();
    result.extend(
        essential
            .iter()
            .filter(|&&b| cells[b].dimension() <= q)
            .map(|&b| (cells[b].dimension(), cells[b].value(), None)),
    );
    result
}

fn print_result(
    times: &[(&str, f64)],
    rss: u64,
    hwm: u64,
    payload: &[(usize, f64, Option<f64>)],
    targets: &[(usize, f64)],
    pairs: &[(usize, usize)],
    essential: &[usize],
    extra: usize,
) {
    print!("{{\"protocol\":\"critical-sets-native-v1\",\"times_ms\":{{");
    for (i, (key, value)) in times.iter().enumerate() {
        if i > 0 {
            print!(",");
        }
        print!("\"{key}\":{value}");
    }
    print!(
        "}},\"rss_input_kib\":{rss},\"hwm_input_kib\":{hwm},\"peak_rss_kib\":{},\"extra_payload_count\":{extra},\"intervals\":[",
        memory("VmHWM:")
    );
    for (i, (p, b, d)) in payload.iter().enumerate() {
        if i > 0 {
            print!(",");
        }
        let end = d.map_or_else(|| "null".into(), |d| d.to_string());
        print!("[{p},{b},{end}]");
    }
    print!("],\"targets\":[");
    for (i, (id, t)) in targets.iter().enumerate() {
        if i > 0 {
            print!(",");
        }
        print!("[{id},{t}]");
    }
    print!("],\"pairs\":[");
    for (i, (b, d)) in pairs.iter().enumerate() {
        if i > 0 {
            print!(",");
        }
        print!("[{b},{d}]");
    }
    println!("],\"essential\":{essential:?}}}");
}

fn ms(start: Instant) -> f64 {
    start.elapsed().as_secs_f64() * 1000.
}

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let mode = &args[1];
    let bytes = std::fs::read_to_string(&args[2])?;
    let mut words = bytes.split_whitespace();
    let mut integer = || words.next().unwrap().parse::<usize>().unwrap();
    let n = integer();
    let q = integer();
    let count = integer();
    let repetitions = integer();
    let scale: f64 = words.next().unwrap().parse()?;
    let mut cells = Vec::new();
    for _ in 0..n {
        let k: usize = words.next().unwrap().parse()?;
        let v: Vec<usize> = (0..k)
            .map(|_| words.next().unwrap().parse().unwrap())
            .collect();
        let value: f64 = words.next().unwrap().parse()?;
        cells.push(Simplex::new(v, value)?);
    }
    let source = SimplicialComplex::new(cells)?;
    let ids: Vec<SimplexId> = source
        .simplices()
        .iter()
        .map(|s| source.find(s.vertices()).unwrap())
        .collect();
    let mut proposals = Vec::new();
    for _ in 0..count {
        let i: usize = words.next().unwrap().parse()?;
        let t: f64 = words.next().unwrap().parse()?;
        proposals.push((ids[i], t));
    }
    let plain: Vec<_> = proposals.iter().map(|&(i, t)| (i.index(), t)).collect();
    match mode.as_str() {
        "critical" => {
            let rss = memory("VmRSS:");
            let hwm = memory("VmHWM:");
            let begin = Instant::now();
            let mut state = CriticalSetWorkspace::new(black_box(&source))?;
            let build = ms(begin);
            let query = Instant::now();
            let output = black_box(state.targets(black_box(&proposals))?);
            let cold = ms(query);
            let pairs: Vec<_> = state
                .finite_pairs()
                .map(|(b, d)| (b.index(), d.index()))
                .collect();
            let essential: Vec<_> = state.essential_births().map(|b| b.index()).collect();
            let payload = normalize(&source, &pairs, &essential, q);
            let cleanup = Instant::now();
            drop(state);
            let cleanup = ms(cleanup);
            let mut state = CriticalSetWorkspace::new(&source)?;
            assert_eq!(state.targets(&proposals)?, output);
            let warm = Instant::now();
            for _ in 0..repetitions {
                black_box(state.targets(black_box(&proposals))?);
            }
            let warm = ms(warm) / repetitions as f64;
            let targets: Vec<_> = output.iter().map(|&(id, t)| (id.index(), t)).collect();
            print_result(
                &[
                    ("build", build),
                    ("cold_query", cold),
                    ("cleanup", cleanup),
                    ("cold_workflow", build + cold + cleanup),
                    ("warm_batch", warm),
                ],
                rss,
                hwm,
                &payload,
                &targets,
                &pairs,
                &essential,
                0,
            );
        }
        "dense" => {
            let input = dense::Input::new(&source);
            let rss = memory("VmRSS:");
            let hwm = memory("VmHWM:");
            let begin = Instant::now();
            let state = dense::State::new(black_box(&input));
            let build = ms(begin);
            let query = Instant::now();
            let output = black_box(state.targets(black_box(&plain)));
            let cold = ms(query);
            let pairs = state.pairs();
            let essential = state.essential();
            let payload = normalize(&source, &pairs, &essential, q);
            let cleanup = Instant::now();
            drop(state);
            let cleanup = ms(cleanup);
            let state = dense::State::new(&input);
            assert_eq!(state.targets(&plain), output);
            let warm = Instant::now();
            for _ in 0..repetitions {
                black_box(state.targets(black_box(&plain)));
            }
            let warm = ms(warm) / repetitions as f64;
            print_result(
                &[
                    ("build", build),
                    ("cold_query", cold),
                    ("cleanup", cleanup),
                    ("cold_workflow", build + cold + cleanup),
                    ("warm_batch", warm),
                ],
                rss,
                hwm,
                &payload,
                &output,
                &pairs,
                &essential,
                0,
            );
        }
        "reference" => {
            let rss = memory("VmRSS:");
            let hwm = memory("VmHWM:");
            let begin = Instant::now();
            let (pairs, essential) = reference::run(black_box(&source), &ids)?;
            let payload = normalize(&source, &pairs, &essential, q);
            let elapsed = ms(begin);
            print_result(
                &[("diagram_workflow", elapsed)],
                rss,
                hwm,
                &payload,
                &[],
                &pairs,
                &essential,
                0,
            );
        }
        "diagram" | "filtered" | "representatives" | "flag" => {
            let requests: Vec<_> = (0..=q)
                .map(|p| RepresentativeRequest::new(p, scale, RepresentativeSelection::Cycles))
                .collect::<Result<_>>()?;
            let edges: Vec<_> = source
                .simplices()
                .iter()
                .filter(|s| s.dimension() == 1)
                .map(|s| WeightedEdge {
                    vertices: [s.vertices()[0], s.vertices()[1]],
                    value: s.value(),
                })
                .collect();
            let flag = FlagFiltration::new(WeightedGraph::new(source.vertex_count(), edges)?);
            let rss = memory("VmRSS:");
            let hwm = memory("VmHWM:");
            let begin = Instant::now();
            let result = match mode.as_str() {
                "filtered" => PersistenceBuilder::from_complex(black_box(&source))
                    .max_homology_dimension(q)
                    .compute()?,
                "flag" => flag.persistence().max_homology_dimension(q).compute()?,
                "representatives" => source
                    .persistence()
                    .max_homology_dimension(q)
                    .representatives(&requests)
                    .compute()?,
                _ => source.persistence().max_homology_dimension(q).compute()?,
            };
            let payload = intervals(result.diagram());
            let extra = result.representatives().map_or(0, |r| r.len());
            let compute = ms(begin);
            // Independent face enumeration validates requested cycle payloads
            // outside timing. It does not reuse the production boundary action.
            if let Some(representatives) = result.representatives() {
                let expected = result
                    .diagram()
                    .intervals()
                    .filter(|bar| {
                        bar.birth() <= scale
                            && match bar.end() {
                                IntervalEnd::Finite(d) => scale < d,
                                IntervalEnd::Essential => true,
                                _ => false,
                            }
                    })
                    .count();
                assert_eq!(extra, expected);
                for representative in representatives {
                    let mut boundary = BTreeSet::new();
                    for term in representative.terms() {
                        assert_eq!(term.coefficient(), 1);
                        let id = source.find(term.vertices()).unwrap();
                        assert!(source.simplex(id).unwrap().value() <= scale);
                        assert_eq!(term.vertices().len(), representative.dimension() + 1);
                        if term.vertices().len() > 1 {
                            for omitted in 0..term.vertices().len() {
                                let mut face = term.vertices().to_vec();
                                face.remove(omitted);
                                if !boundary.insert(face.clone()) {
                                    boundary.remove(&face);
                                }
                            }
                        }
                    }
                    assert!(boundary.is_empty());
                }
            }
            let cleanup = Instant::now();
            drop(result);
            let elapsed = compute + ms(cleanup);
            print_result(
                &[("diagram_workflow", elapsed)],
                rss,
                hwm,
                &payload,
                &[],
                &[],
                &[],
                extra,
            );
        }
        _ => panic!("unknown benchmark path"),
    }
    Ok(())
}
