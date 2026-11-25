//! Independent dense row algebra and public field/resource contracts.
use cocycle::Error;
use cocycle::algebra::PrimeField;
use cocycle::complex::{DirectedGraph, PathBoundaryMode as Mode, PathComplex};
use cocycle::diagram::{Coverage, IntervalEnd, PersistenceData};
use cocycle::execution::Execution;
use cocycle::filtration::{FiltrationKind, FiltrationScale, FiltrationSource};
use cocycle::persistence::{
    PathHomology, PathVector, compute_path_homology, compute_path_homology_with,
};
use cocycle::persistence::{
    PathPersistence, compute_path_persistence, compute_path_persistence_with,
};
use cocycle::persistence::{
    compute_low_dimensional_path_persistence, compute_low_dimensional_path_persistence_with,
};
use cocycle::persistence::{compute_minimum_path_basis, compute_minimum_path_basis_with};
use std::collections::BTreeMap;
use std::sync::atomic::AtomicBool;

// The oracle uses dense row operations, direct word deletion and u128 products;
// it shares neither production path enumeration nor sparse column elimination.
fn product(a: u32, b: u32, p: u32) -> u32 {
    (u128::from(a) * u128::from(b) % u128::from(p)) as u32
}
fn add(a: u32, b: u32, p: u32) -> u32 {
    ((u64::from(a) + u64::from(b)) % u64::from(p)) as u32
}
fn inverse(a: u32, p: u32) -> u32 {
    let (mut b, mut e, mut r) = (a, p - 2, 1);
    while e > 0 {
        if e & 1 != 0 {
            r = product(r, b, p);
        }
        b = product(b, b, p);
        e >>= 1;
    }
    r
}
fn kernel(columns: &[Vec<u32>], rows: usize, p: u32) -> (usize, Vec<Vec<u32>>) {
    let mut a = vec![vec![0; columns.len()]; rows];
    for (j, column) in columns.iter().enumerate() {
        for (i, &v) in column.iter().enumerate() {
            a[i][j] = v;
        }
    }
    let mut pivots = Vec::new();
    for j in 0..columns.len() {
        let r = pivots.len();
        let Some(k) = (r..rows).find(|&k| a[k][j] != 0) else {
            continue;
        };
        a.swap(r, k);
        let scale = inverse(a[r][j], p);
        for v in &mut a[r] {
            *v = product(*v, scale, p);
        }
        let pivot = a[r].clone();
        for (i, row) in a.iter_mut().enumerate() {
            if i == r {
                continue;
            }
            let factor = if row[j] == 0 { 0 } else { p - row[j] };
            for (x, y) in row.iter_mut().zip(&pivot) {
                *x = add(*x, product(factor, *y, p), p);
            }
        }
        pivots.push(j);
    }
    let mut basis = Vec::new();
    for j in 0..columns.len() {
        if pivots.contains(&j) {
            continue;
        }
        let mut v = vec![0; columns.len()];
        v[j] = 1;
        for (i, &k) in pivots.iter().enumerate() {
            v[k] = if a[i][j] == 0 { 0 } else { p - a[i][j] };
        }
        basis.push(v);
    }
    (pivots.len(), basis)
}
fn apply(columns: &[Vec<u32>], vector: &[u32], rows: usize, p: u32) -> Vec<u32> {
    let mut result = vec![0; rows];
    for (column, &scalar) in columns.iter().zip(vector) {
        for (out, &value) in result.iter_mut().zip(column) {
            *out = add(*out, product(scalar, value, p), p);
        }
    }
    result
}
fn dense(columns: &[PathVector], rows: usize, p: u32) -> Vec<Vec<u32>> {
    columns
        .iter()
        .map(|c| {
            let mut v = vec![0; rows];
            for (&i, &value) in c {
                assert!(i < rows && value > 0 && value < p);
                v[i] = value;
            }
            v
        })
        .collect()
}
fn ambient(path: &[usize], mode: Mode, p: u32) -> BTreeMap<Vec<usize>, u32> {
    let mut result = BTreeMap::new();
    if path.len() == 1 {
        return result;
    }
    for i in 0..path.len() {
        let mut face = path.to_vec();
        face.remove(i);
        if mode == Mode::Regular && face.windows(2).any(|v| v[0] == v[1]) {
            continue;
        }
        let v = result.entry(face).or_insert(0);
        *v = add(*v, if i % 2 == 0 { 1 } else { p - 1 }, p);
    }
    result.retain(|_, v| *v != 0);
    result
}
fn certificate(h: &PathHomology) {
    let p = h.field().characteristic();
    let d = h.max_homology_dimension();
    let mut bases: Vec<Vec<Vec<u32>>> = Vec::new();
    let mut boundaries: Vec<Vec<Vec<u32>>> = Vec::new();
    for q in 0..=d + 1 {
        let paths = h.paths(q).unwrap();
        let lower = if q == 0 {
            &[][..]
        } else {
            h.paths(q - 1).unwrap()
        };
        let words: Vec<_> = paths
            .iter()
            .map(|path| ambient(path, h.mode(), p))
            .collect();
        let mut forbidden = BTreeMap::new();
        for word in &words {
            for face in word.keys() {
                if !lower.contains(face) {
                    let n = forbidden.len();
                    forbidden.entry(face).or_insert(n);
                }
            }
        }
        let mut a = vec![vec![0; lower.len()]; paths.len()];
        let mut f = vec![vec![0; forbidden.len()]; paths.len()];
        for (j, word) in words.iter().enumerate() {
            for (face, &v) in word {
                if let Some(i) = lower.iter().position(|p| p == face) {
                    a[j][i] = v;
                } else {
                    f[j][forbidden[face]] = v;
                }
            }
        }
        let (rank, _) = kernel(&f, forbidden.len(), p);
        let s = dense(h.omega_bases(q).unwrap(), paths.len(), p);
        assert_eq!(s.len(), paths.len() - rank);
        assert_eq!(kernel(&s, paths.len(), p).0, s.len());
        for v in &s {
            assert!(apply(&f, v, forbidden.len(), p).iter().all(|&v| v == 0));
        }
        let b = dense(
            h.boundaries(q).unwrap(),
            if q == 0 {
                0
            } else {
                h.omega_bases(q - 1).unwrap().len()
            },
            p,
        );
        assert_eq!(b.len(), s.len());
        for (v, column) in s.iter().zip(&b) {
            let image = if q == 0 {
                vec![]
            } else {
                apply(&bases[q - 1], column, lower.len(), p)
            };
            assert_eq!(image, apply(&a, v, lower.len(), p));
            if q > 0 {
                assert!(
                    apply(
                        &boundaries[q - 1],
                        column,
                        if q == 1 {
                            0
                        } else {
                            h.omega_bases(q - 2).unwrap().len()
                        },
                        p
                    )
                    .iter()
                    .all(|&x| x == 0)
                );
            }
        }
        assert_eq!(
            kernel(
                &b,
                if q == 0 {
                    0
                } else {
                    h.omega_bases(q - 1).unwrap().len()
                },
                p
            )
            .0,
            h.boundary_ranks()[q]
        );
        bases.push(s);
        boundaries.push(b);
    }
    for q in 0..=d {
        let rows = h.paths(q).unwrap().len();
        let images: Vec<_> = boundaries[q + 1]
            .iter()
            .map(|v| apply(&bases[q], v, rows, p))
            .collect();
        let rank = kernel(&images, rows, p).0;
        let reps = dense(h.representatives(q).unwrap(), rows, p);
        let mut span = images;
        span.extend(reps.clone());
        assert_eq!(kernel(&span, rows, p).0, rank + reps.len());
        assert_eq!(reps.len(), bases[q].len() - h.boundary_ranks()[q] - rank);
        for v in reps {
            let paths = h.paths(q).unwrap();
            let mut image = BTreeMap::new();
            for (path, scalar) in paths.iter().zip(v) {
                for (face, value) in ambient(path, h.mode(), p) {
                    let x = image.entry(face).or_insert(0);
                    *x = add(*x, product(scalar, value, p), p);
                }
            }
            assert!(image.values().all(|&v| v == 0));
        }
    }
}

#[test]
fn all_three_vertex_graphs_have_complete_independent_certificates_in_large_fields() {
    let pairs = [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)];
    for mask in 0..64 {
        let edges: Vec<_> = pairs
            .iter()
            .enumerate()
            .filter_map(|(i, &e)| (mask & (1 << i) != 0).then_some(e))
            .collect();
        let source = DirectedGraph::new(3, edges.clone())
            .unwrap()
            .allowed_paths(3, 1000)
            .unwrap();
        for q in 0..=3u32 {
            let mut expected = Vec::new();
            for mut code in 0..3usize.pow(q + 1) {
                let mut word = vec![0; q as usize + 1];
                for v in word.iter_mut().rev() {
                    *v = code % 3;
                    code /= 3;
                }
                if word.windows(2).all(|v| edges.contains(&(v[0], v[1]))) {
                    expected.push(word);
                }
            }
            assert_eq!(source.paths(q as usize), expected);
        }
        for p in [2, 3, 5, 4_294_967_291] {
            for mode in [Mode::Regular, Mode::NonRegular] {
                certificate(
                    &compute_path_homology(&source, 2, PrimeField::new(p).unwrap(), mode, 100000)
                        .unwrap(),
                );
            }
        }
    }
}

#[test]
fn model_truncation_explicit_sources_and_owned_coordinates_are_distinct() {
    let source = DirectedGraph::new(2, vec![(0, 1), (1, 0)])
        .unwrap()
        .allowed_paths(2, 6)
        .unwrap();
    let field = PrimeField::new(3).unwrap();
    let regular = compute_path_homology(&source, 1, field, Mode::Regular, 1000).unwrap();
    let full = compute_path_homology(&source, 1, field, Mode::NonRegular, 1000).unwrap();
    assert_eq!(regular.betti(), [1, 0]);
    assert_eq!(full.betti(), [1, 1]);
    drop(source);
    certificate(&regular);
    certificate(&full);
    let finite = PathComplex::new(2, &[vec![0], vec![1], vec![0, 1], vec![1, 0]], 4).unwrap();
    assert_eq!(
        compute_path_homology(&finite, 1, field, Mode::Regular, 1000)
            .unwrap()
            .betti(),
        [1, 1]
    );
    let short = DirectedGraph::new(2, vec![(0, 1), (1, 0)])
        .unwrap()
        .allowed_paths(1, 4)
        .unwrap();
    assert!(compute_path_homology(&short, 1, field, Mode::Regular, 1000).is_err());
    let empty = PathComplex::new(0, &[], 1).unwrap();
    assert_eq!(
        compute_path_homology(&empty, 3, field, Mode::Regular, 1000)
            .unwrap()
            .betti(),
        [0; 4]
    );
    assert!(compute_path_homology(&empty, usize::MAX, field, Mode::Regular, 1000).is_err());
    assert!(compute_path_homology(&empty, 0, field, Mode::Regular, 0).is_err());
}

#[test]
fn cooperative_controls_fail_inside_algebra_and_allow_recovery() {
    let source = DirectedGraph::new(3, vec![(0, 1), (1, 2), (2, 0)])
        .unwrap()
        .allowed_paths(3, 100)
        .unwrap();
    let field = PrimeField::new(4_294_967_291).unwrap();
    for work in [0, 10, 50, 200] {
        assert!(matches!(
            compute_path_homology_with(
                &source,
                2,
                field,
                Mode::Regular,
                10000,
                &Execution::default().max_work(work)
            ),
            Err(Error::WorkLimitExceeded { .. })
        ));
    }
    let flag = AtomicBool::new(true);
    assert_eq!(
        compute_path_homology_with(
            &source,
            1,
            field,
            Mode::Regular,
            1000,
            &Execution::default().cancellation(&flag)
        ),
        Err(Error::Cancelled)
    );
    assert!(compute_path_homology(&source, 1, field, Mode::Regular, 2).is_err());
    assert_eq!(
        compute_path_homology(&source, 1, field, Mode::Regular, 10000)
            .unwrap()
            .betti(),
        [1, 1]
    );
}

#[test]
fn higher_join_and_characteristic_sensitive_omega_spaces_are_recomputed() {
    let edges: Vec<_> = (0..6)
        .flat_map(|a| (0..6).filter(move |&b| a / 2 < b / 2).map(move |b| (a, b)))
        .collect();
    let source = DirectedGraph::new(6, edges)
        .unwrap()
        .allowed_paths(3, 1000)
        .unwrap();
    for p in [2, 3, 4_294_967_291] {
        let h = compute_path_homology(
            &source,
            2,
            PrimeField::new(p).unwrap(),
            Mode::Regular,
            100000,
        )
        .unwrap();
        assert_eq!(h.betti(), [1, 0, 1]);
        certificate(&h);
    }
    // Fu-Ivanov 2407.17001v1 section 5.4, encoded as directed edge data.
    let edges = vec![
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 4),
        (2, 5),
        (3, 6),
        (1, 5),
        (2, 6),
        (3, 4),
        (4, 7),
        (5, 8),
        (6, 9),
        (4, 8),
        (5, 9),
        (6, 7),
        (7, 10),
        (8, 10),
        (9, 10),
        (1, 7),
        (2, 8),
        (3, 9),
        (1, 9),
        (2, 7),
        (3, 8),
    ];
    let source = DirectedGraph::new(11, edges)
        .unwrap()
        .allowed_paths(4, 10000)
        .unwrap();
    for (p, expected) in [(2, 1), (3, 0), (4_294_967_291, 0)] {
        let h = compute_path_homology(
            &source,
            3,
            PrimeField::new(p).unwrap(),
            Mode::Regular,
            100000,
        )
        .unwrap();
        assert_eq!(h.omega_bases(4).unwrap().len(), expected);
        certificate(&h);
    }
}

fn persistent_ranks(h: &PathPersistence, edges: &[(usize, usize, f64)], mode: Mode) {
    let p = h.context().characteristic();
    for q in 0..=h.diagram().max_dimension() {
        let paths = h.paths(q).unwrap();
        let upper = h.paths(q + 1).unwrap();
        let birth = |path: &[usize]| {
            path.windows(2)
                .map(|v| edges.iter().find(|e| (e.0, e.1) == (v[0], v[1])).unwrap().2)
                .fold(0.0f64, f64::max)
        };
        let mut stages = Vec::new();
        for (event, &time) in h.times().iter().enumerate() {
            let active: Vec<_> = paths
                .iter()
                .enumerate()
                .filter(|(_, w)| birth(w) <= time)
                .map(|(i, _)| i)
                .collect();
            let mut rows = BTreeMap::new();
            let full: Vec<_> = active
                .iter()
                .map(|&i| ambient(&paths[i], mode, p))
                .collect();
            for c in &full {
                for row in c.keys() {
                    let n = rows.len();
                    rows.entry(row.clone()).or_insert(n);
                }
            }
            let columns: Vec<_> = full
                .iter()
                .map(|c| {
                    let mut v = vec![0; rows.len()];
                    for (r, &x) in c {
                        v[rows[r]] = x;
                    }
                    v
                })
                .collect();
            let cycles: Vec<_> = kernel(&columns, rows.len(), p)
                .1
                .iter()
                .map(|v| {
                    let mut out = vec![0; paths.len()];
                    for (&i, &x) in active.iter().zip(v) {
                        out[i] = x;
                    }
                    out
                })
                .collect();
            let upper_words: Vec<_> = upper
                .iter()
                .filter(|w| birth(w) <= time)
                .map(|w| ambient(w, mode, p))
                .collect();
            let mut forbidden = BTreeMap::new();
            for c in &upper_words {
                for row in c.keys() {
                    if !active.iter().any(|&i| paths[i] == *row) {
                        let n = forbidden.len();
                        forbidden.entry(row.clone()).or_insert(n);
                    }
                }
            }
            let mut a = vec![vec![0; paths.len()]; upper_words.len()];
            let mut f = vec![vec![0; forbidden.len()]; upper_words.len()];
            for (j, c) in upper_words.iter().enumerate() {
                for (row, &x) in c {
                    if let Some(&i) = active.iter().find(|&&i| paths[i] == *row) {
                        a[j][i] = x;
                    } else {
                        f[j][forbidden[row]] = x;
                    }
                }
            }
            let images: Vec<_> = kernel(&f, forbidden.len(), p)
                .1
                .iter()
                .map(|v| apply(&a, v, paths.len(), p))
                .collect();
            let active_s: Vec<_> = h
                .omega_bases(q)
                .unwrap()
                .iter()
                .zip(h.births(q).unwrap())
                .filter(|(_, b)| **b <= event)
                .map(|(v, _)| v.clone())
                .collect();
            let dense_s = dense(&active_s, paths.len(), p);
            assert_eq!(kernel(&dense_s, paths.len(), p).0, dense_s.len());
            let cycle_rank = cycles.len();
            let (_, forbidden_kernel) = kernel(&f, forbidden.len(), p);
            assert_eq!(
                forbidden_kernel.len(),
                h.births(q + 1)
                    .unwrap()
                    .iter()
                    .filter(|&&b| b <= event)
                    .count()
            );
            assert!(cycle_rank <= dense_s.len());
            stages.push((cycles, images));
        }
        for (s, (z, _)) in stages.iter().enumerate() {
            for (t, (_, b)) in stages.iter().enumerate().skip(s) {
                let mut combined = z.clone();
                combined.extend(b.clone());
                let expected = kernel(&combined, paths.len(), p).0 - kernel(b, paths.len(), p).0;
                let bars = h
                    .diagram()
                    .dimension(q)
                    .unwrap()
                    .iter()
                    .filter(|bar| {
                        bar.birth() <= h.times()[s]
                            && match bar.end() {
                                IntervalEnd::Essential => true,
                                IntervalEnd::Finite(death) => death > h.times()[t],
                                IntervalEnd::RightCensored { .. } => {
                                    panic!("complete graph has no censored bars")
                                }
                            }
                    })
                    .count();
                assert_eq!(
                    bars, expected,
                    "q={q} s={s} t={t} p={p} mode={mode:?} edges={edges:?}"
                );
            }
        }
    }
}

#[test]
fn complete_path_intervals_match_independent_all_stage_pair_ranks() {
    let pairs = [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)];
    for mask in 0..64 {
        let edges: Vec<_> = pairs
            .iter()
            .enumerate()
            .filter_map(|(i, &(a, b))| {
                (mask & (1 << i) != 0).then_some((a, b, ((i + mask) % 3) as f64))
            })
            .collect();
        for p in [2, 3, 4_294_967_291] {
            for mode in [Mode::Regular, Mode::NonRegular] {
                let h = compute_path_persistence(
                    3,
                    &edges,
                    2,
                    PrimeField::new(p).unwrap(),
                    mode,
                    1000,
                    100000,
                )
                .unwrap();
                assert_eq!(h.diagram().coverage(), Coverage::Complete);
                assert_eq!(
                    h.context().filtration().scale(),
                    FiltrationScale::Unspecified
                );
                persistent_ranks(&h, &edges, mode);
                let low = compute_low_dimensional_path_persistence(
                    3,
                    &edges,
                    PrimeField::new(p).unwrap(),
                    mode,
                    1000,
                    100000,
                )
                .unwrap();
                for q in 0..=1 {
                    assert_eq!(
                        low.diagram()
                            .dimension(q)
                            .unwrap()
                            .iter()
                            .collect::<Vec<_>>(),
                        h.diagram().dimension(q).unwrap().iter().collect::<Vec<_>>()
                    );
                }
                assert_eq!(low.context(), h.context());
                let mut reversed = edges.clone();
                reversed.reverse();
                assert_eq!(
                    h,
                    compute_path_persistence(
                        3,
                        &reversed,
                        2,
                        PrimeField::new(p).unwrap(),
                        mode,
                        1000,
                        100000
                    )
                    .unwrap()
                );
            }
        }
    }
}

#[test]
fn low_dimensional_presentation_caps_and_model_contracts_are_explicit() {
    let edges = [(0, 1, 0.0), (1, 0, 0.0)];
    let field = PrimeField::new(4_294_967_291).unwrap();
    for mode in [Mode::Regular, Mode::NonRegular] {
        let h = compute_low_dimensional_path_persistence(2, &edges, field, mode, 6, 1000).unwrap();
        assert_eq!(h.two_path_count(), 2);
        assert_eq!(h.edges(), [(0, 1), (1, 0)]);
        assert_eq!(
            h.diagram().dimension(1).unwrap().len(),
            usize::from(mode == Mode::NonRegular)
        );
        assert_eq!(h.clone().into_data(), *h.as_ref());
        assert!(compute_low_dimensional_path_persistence(2, &edges, field, mode, 5, 1000).is_err());
        assert!(compute_low_dimensional_path_persistence(2, &edges, field, mode, 6, 1).is_err());
        assert!(matches!(
            compute_low_dimensional_path_persistence_with(
                2,
                &edges,
                field,
                mode,
                6,
                1000,
                &Execution::default().max_work(10)
            ),
            Err(Error::WorkLimitExceeded { .. })
        ));
    }
    assert!(
        compute_low_dimensional_path_persistence(
            2,
            &[(0, 1, f64::NAN)],
            field,
            Mode::Regular,
            10,
            100
        )
        .is_err()
    );
    let empty =
        compute_low_dimensional_path_persistence(0, &[], field, Mode::Regular, 1, 10).unwrap();
    assert!(empty.diagram().is_empty());
    assert_eq!(empty.diagram().max_dimension(), 1);
}

#[test]
fn path_result_context_composition_and_whole_operation_limits_are_preserved() {
    let edges = [(0, 1, 1.0), (1, 2, 1.0), (2, 0, 1.0), (0, 2, 2.0)];
    let field = PrimeField::new(3).unwrap();
    let h = compute_path_persistence(3, &edges, 1, field, Mode::Regular, 100, 10000).unwrap();
    let data: &PersistenceData = h.as_ref();
    assert_eq!(
        data.context().filtration_kind(),
        FiltrationKind::RegularDirectedPaths
    );
    assert!(matches!(
        data.context().filtration().source(),
        FiltrationSource::DirectedPaths {
            mode: Mode::Regular
        }
    ));
    assert_eq!(
        cocycle::descriptors::betti_curve(data.diagram(), 1, &[0.0, 1.0, 2.0]).unwrap(),
        [0, 1, 0]
    );
    assert_eq!(h.clone().into_data(), *data);
    for work in [0, 100, 1000] {
        assert!(matches!(
            compute_path_persistence_with(
                3,
                &edges,
                1,
                field,
                Mode::Regular,
                100,
                10000,
                &Execution::default().max_work(work)
            ),
            Err(Error::WorkLimitExceeded { .. })
        ));
    }
    let flag = AtomicBool::new(true);
    assert!(matches!(
        compute_path_persistence_with(
            3,
            &edges,
            1,
            field,
            Mode::Regular,
            100,
            10000,
            &Execution::default().cancellation(&flag)
        ),
        Err(Error::Cancelled)
    ));
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            compute_path_persistence(2, &[(0, 1, bad)], 1, field, Mode::Regular, 100, 1000)
                .is_err()
        );
    }
    assert!(
        compute_path_persistence(
            2,
            &[(0, 1, 1.0), (0, 1, 2.0)],
            1,
            field,
            Mode::Regular,
            100,
            1000
        )
        .is_err()
    );
    assert!(compute_path_persistence(0, &[], usize::MAX, field, Mode::Regular, 100, 1000).is_err());
    let empty = compute_path_persistence(0, &[], 3, field, Mode::Regular, 1, 1000).unwrap();
    assert_eq!(empty.diagram().max_dimension(), 3);
    assert!(empty.diagram().is_empty());
    let zero =
        compute_path_persistence(2, &[(0, 1, -0.0)], 1, field, Mode::Regular, 100, 1000).unwrap();
    assert_eq!(zero.times(), [0.0]);
    assert_eq!(zero.times()[0].to_bits(), 0.0f64.to_bits());
}

#[test]
fn projections_homotopies_and_dual_pairings_have_independent_dense_checks() {
    let edges = vec![(0, 1), (1, 2), (2, 0), (0, 3), (3, 2)];
    let source = DirectedGraph::new(4, edges)
        .unwrap()
        .allowed_paths(3, 1000)
        .unwrap();
    for p in [2, 3, 4_294_967_291] {
        for mode in [Mode::Regular, Mode::NonRegular] {
            let h = compute_path_homology(&source, 2, PrimeField::new(p).unwrap(), mode, 100000)
                .unwrap();
            let projection = h.projection(100000).unwrap();
            let dual = h.cohomology(100000).unwrap();
            assert_eq!(projection.field(), h.field());
            assert_eq!(dual.mode(), mode);
            for n in 0..=2 {
                let width = h.omega_bases(n).unwrap().len();
                let upper = h.omega_bases(n + 1).unwrap().len();
                let lower = if n == 0 {
                    0
                } else {
                    h.omega_bases(n - 1).unwrap().len()
                };
                let b = dense(h.boundaries(n).unwrap(), lower, p);
                let next = dense(h.boundaries(n + 1).unwrap(), width, p);
                let homotopy = dense(projection.homotopies(n).unwrap(), upper, p);
                let proj = dense(projection.projections(n).unwrap(), width, p);
                let a = dense(
                    projection.annotations(n).unwrap(),
                    h.representatives(n).unwrap().len(),
                    p,
                );
                let cochains = dense(dual.representatives(n).unwrap(), width, p);
                for j in 0..width {
                    let mut rhs = apply(&next, &homotopy[j], width, p);
                    if n > 0 {
                        let lower_h = dense(projection.homotopies(n - 1).unwrap(), width, p);
                        let extra = apply(&lower_h, &b[j], width, p);
                        for (x, y) in rhs.iter_mut().zip(extra) {
                            *x = add(*x, y, p);
                        }
                    }
                    let lhs: Vec<_> = proj[j]
                        .iter()
                        .enumerate()
                        .map(|(i, &x)| add(u32::from(i == j), if x == 0 { 0 } else { p - x }, p))
                        .collect();
                    assert_eq!(lhs, rhs);
                    assert_eq!(apply(&proj, &proj[j], width, p), proj[j]);
                    for (i, c) in cochains.iter().enumerate() {
                        assert_eq!(c[j], a[j][i]);
                    }
                }
                let s = dense(h.omega_bases(n).unwrap(), h.paths(n).unwrap().len(), p);
                for (i, c) in cochains.iter().enumerate() {
                    for column in &next {
                        assert_eq!(
                            column
                                .iter()
                                .zip(c)
                                .fold(0, |z, (&x, &y)| add(z, product(x, y, p), p)),
                            0
                        );
                    }
                    // Evaluations on the original path representatives are checked
                    // by expressing Z through the dense augmented system S x=Z.
                    for (j, r) in h.representatives(n).unwrap().iter().enumerate() {
                        let mut augmented = s.clone();
                        augmented.extend(dense(
                            std::slice::from_ref(r),
                            h.paths(n).unwrap().len(),
                            p,
                        ));
                        let dependencies = kernel(&augmented, h.paths(n).unwrap().len(), p).1;
                        let witness = dependencies.iter().find(|v| v[width] != 0).unwrap();
                        let scale = if witness[width] == 0 {
                            0
                        } else {
                            p - witness[width]
                        };
                        let coords: Vec<_> = witness[..width]
                            .iter()
                            .map(|&x| product(x, inverse(scale, p), p))
                            .collect();
                        let value = coords
                            .iter()
                            .zip(c)
                            .fold(0, |z, (&x, &y)| add(z, product(x, y, p), p));
                        assert_eq!(value, u32::from(i == j));
                    }
                }
            }
            assert!(h.projection(1).is_err());
            assert!(matches!(
                h.projection_with(100000, &Execution::default().max_work(0)),
                Err(Error::WorkLimitExceeded { .. })
            ));
        }
    }
}

#[test]
fn minimum_integer_basis_matches_all_cycle_subsets_and_preserves_large_costs() {
    let pairs = [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)];
    for mask in 0..64 {
        let edges: Vec<_> = pairs
            .iter()
            .enumerate()
            .filter_map(|(i, &e)| (mask & (1 << i) != 0).then_some(e))
            .collect();
        let graph = DirectedGraph::new(3, edges.clone()).unwrap();
        let costs: Vec<_> = (0..edges.len()).map(|i| ((i + mask) % 4) as u64).collect();
        for mode in [Mode::Regular, Mode::NonRegular] {
            let result =
                compute_minimum_path_basis(&graph, &costs, mode, 1000, 100000, 1000, 10000)
                    .unwrap();
            let source = graph.allowed_paths(2, 1000).unwrap();
            let h = compute_path_homology(&source, 1, PrimeField::default(), mode, 100000).unwrap();
            certificate(&h);
            let s = dense(h.omega_bases(1).unwrap(), edges.len(), 2);
            let boundaries = dense(h.boundaries(2).unwrap(), s.len(), 2);
            let mut span: Vec<_> = boundaries
                .iter()
                .map(|v| apply(&s, v, edges.len(), 2))
                .collect();
            let mut rank = kernel(&span, edges.len(), 2).0;
            let mut cycles = Vec::new();
            for subset in 1usize..1 << edges.len() {
                let vector: Vec<_> = (0..edges.len())
                    .map(|i| u32::from(subset & (1 << i) != 0))
                    .collect();
                let mut parity = [0; 3];
                let mut cost = 0u128;
                for (i, &(a, b)) in edges.iter().enumerate() {
                    if vector[i] != 0 {
                        parity[a] ^= 1;
                        parity[b] ^= 1;
                        cost += u128::from(costs[i]);
                    }
                }
                if parity == [0; 3] {
                    cycles.push((cost, vector));
                }
            }
            cycles.sort();
            let mut expected = 0u128;
            let mut selected = 0;
            for (cost, cycle) in cycles {
                let mut candidate = span.clone();
                candidate.push(cycle.clone());
                let next = kernel(&candidate, edges.len(), 2).0;
                if next > rank {
                    span.push(cycle);
                    rank = next;
                    expected += cost;
                    selected += 1;
                }
            }
            assert_eq!(
                result.total_cost(),
                expected,
                "edges={edges:?} mode={mode:?}"
            );
            assert_eq!(result.representatives().len(), selected);
        }
    }
    let graph = DirectedGraph::new(3, vec![(0, 1), (1, 2), (2, 0)]).unwrap();
    assert_eq!(
        compute_minimum_path_basis(&graph, &[u64::MAX; 3], Mode::Regular, 100, 10000, 100, 1000)
            .unwrap()
            .total_cost(),
        3 * u128::from(u64::MAX)
    );
    assert!(
        compute_minimum_path_basis(&graph, &[1; 3], Mode::Regular, 100, 10000, 100, 1).is_err()
    );
    assert!(
        compute_minimum_path_basis(&graph, &[1; 2], Mode::Regular, 100, 10000, 100, 1000).is_err()
    );
    assert!(matches!(
        compute_minimum_path_basis_with(
            &graph,
            &[1; 3],
            Mode::Regular,
            100,
            10000,
            100,
            1000,
            &Execution::default().max_work(100)
        ),
        Err(Error::WorkLimitExceeded { .. })
    ));
}
