//! Public directed-source, ambient-boundary and resource contracts.
use cocycle::Error;
use cocycle::complex::{
    DirectedGraph, PathBoundaryMode as Mode, PathComplex, path_boundary, path_boundary_with,
};
use cocycle::execution::Execution;
use std::collections::BTreeMap;
use std::sync::atomic::AtomicBool;

#[test]
fn direction_return_paths_isolates_and_explicit_truncation_are_preserved() {
    let graph = DirectedGraph::new(3, vec![(1, 0), (0, 1)]).unwrap();
    assert_eq!(graph.edges(), [(0, 1), (1, 0)]);
    assert_eq!(graph.vertex_count(), 3);
    let paths = graph.allowed_paths(2, 7).unwrap();
    assert_eq!(paths.paths(0), [vec![0], vec![1], vec![2]]);
    assert_eq!(paths.paths(2), [vec![0, 1, 0], vec![1, 0, 1]]);
    assert_eq!(paths.truncated_at(), Some(2));
    let explicit =
        PathComplex::new(3, &[vec![0], vec![1], vec![2], vec![0, 1], vec![1, 0]], 5).unwrap();
    assert!(explicit.paths(2).is_empty());
    assert_eq!(explicit.truncated_at(), None);
    assert!(PathComplex::new(3, &[vec![0], vec![1], vec![2], vec![0, 1, 0]], 4).is_err());
    assert!(DirectedGraph::new(2, vec![(0, 0)]).is_err());
    assert!(DirectedGraph::new(2, vec![(0, 1), (0, 1)]).is_err());
}

#[test]
fn full_and_regular_boundaries_have_independent_integer_square_checks() {
    assert!(
        path_boundary(&[0, 1, 0], Mode::NonRegular)
            .unwrap()
            .contains(&(vec![0, 0], -1))
    );
    assert!(
        !path_boundary(&[0, 1, 0], Mode::Regular)
            .unwrap()
            .iter()
            .any(|(w, _)| w == &[0, 0])
    );
    for length in 1..=6 {
        for mask in 0..2usize.pow(length) {
            let word: Vec<_> = (0..length).map(|i| (mask >> i) & 1).collect();
            for mode in [Mode::Regular, Mode::NonRegular] {
                if mode == Mode::Regular && word.windows(2).any(|p| p[0] == p[1]) {
                    continue;
                }
                let mut square = BTreeMap::new();
                for (face, c) in path_boundary(&word, mode).unwrap() {
                    for (term, v) in path_boundary(&face, mode).unwrap() {
                        *square.entry(term).or_insert(0) += c * v;
                    }
                }
                assert!(square.values().all(|&v| v == 0));
            }
        }
    }
}

#[test]
fn cancellation_and_path_work_caps_return_errors_without_partial_sources() {
    let graph = DirectedGraph::new(2, vec![(0, 1), (1, 0)]).unwrap();
    assert!(graph.allowed_paths(2, 5).is_err());
    assert!(matches!(
        graph.allowed_paths_with(2, 100, &Execution::default().max_work(0)),
        Err(Error::WorkLimitExceeded { .. })
    ));
    let flag = AtomicBool::new(true);
    assert!(matches!(
        graph.allowed_paths_with(2, 100, &Execution::default().cancellation(&flag)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        DirectedGraph::new_with(2, vec![(0, 1)], &Execution::default().cancellation(&flag)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        PathComplex::new_with(1, &[vec![0]], 1, &Execution::default().cancellation(&flag)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        path_boundary_with(&[0, 1, 0], Mode::Regular, &Execution::default().max_work(3)),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(matches!(
        path_boundary_with(
            &[0],
            Mode::Regular,
            &Execution::default().cancellation(&flag)
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(
        DirectedGraph::new(0, vec![])
            .unwrap()
            .allowed_paths(usize::MAX, 1)
            .unwrap()
            .path_count(),
        0
    );
}

#[test]
fn explicit_order_empty_extreme_dimensions_and_invalid_words_are_checked() {
    let paths = [vec![1, 0], vec![1], vec![0, 1], vec![0]];
    let p = PathComplex::new(2, &paths, 4).unwrap();
    assert_eq!(p.paths(1), [vec![0, 1], vec![1, 0]]);
    assert!(p.paths(usize::MAX).is_empty());
    assert_eq!(p.max_dimension(), Some(1));
    assert_eq!(PathComplex::new(0, &[], 1).unwrap().max_dimension(), None);
    assert!(PathComplex::new(1, &[vec![0], vec![0]], 2).is_err());
    assert!(PathComplex::new(1, &[vec![0], vec![0, 0]], 2).is_err());
    assert!(PathComplex::new(1, &[vec![1]], 1).is_err());
    assert!(PathComplex::new(1, &[vec![]], 1).is_err());
    assert!(PathComplex::new(usize::MAX, &[], 1).is_err());
    assert!(PathComplex::new(0, &[], 0).is_err());
    assert!(path_boundary(&[], Mode::NonRegular).is_err());
    assert!(path_boundary(&[0, 0], Mode::Regular).is_err());
    assert!(path_boundary(&[0, 0], Mode::NonRegular).unwrap().is_empty());
}

#[test]
fn all_three_vertex_graph_paths_match_independent_vertex_words() {
    let pairs: Vec<_> = (0..3)
        .flat_map(|a| (0..3).filter(move |&b| a != b).map(move |b| (a, b)))
        .collect();
    for mask in 0..64 {
        let edges: Vec<_> = pairs
            .iter()
            .enumerate()
            .filter_map(|(i, &e)| (mask & (1 << i) != 0).then_some(e))
            .collect();
        let graph = DirectedGraph::new(3, edges.clone()).unwrap();
        let complex = graph.allowed_paths(3, 10000).unwrap();
        for length in 1..=4u32 {
            let mut expected = Vec::new();
            for mut code in 0..3usize.pow(length) {
                let mut word = vec![0; length as usize];
                for v in word.iter_mut().rev() {
                    *v = code % 3;
                    code /= 3;
                }
                if word.windows(2).all(|p| edges.contains(&(p[0], p[1]))) {
                    expected.push(word);
                }
            }
            expected.sort();
            assert_eq!(complex.paths(length as usize - 1), expected);
        }
    }
}
