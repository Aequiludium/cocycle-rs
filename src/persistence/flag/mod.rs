//! Shared exact flag engines and supplied-graph computation.
mod clearing;
mod cohomology;
mod dispatch;
mod edges;
mod h0;
mod ordered;
mod selection;
use super::RepresentativeRequest;
use super::simplicial::finish_zero_born;
pub(super) use dispatch::{compute_dense, compute_graph};

use super::execution::WorkBudget;
use super::{ExecutionLimits, PersistenceOptions};
use crate::Result;
use crate::diagram::{ComputationContext, Coverage, FiltrationKind, PersistenceResult};
use crate::filtration::{FlagFiltration, flag::CliqueAccess};

/// Compute ordinary dimension-generic prime-field persistence of a supplied graph's clique filtration.
///
/// Missing edges never enter this filtration. Without a cutoff (or at/above its
/// largest edge), unpaired classes are essential for this supplied graph. A lower
/// cutoff conservatively censors surviving classes. Input is borrowed and never
/// densified. Work/storage include O(n+m) input access plus reduction fill-in.
/// F2 diagram-only requests retain specialized H0/H1 and continue H2+ using
/// complete owned death-triangle clearing keys. H1-only combinatorial IDs must
/// fit `usize`, even for sparse graphs. Higher-dimensional requests fall back
/// to full tuple-based cohomology if that index overflows; odd primes use it directly.
///
/// # Errors
/// Returns size/allocation errors, cancellation, work-limit exhaustion, or an
/// invariant error; never a successful partial result.
pub fn compute_flag(
    input: &FlagFiltration,
    options: &PersistenceOptions,
    limits: &ExecutionLimits<'_>,
) -> Result<PersistenceResult> {
    compute_flag_with_representatives(input, options, &[], limits)
}

/// Compute supplied-graph persistence and requested basis representatives.
///
/// Unlike diagram-only computation, nonempty requests materialize the skeleton
/// through q+1 and retain boundary transformations. Cost can be exponential.
/// The execution budget includes this preparation and representative extraction.
///
/// # Errors
/// Includes [`compute_flag`] errors, uncomputed dimensions and query scales
/// outside coverage. No partial diagram or representative payload is returned.
pub fn compute_flag_with_representatives(
    input: &FlagFiltration,
    options: &PersistenceOptions,
    requests: &[RepresentativeRequest],
    limits: &ExecutionLimits<'_>,
) -> Result<PersistenceResult> {
    let mut budget = WorkBudget::new(limits)?;
    compute_flag_budget(input, options, requests, &mut budget)
}

pub(in crate::persistence) fn compute_flag_budget(
    input: &FlagFiltration,
    options: &PersistenceOptions,
    requests: &[RepresentativeRequest],
    budget: &mut WorkBudget<'_>,
) -> Result<PersistenceResult> {
    let graph = input.graph();
    let (cutoff, coverage) = match options.max_edge() {
        Some(t) if t < graph.max_edge() => (t, Coverage::Through(t)),
        _ => (graph.max_edge(), Coverage::Complete),
    };
    let (diagram, representatives) = finish_zero_born(
        &CliqueAccess::Sparse(graph, cutoff),
        options,
        requests,
        coverage,
        budget,
        |budget| {
            compute_graph(
                graph,
                options.max_homology_dimension(),
                cutoff,
                options.field(),
                budget,
            )
        },
    )?;
    Ok(PersistenceResult::new(
        diagram,
        ComputationContext::new(
            options.field(),
            FiltrationKind::SuppliedFlag,
            graph.vertex_count(),
            options.max_edge(),
            None,
            None,
        ),
        representatives,
    ))
}
