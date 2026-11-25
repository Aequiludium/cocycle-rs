//! Exact non-reduced path chains and homology over a bound prime field.
//!
//! Allowed paths are data, not integer filtered cells. Omega is rebuilt in the
//! selected field; its modular coordinate boundaries never implement that trait.
mod linear;
use crate::algebra::PrimeField;
use crate::complex::{Path, PathBoundaryMode, PathComplex, path_boundary_budget};
use crate::execution::{Execution, WorkBudget};
use crate::{Error, Result};
use linear::*;
use std::collections::BTreeMap;

/// Sparse nonzero canonical residues indexed in a documented path/Omega basis.
pub type PathVector = BTreeMap<usize, u32>;

/// Owned exact path homology and its field-bound coordinate chain through d+1.
///
/// Bases are deterministic, without a canonical/minimal/orthogonal claim. This
/// result owns path identities and stays interpretable after releasing its source.
/// Read-only access prevents changing the field without recomputing Omega.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathHomology {
    field: PrimeField,
    mode: PathBoundaryMode,
    paths: Vec<Vec<Path>>,
    omega_bases: Vec<Vec<PathVector>>,
    boundaries: Vec<Vec<PathVector>>,
    boundary_ranks: Vec<usize>,
    representatives: Vec<Vec<PathVector>>,
}

impl PathHomology {
    /// Field used for every basis and boundary, including the upper layer.
    pub fn field(&self) -> PrimeField {
        self.field
    }
    /// Ambient boundary convention; allowed inputs remain regular.
    pub fn mode(&self) -> PathBoundaryMode {
        self.mode
    }
    /// Largest computed homology degree; the chain retains one extra degree.
    pub fn max_homology_dimension(&self) -> usize {
        self.representatives.len() - 1
    }
    /// Allowed paths in one retained degree; None means outside the stored window.
    pub fn paths(&self, dimension: usize) -> Option<&[Path]> {
        self.paths.get(dimension).map(Vec::as_slice)
    }
    /// Omega columns in the corresponding allowed-path coordinates, through d+1.
    pub fn omega_bases(&self, dimension: usize) -> Option<&[PathVector]> {
        self.omega_bases.get(dimension).map(Vec::as_slice)
    }
    /// Boundary columns in the preceding Omega basis; degree-zero columns are empty.
    pub fn boundaries(&self, dimension: usize) -> Option<&[PathVector]> {
        self.boundaries.get(dimension).map(Vec::as_slice)
    }
    /// Boundary ranks through the upper degree d+1.
    pub fn boundary_ranks(&self) -> &[usize] {
        &self.boundary_ranks
    }
    /// Homology basis cycles in allowed-path coordinates, through degree d.
    pub fn representatives(&self, dimension: usize) -> Option<&[PathVector]> {
        self.representatives.get(dimension).map(Vec::as_slice)
    }
    /// Betti numbers in the contiguous computed domain 0..=d.
    pub fn betti(&self) -> Vec<usize> {
        self.representatives.iter().map(Vec::len).collect()
    }
}

/// Compute exact non-reduced path homology using unlimited cooperative controls.
///
/// `max_entries` is positive and bounds output layer/column counts, each sparse
/// matrix and each elimination's retained plus temporary entries/transforms.
/// It is not a byte or aggregate live-memory cap. H_d retains Omega_(d+1).
///
/// # Errors
/// Rejects insufficient graph truncation, overflowing dimensions, exhausted
/// entry limits, failed reservations or broken internal chain identities.
pub fn compute_path_homology(
    source: &PathComplex,
    dimension: usize,
    field: PrimeField,
    mode: PathBoundaryMode,
    max_entries: usize,
) -> Result<PathHomology> {
    compute_path_homology_with(
        source,
        dimension,
        field,
        mode,
        max_entries,
        &Execution::default(),
    )
}

/// Compute with one work budget across boundary construction, kernels and witnesses.
///
/// Coordinate, face, column, elimination and arithmetic visits count repeatedly.
/// Already validated source construction is outside this terminal's budget.
/// B-tree node allocations and individual container operations cannot be made
/// fallible/interruptible; no unsafe code or partial successful output is used.
///
/// # Errors
/// Returns the errors of [`compute_path_homology`], cancellation or work exhaustion.
pub fn compute_path_homology_with(
    source: &PathComplex,
    dimension: usize,
    field: PrimeField,
    mode: PathBoundaryMode,
    max_entries: usize,
    execution: &Execution<'_>,
) -> Result<PathHomology> {
    let mut budget = WorkBudget::new(execution)?;
    homology(source, dimension, field, mode, max_entries, &mut budget)
}

fn layers(source: &PathComplex, dimension: usize, limit: usize) -> Result<usize> {
    if limit == 0 {
        return Err(Error::InvalidParameter {
            parameter: "max_entries",
            reason: "must be positive",
        });
    }
    let layers = dimension.checked_add(2).ok_or(Error::SizeOverflow {
        operation: "path chain degrees",
    })?;
    count(layers, limit)?;
    if source
        .truncated_at()
        .is_some_and(|available| available < layers - 1)
    {
        return Err(Error::InvalidParameter {
            parameter: "path dimension",
            reason: "H_d requires allowed paths through d+1 unless the source is complete",
        });
    }
    Ok(layers)
}

fn copy_paths(
    source: &PathComplex,
    layers: usize,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<Vec<Path>>> {
    let mut result = Vec::new();
    for q in 0..layers {
        budget.step()?;
        let mut layer = Vec::new();
        for path in source.paths(q) {
            budget.step_by(path.len())?;
            let mut copied = Vec::new();
            copied
                .try_reserve_exact(path.len())
                .map_err(|_| Error::AllocationFailed {
                    context: "path coordinates",
                })?;
            copied.extend_from_slice(path);
            push(&mut layer, copied)?;
        }
        push(&mut result, layer)?;
    }
    Ok(result)
}

fn homology(
    source: &PathComplex,
    dimension: usize,
    field: PrimeField,
    mode: PathBoundaryMode,
    limit: usize,
    budget: &mut WorkBudget<'_>,
) -> Result<PathHomology> {
    let layers = layers(source, dimension, limit)?;
    count(source.vertex_count(), limit)?;
    let mut vertices = Vec::new();
    let mut zero = Vec::new();
    for i in 0..source.vertex_count() {
        budget.step()?;
        push(&mut vertices, BTreeMap::from([(i, 1)]))?;
        push(&mut zero, PathVector::new())?;
    }
    let mut previous = Echelon::from_columns(&vertices, field, limit, budget)?;
    let mut omega_bases = Vec::new();
    let mut boundaries = Vec::new();
    push(&mut omega_bases, vertices)?;
    push(&mut boundaries, zero)?;
    for q in 1..layers {
        budget.step()?;
        count(source.paths(q).len(), limit)?;
        let mut lower = BTreeMap::new();
        for (i, path) in source.paths(q - 1).iter().enumerate() {
            budget.step_by(path.len())?;
            lower.insert(path.as_slice(), i);
        }
        let mut forbidden_index = BTreeMap::new();
        let (mut allowed, mut forbidden) = (Vec::new(), Vec::new());
        let mut entries = 0;
        for path in source.paths(q) {
            budget.step()?;
            let (mut a, mut f) = (PathVector::new(), PathVector::new());
            for (face, coefficient) in path_boundary_budget(path, mode, budget)? {
                budget.step_by(face.len())?;
                let value = coefficient.rem_euclid(i64::from(field.characteristic())) as u32;
                if value == 0 {
                    continue;
                }
                entries = sum(entries, 1)?;
                count(entries, limit)?;
                if let Some(&row) = lower.get(face.as_slice()) {
                    a.insert(row, value);
                } else {
                    let next = forbidden_index.len();
                    let row = *forbidden_index.entry(face).or_insert(next);
                    f.insert(row, value);
                }
            }
            push(&mut allowed, a)?;
            push(&mut forbidden, f)?;
        }
        let basis = kernel(&forbidden, field, limit, budget)?.1;
        if apply(&forbidden, &basis, field, limit, budget)?
            .iter()
            .any(|v| !v.is_empty())
        {
            return Err(invariant("forbidden boundary times Omega is zero"));
        }
        let images = apply(&allowed, &basis, field, limit, budget)?;
        let mut boundary = Vec::new();
        let mut entries = 0;
        for image in &images {
            budget.step()?;
            let column = previous.coordinates(image, budget)?;
            entries = sum(entries, column.len())?;
            count(entries, limit)?;
            push(&mut boundary, column)?;
        }
        matrix(&boundary, limit)?;
        if apply(&omega_bases[q - 1], &boundary, field, limit, budget)? != images
            || apply(&boundaries[q - 1], &boundary, field, limit, budget)?
                .iter()
                .any(|v| !v.is_empty())
        {
            return Err(invariant("Omega boundary coordinates and square"));
        }
        previous = Echelon::from_columns(&basis, field, limit, budget)?;
        if previous.rank() != basis.len() {
            return Err(invariant("Omega basis is independent"));
        }
        push(&mut omega_bases, basis)?;
        push(&mut boundaries, boundary)?;
    }
    let mut boundary_ranks = Vec::new();
    let mut representatives = Vec::new();
    for q in 0..layers {
        budget.step()?;
        let (rank, cycles) = kernel(&boundaries[q], field, limit, budget)?;
        push(&mut boundary_ranks, rank)?;
        if q + 1 == layers {
            break;
        }
        let images = &boundaries[q + 1];
        let mut quotient = Echelon::from_columns(images, field, limit, budget)?;
        let expected = omega_bases[q]
            .len()
            .checked_sub(sum(rank, quotient.rank())?)
            .ok_or(invariant("nonnegative path Betti number"))?;
        let mut reps = Vec::new();
        for (i, cycle) in cycles.iter().enumerate() {
            if quotient
                .insert(cycle, sum(images.len(), i)?, 0, budget)?
                .is_none()
            {
                let mut image = apply(
                    &omega_bases[q],
                    std::slice::from_ref(cycle),
                    field,
                    limit,
                    budget,
                )?;
                push(
                    &mut reps,
                    image.pop().ok_or(invariant("one image per cycle"))?,
                )?;
                matrix(&reps, limit)?;
            }
        }
        if reps.len() != expected {
            return Err(invariant("complete path homology representatives"));
        }
        push(&mut representatives, reps)?;
    }
    let paths = copy_paths(source, layers, budget)?;
    budget.check()?;
    Ok(PathHomology {
        field,
        mode,
        paths,
        omega_bases,
        boundaries,
        boundary_ranks,
        representatives,
    })
}
