//! Big-steps critical sets on a frozen simplicial filtration over F2.
//!
//! This opt-in workspace reuses the serial lazy boundary reducer. It retains
//! sparse R/V, builds the dual decomposition only when requested, and solves
//! bounded U rows from V instead of materializing the inverse. Ordinary
//! persistence and representative requests keep their existing implementations.
//! Target combination uses maximum absolute displacement, a heuristic rather
//! than an exact diagram-loss derivative or a general descent guarantee.
//!
//! See Nigmetov and Morozov, *Topological Optimization with Big Steps*,
//! <https://arxiv.org/abs/2203.16748v2>, and the
//! [usage guide](https://github.com/Aequiludium/cocycle-rs/blob/main/docs/guides/critical-sets.md).

use crate::algebra::{PrimeField, column::Column, reduction};
use crate::complex::{SimplexId, SimplicialComplex};
use crate::execution::{Execution, WorkBudget};
use crate::{Error, Result};
use std::collections::BTreeMap;

/// Borrowed-source working state for singleton critical sets over F2.
///
/// The source is the entire supplied, immutable complex; no cutoff or implicit
/// Rips expansion is inferred. IDs and units belong to that source. Sparse
/// reduction may still have quadratic fill-in and cubic worst-case work.
/// Primal R/V are retained; dual R/V and per-dimension V transposes are lazy.
/// Each U query solves only the requested row through its target-value bound;
/// U rows are not retained. This type never changes simplex or vertex values.
///
/// Supports finite birth/death endpoints and essential births. Same-dimensional
/// sets omit the faces/cofaces needed for arbitrary explicit value updates;
/// reconstruct a lower-star filtration or apply closure before updating data.
/// Parallel/cleared decompositions cannot be imported: the serial F2 reduction
/// guarantees the ELZ/lazy support property required by the critical-set formulas.
pub struct CriticalSetWorkspace<'a> {
    complex: &'a SimplicialComplex,
    primal: Decomposition,
    dual: Option<Decomposition>,
}

struct Decomposition {
    reduction: reduction::BoundaryReduction,
    // Only requested dimensions are transposed. A sparse map avoids an n-wide
    // row array per dimension. Complete caches are committed after successful work.
    rows: BTreeMap<usize, BTreeMap<usize, Column<usize>>>,
}

impl<'a> CriticalSetWorkspace<'a> {
    /// Retain a sparse primal R/V decomposition of the complete source over F2.
    ///
    /// # Errors
    /// Returns allocation or internal-invariant errors during reduction.
    pub fn new(complex: &'a SimplicialComplex) -> Result<Self> {
        Self::new_with(complex, &Execution::default())
    }

    /// Construct under one cooperative reduction budget.
    ///
    /// # Errors
    /// Includes [`Self::new`] errors, cancellation and work exhaustion.
    pub fn new_with(complex: &'a SimplicialComplex, execution: &Execution<'_>) -> Result<Self> {
        let mut budget = WorkBudget::new(execution)?;
        let primal = decompose(complex, false, &mut budget)?;
        budget.check()?;
        Ok(Self {
            complex,
            primal,
            dual: None,
        })
    }

    /// Finite (birth, death) IDs in increasing birth-ID order, including zero bars.
    ///
    /// These pairs come from the full supplied complex over F2. They do not
    /// imply coverage of a larger geometric source or an omitted skeleton.
    pub fn finite_pairs(&self) -> impl Iterator<Item = (SimplexId, SimplexId)> + '_ {
        self.primal
            .reduction
            .deaths
            .iter()
            .enumerate()
            .filter_map(|(b, &d)| d.map(|d| (SimplexId(b), SimplexId(d))))
    }

    /// Essential birth IDs in increasing filtration order, relative to this source.
    pub fn essential_births(&self) -> impl Iterator<Item = SimplexId> + '_ {
        self.primal
            .reduction
            .reduced
            .iter()
            .enumerate()
            .filter_map(|(i, r)| {
                (r.is_empty() && self.primal.reduction.deaths[i].is_none()).then_some(SimplexId(i))
            })
    }

    /// Find the same-dimensional critical set for moving an endpoint to `target`.
    ///
    /// Role is determined from the primal pairing. Increase death / decrease
    /// finite birth solve bounded primal / dual U rows, respectively. Decrease
    /// death / increase birth inspect the corresponding V column. Decreasing an
    /// essential birth uses primal V; increasing it uses dual V. Returned IDs are
    /// sorted in the source's original filtration order. Ties are included in
    /// the closed value window; a zero displacement returns an empty set and
    /// performs no lazy decomposition or U work.
    ///
    /// The set is a proposal, not a modified filtration or tracked-point guarantee
    /// after arbitrary simultaneous updates. Multiplicity and nondifferentiable
    /// ties need the paper's qualifications; faces/cofaces are not included.
    ///
    /// # Errors
    /// Rejects out-of-range source IDs and nonfinite targets. Includes reduction,
    /// allocation and internal-invariant errors. An ID from another complex with
    /// the same numeric index cannot be distinguished; source identity is required.
    pub fn critical_set(&mut self, endpoint: SimplexId, target: f64) -> Result<Vec<SimplexId>> {
        self.critical_set_with(endpoint, target, &Execution::default())
    }

    /// Query under a fresh cooperative budget, including lazy preparation.
    ///
    /// Complete cached V transposes may survive interruption; incomplete ones
    /// are discarded. No partial critical set is returned. A retry creates a
    /// fresh budget and preserves the source and already-valid decompositions.
    ///
    /// # Errors
    /// Includes [`Self::critical_set`] errors, cancellation and work exhaustion.
    pub fn critical_set_with(
        &mut self,
        endpoint: SimplexId,
        target: f64,
        execution: &Execution<'_>,
    ) -> Result<Vec<SimplexId>> {
        let mut budget = WorkBudget::new(execution)?;
        self.critical(endpoint, target, &mut budget)
    }

    /// Combine singleton endpoint proposals by maximum absolute displacement.
    ///
    /// For each `(endpoint, target)` compute its critical set. A simplex affected
    /// by multiple sets receives the target farthest from its original value,
    /// with the first proposal winning equal-displacement ties. Output contains
    /// only affected IDs, sorted by original filtration order. Empty input and
    /// zero-displacement proposals produce no targets. Duplicates are allowed.
    ///
    /// This is the paper's maximum-displacement heuristic and Oineus's `Max`
    /// rule at Oineus revision `e52814a1ffb5b8a81e71ff1e93b4c14194673f0f`,
    /// not a numeric maximum target or sum
    /// of gradients. The frozen-target surrogate gradient is `2*(f-target)`;
    /// callers must separately handle data derivatives and filtration validity.
    /// No generic descent or convergence guarantee is provided.
    ///
    /// # Errors
    /// Includes [`Self::critical_set`] errors. Rejects an unrepresentable target
    /// displacement rather than using overflowed infinity to compare proposals.
    pub fn targets(&mut self, endpoints: &[(SimplexId, f64)]) -> Result<Vec<(SimplexId, f64)>> {
        self.targets_with(endpoints, &Execution::default())
    }

    /// Combine under one budget spanning all queries and conflict handling.
    ///
    /// # Errors
    /// Includes [`Self::targets`] errors, cancellation and work exhaustion.
    pub fn targets_with(
        &mut self,
        endpoints: &[(SimplexId, f64)],
        execution: &Execution<'_>,
    ) -> Result<Vec<(SimplexId, f64)>> {
        let mut budget = WorkBudget::new(execution)?;
        let mut targets: BTreeMap<SimplexId, f64> = BTreeMap::new();
        for &(endpoint, target) in endpoints {
            budget.step()?;
            for id in self.critical(endpoint, target, &mut budget)? {
                budget.step()?;
                let value = self.complex.simplices()[id.0].value();
                let distance = displacement(value, target)?;
                if let Some(previous) = targets.get_mut(&id) {
                    if distance > displacement(value, *previous)? {
                        *previous = target;
                    }
                } else {
                    targets.insert(id, target);
                }
            }
        }
        let mut result = Vec::new();
        result
            .try_reserve_exact(targets.len())
            .map_err(|_| allocation())?;
        for item in targets {
            budget.step()?;
            result.push(item);
        }
        budget.check()?;
        Ok(result)
    }

    fn critical(
        &mut self,
        endpoint: SimplexId,
        target: f64,
        budget: &mut WorkBudget<'_>,
    ) -> Result<Vec<SimplexId>> {
        budget.step()?;
        let cell = self
            .complex
            .simplex(endpoint)
            .ok_or(Error::InvalidComplex {
                cell: Some(endpoint.0),
                reason: "critical-set ID is outside the source",
            })?;
        if !target.is_finite() {
            return Err(Error::NonFiniteValue {
                field: "critical-set target",
                index: None,
            });
        }
        let value = cell.value();
        if target == value {
            return Ok(Vec::new());
        }
        let negative = !self.primal.reduction.reduced[endpoint.0].is_empty();
        let unpaired = !negative && self.primal.reduction.deaths[endpoint.0].is_none();
        let increase = target > value;
        let dualize = !negative && (!unpaired || increase);
        // Unpaired births use a V column in both directions (paper's remarks).
        let use_u = !unpaired && (negative == increase);
        if dualize && self.dual.is_none() {
            self.dual = Some(decompose(self.complex, true, budget)?);
        }
        let state = if dualize {
            self.dual.as_mut().unwrap()
        } else {
            &mut self.primal
        };
        let n = self.complex.len();
        let index = if dualize {
            n - 1 - endpoint.0
        } else {
            endpoint.0
        };
        let support = if use_u {
            state.u_row(
                self.complex,
                dualize,
                index,
                cell.dimension(),
                target,
                budget,
            )?
        } else {
            // Cloning a sparse V column keeps the query's ownership simple; no
            // full inverse or dense matrix is constructed.
            let mut result = Column::new();
            for (&j, &coefficient) in state.reduction.transforms[index].entries() {
                budget.step()?;
                let original = if dualize { n - 1 - j } else { j };
                let v = self.complex.simplices()[original].value();
                if (increase && v <= target) || (!increase && v >= target) {
                    result.add_term(j, coefficient, PrimeField::default());
                }
            }
            result
        };
        let mut ids = Vec::new();
        for (&j, _) in support.entries() {
            budget.step()?;
            ids.try_reserve(1).map_err(|_| allocation())?;
            ids.push(SimplexId(if dualize { n - 1 - j } else { j }));
        }
        if dualize {
            ids.reverse();
        }
        budget.check()?;
        Ok(ids)
    }
}

impl Decomposition {
    fn u_row(
        &mut self,
        complex: &SimplicialComplex,
        dualize: bool,
        row: usize,
        dimension: usize,
        target: f64,
        budget: &mut WorkBudget<'_>,
    ) -> Result<Column<usize>> {
        let n = complex.len();
        if !self.rows.contains_key(&dimension) {
            let mut rows: BTreeMap<usize, Column<usize>> = BTreeMap::new();
            for (j, column) in self.reduction.transforms.iter().enumerate() {
                budget.step()?;
                let original = if dualize { n - 1 - j } else { j };
                if complex.simplices()[original].dimension() == dimension {
                    for (&i, &coefficient) in column.entries() {
                        budget.step()?;
                        rows.entry(i).or_insert_with(Column::new).add_term(
                            j,
                            coefficient,
                            PrimeField::default(),
                        );
                    }
                }
            }
            budget.check()?;
            self.rows.insert(dimension, rows);
        }
        let rows = &self.rows[&dimension];
        let mut residual = Column::unit(row);
        let mut result = Column::new();
        while let Some((&pivot, &coefficient)) = residual.first() {
            budget.step()?;
            let original = if dualize { n - 1 - pivot } else { pivot };
            let value = complex.simplices()[original].value();
            if (!dualize && value > target) || (dualize && value < target) {
                break;
            }
            let v_row = rows.get(&pivot).ok_or(Error::InternalInvariant {
                reason: "critical-set V transpose is missing a diagonal row",
            })?;
            if v_row.first() != Some((&pivot, &1)) {
                return Err(Error::InternalInvariant {
                    reason: "critical-set V must be unit upper triangular over F2",
                });
            }
            result.add_term(pivot, coefficient, PrimeField::default());
            residual.add_scaled(v_row, coefficient, PrimeField::default(), &mut || {
                budget.step()
            })?;
        }
        Ok(result)
    }
}

fn decompose(
    complex: &SimplicialComplex,
    dualize: bool,
    budget: &mut WorkBudget<'_>,
) -> Result<Decomposition> {
    let n = complex.len();
    let mut columns = Vec::new();
    columns.try_reserve_exact(n).map_err(|_| allocation())?;
    columns.resize_with(n, Column::new);
    for i in 0..n {
        budget.step()?;
        for term in complex.boundary(SimplexId(i)).unwrap() {
            budget.step()?;
            let (column, row) = if dualize {
                (n - 1 - term.face.0, n - 1 - i)
            } else {
                (i, term.face.0)
            };
            columns[column].add_term(row, 1, PrimeField::default());
        }
    }
    Ok(Decomposition {
        reduction: reduction::reduce(columns, PrimeField::default(), &mut || budget.step())?,
        rows: BTreeMap::new(),
    })
}

fn displacement(value: f64, target: f64) -> Result<f64> {
    let distance = (target - value).abs();
    if !distance.is_finite() {
        return Err(Error::NumericalFailure {
            context: "critical-set target displacement",
        });
    }
    Ok(distance)
}

fn allocation() -> Error {
    Error::AllocationFailed {
        context: "critical-set workspace",
    }
}

#[cfg(test)]
mod tests;
