//! Ordinary interval normalization shared by the existing PH algorithms.
use crate::diagram::{Coverage, IntervalEnd, PersistenceDiagram, PersistenceInterval};
use crate::{Error, Result};

/// Raw births and optional deaths in the computed range.
pub(super) type RawIntervals = Vec<(usize, f64, Option<f64>)>;

/// Shared by algorithms; None means unpaired in the computed range, not necessarily essential.
pub(super) fn assemble_diagram(
    max_dimension: usize,
    coverage: Coverage,
    raw: impl IntoIterator<Item = (usize, f64, Option<f64>)>,
) -> Result<PersistenceDiagram> {
    let mut intervals = Vec::new();
    for (dimension, birth, death) in raw {
        if dimension > max_dimension || death == Some(birth) {
            continue;
        }
        let end = match death {
            Some(value) => IntervalEnd::Finite(value),
            None => match coverage {
                Coverage::Complete => IntervalEnd::Essential,
                Coverage::Through(through) => IntervalEnd::RightCensored { through },
            },
        };
        let interval = PersistenceInterval::new(dimension, birth, end)?;
        intervals
            .try_reserve(1)
            .map_err(|_| Error::AllocationFailed {
                context: "persistence diagram",
            })?;
        intervals.push(interval);
    }
    PersistenceDiagram::new(max_dimension, coverage, intervals)
}
