//! Bounded sparse kernel/coordinate solving for field-specific Omega bases.
//!
//! Unlike the ordinary boundary reducer, this also solves rectangular systems
//! and retains dependency witnesses. All callers share their terminal's budget.
use super::*;

pub(super) fn push<T>(output: &mut Vec<T>, value: T) -> Result<()> {
    output.try_reserve(1).map_err(|_| Error::AllocationFailed {
        context: "path homology",
    })?;
    output.push(value);
    Ok(())
}

pub(super) fn count(count: usize, limit: usize) -> Result<()> {
    if count > limit {
        return Err(Error::InvalidParameter {
            parameter: "max_entries",
            reason: "path algebra entry or column limit exceeded",
        });
    }
    Ok(())
}

pub(super) fn sum(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or(Error::SizeOverflow {
        operation: "path algebra entries",
    })
}

pub(super) fn matrix(columns: &[SparseVector], limit: usize) -> Result<()> {
    count(columns.len(), limit)?;
    let mut entries = 0;
    for column in columns {
        entries = sum(entries, column.len())?;
        count(entries, limit)?;
    }
    Ok(())
}

pub(super) fn copy_vector(v: &SparseVector, budget: &mut WorkBudget<'_>) -> Result<SparseVector> {
    let mut output = SparseVector::new();
    budget.step()?;
    for (&row, &value) in v {
        budget.step()?;
        output.insert(row, value);
    }
    Ok(output)
}

pub(super) fn add_scaled(
    target: &mut SparseVector,
    source: &SparseVector,
    factor: u32,
    field: PrimeField,
    retained: usize,
    limit: usize,
    budget: &mut WorkBudget<'_>,
) -> Result<()> {
    count(sum(retained, target.len())?, limit)?;
    for (&row, &value) in source {
        budget.step()?;
        let value = field.add(
            target.get(&row).copied().unwrap_or(0),
            field.multiply(factor, value),
        );
        if value == 0 {
            target.remove(&row);
        } else {
            if !target.contains_key(&row) {
                count(sum(sum(retained, target.len())?, 1)?, limit)?;
            }
            target.insert(row, value);
        }
    }
    Ok(())
}

pub(super) fn apply(
    columns: &[SparseVector],
    vectors: &[SparseVector],
    field: PrimeField,
    limit: usize,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<SparseVector>> {
    count(vectors.len(), limit)?;
    let mut result = Vec::new();
    let mut retained = 0;
    for vector in vectors {
        budget.step()?;
        let mut image = SparseVector::new();
        for (&index, &value) in vector {
            budget.step()?;
            let column = columns
                .get(index)
                .ok_or(invariant("matrix coordinate range"))?;
            add_scaled(&mut image, column, value, field, retained, limit, budget)?;
        }
        retained = sum(retained, image.len())?;
        push(&mut result, image)?;
    }
    Ok(result)
}

pub(super) struct Echelon {
    pivots: BTreeMap<usize, (SparseVector, SparseVector)>,
    entries: usize,
    field: PrimeField,
    limit: usize,
}

impl Echelon {
    pub(super) fn new(field: PrimeField, limit: usize) -> Self {
        Self {
            pivots: BTreeMap::new(),
            entries: 0,
            field,
            limit,
        }
    }

    pub(super) fn from_columns(
        columns: &[SparseVector],
        field: PrimeField,
        limit: usize,
        budget: &mut WorkBudget<'_>,
    ) -> Result<Self> {
        matrix(columns, limit)?;
        let mut solver = Self::new(field, limit);
        for (index, column) in columns.iter().enumerate() {
            solver.insert(column, index, 0, budget)?;
        }
        Ok(solver)
    }

    pub(super) fn rank(&self) -> usize {
        self.pivots.len()
    }

    fn reduce(
        &self,
        vector: &mut SparseVector,
        witness: &mut SparseVector,
        reserved: usize,
        budget: &mut WorkBudget<'_>,
    ) -> Result<()> {
        let retained = sum(self.entries, reserved)?;
        count(
            sum(sum(retained, vector.len())?, witness.len())?,
            self.limit,
        )?;
        while let Some((&row, &value)) = vector.last_key_value() {
            budget.step()?;
            let Some((other, transform)) = self.pivots.get(&row) else {
                break;
            };
            let factor = self.field.negate(value);
            add_scaled(
                vector,
                other,
                factor,
                self.field,
                sum(retained, witness.len())?,
                self.limit,
                budget,
            )?;
            add_scaled(
                witness,
                transform,
                factor,
                self.field,
                sum(retained, vector.len())?,
                self.limit,
                budget,
            )?;
        }
        Ok(())
    }

    // None adds an independent column; Some is a dependency in input coordinates.
    pub(super) fn insert(
        &mut self,
        column: &SparseVector,
        index: usize,
        reserved: usize,
        budget: &mut WorkBudget<'_>,
    ) -> Result<Option<SparseVector>> {
        budget.step()?;
        count(
            sum(sum(sum(self.entries, reserved)?, column.len())?, 1)?,
            self.limit,
        )?;
        let mut vector = copy_vector(column, budget)?;
        let mut witness = BTreeMap::from([(index, 1)]);
        self.reduce(&mut vector, &mut witness, reserved, budget)?;
        if let Some((&row, &value)) = vector.last_key_value() {
            let scale = self.field.inverse(value)?;
            for value in vector.values_mut().chain(witness.values_mut()) {
                budget.step()?;
                *value = self.field.multiply(*value, scale);
            }
            self.entries = sum(sum(self.entries, vector.len())?, witness.len())?;
            self.pivots.insert(row, (vector, witness));
            Ok(None)
        } else {
            Ok(Some(witness))
        }
    }

    pub(super) fn coordinates(
        &self,
        column: &SparseVector,
        budget: &mut WorkBudget<'_>,
    ) -> Result<SparseVector> {
        count(sum(self.entries, column.len())?, self.limit)?;
        let mut vector = copy_vector(column, budget)?;
        let mut witness = SparseVector::new();
        self.reduce(&mut vector, &mut witness, 0, budget)?;
        if !vector.is_empty() {
            return Err(invariant("boundary belongs to lower Omega"));
        }
        for value in witness.values_mut() {
            budget.step()?;
            *value = self.field.negate(*value);
        }
        Ok(witness)
    }
}

pub(super) fn kernel(
    columns: &[SparseVector],
    field: PrimeField,
    limit: usize,
    budget: &mut WorkBudget<'_>,
) -> Result<(usize, Vec<SparseVector>)> {
    matrix(columns, limit)?;
    let mut solver = Echelon::new(field, limit);
    let mut basis = Vec::new();
    let mut retained = 0;
    for (index, column) in columns.iter().enumerate() {
        if let Some(witness) = solver.insert(column, index, retained, budget)? {
            retained = sum(retained, witness.len())?;
            push(&mut basis, witness)?;
        }
    }
    Ok((solver.rank(), basis))
}

pub(super) fn invariant(reason: &'static str) -> Error {
    Error::InternalInvariant { reason }
}
