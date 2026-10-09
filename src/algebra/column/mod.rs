//! Ordered sparse coefficient columns, independent of filtration and execution.
use super::PrimeField;
use crate::Result;
use std::collections::BTreeMap;
use std::collections::btree_map::Entry;

#[derive(Clone, Debug)]
pub(crate) struct Column<K>(BTreeMap<K, u32>);
impl<K: Ord + Clone> Column<K> {
    pub(crate) fn new() -> Self {
        Self(BTreeMap::new())
    }
    pub(crate) fn unit(key: K) -> Self {
        Self(BTreeMap::from([(key, 1)]))
    }
    pub(crate) fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    pub(crate) fn entries(&self) -> impl DoubleEndedIterator<Item = (&K, &u32)> {
        self.0.iter()
    }
    pub(crate) fn first(&self) -> Option<(&K, &u32)> {
        self.0.first_key_value()
    }
    pub(crate) fn last(&self) -> Option<(&K, &u32)> {
        self.0.last_key_value()
    }
    pub(crate) fn get(&self, key: &K) -> u32 {
        self.0.get(key).copied().unwrap_or(0)
    }
    pub(crate) fn add_term(&mut self, key: K, coefficient: u32, field: PrimeField) {
        match self.0.entry(key) {
            Entry::Occupied(mut entry) => {
                let value = field.add(*entry.get(), coefficient);
                if value == 0 {
                    entry.remove();
                } else {
                    *entry.get_mut() = value;
                }
            }
            Entry::Vacant(entry) => {
                let value = field.add(0, coefficient);
                if value != 0 {
                    entry.insert(value);
                }
            }
        }
    }
    pub(crate) fn add_scaled(
        &mut self,
        other: &Self,
        factor: u32,
        field: PrimeField,
        checkpoint: &mut impl FnMut() -> Result<()>,
    ) -> Result<()> {
        for (key, &value) in other.entries() {
            checkpoint()?;
            let value = if factor == 1 {
                value
            } else {
                field.multiply(factor, value)
            };
            self.add_term(key.clone(), value, field);
        }
        Ok(())
    }
    pub(crate) fn scale(
        &mut self,
        factor: u32,
        field: PrimeField,
        checkpoint: &mut impl FnMut() -> Result<()>,
    ) -> Result<()> {
        for value in self.0.values_mut() {
            checkpoint()?;
            *value = field.multiply(*value, factor);
        }
        self.0.retain(|_, v| *v != 0);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn updates_match_integer_modular_arithmetic_with_noncanonical_inputs() {
        for prime in [2, 3, 7] {
            let field = PrimeField::new(prime).unwrap();
            for a in 0..2 * prime {
                for b in 0..2 * prime {
                    let mut column = Column::new();
                    column.add_term(0, a, field);
                    column.add_term(0, b, field);
                    assert_eq!(column.get(&0), (a + b) % prime);
                    assert_eq!(column.is_empty(), (a + b) % prime == 0);
                    for factor in 0..2 * prime {
                        let mut left = Column::new();
                        left.add_term(0, a, field);
                        let mut right = Column::new();
                        right.add_term(0, b, field);
                        right.add_term(1, b, field);
                        left.add_scaled(&right, factor, field, &mut || Ok(()))
                            .unwrap();
                        assert_eq!(left.get(&0), (a + factor * b) % prime);
                        assert_eq!(left.get(&1), factor * b % prime);
                        assert!(left.entries().all(|(_, &v)| v > 0 && v < prime));
                    }
                }
            }
        }
    }
}
