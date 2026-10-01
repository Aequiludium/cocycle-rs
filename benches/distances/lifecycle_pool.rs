//! Capacity-only experiment state, compiled solely by the lifecycle worker.
use crate::{Error, Result};
use std::{
    any::Any,
    cell::RefCell,
    collections::HashMap,
    ops::{Deref, DerefMut},
    time::Instant,
};

type Entry = (u64, usize, Box<dyn Any>);
#[derive(Default)]
struct State {
    retain: bool,
    detailed: bool,
    generation: u64,
    buffers: HashMap<&'static str, Entry>,
    phases_ns: [u64; 10],
    hits: usize,
    growths: usize,
    limit: usize,
}
thread_local! { static STATE: RefCell<State> = RefCell::new(State::default()); }

pub fn configure(retain: bool, detailed: bool, limit: usize) {
    STATE.with(|s| {
        let mut s = s.borrow_mut();
        s.retain = retain;
        s.detailed = detailed;
        s.limit = limit;
    });
}
pub fn begin() {
    STATE.with(|s| s.borrow_mut().generation += 1);
}
pub fn clear() {
    STATE.with(|s| s.borrow_mut().buffers.clear());
}
pub fn reset_stats() {
    STATE.with(|s| {
        let mut s = s.borrow_mut();
        s.phases_ns = [0; 10];
        s.hits = 0;
        s.growths = 0;
    });
}
pub fn stats() -> ([u64; 10], usize, usize, usize) {
    STATE.with(|s| {
        let s = s.borrow();
        (
            s.phases_ns,
            s.hits,
            s.growths,
            s.buffers.values().map(|x| x.1).sum(),
        )
    })
}
pub struct Clock(Option<(Instant, usize)>);
impl Clock {
    pub fn new(phase: usize) -> Self {
        Self(STATE.with(|s| s.borrow().detailed.then(|| (Instant::now(), phase))))
    }
}
impl Drop for Clock {
    fn drop(&mut self) {
        if let Some((start, phase)) = self.0 {
            STATE.with(|s| s.borrow_mut().phases_ns[phase] += start.elapsed().as_nanos() as u64);
        }
    }
}
pub fn empty<T: 'static>(key: &'static str) -> Vec<T> {
    STATE.with(|s| {
        let mut s = s.borrow_mut();
        // Never warm-start semantic state or reuse an allocation from this call.
        let usable = s.retain && s.buffers.get(key).is_some_and(|x| x.0 != s.generation);
        if usable {
            let (_, _, value) = s.buffers.remove(key).unwrap();
            s.hits += 1;
            *value
                .downcast::<Vec<T>>()
                .expect("fixed typed benchmark slot")
        } else {
            Vec::new()
        }
    })
}
pub fn filled<T: Clone + 'static>(key: &'static str, len: usize, value: T) -> Result<Vec<T>> {
    let _clock = Clock::new(8);
    let mut result = empty(key);
    let before = result.capacity();
    result.clear();
    result
        .try_reserve_exact(len)
        .map_err(|_| Error::AllocationFailed {
            context: "lifecycle scratch",
        })?;
    result.resize(len, value);
    if result.capacity() != before {
        STATE.with(|s| s.borrow_mut().growths += 1);
    }
    Ok(result)
}
pub fn put<T: 'static>(key: &'static str, mut value: Vec<T>) {
    value.clear();
    let bytes = value.capacity() * std::mem::size_of::<T>();
    STATE.with(|s| {
        let mut s = s.borrow_mut();
        if s.retain && bytes <= s.limit {
            let generation = s.generation;
            s.buffers.insert(key, (generation, bytes, Box::new(value)));
        }
    });
}
pub struct Vector<T: 'static>(&'static str, Vec<T>);
impl<T: 'static> Deref for Vector<T> {
    type Target = Vec<T>;
    fn deref(&self) -> &Vec<T> {
        &self.1
    }
}
impl<T: 'static> DerefMut for Vector<T> {
    fn deref_mut(&mut self) -> &mut Vec<T> {
        &mut self.1
    }
}
impl<T: 'static> Drop for Vector<T> {
    fn drop(&mut self) {
        put(self.0, std::mem::take(&mut self.1));
    }
}
pub fn vector<T: Clone + 'static>(key: &'static str, len: usize, value: T) -> Result<Vector<T>> {
    Ok(Vector(key, filled(key, len, value)?))
}
