use crate::cache::MmapTree;
use crate::tree::PyEntry;
use pyo3::prelude::*;
use std::collections::BTreeSet;
use std::sync::Mutex;

#[pyclass(frozen)]
pub struct ChildrenIterator {
    tree: MmapTree,
    stack: Mutex<Vec<(usize, Vec<PyEntry>)>>,
    skip: BTreeSet<usize>,
    stop_rank: Option<String>,
}

impl ChildrenIterator {
    pub fn new(
        tree: MmapTree,
        start: usize,
        skip: Option<BTreeSet<usize>>,
        stop_rank: Option<String>,
    ) -> Self {
        Self {
            tree,
            stack: Mutex::new(vec![(start, vec![])]),
            skip: skip.unwrap_or_default(),
            stop_rank,
        }
    }
}

#[pymethods]
impl ChildrenIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self) -> Option<Vec<PyEntry>> {
        let tree = self.tree.as_ref();
        let mut stack = self.stack.lock().unwrap();
        while let Some((index, lineage)) = stack.pop() {
            if self.skip.contains(&index) {
                continue;
            }
            let entry = &tree.entries[index];
            let lineage = std::iter::once(PyEntry::from(entry))
                .chain(lineage.into_iter())
                .collect::<Vec<_>>();
            let mut should_continue = true;
            if let Some(stop_rank) = &self.stop_rank {
                if let Some(rank) = &entry.rank.as_ref() {
                    if rank.as_ref() == stop_rank {
                        should_continue = false;
                    }
                }
            }
            if should_continue {
                stack.extend(
                    entry
                        .children
                        .iter()
                        .rev()
                        .map(|i| (i.to_native() as usize, lineage.clone())),
                );
            }
            return Some(lineage);
        }
        None
    }
}

#[pyclass(frozen)]
pub struct ParentIterator {
    tree: MmapTree,
    current: Mutex<Option<usize>>,
    stop_rank: Option<String>,
}

impl ParentIterator {
    pub fn new(tree: MmapTree, start: usize, stop_rank: Option<String>) -> Self {
        Self {
            tree,
            current: Mutex::new(Some(start)),
            stop_rank,
        }
    }
}

#[pymethods]
impl ParentIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self) -> Option<PyEntry> {
        let tree = self.tree.as_ref();
        let mut current = self.current.lock().unwrap();
        if let Some(index) = *current {
            let entry = &tree.entries[index];
            if let Some(stop_rank) = &self.stop_rank {
                if let Some(rank) = &entry.rank.as_ref() {
                    if rank.as_ref() == stop_rank {
                        *current = None;
                        return Some(PyEntry::from(entry));
                    }
                }
            }
            *current = entry.parent.as_ref().map(|p| p.to_native() as usize);
            Some(PyEntry::from(entry))
        } else {
            None
        }
    }
}

#[pyclass(frozen)]
pub struct SearchIterator {
    tree: MmapTree,
    state: Mutex<(
        Option<usize>,            // current index
        Option<ChildrenIterator>, // current children iterator
    )>,
    stop_rank_top: Option<String>,
    stop_rank_bottom: Option<String>,
}

impl SearchIterator {
    pub fn new(
        tree: MmapTree,
        start: usize,
        stop_rank_top: Option<String>,
        stop_rank_bottom: Option<String>,
    ) -> Self {
        Self {
            tree: tree.clone(),
            state: Mutex::new((
                Some(start),
                Some(ChildrenIterator::new(
                    tree,
                    start,
                    None,
                    stop_rank_bottom.clone(),
                )),
            )),
            stop_rank_top,
            stop_rank_bottom,
        }
    }
}

#[pymethods]
impl SearchIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self) -> Option<PyEntry> {
        let tree = self.tree.as_ref();
        let (current, current_children_iterator) = &mut *self.state.lock().unwrap();
        loop {
            if let Some(children_iterator) = current_children_iterator {
                if let Some(entry) = children_iterator.__next__() {
                    return Some(entry.into_iter().next().unwrap());
                } else {
                    *current_children_iterator = None;
                }
            } else if let Some(old) = *current {
                let entry = &tree.entries[old];
                if let Some(stop_rank) = &self.stop_rank_top {
                    if let Some(rank) = &entry.rank.as_ref() {
                        if rank.as_ref() == stop_rank {
                            *current = None;
                        }
                    }
                } else if let Some(parent) = entry.parent.as_ref() {
                    let parent = parent.to_native() as usize;
                    *current = Some(parent);
                    *current_children_iterator = Some(ChildrenIterator::new(
                        self.tree.clone(),
                        parent,
                        Some(vec![old].into_iter().collect()),
                        self.stop_rank_bottom.clone(),
                    ));
                } else {
                    *current = None;
                }
            } else {
                return None;
            }
        }
    }
}
