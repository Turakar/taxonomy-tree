use anyhow::anyhow;
use anyhow::Result;
use pyo3::types::PyString;
use pyo3::{exceptions::PyRuntimeError, prelude::*};
use regex::Regex;
use rkyv::{Archive, Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

use crate::cache::MmapTree;
use crate::iterators::{ChildrenIterator, ParentIterator, SearchIterator};

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq, Clone)]
pub struct Entry {
    pub identifier: String,
    pub name: String,
    pub rank: Option<String>,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
}

#[derive(Clone)]
#[pyclass(frozen, name = "TaxonomyEntry")]
pub struct PyEntry {
    #[pyo3(get)]
    pub identifier: String,
    #[pyo3(get)]
    pub name: String,
    #[pyo3(get)]
    pub rank: Option<String>,
}

#[pymethods]
impl PyEntry {
    fn __repr__(slf: &Bound<'_, Self>) -> PyResult<String> {
        // This is the equivalent of `self.__class__.__name__` in Python.
        let class_name: Bound<'_, PyString> = slf.get_type().qualname()?;
        fn format_option_field(field: &Option<String>) -> String {
            match field {
                Some(value) => format!("'{}'", value),
                None => "None".to_string(),
            }
        }
        let slf_ref = slf.borrow();
        Ok(format!(
            "{}({}, {}, {})",
            class_name,
            slf_ref.identifier,
            slf_ref.name,
            format_option_field(&slf_ref.rank)
        ))
    }
}

impl From<&ArchivedEntry> for PyEntry {
    fn from(entry: &ArchivedEntry) -> Self {
        Self {
            identifier: entry.identifier.to_string(),
            name: entry.name.to_string(),
            rank: entry.rank.as_deref().map(|s| s.to_string()),
        }
    }
}

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq, Clone)]
pub struct Tree {
    pub entries: Vec<Entry>,
    pub id_to_index: BTreeMap<String, usize>,
}

#[pyclass(frozen, name = "TaxonomyTree")]
pub struct PyTree {
    tree: MmapTree,
}

impl PyTree {
    pub fn contains(&self, id: &str) -> bool {
        self.tree.as_ref().id_to_index.contains_key(id)
    }

    pub fn find_scientific_name(&self, id: &str) -> Result<&str> {
        let tree = self.tree.as_ref();
        match tree.id_to_index.get(id) {
            Some(&index) => Ok(tree.entries[index.to_native() as usize].name.as_str()),
            None => Result::Err(anyhow!("ID not found: {}", id)),
        }
    }

    pub fn search_scientific_name<'a>(
        &'a self,
        search: &str,
        lower_case: bool,
        regex: bool,
    ) -> Result<Vec<&'a ArchivedEntry>> {
        let re = if regex {
            if lower_case {
                Regex::new(&format!("(?i){}", search))
                    .map_err(|e| anyhow!("Invalid regex: {}", e))?
            } else {
                Regex::new(search).map_err(|e| anyhow!("Invalid regex: {}", e))?
            }
        } else {
            let pattern = if lower_case {
                search.to_lowercase()
            } else {
                search.to_string()
            };
            Regex::new(&format!("^{}$", regex::escape(&pattern)))
                .map_err(|e| anyhow!("Invalid regex: {}", e))?
        };
        let results = self
            .tree
            .as_ref()
            .entries
            .iter()
            .filter(|entry| {
                let name = &entry.name;
                let name_to_check = if lower_case {
                    name.to_lowercase()
                } else {
                    name.to_string()
                };
                re.is_match(&name_to_check)
            })
            .collect();
        Ok(results)
    }

    pub fn find_lineage(&self, id: &str, stop_rank: Option<&str>) -> Result<ParentIterator> {
        match self.tree.as_ref().id_to_index.get(id) {
            Some(&index) => Ok(ParentIterator::new(
                self.tree.clone(),
                index.to_native() as usize,
                stop_rank.map(|s| s.to_string()),
            )),
            None => Result::Err(anyhow!("ID not found: {}", id)),
        }
    }

    pub fn find_children(
        &self,
        id: &str,
        skip: Option<BTreeSet<String>>,
        stop_rank: Option<&str>,
    ) -> Result<ChildrenIterator> {
        let tree = self.tree.as_ref();
        let skip = match skip {
            Some(skip) => skip
                .into_iter()
                .map(|id| {
                    tree.id_to_index
                        .get(id.as_str())
                        .map(|&i| i.to_native() as usize)
                        .ok_or(anyhow!("ID in skip set not found: {}", id))
                })
                .collect::<Result<BTreeSet<_>>>()?,
            None => BTreeSet::new(),
        };
        match tree.id_to_index.get(id) {
            Some(&index) => Ok(ChildrenIterator::new(
                self.tree.clone(),
                index.to_native() as usize,
                Some(skip),
                stop_rank.map(|s| s.to_string()),
            )),
            None => Result::Err(anyhow!("ID not found: {}", id)),
        }
    }

    pub fn depth_first_search(
        &self,
        id: &str,
        stop_rank_top: Option<&str>,
        stop_rank_bottom: Option<&str>,
    ) -> Result<SearchIterator> {
        match self.tree.as_ref().id_to_index.get(id) {
            Some(&index) => Ok(SearchIterator::new(
                self.tree.clone(),
                index.to_native() as usize,
                stop_rank_top.map(|s| s.to_string()),
                stop_rank_bottom.map(|s| s.to_string()),
            )),
            None => Result::Err(anyhow!("ID not found: {}", id)),
        }
    }
}

#[pymethods]
impl PyTree {
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let mmap_tree = MmapTree::open(std::path::Path::new(path))
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(PyTree { tree: mmap_tree })
    }

    #[pyo3(name = "contains")]
    fn py_contains(&self, id: &str) -> bool {
        self.contains(id)
    }

    #[pyo3(name = "find_scientific_name")]
    fn py_find_scientific_name(&self, id: &str) -> PyResult<String> {
        self.find_scientific_name(id)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))
            .map(|s| s.to_string())
    }

    #[pyo3(name = "search_scientific_name")]
    fn py_search_scientific_name(
        &self,
        search: &str,
        lower_case: bool,
        regex: bool,
    ) -> PyResult<Vec<PyEntry>> {
        self.search_scientific_name(search, lower_case, regex)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))
            .map(|entries| entries.into_iter().map(PyEntry::from).collect())
    }

    #[pyo3(name = "find_lineage")]
    fn py_find_lineage(&self, id: String, stop_rank: Option<&str>) -> PyResult<ParentIterator> {
        self.find_lineage(&id, stop_rank)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    #[pyo3(name = "find_children")]
    fn py_find_children(
        &self,
        id: String,
        skip: Option<BTreeSet<String>>,
        stop_rank: Option<&str>,
    ) -> PyResult<ChildrenIterator> {
        self.find_children(&id, skip, stop_rank)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }

    #[pyo3(name = "depth_first_search")]
    fn py_depth_first_search(
        &self,
        id: String,
        stop_rank_top: Option<&str>,
        stop_rank_bottom: Option<&str>,
    ) -> PyResult<SearchIterator> {
        self.depth_first_search(&id, stop_rank_top, stop_rank_bottom)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))
    }
}
