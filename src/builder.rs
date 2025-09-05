use crate::cache::MmapTree;
use crate::tree::{Entry, Tree};
use anyhow::{anyhow, Result};
use std::collections::BTreeMap;
use std::sync::Mutex;

use pyo3::{exceptions::PyRuntimeError, prelude::*};

struct TreeBuilderEntry {
    identifier: String,
    name: String,
    rank: Option<String>,
    parent: Option<String>,
    parent_idx: Option<usize>,
    children_idx: Vec<usize>,
}

pub struct TreeBuilder {
    has_root: bool,
    entries: Vec<TreeBuilderEntry>,
    id_to_index: BTreeMap<String, usize>,
}

impl TreeBuilder {
    pub fn new() -> Self {
        Self {
            has_root: false,
            entries: Vec::new(),
            id_to_index: BTreeMap::new(),
        }
    }

    pub fn insert(
        &mut self,
        identifier: String,
        name: String,
        rank: Option<String>,
        parent: Option<String>,
        overwrite: bool,
    ) -> Result<()> {
        if parent.is_none() {
            if self.has_root {
                return Err(anyhow!("Multiple root entries found"));
            }
            self.has_root = true;
        }
        if self.id_to_index.contains_key(&identifier) {
            if overwrite {
                self.id_to_index.remove(&identifier);
            } else {
                return Err(anyhow!("Duplicate identifier: {}", identifier));
            }
        }
        let index = self.entries.len();
        self.id_to_index.insert(identifier.clone(), index);
        self.entries.push(TreeBuilderEntry {
            identifier,
            name,
            rank,
            parent,
            parent_idx: None,
            children_idx: Vec::new(),
        });
        Ok(())
    }

    pub fn build(mut self) -> Result<Tree> {
        if !self.has_root {
            return Err(anyhow!("No root entry found"));
        }
        // Create parent and children links
        for i in 0..self.entries.len() {
            if let Some(parent) = &self.entries[i].parent {
                if let Some(&parent_idx) = self.id_to_index.get(parent) {
                    self.entries[i].parent_idx = Some(parent_idx);
                    self.entries[parent_idx].children_idx.push(i);
                } else {
                    return Err(anyhow!("Parent identifier not found: {}", parent));
                }
            }
        }
        // Convert to final Entry structure
        let entries = self
            .entries
            .into_iter()
            .map(|e| Entry {
                identifier: e.identifier,
                name: e.name,
                rank: e.rank,
                parent: e.parent_idx,
                children: e.children_idx,
            })
            .collect();
        Ok(Tree {
            entries,
            id_to_index: self.id_to_index,
        })
    }
}

#[pyclass(name = "TaxonomyTreeBuilder", frozen)]
pub struct PyTreeBuilder {
    builder: Mutex<Option<TreeBuilder>>,
}

#[pymethods]
impl PyTreeBuilder {
    #[new]
    fn new() -> Self {
        PyTreeBuilder {
            builder: Mutex::new(Some(TreeBuilder::new())),
        }
    }

    fn insert(
        &self,
        identifier: String,
        name: String,
        rank: Option<String>,
        parent: Option<String>,
        overwrite: bool,
    ) -> PyResult<()> {
        let mut builder = self.builder.lock().unwrap();
        if let Some(builder) = &mut *builder {
            builder
                .insert(identifier, name, rank, parent, overwrite)
                .map_err(|e| PyRuntimeError::new_err(e.to_string()))
        } else {
            Err(PyRuntimeError::new_err(
                "Builder has already been finalized.",
            ))
        }
    }

    fn write(&self, path: &str) -> PyResult<()> {
        let builder = self
            .builder
            .lock()
            .unwrap()
            .take()
            .ok_or_else(|| PyRuntimeError::new_err("Builder has already been finalized."))?;
        let tree = builder
            .build()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        MmapTree::write(std::path::Path::new(path), &tree)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(())
    }
}
