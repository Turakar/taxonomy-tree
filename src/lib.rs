mod builder;
mod cache;
mod iterators;
mod tree;

use pyo3::prelude::*;

#[pymodule]
fn taxonomy_tree(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<crate::tree::PyEntry>()?;
    m.add_class::<crate::tree::PyTree>()?;
    m.add_class::<crate::builder::PyTreeBuilder>()?;
    m.add_class::<crate::iterators::ChildrenIterator>()?;
    m.add_class::<crate::iterators::ParentIterator>()?;
    m.add_class::<crate::iterators::SearchIterator>()?;
    Ok(())
}
