use std::sync::Arc;
use std::{fs::File, path::Path};

use memmap2::Mmap;

use crate::tree::{ArchivedTree, Tree};
use anyhow::Context;
use anyhow::Result;
use rkyv::{self, ser::writer::IoWriter};

#[derive(Clone)]
pub struct MmapTree {
    pub mmap: Arc<Mmap>,
}

impl MmapTree {
    pub fn as_ref(&self) -> &ArchivedTree {
        let bytes = self.mmap.as_ref().as_ref();
        unsafe { rkyv::access_unchecked::<ArchivedTree>(bytes) }
    }

    pub fn write(path: &Path, tree: &Tree) -> Result<()> {
        // Write tree to file directly without buffering the serialized data in memory
        let file = File::create(path)?;
        let mut writer = IoWriter::new(file);
        rkyv::api::high::to_bytes_in::<_, rkyv::rancor::Error>(tree, &mut writer)?;
        Ok(())
    }

    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        // Check that the data is valid
        rkyv::access::<ArchivedTree, rkyv::rancor::Error>(&mmap)
            .context("rkyv validation failed")?;
        Ok(MmapTree {
            mmap: Arc::new(mmap),
        })
    }
}
