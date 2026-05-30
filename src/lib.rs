//! FASTQ → unaligned BAM conversion, porting `samtools import` (MIT).
//!
//! The hot path avoids noodles' record model entirely: each FASTQ record is
//! serialised directly to a BAM block in memory and flushed through rsomics-bamio's
//! parallel BGZF writer. This skips the decode→struct→encode round-trip that
//! noodles' `write_alignment_record` performs and means the only non-trivial work
//! per read is the 4-bit SEQ packing and one copy of the quality bytes.
//!
//! Semantics match samtools import (`bam_import.c`, MIT): FLAG, RNEXT/PNEXT/TLEN,
//! MAPQ, RNAME, CIGAR, and read-name /1 /2 stripping are all faithfully ported
//! from reading the 1.23.1 source.

mod encode;
mod flags;
mod header;
mod writer;

use std::num::NonZero;
use std::path::PathBuf;

use noodles::bam;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

/// Which FASTQ input mode the user selected.
#[derive(Debug, Clone)]
pub enum ImportMode {
    /// Single-ended reads: no pairing flags set (samtools `-0`).
    Single(PathBuf),
    /// Interleaved: alternating R1/R2 in one file (samtools `-s`).
    Interleaved(PathBuf),
    /// Paired-end from two separate files (samtools `-1`/`-2`).
    Paired(PathBuf, PathBuf),
}

/// Options passed through from the CLI layer.
pub struct ImportOpts {
    pub mode: ImportMode,
    pub output: Option<PathBuf>,
    pub uncompressed: bool,
    pub rg_line: Option<String>,
    pub workers: NonZero<usize>,
}

#[derive(Debug, Default, Clone, Serialize)]
pub struct ImportStats {
    pub records_written: u64,
}

/// Convert FASTQ reads to unaligned BAM.
pub fn import(opts: &ImportOpts) -> Result<ImportStats> {
    let hdr = header::build_header(opts.rg_line.as_deref())?;
    let rg_id: Option<String> = opts
        .rg_line
        .as_deref()
        .and_then(header::extract_rg_id)
        .map(str::to_owned);

    if let Some(out_path) = &opts.output {
        let records_written = if opts.uncompressed {
            let file = std::fs::File::create(out_path).map_err(|e| {
                RsomicsError::InvalidInput(format!("creating {}: {e}", out_path.display()))
            })?;
            let mut w = bam::io::Writer::new(file);
            w.write_header(&hdr).map_err(RsomicsError::Io)?;
            writer::dispatch_raw(opts, w.get_mut(), rg_id.as_deref())?
        } else {
            let mut w = rsomics_bamio::create_with_workers(out_path, opts.workers)?;
            w.write_header(&hdr).map_err(RsomicsError::Io)?;
            writer::dispatch_raw(opts, w.get_mut(), rg_id.as_deref())?
        };
        Ok(ImportStats { records_written })
    } else {
        // stdout: SAM text
        let records_written = writer::dispatch_sam(opts, &hdr, rg_id.as_deref())?;
        Ok(ImportStats { records_written })
    }
}
