use std::io::Write;

use noodles::sam;
use noodles::sam::alignment::io::Write as AlnWrite;
use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::FastqSource;

use crate::encode::{encode_bam_record, strip_pair_suffix, write_block};
use crate::flags::{PE_R1_FLAGS, PE_R2_FLAGS, SE_FLAGS};
use crate::{ImportMode, ImportOpts};

pub(crate) fn dispatch_raw<W: Write>(
    opts: &ImportOpts,
    out: &mut W,
    rg_id: Option<&str>,
) -> Result<u64> {
    match &opts.mode {
        ImportMode::Single(path) => {
            let src = rsomics_seqio::open_fastq(path)?;
            write_single_raw(out, src, rg_id)
        }
        ImportMode::Interleaved(path) => {
            let src = rsomics_seqio::open_fastq(path)?;
            write_interleaved_raw(out, src, rg_id)
        }
        ImportMode::Paired(r1, r2) => {
            let src1 = rsomics_seqio::open_fastq(r1)?;
            let src2 = rsomics_seqio::open_fastq(r2)?;
            write_paired_raw(out, src1, src2, rg_id)
        }
    }
}

pub(crate) fn dispatch_sam(
    opts: &ImportOpts,
    header: &sam::Header,
    rg_id: Option<&str>,
) -> Result<u64> {
    let stdout = std::io::stdout();
    let mut buf = std::io::BufWriter::with_capacity(256 * 1024, stdout.lock());
    let mut writer = sam::io::Writer::new(&mut buf);
    writer.write_header(header).map_err(RsomicsError::Io)?;
    let count = match &opts.mode {
        ImportMode::Single(path) => {
            let src = rsomics_seqio::open_fastq(path)?;
            write_single_sam(&mut writer, header, src, rg_id)?
        }
        ImportMode::Interleaved(path) => {
            let src = rsomics_seqio::open_fastq(path)?;
            write_interleaved_sam(&mut writer, header, src, rg_id)?
        }
        ImportMode::Paired(r1, r2) => {
            let src1 = rsomics_seqio::open_fastq(r1)?;
            let src2 = rsomics_seqio::open_fastq(r2)?;
            write_paired_sam(&mut writer, header, src1, src2, rg_id)?
        }
    };
    buf.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}

fn write_single_raw<W: Write>(out: &mut W, src: FastqSource, rg_id: Option<&str>) -> Result<u64> {
    let mut buf = Vec::with_capacity(512);
    let mut count = 0u64;
    for result in src {
        let rec = result?;
        let name = strip_pair_suffix(&rec.id);
        encode_bam_record(&mut buf, name, &rec.seq, &rec.qual, SE_FLAGS, rg_id)?;
        write_block(out, &buf)?;
        count += 1;
    }
    Ok(count)
}

fn write_interleaved_raw<W: Write>(
    out: &mut W,
    src: FastqSource,
    rg_id: Option<&str>,
) -> Result<u64> {
    let mut src = src;
    let mut buf = Vec::with_capacity(512);
    let mut count = 0u64;
    while let Some(r1_res) = src.next() {
        let r1 = r1_res?;
        let r2 = src.next().ok_or_else(|| {
            RsomicsError::InvalidInput("interleaved FASTQ has odd number of records".into())
        })??;
        let name = strip_pair_suffix(&r1.id);
        encode_bam_record(&mut buf, name, &r1.seq, &r1.qual, PE_R1_FLAGS, rg_id)?;
        write_block(out, &buf)?;
        encode_bam_record(&mut buf, name, &r2.seq, &r2.qual, PE_R2_FLAGS, rg_id)?;
        write_block(out, &buf)?;
        count += 2;
    }
    Ok(count)
}

fn write_paired_raw<W: Write>(
    out: &mut W,
    src1: FastqSource,
    src2: FastqSource,
    rg_id: Option<&str>,
) -> Result<u64> {
    let mut src2 = src2.peekable();
    let mut buf = Vec::with_capacity(512);
    let mut count = 0u64;
    for r1_res in src1 {
        let r1 = r1_res?;
        let r2 = src2.next().ok_or_else(|| {
            RsomicsError::InvalidInput("read-2 file has fewer records than read-1".into())
        })??;
        let name = strip_pair_suffix(&r1.id);
        encode_bam_record(&mut buf, name, &r1.seq, &r1.qual, PE_R1_FLAGS, rg_id)?;
        write_block(out, &buf)?;
        encode_bam_record(&mut buf, name, &r2.seq, &r2.qual, PE_R2_FLAGS, rg_id)?;
        write_block(out, &buf)?;
        count += 2;
    }
    if src2.next().is_some() {
        return Err(RsomicsError::InvalidInput(
            "read-2 file has more records than read-1".into(),
        ));
    }
    Ok(count)
}

fn build_sam_record(
    name: &[u8],
    seq: &[u8],
    qual: &[u8],
    flags: u16,
    rg_id: Option<&str>,
) -> sam::alignment::RecordBuf {
    use sam::alignment::record::data::field::Tag;
    use sam::alignment::record::{Flags, MappingQuality};
    use sam::alignment::record_buf::data::field::Value;
    use sam::alignment::record_buf::{Data, QualityScores, Sequence};

    let sequence = Sequence::from(seq.to_vec());
    // FASTQ qual is ASCII (Phred+33); SAM stores raw Phred values
    let qual_phred: Vec<u8> = qual.iter().map(|&q| q - 33).collect();
    let quality_scores = QualityScores::from(qual_phred);

    let mut data = Data::default();
    if let Some(id) = rg_id {
        data.insert(Tag::READ_GROUP, Value::String(id.into()));
    }

    // MAPQ=0 for all unaligned reads — MappingQuality::new(255) is None (missing);
    // samtools import uses 0, not the missing sentinel.
    let mapq = MappingQuality::new(0).expect("0 is a valid MAPQ");

    let mut builder = sam::alignment::RecordBuf::builder()
        .set_flags(Flags::from(flags))
        .set_mapping_quality(mapq)
        .set_sequence(sequence)
        .set_quality_scores(quality_scores)
        .set_data(data);

    if !name.is_empty() {
        builder = builder.set_name(bstr::BString::from(name));
    }

    builder.build()
}

fn write_single_sam<W: Write>(
    writer: &mut sam::io::Writer<W>,
    header: &sam::Header,
    src: FastqSource,
    rg_id: Option<&str>,
) -> Result<u64> {
    let mut count = 0u64;
    for result in src {
        let rec = result?;
        let name = strip_pair_suffix(&rec.id);
        let record = build_sam_record(name, &rec.seq, &rec.qual, SE_FLAGS, rg_id);
        writer
            .write_alignment_record(header, &record)
            .map_err(RsomicsError::Io)?;
        count += 1;
    }
    Ok(count)
}

fn write_interleaved_sam<W: Write>(
    writer: &mut sam::io::Writer<W>,
    header: &sam::Header,
    src: FastqSource,
    rg_id: Option<&str>,
) -> Result<u64> {
    let mut src = src;
    let mut count = 0u64;
    while let Some(r1_res) = src.next() {
        let r1 = r1_res?;
        let r2 = src.next().ok_or_else(|| {
            RsomicsError::InvalidInput("interleaved FASTQ has odd number of records".into())
        })??;
        let name = strip_pair_suffix(&r1.id);
        let rec1 = build_sam_record(name, &r1.seq, &r1.qual, PE_R1_FLAGS, rg_id);
        let rec2 = build_sam_record(name, &r2.seq, &r2.qual, PE_R2_FLAGS, rg_id);
        writer
            .write_alignment_record(header, &rec1)
            .map_err(RsomicsError::Io)?;
        writer
            .write_alignment_record(header, &rec2)
            .map_err(RsomicsError::Io)?;
        count += 2;
    }
    Ok(count)
}

fn write_paired_sam<W: Write>(
    writer: &mut sam::io::Writer<W>,
    header: &sam::Header,
    src1: FastqSource,
    src2: FastqSource,
    rg_id: Option<&str>,
) -> Result<u64> {
    let mut src2 = src2.peekable();
    let mut count = 0u64;
    for r1_res in src1 {
        let r1 = r1_res?;
        let r2 = src2.next().ok_or_else(|| {
            RsomicsError::InvalidInput("read-2 file has fewer records than read-1".into())
        })??;
        let name = strip_pair_suffix(&r1.id);
        let rec1 = build_sam_record(name, &r1.seq, &r1.qual, PE_R1_FLAGS, rg_id);
        let rec2 = build_sam_record(name, &r2.seq, &r2.qual, PE_R2_FLAGS, rg_id);
        writer
            .write_alignment_record(header, &rec1)
            .map_err(RsomicsError::Io)?;
        writer
            .write_alignment_record(header, &rec2)
            .map_err(RsomicsError::Io)?;
        count += 2;
    }
    if src2.next().is_some() {
        return Err(RsomicsError::InvalidInput(
            "read-2 file has more records than read-1".into(),
        ));
    }
    Ok(count)
}
