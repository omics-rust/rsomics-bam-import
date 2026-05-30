use std::io::Write;

use rsomics_common::{Result, RsomicsError};

// SAMv1 §4.2 fixed record header size
pub(crate) const FIXED_HEADER: usize = 32;

/// Serialise one FASTQ record to BAM bytes in `buf` (cleared first).
///
/// Unmapped unaligned defaults: refID=-1, pos=-1 (BAM 0-based; SAM POS=0),
/// MAPQ=0, CIGAR empty, RNEXT=-1, PNEXT=-1, TLEN=0 — exactly what samtools
/// import emits (black-box verified against 1.23.1).
pub(crate) fn encode_bam_record(
    buf: &mut Vec<u8>,
    name: &[u8],
    seq: &[u8],
    qual: &[u8],
    flags: u16,
    rg_id: Option<&str>,
) -> Result<()> {
    buf.clear();

    let l_read_name: u8 = (name.len() + 1)
        .try_into()
        .map_err(|_| RsomicsError::InvalidInput("read name too long for BAM".into()))?;
    let l_seq = u32::try_from(seq.len())
        .map_err(|_| RsomicsError::InvalidInput("sequence too long for BAM l_seq field".into()))?;
    let seq_bytes = seq.len().div_ceil(2);

    buf.extend_from_slice(&(-1i32).to_le_bytes()); // refID = -1 (unmapped)
    buf.extend_from_slice(&(-1i32).to_le_bytes()); // pos = -1 (none)
    buf.push(l_read_name);
    buf.push(0); // MAPQ = 0 (unmapped, matches samtools import)
    buf.extend_from_slice(&4680u16.to_le_bytes()); // bin = 4680 (unmapped)
    buf.extend_from_slice(&0u16.to_le_bytes()); // n_cigar_op = 0
    buf.extend_from_slice(&flags.to_le_bytes());
    buf.extend_from_slice(&l_seq.to_le_bytes());
    buf.extend_from_slice(&(-1i32).to_le_bytes()); // next_refID = -1
    buf.extend_from_slice(&(-1i32).to_le_bytes()); // next_pos = -1
    buf.extend_from_slice(&0i32.to_le_bytes()); // tlen = 0

    debug_assert_eq!(buf.len(), FIXED_HEADER);

    buf.extend_from_slice(name);
    buf.push(0); // NUL-terminate read_name

    // SEQ: 4-bit encoding, two bases per byte (SAM spec Table 3)
    let full_pairs = seq.len() / 2;
    for i in 0..full_pairs {
        buf.push((nt16(seq[i * 2]) << 4) | nt16(seq[i * 2 + 1]));
    }
    if seq.len() % 2 == 1 {
        buf.push(nt16(seq[seq.len() - 1]) << 4);
    }
    debug_assert_eq!(buf.len(), FIXED_HEADER + name.len() + 1 + seq_bytes);

    // QUAL: FASTQ ASCII (Phred+33) → raw Phred
    for &q in qual {
        buf.push(q - 33);
    }

    if let Some(id) = rg_id {
        buf.push(b'R');
        buf.push(b'G');
        buf.push(b'Z');
        buf.extend_from_slice(id.as_bytes());
        buf.push(0); // NUL-terminate Z field
    }

    Ok(())
}

/// Write one BAM record block: 4-byte LE `block_size` then payload.
#[inline]
pub(crate) fn write_block<W: Write>(out: &mut W, payload: &[u8]) -> Result<()> {
    let size = u32::try_from(payload.len())
        .map_err(|e| RsomicsError::InvalidInput(format!("record too large: {e}")))?;
    out.write_all(&size.to_le_bytes())
        .map_err(RsomicsError::Io)?;
    out.write_all(payload).map_err(RsomicsError::Io)
}

/// `seq_nt16` code (BAM spec Table 3), matching htslib's `seq_nt16_table`.
#[inline]
pub(crate) fn nt16(base: u8) -> u8 {
    match base {
        b'=' => 0,
        b'A' | b'a' => 1,
        b'C' | b'c' => 2,
        b'M' | b'm' => 3,
        b'G' | b'g' => 4,
        b'R' | b'r' => 5,
        b'S' | b's' => 6,
        b'V' | b'v' => 7,
        b'T' | b't' => 8,
        b'W' | b'w' => 9,
        b'Y' | b'y' => 10,
        b'H' | b'h' => 11,
        b'K' | b'k' => 12,
        b'D' | b'd' => 13,
        b'B' | b'b' => 14,
        _ => 15, // N and anything unrecognised
    }
}

/// Strip trailing `/1` or `/2` and drop the comment (space-delimited suffix),
/// matching htslib fastq.c so both mates share the same QNAME in BAM.
#[inline]
pub(crate) fn strip_pair_suffix(id: &[u8]) -> &[u8] {
    let name = match id.iter().position(|&b| b == b' ') {
        Some(sp) => &id[..sp],
        None => id,
    };
    if name.ends_with(b"/1") || name.ends_with(b"/2") {
        &name[..name.len() - 2]
    } else {
        name
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flags::SE_FLAGS;

    #[test]
    fn strip_suffix_removes_1_and_2() {
        assert_eq!(strip_pair_suffix(b"read1/1"), b"read1");
        assert_eq!(strip_pair_suffix(b"read1/2"), b"read1");
        assert_eq!(strip_pair_suffix(b"read1"), b"read1");
        assert_eq!(strip_pair_suffix(b"read1 comment"), b"read1");
        assert_eq!(strip_pair_suffix(b"read1/1 comment"), b"read1");
    }

    #[test]
    fn nt16_standard_bases() {
        assert_eq!(nt16(b'A'), 1);
        assert_eq!(nt16(b'C'), 2);
        assert_eq!(nt16(b'G'), 4);
        assert_eq!(nt16(b'T'), 8);
        assert_eq!(nt16(b'N'), 15);
        assert_eq!(nt16(b'n'), 15);
        assert_eq!(nt16(b'a'), 1);
    }

    #[test]
    fn bam_record_layout_single() {
        let mut buf = Vec::new();
        encode_bam_record(&mut buf, b"r1", b"ACGT", b"IIII", SE_FLAGS, None).unwrap();
        assert_eq!(&buf[0..4], &(-1i32).to_le_bytes()); // refID = -1
        assert_eq!(&buf[4..8], &(-1i32).to_le_bytes()); // pos = -1
        assert_eq!(buf[8], 3); // l_read_name = len("r1") + 1
        assert_eq!(buf[9], 0); // MAPQ = 0
        assert_eq!(
            u16::from_le_bytes(buf[14..16].try_into().unwrap()),
            SE_FLAGS
        );
        assert_eq!(u32::from_le_bytes(buf[16..20].try_into().unwrap()), 4); // l_seq
        assert_eq!(&buf[32..35], b"r1\0");
        assert_eq!(buf[35], 0x12); // A=1,C=2 → 0x12
        assert_eq!(buf[36], 0x48); // G=4,T=8 → 0x48
        assert_eq!(&buf[37..41], &[40u8, 40, 40, 40]); // 'I' - 33 = 40
    }

    #[test]
    fn bam_record_rg_aux() {
        let mut buf = Vec::new();
        encode_bam_record(&mut buf, b"r1", b"ACGT", b"IIII", SE_FLAGS, Some("lib1")).unwrap();
        let expected_aux = b"RGZlib1\0";
        let len = buf.len();
        assert_eq!(&buf[len - expected_aux.len()..], expected_aux);
    }
}
