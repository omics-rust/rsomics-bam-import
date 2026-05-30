use bstr::BString;
use noodles::sam;
use noodles::sam::header::record::value::{
    Map,
    map::{Header as MapHeader, ReadGroup},
};
use rsomics_common::{Result, RsomicsError};

pub(crate) fn build_header(rg_line: Option<&str>) -> Result<sam::Header> {
    use noodles::sam::header::record::value::map::header::tag;

    let mut hd =
        Map::<MapHeader>::new(noodles::sam::header::record::value::map::header::Version::new(1, 6));
    // SO:unsorted and GO:query — stored as other_fields (noodles has no enum for these)
    hd.other_fields_mut()
        .insert(tag::SORT_ORDER, BString::from("unsorted"));
    hd.other_fields_mut()
        .insert(tag::GROUP_ORDER, BString::from("query"));

    let mut builder = sam::Header::builder().set_header(hd);

    if let Some(line) = rg_line {
        let (rg_id, rg_map) = parse_rg_line(line)?;
        builder = builder.add_read_group(rg_id, rg_map);
    }

    Ok(builder.build())
}

/// Parse "@RG\tID:foo\tSM:bar" into (id, Map<ReadGroup>). ID is returned
/// separately because it is both the map key and the RG:Z aux value.
pub(crate) fn parse_rg_line(line: &str) -> Result<(BString, Map<ReadGroup>)> {
    let body = line.strip_prefix("@RG\t").unwrap_or(line);
    let mut id: Option<String> = None;
    let mut map = Map::<ReadGroup>::default();
    for field in body.split('\t') {
        let (tag, value) = field.split_once(':').ok_or_else(|| {
            RsomicsError::InvalidInput(format!("malformed @RG field (no ':'): {field:?}"))
        })?;
        if tag == "ID" {
            id = Some(value.to_string());
        } else {
            use noodles::sam::header::record::value::map::read_group::tag::Standard;
            use noodles::sam::header::record::value::map::tag::Other;
            let tag_bytes = tag.as_bytes();
            if tag_bytes.len() != 2 {
                return Err(RsomicsError::InvalidInput(format!(
                    "RG tag must be 2 chars: {tag:?}"
                )));
            }
            let buf: [u8; 2] = [tag_bytes[0], tag_bytes[1]];
            let other: Other<Standard> = Other::try_from(buf).map_err(|_| {
                RsomicsError::InvalidInput(format!("RG tag is a reserved standard tag: {tag:?}"))
            })?;
            map.other_fields_mut().insert(other, BString::from(value));
        }
    }
    let rg_id = id.ok_or_else(|| RsomicsError::InvalidInput("@RG line missing ID field".into()))?;
    Ok((BString::from(rg_id), map))
}

/// Extract the ID value from a raw @RG line for embedding in RG:Z aux tags.
pub(crate) fn extract_rg_id(rg_line: &str) -> Option<&str> {
    let body = rg_line.strip_prefix("@RG\t").unwrap_or(rg_line);
    for field in body.split('\t') {
        if let Some(v) = field.strip_prefix("ID:") {
            return Some(v);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rg_id_extracted() {
        assert_eq!(extract_rg_id("@RG\tID:lib1\tSM:s1"), Some("lib1"));
        assert_eq!(extract_rg_id("ID:lib1\tSM:s1"), Some("lib1"));
        assert_eq!(extract_rg_id("@RG\tSM:s1"), None);
    }
}
