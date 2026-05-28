use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bam_import(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bam-import");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let se = manifest.join("tests/golden/se.fastq");
    c.bench_function("rsomics-bam-import golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args(["--single", se.to_str().unwrap(), "-o", "/dev/null"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_bam_import);
criterion_main!(benches);
