/// Produces the BGZF-compressed files that samtools/tabix need to generate
/// the committed reference index files used by integration tests.
///
/// Usage:
///   cargo run --example gen_references -- <fasta_in> <fasta_bgz_out> <gff_in> <gff_bgz_out>
///
/// Input files may be plain or gzip-compressed (detected by magic bytes).
///
/// After running this, use tests/generate_references.sh to invoke samtools/tabix
/// on the outputs and commit the resulting .fai, .gzi, and .csi files.

use std::fs;
use std::io::{Cursor, Read};

use flate2::read::MultiGzDecoder;
use mgnify_wasm::htslib::bgzf_compress;
use mgnify_wasm::gff_preprocess;

/// Read a file, transparently decompressing if it begins with the gzip magic bytes.
fn read_maybe_gz(path: &str) -> Vec<u8> {
    let raw = fs::read(path)
        .unwrap_or_else(|e| panic!("cannot read {}: {}", path, e));
    if raw.starts_with(&[0x1F, 0x8B]) {
        let mut decoder = MultiGzDecoder::new(raw.as_slice());
        let mut out = Vec::new();
        decoder.read_to_end(&mut out)
            .unwrap_or_else(|e| panic!("gzip decode failed for {}: {}", path, e));
        out
    } else {
        raw
    }
}

/// Read a (possibly gzip-compressed) file as a UTF-8 string.
fn read_text_maybe_gz(path: &str) -> String {
    let bytes = read_maybe_gz(path);
    String::from_utf8(bytes)
        .unwrap_or_else(|e| panic!("non-UTF8 content in {}: {}", path, e))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 5 {
        eprintln!("Usage: gen_references <fasta_in> <fasta_bgz_out> <gff_in> <gff_bgz_out>");
        std::process::exit(1);
    }

    let (fasta_in, fasta_out, gff_in, gff_out) = (&args[1], &args[2], &args[3], &args[4]);

    // FASTA: read (decompressing if needed) and BGZF-compress
    let fa_bytes = read_maybe_gz(fasta_in);
    let mut fasta_bgz = Vec::new();
    bgzf_compress(Cursor::new(&fa_bytes), &mut fasta_bgz).expect("FASTA bgzf_compress failed");
    fs::write(fasta_out, &fasta_bgz)
        .unwrap_or_else(|e| panic!("cannot write {}: {}", fasta_out, e));
    eprintln!("Wrote {} bytes → {}", fasta_bgz.len(), fasta_out);

    // GFF3: read (decompressing if needed), preprocess (sort + strip ##FASTA), then BGZF-compress
    let gff_raw = read_text_maybe_gz(gff_in);
    let gff_preprocessed = gff_preprocess(&gff_raw);
    let mut gff_bgz = Vec::new();
    bgzf_compress(Cursor::new(gff_preprocessed.as_bytes()), &mut gff_bgz)
        .expect("GFF bgzf_compress failed");
    fs::write(gff_out, &gff_bgz)
        .unwrap_or_else(|e| panic!("cannot write {}: {}", gff_out, e));
    eprintln!("Wrote {} bytes → {}", gff_bgz.len(), gff_out);
}
