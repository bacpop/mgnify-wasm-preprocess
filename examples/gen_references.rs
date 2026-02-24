/// Produces the BGZF-compressed files that samtools/tabix need to generate
/// the committed reference index files used by integration tests.
///
/// Usage:
///   cargo run --example gen_references -- <fasta_in> <fasta_bgz_out> <gff_in> <gff_bgz_out>
///
/// After running this, use tests/generate_references.sh to invoke samtools/tabix
/// on the outputs and commit the resulting .fai, .gzi, and .tbi files.

use std::fs;
use std::io::Cursor;

use mgnify_wasm::htslib::bgzf_compress;
use mgnify_wasm::gff_preprocess;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 5 {
        eprintln!("Usage: gen_references <fasta_in> <fasta_bgz_out> <gff_in> <gff_bgz_out>");
        std::process::exit(1);
    }

    let (fasta_in, fasta_out, gff_in, gff_out) = (&args[1], &args[2], &args[3], &args[4]);

    // FASTA: read and BGZF-compress
    let fa_bytes = fs::read(fasta_in)
        .unwrap_or_else(|e| panic!("cannot read {}: {}", fasta_in, e));
    let mut fasta_bgz = Vec::new();
    bgzf_compress(Cursor::new(&fa_bytes), &mut fasta_bgz).expect("FASTA bgzf_compress failed");
    fs::write(fasta_out, &fasta_bgz)
        .unwrap_or_else(|e| panic!("cannot write {}: {}", fasta_out, e));
    eprintln!("Wrote {} bytes → {}", fasta_bgz.len(), fasta_out);

    // GFF3: read, preprocess (sort + strip ##FASTA), then BGZF-compress
    let gff_raw = fs::read_to_string(gff_in)
        .unwrap_or_else(|e| panic!("cannot read {}: {}", gff_in, e));
    let gff_preprocessed = gff_preprocess(&gff_raw);
    let mut gff_bgz = Vec::new();
    bgzf_compress(Cursor::new(gff_preprocessed.as_bytes()), &mut gff_bgz)
        .expect("GFF bgzf_compress failed");
    fs::write(gff_out, &gff_bgz)
        .unwrap_or_else(|e| panic!("cannot write {}: {}", gff_out, e));
    eprintln!("Wrote {} bytes → {}", gff_bgz.len(), gff_out);
}
