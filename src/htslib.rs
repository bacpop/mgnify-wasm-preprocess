use wasm_bindgen::prelude::*;
use std::io::Cursor;

mod bgzf;
mod tabix;
mod faidx;

pub use bgzf::{BgzfWriter, BgzfReader, bgzf_compress};
pub use tabix::tabix_index_gff;
pub use faidx::faidx_index_fasta;

// ---------------------------------------------------------------------------
// WASM-bindgen exports
// ---------------------------------------------------------------------------

/// Compress raw bytes into BGZF format.
#[wasm_bindgen]
pub fn compress_bgzf(input: &[u8]) -> Vec<u8> {
    let mut output = Vec::new();
    bgzf_compress(Cursor::new(input), &mut output)
        .expect("bgzf_compress failed");
    output
}

/// Build a tabix `.tbi` index from a BGZF-compressed GFF3 byte slice.
#[wasm_bindgen]
pub fn index_gff_tbi(bgzf_input: &[u8]) -> Vec<u8> {
    let mut tbi = Vec::new();
    tabix_index_gff(Cursor::new(bgzf_input), &mut tbi)
        .expect("tabix_index_gff failed");
    tbi
}

/// Result of indexing a BGZF-compressed FASTA file.
/// Fields are exposed via getter methods because `Vec<u8>` is not `Copy`.
#[wasm_bindgen]
pub struct FaidxResult {
    fai: Vec<u8>,
    gzi: Vec<u8>,
}

#[wasm_bindgen]
impl FaidxResult {
    /// Returns the `.fai` index bytes.
    pub fn fai(&self) -> Vec<u8> {
        self.fai.clone()
    }
    /// Returns the `.gzi` index bytes.
    pub fn gzi(&self) -> Vec<u8> {
        self.gzi.clone()
    }
}

/// Build `.fai` and `.gzi` indexes from a BGZF-compressed FASTA byte slice.
#[wasm_bindgen]
pub fn index_fasta_fai(bgzf_input: &[u8]) -> FaidxResult {
    let mut fai = Vec::new();
    let mut gzi = Vec::new();
    faidx_index_fasta(Cursor::new(bgzf_input), &mut fai, &mut gzi)
        .expect("faidx_index_fasta failed");
    FaidxResult { fai, gzi }
}
