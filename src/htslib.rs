use wasm_bindgen::prelude::*;
use std::io::Cursor;

mod bgzf;
mod tabix;
mod faidx;

pub use bgzf::{BgzfWriter, BgzfReader, bgzf_compress};
pub use tabix::csi_index_gff;
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

/// Build a tabix `.csi` index from a BGZF-compressed GFF3 byte slice.
#[wasm_bindgen]
pub fn index_gff_csi(bgzf_input: &[u8]) -> Vec<u8> {
    let mut csi = Vec::new();
    csi_index_gff(Cursor::new(bgzf_input), &mut csi)
        .expect("csi_index_gff failed");
    csi
}

/// Result of indexing a BGZF-compressed FASTA file.
#[wasm_bindgen]
pub struct FaidxResult {
    pub(crate) fai: Vec<u8>,
    pub(crate) gzi: Vec<u8>,
}

#[wasm_bindgen]
impl FaidxResult {
    /// Moves the `.fai` index bytes out. May only be called once meaningfully.
    pub fn fai(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.fai)
    }
    /// Moves the `.gzi` index bytes out. May only be called once meaningfully.
    pub fn gzi(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.gzi)
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
