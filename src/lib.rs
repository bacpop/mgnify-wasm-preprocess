use std::io::Read;

use wasm_bindgen::prelude::*;
use wasm_bindgen_file_reader::WebSysFile;

use crate::decompress::open_file_maybe_gz;
extern crate console_error_panic_hook;
mod decompress;

pub mod htslib;
use crate::htslib::{compress_bgzf, index_gff_tbi, index_fasta_fai, FaidxResult};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);

    #[wasm_bindgen(js_name = postMessage)]
    fn post_message(data: &JsValue);
}

/// Logging wrapper function
pub fn logw(text : &str, typ : Option<&str>) {
    if typ.is_some() {
        log((String::from("mgnify_preprocess::") + typ.unwrap() + "::" + text).as_str());
    } else {
        log(text);
    }
}

#[wasm_bindgen]
/// Function that allows to propagate panic error messages when compiling to wasm, see https://github.com/rustwasm/console_error_panic_hook
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

/// Convert an owned `Vec<u8>` into a JS `Blob` with one copy (Rust heap → JS heap).
fn vec_to_blob(data: Vec<u8>) -> Result<web_sys::Blob, JsValue> {
    let arr = js_sys::Uint8Array::from(data.as_slice());
    let seq = js_sys::Array::of1(&arr);
    web_sys::Blob::new_with_u8_array_sequence(&seq)
}

#[wasm_bindgen]
/// Main struct that acts as wrapper of the assembler when compiling to wasm
pub struct IndexGen {
    fasta_bgz: Vec<u8>,
    fasta_fai: Vec<u8>,
    fasta_gzi: Vec<u8>,
    gff_bgz: Vec<u8>,
    gff_idx: Vec<u8>,
}


#[wasm_bindgen]
impl IndexGen {
    /// Constructor/initialiser of the wasm assembler. It also performs the preprocessing.
    pub fn new(fa_file : web_sys::File, gff_file : web_sys::File) -> Self {
        if cfg!(debug_assertions) {
            init_panic_hook();
        }

        // Read in files and preprocess
        logw("Reading fasta and gff into memory", None);
        let mut wf_fa = WebSysFile::new(fa_file);
        let mut wf_gff = WebSysFile::new(gff_file);

        let mut fa_reader = open_file_maybe_gz(&mut wf_fa);
        let mut gff_reader = open_file_maybe_gz(&mut wf_gff);

        let mut fa_bytes = Vec::new();
        fa_reader.read_to_end(&mut fa_bytes).expect_throw("fasta read failed");

        let mut gff_string = String::new();
        gff_reader.read_to_string(&mut gff_string).expect_throw("GFF read failed");
        gff_string = gff_preprocess(&gff_string);

        // Output fasta files
        logw("Compressing and indexing fasta", None);
        // bgzip
        let fasta_bgz = compress_bgzf(&fa_bytes);
        // faidx
        let FaidxResult { fai: fasta_fai, gzi: fasta_gzi } = index_fasta_fai(&fasta_bgz);

        // Output gff files
        logw("Compressing and indexing gff", None);
        // bgzip
        let gff_bgz = compress_bgzf(gff_string.as_bytes());
        let gff_idx = index_gff_tbi(&gff_bgz);

        Self {
            fasta_bgz,
            fasta_fai,
            fasta_gzi,
            gff_bgz,
            gff_idx,
        }
    }

    /// Returns the BGZF-compressed FASTA as a Blob. Drains the field; call once.
    pub fn fasta_bgz_blob(&mut self) -> Result<web_sys::Blob, JsValue> {
        vec_to_blob(std::mem::take(&mut self.fasta_bgz))
    }

    /// Returns the FASTA `.fai` index as a Blob. Drains the field; call once.
    pub fn fasta_fai_blob(&mut self) -> Result<web_sys::Blob, JsValue> {
        vec_to_blob(std::mem::take(&mut self.fasta_fai))
    }

    /// Returns the FASTA `.gzi` block index as a Blob. Drains the field; call once.
    pub fn fasta_gzi_blob(&mut self) -> Result<web_sys::Blob, JsValue> {
        vec_to_blob(std::mem::take(&mut self.fasta_gzi))
    }

    /// Returns the BGZF-compressed GFF3 as a Blob. Drains the field; call once.
    pub fn gff_bgz_blob(&mut self) -> Result<web_sys::Blob, JsValue> {
        vec_to_blob(std::mem::take(&mut self.gff_bgz))
    }

    /// Returns the GFF3 `.tbi` tabix index as a Blob. Drains the field; call once.
    pub fn gff_tbi_blob(&mut self) -> Result<web_sys::Blob, JsValue> {
        vec_to_blob(std::mem::take(&mut self.gff_idx))
    }
}

// Reorders start for indexing and removes sequence if present
pub fn gff_preprocess(gff_string: &str) -> String {
    let mut outbuf = String::new();
    let mut records: Vec<&str> = Vec::new();

    for line in gff_string.split('\n') {
        if line.starts_with("##FASTA") {
            break;
        }
        if line.starts_with('#') {
            outbuf.push_str(line);
            outbuf.push('\n');
        } else if !line.is_empty() {
            records.push(line);
        }
    }

    // Emulating `sort -k1,1d -k4,4n -k5,5n`
    records.sort_by(|a, b| {
        let a_fields: Vec<&str> = a.split('\t').collect();
        let b_fields: Vec<&str> = b.split('\t').collect();

        // k1,1d - dictionary order on field 1 (index 0)
        a_fields[0].cmp(&b_fields[0])
            // k4,4n - numeric on field 4 (index 3)
            .then_with(|| {
                let a4: i64 = a_fields[3].parse().unwrap_or(0);
                let b4: i64 = b_fields[3].parse().unwrap_or(0);
                a4.cmp(&b4)
            })
            // k5,5n - numeric on field 5 (index 4)
            .then_with(|| {
                let a5: i64 = a_fields[4].parse().unwrap_or(0);
                let b5: i64 = b_fields[4].parse().unwrap_or(0);
                a5.cmp(&b5)
            })
    });

    for rec in &records {
        outbuf.push_str(rec);
        outbuf.push('\n');
    }

    outbuf
}
