use std::collections::HashMap;
use std::io::{self, Read, Write};
use super::bgzf::BgzfReader;

// ---------------------------------------------------------------------------
// GFF3 preset constants (hts_idx_t TBX_GENERIC/TBX_GFF)
// ---------------------------------------------------------------------------
const COL_SEQ: i32 = 1;
const COL_BEG: i32 = 4;
const COL_END: i32 = 5;
const META_CHAR: i32 = b'#' as i32;
const SKIP: i32 = 0;
const MIN_SHIFT: u32 = 14;
const N_LVLS: u32 = 5;
const FORMAT: i32 = 0; // TBX_GENERIC

// ---------------------------------------------------------------------------
// Binning scheme (hts_reg2bin equivalent)
// ---------------------------------------------------------------------------

/// Compute the bin number for a 0-based half-open interval [beg, end).
/// Uses the same scheme as hts_reg2bin(beg, end, min_shift=14, n_lvls=5).
fn reg2bin(beg: u64, end: u64) -> u32 {
    // t = offset of first bin at each level:
    // level k: offset = sum_{i=1}^{k} 8^i = (8^(k+1) - 8) / 7
    // For min_shift=14, n_lvls=5:
    //   level 5 (finest): bin offset starts at (8^5 - 8)/7 + ... but we use
    //   the htslib formula directly.
    //
    // hts_reg2bin counts from the bottom (finest) up:
    //   start with s = min_shift, t = accumulated offset at finest level
    //   t for 5 levels: (8^5 - 1)/7 = 4681 - 0 = 4681? Let's compute:
    //   total bins = sum_{k=0}^{5} 8^k = (8^6 - 1)/7 = 37449
    //   finest-level offset = 37449 - 8^5 = 37449 - 32768 = 4681
    let mut e = end.saturating_sub(1);
    let mut s: u32 = MIN_SHIFT;
    // offset of the finest level: (8^(n_lvls+1) - 1)/7 - 8^n_lvls
    // For n_lvls=5: (8^6-1)/7 - 8^5 = 37449 - 32768 = 4681
    // But htslib accumulates from finest down; we replicate the loop:
    let mut t: u64 = ((1u64 << (3 * N_LVLS + 3)) - 1) / 7; // = 37449 for n_lvls=5
    // The loop goes from finest level (n_lvls) down to 1
    for l in (1..=N_LVLS).rev() {
        t -= 1u64 << (3 * l);
        if (beg >> s) == (e >> s) {
            return (t + (beg >> s)) as u32;
        }
        s += 3;
        e >>= 0; // already computed above
    }
    // Level 0 (whole sequence): bin 0
    0
}

// ---------------------------------------------------------------------------
// Index data structures
// ---------------------------------------------------------------------------

#[derive(Clone)]
struct Chunk {
    start: u64,
    end: u64,
}

struct SeqIdx {
    name: String,
    bins: HashMap<u32, Vec<Chunk>>,
    lidx: Vec<u64>, // one slot per 16384 bp window
}

impl SeqIdx {
    fn new(name: String) -> Self {
        SeqIdx {
            name,
            bins: HashMap::new(),
            lidx: Vec::new(),
        }
    }

    fn add_chunk(&mut self, bin: u32, chunk: Chunk) {
        self.bins.entry(bin).or_default().push(chunk);
    }

    fn update_lidx(&mut self, beg: u64, end: u64, voff: u64) {
        if end == 0 {
            return;
        }
        let win_beg = (beg >> MIN_SHIFT) as usize;
        let win_end = ((end - 1) >> MIN_SHIFT) as usize;
        if win_end >= self.lidx.len() {
            self.lidx.resize(win_end + 1, 0);
        }
        for i in win_beg..=win_end {
            if self.lidx[i] == 0 {
                self.lidx[i] = voff;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Build a tabix index for a BGZF-compressed GFF3 file.
///
/// Reads from `bgzf_input` (a BGZF-compressed byte stream) and writes the
/// binary `.tbi` index to `tbi_output`.
pub fn tabix_index_gff<R: Read, W: Write>(bgzf_input: R, mut tbi_output: W) -> io::Result<()> {
    let mut reader = BgzfReader::new(bgzf_input);

    let mut seqs: Vec<SeqIdx> = Vec::new();
    let mut seq_map: HashMap<String, usize> = HashMap::new();

    let mut line_buf = Vec::with_capacity(4096);

    loop {
        line_buf.clear();
        let (n, voff_start) = reader.read_line(&mut line_buf)?;
        if n == 0 {
            break;
        }

        // Strip trailing newline/CR for parsing, but keep voff_start
        let line = strip_newline(&line_buf);

        // Skip empty lines and comment/meta lines
        if line.is_empty() || line[0] == b'#' {
            continue;
        }

        // Split on tabs
        let fields: Vec<&[u8]> = line.splitn(6, |&b| b == b'\t').collect();
        if fields.len() < 5 {
            // Not enough fields — skip
            continue;
        }

        let seqname = std::str::from_utf8(fields[0])
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "non-UTF8 sequence name"))?
            .to_owned();

        let start_1: u64 = parse_u64(fields[3])?;
        let end_1: u64 = parse_u64(fields[4])?;

        // GFF3 columns are 1-based, inclusive → convert to 0-based half-open
        let beg = start_1.saturating_sub(1);
        let end = end_1; // end_col is already 1-based inclusive = 0-based exclusive

        // Virtual offset after the line
        let voff_end = reader.virtual_offset();

        let tid = match seq_map.get(&seqname) {
            Some(&id) => id,
            None => {
                let id = seqs.len();
                seqs.push(SeqIdx::new(seqname.clone()));
                seq_map.insert(seqname, id);
                id
            }
        };

        let bin = reg2bin(beg, end);
        let chunk = Chunk { start: voff_start, end: voff_end };
        seqs[tid].add_chunk(bin, chunk);
        seqs[tid].update_lidx(beg, end, voff_start);
    }

    // Fill trailing zeros in lidx: any 0 slot after the last populated entry
    // should be set to the current EOF virtual offset.
    let eof_voff = reader.virtual_offset();
    for seq in &mut seqs {
        // Walk forward: set zero slots after last non-zero to eof_voff
        let mut seen_nonzero = false;
        for slot in seq.lidx.iter_mut() {
            if *slot != 0 {
                seen_nonzero = true;
            } else if seen_nonzero {
                *slot = eof_voff;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Write the .tbi binary format (all little-endian)
    // -----------------------------------------------------------------------

    // Magic
    tbi_output.write_all(b"TBI\x01")?;

    // n_ref
    write_i32(&mut tbi_output, seqs.len() as i32)?;

    // Header: format, col_seq, col_beg, col_end, meta_char, skip
    write_i32(&mut tbi_output, FORMAT)?;
    write_i32(&mut tbi_output, COL_SEQ)?;
    write_i32(&mut tbi_output, COL_BEG)?;
    write_i32(&mut tbi_output, COL_END)?;
    write_i32(&mut tbi_output, META_CHAR)?;
    write_i32(&mut tbi_output, SKIP)?;

    // l_nm + names (null-terminated, concatenated)
    let mut names_buf: Vec<u8> = Vec::new();
    for seq in &seqs {
        names_buf.extend_from_slice(seq.name.as_bytes());
        names_buf.push(0);
    }
    write_i32(&mut tbi_output, names_buf.len() as i32)?;
    tbi_output.write_all(&names_buf)?;

    // Per-sequence index data
    for seq in &seqs {
        // Collect bins in a stable order (sort by bin number)
        let mut bin_ids: Vec<u32> = seq.bins.keys().cloned().collect();
        bin_ids.sort_unstable();

        write_i32(&mut tbi_output, bin_ids.len() as i32)?;
        for bin in &bin_ids {
            let chunks = &seq.bins[bin];
            tbi_output.write_all(&bin.to_le_bytes())?;
            write_i32(&mut tbi_output, chunks.len() as i32)?;
            for chunk in chunks {
                tbi_output.write_all(&chunk.start.to_le_bytes())?;
                tbi_output.write_all(&chunk.end.to_le_bytes())?;
            }
        }

        // Linear index
        write_i32(&mut tbi_output, seq.lidx.len() as i32)?;
        for &slot in &seq.lidx {
            tbi_output.write_all(&slot.to_le_bytes())?;
        }
    }

    // n_no_coor (unaligned read count) = 0
    tbi_output.write_all(&0u64.to_le_bytes())?;

    Ok(())
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn strip_newline(buf: &[u8]) -> &[u8] {
    let mut end = buf.len();
    while end > 0 && (buf[end - 1] == b'\n' || buf[end - 1] == b'\r') {
        end -= 1;
    }
    &buf[..end]
}

fn parse_u64(bytes: &[u8]) -> io::Result<u64> {
    let s = std::str::from_utf8(bytes)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "non-UTF8 field"))?
        .trim();
    s.parse::<u64>()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("cannot parse integer: {:?}", s)))
}

fn write_i32<W: Write>(w: &mut W, v: i32) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}
