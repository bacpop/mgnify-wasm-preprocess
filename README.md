# mgnify-wasm

A Rust library compiled to WebAssembly that preprocesses genomic files in the
browser.  Given a FASTA and a GFF3 file it:

1. Transparently decompresses either file if it arrives as gzip.
2. Preprocesses the GFF3 (strips any embedded `##FASTA` section, sorts records
   by seqname / start / end).
3. BGZF-compresses both files.
4. Generates the indexes needed for random-access queries: `.fai` + `.gzi` for
   the FASTA, `.csi` (tabix CSI format) for the GFF3.

All index bytes are produced entirely in-browser; no server-side preprocessing
is required.

---

## Prerequisites

### Rust tests (no external tools needed)

```
cargo  (any recent stable Rust toolchain)
```

### Reference file generation (requires external bioinformatics tools)

```
bgzip      (from htslib — used to validate BGZF output)
samtools   (≥ 1.10 recommended)
tabix      (from htslib)
```

All three ship together in the
[htslib](https://github.com/samtools/htslib) distribution and are commonly
available via conda/mamba:

```bash
mamba install -c bioconda samtools htslib
```

### WebAssembly build (optional)

```
wasm-pack  (https://rustwasm.github.io/wasm-pack/)
```

---

## Running Tests

Tests run against committed reference files and have no external tool
dependencies:

```bash
cargo test
```

The test suite contains 10 integration tests across two fixture pairs:

| Test | What it checks |
|------|---------------|
| `bgzf_roundtrip_fasta` / `bgzf_roundtrip_gff` | BGZF compress → decompress round-trips |
| `fai_matches_samtools` / `gzi_matches_samtools` | `.fai` and `.gzi` match `samtools faidx` output |
| `csi_matches_tabix` | `.csi` matches `tabix -C -p gff` output |
| `bgzf_roundtrip_bu_fasta` / `bgzf_roundtrip_bu_gff` | Same round-trip tests for the larger BU fixture |
| `bu_fai_matches_samtools` / `bu_gzi_matches_samtools` | `.fai` and `.gzi` for the BU fixture |
| `bu_csi_matches_tabix` | `.csi` for the BU fixture |

The two fixtures used are:

* `tests/fixtures/test.fasta` + `tests/fixtures/test.gff3` — small synthetic
  files for fast iteration.
* `tests/fixtures/BU_ATCC8492VPI0062_NT5002.1.fa.gz` +
  `tests/fixtures/BU_ATCC8492_annotations.gff.gz` — a real gzip-compressed
  genome/annotation pair that exercises multi-block BGZF paths.

---

## Generating Reference Files

Reference index files live in `tests/fixtures/reference/` and are committed to
the repository.  They must be regenerated whenever:

* A fixture file changes, or
* The BGZF writer or indexing logic changes in a way that affects virtual
  offsets (which would make samtools/tabix produce different index bytes).

### Step-by-step

```bash
# Activate an environment that has bgzip, samtools, and tabix on $PATH.
# Example with mamba:
mamba activate mgnify_viewer

# From the repository root:
tests/generate_references.sh
```

The script:

1. Builds the `gen_references` example binary (`cargo build --example
   gen_references`).
2. Uses **our own BGZF implementation** to compress each fixture into a
   temporary `.bgz` file.  This is critical: samtools and tabix must index the
   same BGZF stream that our code produces so that virtual offsets in the
   reference indexes correspond to the offsets our reader observes at test time.
3. Validates each `.bgz` with `bgzip -t`.
4. Runs `samtools faidx` on the FASTA `.bgz` to produce `.fai` and `.gzi`.
5. Runs `tabix -C -p gff` on the GFF3 `.bgz` to produce `.csi`.
6. Deletes the temporary `.bgz` files and leaves only the index files in
   `tests/fixtures/reference/`.

After running the script, commit the updated reference files.

---

## Building for WebAssembly

```bash
# Browser (ES module, outputs to pkg/)
wasm-pack build --target web

# Node.js / bundler
wasm-pack build --target bundler
wasm-pack build --target nodejs
```

The main WASM entry point is `IndexGen::new(fa_file, gff_file)` in `src/lib.rs`.
It returns a struct with getter methods for each output blob (`.bgz`, `.fai`,
`.gzi`, `.csi`).

The lower-level functions are also exported directly via `wasm-bindgen`:

| Function | Description |
|----------|-------------|
| `compress_bgzf(input)` | Compress raw bytes to BGZF |
| `index_fasta_fai(bgzf_input)` | Build `.fai` + `.gzi` from a BGZF FASTA |
| `index_gff_csi(bgzf_input)` | Build `.csi` from a BGZF GFF3 |

---

## Differences from htslib

The implementation closely follows htslib's algorithms but differs in a few
places.  Understanding these differences matters when comparing output
byte-for-byte.

### BGZF block contents

htslib's `bgzip` uses its own zlib invocation, which may choose different
Deflate literal/length/distance trees than the `flate2` crate used here.
Consequently, **BGZF blocks produced by this library are not byte-for-byte
identical to those produced by `bgzip`**, even though they decompress to the
same data and pass `bgzip -t` validation.

This is why the reference generation script compresses fixtures with *our*
implementation before running samtools/tabix: the virtual offsets encoded in
`.fai`, `.gzi`, and `.csi` depend on exactly where BGZF block boundaries fall.
If the blocks were placed differently, all virtual offsets in every index would
shift, and the tests would fail.

### CSI bin ordering

htslib stores bins in a hash table (khash) and iterates them in bucket order,
which is neither ascending by bin number nor insertion order.  Our
implementation sorts bins ascending by bin number before writing.

The binary `.csi` files produced by our code and by `tabix -C` are therefore
not byte-for-byte identical.  The integration tests use a `normalize_csi()`
helper that sorts bins and chunks in both files before comparing, so the tests
pass despite this ordering difference.

If you need to compare the raw `.csi` output against tabix manually, decompress
both and parse the bin entries rather than doing a plain `diff`:

```bash
# Decompress our output for inspection
bgzip -d < our_output.csi | xxd | head -80

# Or compare structurally with htslib's htsfile utility
htsfile our_output.csi
```

### CSI parameters

We use tabix's default CSI parameters for GFF files:

| Parameter | Value |
|-----------|-------|
| `MIN_SHIFT` | 14 |
| `N_LVLS` | 8 |
| Binning scheme | Finest-first (old BAM / htslib `hts_reg2bin`) |
| Number of bins | 19,173,961 |
| Meta-bin | 19,173,962 |

These are the values tabix uses for `tabix -C -p gff`.

### Chunk merging

htslib merges consecutive index chunks whose virtual-offset gap is ≤
`HTS_MIN_MARKER_DIST` (0x10000 = 65,536 bytes, i.e., one BGZF block).  The
comparison is on full 64-bit virtual offsets:

```
merge if: next_chunk.start ≤ prev_chunk.end + 0x10000
```

We apply the same criterion.

### n_no_coor

The `.csi` `n_no_coor` field (count of records with no assigned coordinates) is
always written as 0.  We only index GFF3 records that have valid seqname, start,
and end fields; records that fail to parse are silently skipped rather than
counted.

### FAI seq_offset field

The `seq_offset` field in `.fai` records is a **plain uncompressed byte offset**,
not a BGZF virtual offset.  This matches samtools' own convention.

### GFF3 preprocessing

Before BGZF-compressing and indexing, the GFF3 is preprocessed:

* Any embedded `##FASTA` section (and everything after it) is stripped.
* Comment and directive lines (starting with `#`) are preserved in their
  original order, before any data records.
* Data records are sorted by `(seqname, start, end)` — equivalent to
  `sort -k1,1d -k4,4n -k5,5n`.

This sort is required because tabix indexing assumes the file is sorted; tabix
itself will refuse to index an unsorted file.

---

## Source layout

```
src/
  lib.rs              — WASM entry point (IndexGen), gff_preprocess()
  decompress.rs       — transparent gzip detection/decompression
  htslib/
    bgzf.rs           — BgzfWriter, BgzfReader, bgzf_compress()
    faidx.rs          — faidx_index_fasta() → .fai + .gzi
    tabix.rs          — csi_index_gff() → .csi
    mod.rs            — wasm-bindgen exports

examples/
  gen_references.rs   — CLI tool used by generate_references.sh

tests/
  integration_test.rs — 10 integration tests
  generate_references.sh
  fixtures/
    test.fasta
    test.gff3
    BU_ATCC8492VPI0062_NT5002.1.fa.gz
    BU_ATCC8492_annotations.gff.gz
    reference/        — committed reference index files
```
