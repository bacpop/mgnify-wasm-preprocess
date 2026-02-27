#!/usr/bin/env bash
# Regenerate the committed reference index files used by integration tests.
#
# Run this script whenever the fixture files change, or when the BGZF/indexing
# implementation changes in a way that would alter output.
#
# Requirements: bgzip, samtools, tabix must be on $PATH.
#
# Usage (from repo root):
#   tests/generate_references.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FIXTURES="$REPO_ROOT/tests/fixtures"
REF="$FIXTURES/reference"

FASTA="$FIXTURES/test.fasta"
GFF="$FIXTURES/test.gff3"
BU_FASTA="$FIXTURES/BU_ATCC8492VPI0062_NT5002.1.fa.gz"
BU_GFF="$FIXTURES/BU_ATCC8492_annotations.gff.gz"

# Check fixtures exist
for f in "$FASTA" "$GFF" "$BU_FASTA" "$BU_GFF"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing fixture $f" >&2
        exit 1
    fi
done

# Check required tools
for tool in bgzip samtools tabix; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: '$tool' not found on \$PATH" >&2
        exit 1
    fi
done

mkdir -p "$REF"

# Step 1: produce BGZF files using our implementation so that the virtual
# offsets in the reference indexes match what our code produces.
echo "Building gen_references example..."
cargo build --example gen_references --manifest-path "$REPO_ROOT/Cargo.toml"

# --- test.fasta + test.gff3 ---
FASTA_BGZ="$REF/test.fasta.bgz"
GFF_BGZ="$REF/test.gff3.bgz"

echo "Compressing test fixtures with our BGZF implementation..."
"$REPO_ROOT/target/debug/examples/gen_references" \
    "$FASTA" "$FASTA_BGZ" \
    "$GFF"   "$GFF_BGZ"

echo "Verifying BGZF output with bgzip -t..."
bgzip -t "$FASTA_BGZ"
bgzip -t "$GFF_BGZ"
echo "  bgzip validation passed."

echo "Running samtools faidx (test.fasta)..."
samtools faidx "$FASTA_BGZ"

echo "Running tabix (test.gff3)..."
tabix -C -p gff "$GFF_BGZ"

rm -f "$FASTA_BGZ" "$GFF_BGZ"

# --- BU_ATCC8492 ---
BU_FASTA_BGZ="$REF/BU_ATCC8492.fasta.bgz"
BU_GFF_BGZ="$REF/BU_ATCC8492.gff3.bgz"

echo "Compressing BU_ATCC8492 fixtures with our BGZF implementation..."
"$REPO_ROOT/target/debug/examples/gen_references" \
    "$BU_FASTA" "$BU_FASTA_BGZ" \
    "$BU_GFF"   "$BU_GFF_BGZ"

echo "Verifying BU BGZF output with bgzip -t..."
bgzip -t "$BU_FASTA_BGZ"
bgzip -t "$BU_GFF_BGZ"
echo "  bgzip validation passed."

echo "Running samtools faidx (BU_ATCC8492)..."
samtools faidx "$BU_FASTA_BGZ"

echo "Running tabix (BU_ATCC8492)..."
tabix -C -p gff "$BU_GFF_BGZ"

rm -f "$BU_FASTA_BGZ" "$BU_GFF_BGZ"

echo ""
echo "Reference files written to $REF:"
ls -lh "$REF"
echo ""
echo "Commit these files alongside the fixtures."
