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

# Check fixtures exist
if [[ ! -f "$FASTA" ]]; then
    echo "ERROR: missing fixture $FASTA" >&2
    exit 1
fi
if [[ ! -f "$GFF" ]]; then
    echo "ERROR: missing fixture $GFF" >&2
    exit 1
fi

# Check required tools
for tool in bgzip samtools tabix; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: '$tool' not found on \$PATH" >&2
        exit 1
    fi
done

mkdir -p "$REF"

FASTA_BGZ="$REF/test.fasta.bgz"
GFF_BGZ="$REF/test.gff3.bgz"

# Step 1: produce BGZF files using our implementation so that the virtual
# offsets in the reference indexes match what our code produces.
echo "Building gen_references example..."
cargo build --example gen_references --manifest-path "$REPO_ROOT/Cargo.toml"

echo "Compressing fixtures with our BGZF implementation..."
"$REPO_ROOT/target/debug/examples/gen_references" \
    "$FASTA" "$FASTA_BGZ" \
    "$GFF"   "$GFF_BGZ"

# Step 2: sanity-check our BGZF with the reference tool
echo "Verifying BGZF output with bgzip -t..."
bgzip -t "$FASTA_BGZ"
bgzip -t "$GFF_BGZ"
echo "  bgzip validation passed."

# Step 3: generate reference index files from our BGZF output
# samtools and tabix write their output alongside the input file, which is
# already in $REF, so no copying is needed.
echo "Running samtools faidx..."
samtools faidx "$FASTA_BGZ"

echo "Running tabix..."
tabix -p gff "$GFF_BGZ"

# Clean up the intermediate .bgz files — tests regenerate them at runtime
rm -f "$FASTA_BGZ" "$GFF_BGZ"

echo ""
echo "Reference files written to $REF:"
ls -lh "$REF"
echo ""
echo "Commit these files alongside the fixtures."
