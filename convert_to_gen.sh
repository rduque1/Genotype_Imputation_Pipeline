#!/usr/bin/env bash
# Convert imputed VCFs to GEN/SAMPLE format using plink2.
# Usage: $0 --input-dir <dir> --out <dir>
set -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/utils.sh"

INPUT_DIR=""
OUTDIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --input-dir) INPUT_DIR="$2"; shift 2 ;;
        --out)       OUTDIR="$2"; shift 2 ;;
        *)           echo "Unknown argument: $1"; exit 1 ;;
    esac
done

[[ -z "$INPUT_DIR" || -z "$OUTDIR" ]] && { echo "Missing required arguments"; exit 1; }

mkdir -p "$OUTDIR"

for vcf in "${INPUT_DIR}"/*.dose.vcf.gz; do
    [[ -f "$vcf" ]] || continue
    base=$(basename "$vcf" .dose.vcf.gz)
    out="${OUTDIR}/${base}"

    if [[ -f "${out}.gen" || -f "${out}.gen.gz" ]]; then
        print_info "GEN output already exists, skipping ${base}"
        continue
    fi

    print_info "Converting ${base}..."
    plink2 --vcf "$vcf" dosage=GP-force \
           --export oxford \
           --out "$out"
    print_success "${base}.gen"
done

print_success "GEN conversion complete"
