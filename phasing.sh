#!/usr/bin/env bash
# Phase per-chromosome VCF files with Eagle.
# Usage: $0 --input-dir <dir> --out <dir> --ref-bcf <dir> [--threads N]
#   --input-dir: directory with *.chr{1..22,X}.vcf.gz files
#   --ref-bcf:   directory with 1000G BCF reference files
set -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/utils.sh"

INPUT_DIR=""
OUTDIR=""
REF_BASE=""
THREADS=4

while [[ $# -gt 0 ]]; do
    case $1 in
        --input-dir) INPUT_DIR="$2"; shift 2 ;;
        --out)       OUTDIR="$2"; shift 2 ;;
        --ref-bcf)   REF_BASE="$2"; shift 2 ;;
        --threads)   THREADS="$2"; shift 2 ;;
        *)           echo "Unknown argument: $1"; exit 1 ;;
    esac
done

[[ -z "$INPUT_DIR" || -z "$OUTDIR" || -z "$REF_BASE" ]] && { echo "Missing required arguments"; exit 1; }

MYMAP="${SCRIPT_DIR}/required_tools/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

phase_chr() {
    local chr=$1
    local vcf=$(ls "${INPUT_DIR}"/*.chr${chr}.vcf.gz 2>/dev/null | head -1)
    if [[ -z "$vcf" ]]; then
        print_warning "No input VCF found for chromosome ${chr}, skipping"
        return
    fi
    local prefix=$(basename "$vcf" .vcf.gz)

    if [[ -f "${prefix}.phased.vcf.gz" ]]; then
        print_info "Chromosome ${chr} already phased, skipping: ${prefix}.phased.vcf.gz"
        return
    fi

    print_info "Phasing chromosome ${chr}..."

    # Reference file
    local myref
    if [[ "$chr" == "X" ]]; then
        myref="${REF_BASE}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.bcf"
    else
        myref="${REF_BASE}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf"
    fi

    # Run Eagle
    eagle --vcfTarget="$vcf" \
        --vcfRef="$myref" \
        --noImpMissing \
        --geneticMapFile="$MYMAP" \
        --Kpbwt=100000 --numThreads="$THREADS" \
        --chrom="$chr" --allowRefAltSwap \
        --outPrefix="${prefix}.phased"

    tabix -p vcf "${prefix}.phased.vcf.gz"
    print_success "Chromosome ${chr} phased"
}

export -f phase_chr print_info print_warning print_success
export MYMAP INPUT_DIR REF_BASE OUTDIR THREADS
export RED GREEN YELLOW BLUE NC

CHRS=$(ls "${INPUT_DIR}"/*.chr*.vcf.gz 2>/dev/null | sed 's/.*\.chr//;s/\.vcf\.gz//' | sort -Vu) # Extract chromosome names from filenames
print_info "Detected chromosomes for phasing: $CHRS"

if command -v parallel &>/dev/null; then
    N_JOBS=$(( THREADS / 4 ))
    [[ $N_JOBS -lt 1 ]] && N_JOBS=1
    parallel -j "$N_JOBS" phase_chr ::: $CHRS
else
    for chr in $CHRS; do
        phase_chr "$chr"
    done
fi

print_success "Phasing complete"
