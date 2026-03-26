#!/usr/bin/env bash
# Impute per-chromosome phased VCFs with Minimac4.
# Usage: $0 --input-dir <dir> --out <dir> --ref-msav <dir> [--threads N]
#   --input-dir: directory with *.phased.vcf.gz files (from Eagle)
#   --ref-msav:  directory with 1000G MSAV reference files
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
        --ref-msav)  REF_BASE="$2"; shift 2 ;;
        --threads)   THREADS="$2"; shift 2 ;;
        *)           echo "Unknown argument: $1"; exit 1 ;;
    esac
done

[[ -z "$INPUT_DIR" || -z "$OUTDIR" || -z "$REF_BASE" ]] && { echo "Missing required arguments"; exit 1; }

mkdir -p "$OUTDIR"
cd "$OUTDIR"

impute_chr() {
    local chr=$1
    local vcf=$(ls "${INPUT_DIR}"/*.chr${chr}.phased.vcf.gz 2>/dev/null | head -1)
    if [[ -z "$vcf" ]]; then
        print_warning "No phased input found for chromosome ${chr}, skipping"
        return
    fi
    local prefix=$(basename "$vcf" .phased.vcf.gz)
    print_info "Imputing chromosome ${chr}..."

    if [[ "$chr" == "X" ]]; then
        # for region in PAR1 nonPAR PAR2; do
        for region in nonPAR; do
            # local out_file="${prefix}.imputed.chrX_${region}.dose.vcf.gz"
            local out_file="${prefix}.imputed.chrX.dose.vcf.gz"
            if [[ -f "$out_file" ]]; then
                print_info "Chromosome X ${region} already imputed, skipping: ${out_file}"
                continue
            fi
            local myref="${REF_BASE}/ALL.chrX_${region}.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.msav"
            print_info "Imputing chromosome X ${region} with reference ${myref}"
            minimac4 "$myref" "$vcf" \
                --output "$out_file" \
                --output-format vcf.gz \
                --format GT,DS,GP \
                --min-ratio 0.001 --chunk 30000000 --overlap 3000000 --threads "$THREADS"
        done
    else
        local out_file="${prefix}.imputed.chr${chr}.dose.vcf.gz"
        if [[ -f "$out_file" ]]; then
            print_info "Chromosome ${chr} already imputed, skipping: ${out_file}"
            return
        fi
        local myref="${REF_BASE}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.msav"
        print_info "Imputing chromosome ${chr} with reference ${myref}"
        minimac4 "$myref" "$vcf" \
            --output "$out_file" \
            --output-format vcf.gz \
            --format GT,DS,GP \
            --min-ratio 0.001 --chunk 30000000 --overlap 3000000 --threads "$THREADS"
    fi

    print_success "Chromosome ${chr} imputed"
}

export -f impute_chr print_info print_warning print_success
export INPUT_DIR REF_BASE OUTDIR THREADS
export RED GREEN YELLOW BLUE NC

# Get list of chromosomes to process based on input files
CHRS=$(ls "${INPUT_DIR}"/*.chr*.phased.vcf.gz 2>/dev/null | sed 's/.*\.chr//;s/\.phased\.vcf\.gz//' | sort -Vu)
print_info "Detected chromosomes for imputation: $CHRS"

if command -v parallel &>/dev/null; then
    N_JOBS=$(( THREADS / 4 ))
    [[ $N_JOBS -lt 1 ]] && N_JOBS=1
    parallel -j "$N_JOBS" impute_chr ::: $CHRS
else
    for chr in $CHRS; do
        impute_chr "$chr"
    done
fi

print_success "Imputation complete"
