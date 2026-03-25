#!/usr/bin/env bash
# check_chr_naming.sh
#
# Validates that target VCF files, reference BCF files, and the genetic map
# all use the same chromosome naming convention:
#   "chr"-prefixed  →  chr1, chrX
#   bare            →  1, X
#
# Style is detected from FILE CONTENT, not filenames.
# MSAV files (binary format) are not inspectable and are excluded from the check.
#
# Usage:
#   check_chr_naming.sh [--target-dir DIR] [--ref-bcf DIR] \
#                       [--ref-msav DIR]   [--genetic-map FILE]
#
# Returns 0 if all supplied inputs agree, 1 on mismatch.

set -eu

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/utils.sh"

TARGET_DIR=""
REF_BCF_DIR=""
REF_MSAV_DIR=""
GENETIC_MAP=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --target-dir)  TARGET_DIR="$2";  shift 2 ;;
        --ref-bcf)     REF_BCF_DIR="$2"; shift 2 ;;
        --ref-msav)    REF_MSAV_DIR="$2";shift 2 ;;  # accepted but not checked (binary)
        --genetic-map) GENETIC_MAP="$2"; shift 2 ;;
        --help|-h)
            echo "Usage: $0 [--target-dir DIR] [--ref-bcf DIR] [--ref-msav DIR] [--genetic-map FILE]"
            exit 0 ;;
        *) print_error "Unknown argument: $1"; exit 1 ;;
    esac
done

###############################################################################
# Helpers
###############################################################################

# Classify a chromosome token: "chr" | "nochr" | "unknown"
chr_style() {
    case "$1" in
        chr*)             echo "chr"   ;;
        [0-9]*|X|Y|MT|M) echo "nochr" ;;
        *)                echo "unknown" ;;
    esac
}

# Detect convention from a VCF or BCF file by inspecting its content
vcf_chr_style() {
    local chrom
    # Prefer ##contig header lines; fall back to first data record
    chrom=$(bcftools view -h "$1" 2>/dev/null \
        | grep "^##contig" | head -1 \
        | sed 's/.*ID=//; s/[,>].*//')
    [[ -z "$chrom" ]] && \
        chrom=$(bcftools view -H "$1" 2>/dev/null | head -1 | cut -f1)
    chr_style "$chrom"
}

# Detect convention from an Eagle genetic map (gzipped or plain text)
# Expected format: header line, then data with chromosome in column 1
map_chr_style() {
    local chrom
    if [[ "$1" == *.gz ]]; then
        chrom=$(zcat "$1" | awk 'NR>1 && !/^#/{ print $1; exit }')
    else
        chrom=$(awk 'NR>1 && !/^#/{ print $1; exit }' "$1")
    fi
    chr_style "$chrom"
}

# Detect convention from an MSAV file (binary) by running minimac4 with a
# throwaway dummy VCF and parsing "Imputing CHROM:..." from its stderr/stdout.
# Requires minimac4 on PATH (or MINIMAC4 env var).
msav_chr_style() {
    local msav="$1"
    local minimac="${MINIMAC4:-minimac4}"
    command -v "$minimac" &>/dev/null || { echo "unknown"; return; }

    local tmpdir
    tmpdir=$(mktemp -d)
    # Create a dummy VCF with a dummy chromosome; minimac4 will fail to impute
    # but will still emit "Imputing CHROM:..." before it errors out.
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n0\t1\t.\tA\tG\t.\t.\t.\tGT\t0/1\n' \
        > "${tmpdir}/dummy.vcf"
    bgzip "${tmpdir}/dummy.vcf"
    tabix -p vcf "${tmpdir}/dummy.vcf.gz"

    local chrom
    chrom=$( ("$minimac" "$msav" "${tmpdir}/dummy.vcf.gz" 2>&1 || true) \
        | grep -aP 'Imputing [^:]+:' | grep -oP 'Imputing \K[^:]+' | head -1 )

    rm -rf "$tmpdir"
    chr_style "$chrom"
}

###############################################################################
# Main
###############################################################################

print_step "Chromosome Naming Convention Check"

MISMATCH=0
REF_STYLE=""

check_style() {
    local label="$1" style="$2"
    [[ "$style" == "unknown" ]] && return
    print_info "${label}: ${style}"
    if [[ -z "$REF_STYLE" ]]; then
        REF_STYLE="$style"
    elif [[ "$style" != "$REF_STYLE" ]]; then
        print_error "Naming mismatch: ${label} uses '${style}' but expected '${REF_STYLE}'"
        MISMATCH=1
    fi
}

if [[ -n "$TARGET_DIR" ]]; then
    f=$(ls "${TARGET_DIR}"/*.vcf.gz 2>/dev/null | head -1 || true)
    [[ -n "$f" ]] && check_style "Target VCF   " "$(vcf_chr_style "$f")"
fi

if [[ -n "$REF_BCF_DIR" ]]; then
    f=$(ls "${REF_BCF_DIR}"/*.bcf 2>/dev/null | head -1 || true)
    [[ -n "$f" ]] && check_style "Reference BCF" "$(vcf_chr_style "$f")"
fi

if [[ -n "$GENETIC_MAP" ]]; then
    check_style "Genetic map  " "$(map_chr_style "$GENETIC_MAP")"
fi

if [[ -n "$REF_MSAV_DIR" ]]; then
    f=$(ls "${REF_MSAV_DIR}"/*.msav 2>/dev/null | head -1 || true)
    if [[ -n "$f" ]]; then
        check_style "Reference MSAV" "$(msav_chr_style "$f")"
    fi
fi

if [[ "$MISMATCH" -eq 1 ]]; then
    print_error "Chromosome naming mismatch — fix before proceeding."
    exit 1
fi

print_success "Chromosome naming OK (convention: ${REF_STYLE:-nothing to check})"
