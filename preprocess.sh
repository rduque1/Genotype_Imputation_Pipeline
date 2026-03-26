#!/usr/bin/env bash

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source "${SCRIPT_DIR}/utils.sh"

usage() {
    cat << EOF
Usage: $0 --vcf <input_file> --out <output_dir> [OPTIONS]

Preprocessing Pipeline:
  1. Clean input file (VCF or 23andMe format)
  2. Detect genome build and lift to GRCh37 if needed
  3. Run Genotype Harmonizer for QC1 (optional)

Required Arguments:
  --vcf, -v PATH        Input file (VCF/VCF.gz or 23andMe .txt/.txt.gz/.zip)
  --out, -o PATH        Output directory for all results

Optional Arguments:
  --ref-path PATH       Reference panel path for Genotype Harmonizer (VCF files)
  --fasta-ref PATH      Reference FASTA file for fixref step
  --ref-bcf PATH        Reference BCF directory for Eagle phasing
  --ref-msav PATH       Reference MSAV directory for Minimac4 imputation
  --threads NUM         Number of threads (default: 4)

Example:
  $0 --vcf /path/to/input.txt --out /path/to/output --ref-path /path/to/reference --ref-bcf /path/to/bcf
EOF
    exit 1
}

# Parse arguments
MYINPUT=""
OUTROOT=""
REF_PATH=""
FASTA_REF=""
REF_BCF=""
REF_MSAV=""
THREADS=4

while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf|-v) MYINPUT="$2"; shift 2 ;;
        --out|-o) OUTROOT="$2"; shift 2 ;;
        --ref-path) REF_PATH="$2"; shift 2 ;;
        --fasta-ref) FASTA_REF="$2"; shift 2 ;;
        --ref-bcf) REF_BCF="$2"; shift 2 ;;
        --ref-msav) REF_MSAV="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --help|-h) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

if [[ -z "$MYINPUT" || -z "$OUTROOT" ]]; then
    usage
fi

if [[ ! -f "$MYINPUT" ]]; then
    print_error "Input file not found: $MYINPUT"
    exit 1
fi

print_step "Preprocessing Pipeline"
print_info "Input: $MYINPUT"
print_info "Output: $OUTROOT"

################################################################################
# Step 1: Clean input file
################################################################################

print_step "Step 1: Cleaning Input File"

bash "${SCRIPT_DIR}/clean_input_file.sh" --vcf "$MYINPUT" --out "$OUTROOT"

CLEANED_DATA="${OUTROOT}/sorted_rawData.txt.dedup.txt"

if [[ ! -f "$CLEANED_DATA" ]]; then
    print_error "Cleaned data file not found: $CLEANED_DATA"
    exit 1
fi

awk '{print $2}' "$CLEANED_DATA" | sort -Vu

print_success "Step 1 complete: $CLEANED_DATA"

################################################################################
# Step 2: Detect build and lift to GRCh37
################################################################################

print_step "Step 2: Build Detection and Lifting"

cd "$OUTROOT"
bash "${SCRIPT_DIR}/lifting.sh" "$CLEANED_DATA"

if [[ ! -f "${OUTROOT}/result.bed" ]]; then
    print_error "Lifting failed - result.bed not found"
    exit 1
fi

print_success "Step 2 complete: ${OUTROOT}/result.bed"

################################################################################
# Step 2b: Split by chromosome (Autosomes + X only)
################################################################################

# Filter the chromosome list to only include 1-22 and X
# This ignores Y (24), XY (25), and MT (26)
target_chroms=$(awk '{print $1}' "${OUTROOT}/result.bim" | sort -u | grep -E '^([1-9]|1[0-9]|2[0-2]|X)$')

print_step "Step 2b: Splitting by Chromosome"

cut -f1 "${OUTROOT}/result.bim" | sort -u

for chr in $target_chroms; do
    plink --bfile "${OUTROOT}/result" --chr "$chr" --make-bed --out "${OUTROOT}/result.chr${chr}.tmp"
    # Remove monomorphic variants (allele "0" in BIM) — Genotype Harmonizer crashes on these
    awk '$5 != "0" && $6 != "0" {print $2}' "${OUTROOT}/result.chr${chr}.tmp.bim" > "${OUTROOT}/result.chr${chr}.polymorphic.snps"
    plink --bfile "${OUTROOT}/result.chr${chr}.tmp" --extract "${OUTROOT}/result.chr${chr}.polymorphic.snps" --make-bed --output-chr M --out "${OUTROOT}/result.chr${chr}"
    rm -f "${OUTROOT}/result.chr${chr}.tmp".{bed,bim,fam,log,nosex} "${OUTROOT}/result.chr${chr}.polymorphic.snps"
    print_info "Split chromosome ${chr}: ${OUTROOT}/result.chr${chr}.bed"
done

print_success "Step 2b complete: split by chromosome"

# ################################################################################
# # Step 2b: Split by chromosome (PLINK 2 Optimized)
# ################################################################################

# print_step "Step 2b: Splitting by Chromosome (PLINK 2)"

# # Extract unique chromosome names from the .bim file
# # Filter the chromosome list to only include 1-22 and X
# # This ignores Y (24), XY (25), and MT (26)
# target_chroms=$(awk '{print $1}' "${OUTROOT}/result.bim" | sort -u | grep -E '^([1-9]|1[0-9]|2[0-2]|X)$')

# for chr in $target_chroms; do
#     print_info "Processing Chromosome ${chr}..."

#     # 1. Initial split using PLINK 2 (Fastest method)
#     # Using --allow-no-sex to keep the directory clean
#     plink2 --bfile "${OUTROOT}/result" \
#            --chr "$chr" \
#            --make-bed \
#            --allow-no-sex \
#            --out "${OUTROOT}/result.chr${chr}.tmp"

#     # 2. Identify Polymorphic SNPs (Alleles 5 and 6 must be A,C,T,G, not 0)
#     # This specifically targets the Genotype Harmonizer crash issue.
#     awk '$5 != "0" && $6 != "0" {print $2}' "${OUTROOT}/result.chr${chr}.tmp.bim" > "${OUTROOT}/result.chr${chr}.polymorphic.snps"

#     # 3. Create the final clean chromosome file
#     # We use --extract to keep only the SNPs identified in the previous step
#     plink2 --bfile "${OUTROOT}/result.chr${chr}.tmp" \
#            --extract "${OUTROOT}/result.chr${chr}.polymorphic.snps" \
#            --make-bed \
#            --allow-no-sex \
#            --out "${OUTROOT}/result.chr${chr}"

#     # 4. Clean up temporary files
#     # Only the final .bed, .bim, and .fam for each chromosome will remain
#     rm -f "${OUTROOT}/result.chr${chr}.tmp".{bed,bim,fam,log}
#     rm -f "${OUTROOT}/result.chr${chr}.polymorphic.snps"

#     print_info "Split complete for Chromosome ${chr}: ${OUTROOT}/result.chr${chr}.bed"
# done

# print_success "Step 2b complete: All chromosomes split and cleaned for Genotype Harmonizer."

################################################################################
# Step 3: Genotype Harmonizer
################################################################################
if [[ -z "$REF_PATH" ]]; then
    print_error "No reference path provided for Genotype Harmonizer."
    exit 1 # We exit with an error here because GH is a critical QC step. If the user doesn't provide a reference, it's safer to stop the pipeline than to continue without this QC.
else
    print_step "Step 3: Genotype Harmonizer QC1"

    GH_OUT="${OUTROOT}/2_GH"

    GH_ARGS=(--input "${OUTROOT}/result" --out "$GH_OUT" --ref-path "$REF_PATH")
    if [[ -n "$FASTA_REF" ]]; then
        GH_ARGS+=(--fasta-ref "$FASTA_REF")
    fi
    bash "${SCRIPT_DIR}/genotype_harmonizer.sh" "${GH_ARGS[@]}"

    print_success "Step 3 complete: $GH_OUT"
fi

################################################################################
# Step 4: Phasing with Eagle
################################################################################

if [[ -z "$REF_BCF" ]]; then
    print_warning "Skipping phasing step (no --ref-bcf provided)"
else
    print_step "Step 4: Phasing with Eagle"

    PHASE_OUT="${OUTROOT}/4_phase"

    # Use GH output if present, otherwise fall back to OUTROOT chromosome VCFs
    if ls "${OUTROOT}/2_GH"/*.chr*.vcf.gz &>/dev/null; then
        PHASE_INPUT="${OUTROOT}/2_GH"
    else
        PHASE_INPUT="$OUTROOT"
    fi

    bash "${SCRIPT_DIR}/phasing.sh" \
        --input-dir "$PHASE_INPUT" \
        --out "$PHASE_OUT" \
        --ref-bcf "$REF_BCF" \
        --threads "$THREADS"

    print_success "Step 4 complete: $PHASE_OUT"
fi

################################################################################
# Step 4b: Chromosome Naming Convention Check
################################################################################

MYMAP="${SCRIPT_DIR}/required_tools/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz"

if [[ -n "$REF_BCF" || -n "$REF_MSAV" ]]; then
    print_step "Step 4b: Chromosome Naming Convention Check"

    CHR_CHECK_ARGS=()
    [[ -d "${OUTROOT}/4_phase" ]] && CHR_CHECK_ARGS+=(--target-dir "${OUTROOT}/4_phase")
    [[ -n "$REF_BCF"           ]] && CHR_CHECK_ARGS+=(--ref-bcf "$REF_BCF")
    [[ -n "$REF_MSAV"          ]] && CHR_CHECK_ARGS+=(--ref-msav "$REF_MSAV")
    [[ -f "$MYMAP"             ]] && CHR_CHECK_ARGS+=(--genetic-map "$MYMAP")

    bash "${SCRIPT_DIR}/check_chr_naming.sh" "${CHR_CHECK_ARGS[@]}"

    print_success "Step 4b complete: chromosome naming consistent"
fi

################################################################################
# Step 5: Imputation with Minimac4
################################################################################

if [[ -z "$REF_MSAV" ]]; then
    print_warning "Skipping imputation step (no --ref-msav provided)"
else
    print_step "Step 5: Imputation with Minimac4"

    IMPUTE_OUT="${OUTROOT}/5_impute"

    bash "${SCRIPT_DIR}/imputation.sh" \
        --input-dir "${OUTROOT}/4_phase" \
        --out "$IMPUTE_OUT" \
        --ref-msav "$REF_MSAV" \
        --threads "$THREADS"

    print_success "Step 5 complete: $IMPUTE_OUT"
fi

################################################################################
# Step 6: Convert to GEN format
################################################################################

if [[ -d "${OUTROOT}/5_impute" ]]; then
    print_step "Step 6: Convert VCF to GEN format"

    GEN_OUT="${OUTROOT}/6_gen_format"

    bash "${SCRIPT_DIR}/convert_to_gen.sh" \
        --input-dir "${OUTROOT}/5_impute" \
        --out "$GEN_OUT"

    print_success "Step 6 complete: $GEN_OUT"
fi

################################################################################
# Done
################################################################################

print_step "Preprocessing Complete!"

print_info "Output files:"
print_info "  Cleaned data: $CLEANED_DATA"
print_info "  PLINK files:  ${OUTROOT}/result.bed/.bim/.fam"
if [[ -n "$REF_PATH" ]]; then
    print_info "  Harmonized:   ${OUTROOT}/2_GH/"
fi

echo
