#!/usr/bin/env bash

################################################################################
# Genotype Imputation Pipeline - Simple Bash Script
#
# A simplified version of the imputation pipeline that runs all steps
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

################################################################################
# Functions
################################################################################

print_header() {
    echo
    echo -e "${BLUE}################################################################################${NC}"
    echo -e "${BLUE}##                                                                            ##${NC}"
    echo -e "${BLUE}##              Genotype Imputation Pipeline - Bash Script                    ##${NC}"
    echo -e "${BLUE}##                          Torkamani Lab                                     ##${NC}"
    echo -e "${BLUE}##                                                                            ##${NC}"
    echo -e "${BLUE}################################################################################${NC}"
    echo
}

print_step() {
    echo
    echo -e "${GREEN}======================================================================${NC}"
    echo -e "${GREEN} $1${NC}"
    echo -e "${GREEN}======================================================================${NC}"
    echo
}

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

usage() {
    cat << EOF
Usage: bash imputation_pipeline.sh [OPTIONS]

Required Arguments:
  --vcf, -v PATH          Input file (VCF/VCF.gz or 23andMe .txt/.txt.gz/.zip)
  --out, -o PATH          Output directory for all results

Optional Arguments:
  --ref, -r TYPE          Reference panel: HRC or 1KG (default: 1KG)
  --start, -s NUM         First step to run (default: 0)
  --end, -e NUM           Last step to run (default: 7)
  --chr, -c NUM           Process single chromosome only (1-22 or X, default: all)
  --threads, -t NUM       Number of threads to use (default: 8)
  --temp PATH             Custom temporary directory
  --ref-path PATH         Custom reference panel path (VCF files)
  --ref-bcf PATH          Custom reference BCF path (for phasing)
  --ref-msav PATH        Custom reference MSAV path (for imputation)
  --wgs, -w               Enable WGS mode (variant down-sampling in step 3)
  --hwe, -h               Enable HWE filtering in step 4 (for non-mixed populations)
  --help                  Show this help message

Pipeline Steps:
  0. Check VCF genome build version
  1. Lift to GRCh37 (hg19)
  2. Quality control 1: LD-based fixes, strand corrections
  3. Quality control 2: missingness and HWE filtering
  4. Phasing with Eagle
  5. Imputation with Minimac4
  6. Convert VCF to GEN format
  7. Create ZIP archives

Note: 23andMe format (.txt) is automatically detected and converted to VCF

Example:
  bash imputation_pipeline.sh \\
    --vcf /path/to/input.vcf.gz \\
    --out /path/to/output \\
    --ref HRC \\
    --threads 16

EOF
}

################################################################################
# Parse Arguments
################################################################################

MYINPUT=""
OUTROOT=""
REF_MODE="1KG"
START_FROM=0
STOP_AFTER=7
CHR_ONLY=""
THREADS=8
TEMPDIR=""
REF_PATH_CUSTOM=""
REF_BCF_CUSTOM=""
REF_MSAV_CUSTOM=""
WGS_MODE="no"
HWE_MODE="no"

while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf|-v)
            MYINPUT="$2"
            shift 2
            ;;
        --out|-o)
            OUTROOT="$2"
            shift 2
            ;;
        --ref|-r)
            REF_MODE="$2"
            shift 2
            ;;
        --start|-s)
            START_FROM="$2"
            shift 2
            ;;
        --end|-e)
            STOP_AFTER="$2"
            shift 2
            ;;
        --chr|-c)
            CHR_ONLY="$2"
            shift 2
            ;;
        --threads|-t)
            THREADS="$2"
            shift 2
            ;;
        --temp)
            TEMPDIR="$2"
            shift 2
            ;;
        --ref-path)
            REF_PATH_CUSTOM="$2"
            shift 2
            ;;
        --ref-bcf)
            REF_BCF_CUSTOM="$2"
            shift 2
            ;;
        --ref-msav)
            REF_MSAV_CUSTOM="$2"
            shift 2
            ;;
        --wgs|-w)
            WGS_MODE="yes"
            shift
            ;;
        --hwe|-h)
            HWE_MODE="yes"
            shift
            ;;
        --help)
            print_header
            usage
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

################################################################################
# Validate Arguments
################################################################################

print_header

print_step "Validating Arguments"

if [[ -z "$MYINPUT" ]]; then
    print_error "Input VCF file is required (--vcf)"
    exit 1
fi

if [[ ! -f "$MYINPUT" ]]; then
    print_error "Input file does not exist: $MYINPUT"
    exit 1
fi

if [[ -z "$OUTROOT" ]]; then
    print_error "Output directory is required (--out)"
    exit 1
fi

if [[ "$REF_MODE" != "1KG" ]]; then
    print_error "Reference panel must be '1KG'"
    exit 1
fi

if [[ $START_FROM -gt $STOP_AFTER ]]; then
    print_error "Start step ($START_FROM) cannot be greater than end step ($STOP_AFTER)"
    exit 1
fi

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Set up paths
export PLINK="${SCRIPT_DIR}/required_tools/plink"
export PLINK2="${SCRIPT_DIR}/required_tools/plink2"
export GH="${SCRIPT_DIR}/required_tools/GenotypeHarmonizer/GenotypeHarmonizer.jar"
export EAGLE="${SCRIPT_DIR}/required_tools/Eagle_v2.4.1/eagle"
export LIFT="${SCRIPT_DIR}/required_tools/lift/LiftMap.py"
export CPATH="${SCRIPT_DIR}/required_tools/chainfiles"
export SPLIT_BY_ANCESTRY="${SCRIPT_DIR}/required_tools/split_by_ancestry/split_by_ancestry.R"

################################################################################
# Detect and convert 23andMe format
################################################################################

# Helper to read file (handles compressed and uncompressed)
read_file() {
    if [[ "$1" =~ \.(gz|zip)$ ]]; then
        zcat "$1" 2>/dev/null || gunzip -c "$1" 2>/dev/null
    else
        cat "$1"
    fi
}

# Check if input is 23andMe format (tab-delimited with rsid/chromosome/position/genotype)
if [[ "$MYINPUT" =~ \.(txt|tsv|gz|zip)$ ]] || read_file "$MYINPUT" | head -1 | grep -q "^#.*rsid"; then
    print_step "Detecting 23andMe format"

    # Verify 23andMe format
    if read_file "$MYINPUT" | head -100 | grep -v "^#" | head -1 | awk '{print NF}' | grep -q "^4$"; then
        print_info "23andMe format detected, converting to VCF..."

        CONVERT_DIR="${OUTROOT}/0_convert_23andme"
        mkdir -p "$CONVERT_DIR"

        INFILE=$(basename "$MYINPUT")
        PREFIX=$(echo "$INFILE" | sed -e 's/\.txt$//g' | sed -e 's/\.tsv$//g' | sed -e 's/\.gz$//g' | sed -e 's/\.zip$//g')

        # Convert to VCF
        (
            echo '##fileformat=VCFv4.2'
            echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"

            read_file "$MYINPUT" | grep -v "^#" | awk -F'\t' 'NF==4 {
                # 23andMe format: rsid, chromosome, position, genotype
                rsid = $1
                chr = $2
                pos = $3
                gt = $4

                # Filter invalid chromosomes and positions
                if (chr !~ /^[0-9]+$/ && chr !~ /^[XY]$/ && chr !~ /^MT$/) next
                if (pos !~ /^[0-9]+$/) next

                # Skip missing/invalid genotypes
                if (gt == "--" || gt == "" || length(gt) > 2) next
                if (gt !~ /^[ACGTDI-]+$/) next
                if (gt ~ /[ID]/) next

                # Convert MT to M
                if (chr == "MT") chr = "M"

                # Determine REF/ALT and genotype
                if (length(gt) == 1) {
                    ref = gt
                    alt = "."
                    genotype = "0/0"
                } else if (length(gt) == 2) {
                    a1 = substr(gt, 1, 1)
                    a2 = substr(gt, 2, 1)
                    if (a1 == a2) {
                        ref = a1
                        alt = "."
                        genotype = "0/0"
                    } else {
                        ref = a1
                        alt = a2
                        genotype = "0/1"
                    }
                } else {
                    next
                }

                print chr "\t" pos "\t" rsid "\t" ref "\t" alt "\t.\t.\t.\tGT\t" genotype
            }'
        ) | bgzip -c > "${CONVERT_DIR}/${PREFIX}.vcf.gz"

        tabix -p vcf "${CONVERT_DIR}/${PREFIX}.vcf.gz"

        # Update input to converted VCF
        MYINPUT="${CONVERT_DIR}/${PREFIX}.vcf.gz"
        print_success "Converted to VCF: $MYINPUT"
    else
        print_error "File appears to be 23andMe format but structure is invalid"
        print_error "Expected: rsid chromosome position genotype"
        exit 1
    fi
fi

# Parse input
INDIR=$(dirname "$MYINPUT")
INFILE=$(basename "$MYINPUT")
PREFIX=$(echo "$INFILE" | sed -e 's/\.vcf.gz$//g' | sed -e 's/\.vcf$//g')

# Create output directory
mkdir -p "$OUTROOT"

print_info "Input VCF: $MYINPUT"
print_info "Output directory: $OUTROOT"
print_info "Reference panel: $REF_MODE"
print_info "Threads: $THREADS"
print_info "Steps: $START_FROM to $STOP_AFTER"
if [[ -n "$CHR_ONLY" ]]; then
    print_info "Chromosome: $CHR_ONLY only"
else
    print_info "Chromosomes: 1-22 and X (if present)"
fi
print_info "WGS mode: $WGS_MODE"
print_info "HWE filtering: $HWE_MODE"

# Helper function to get chromosome range
# Detect which chromosomes are present in the VCF
detect_chromosomes() {
    local vcf_file="$1"
    # Get unique chromosome names from VCF, filter for autosomes (1-22) and X
    bcftools query -f '%CHROM\n' "$vcf_file" 2>/dev/null | \
        sort -u | \
        grep -E '^(chr)?(([1-9]|1[0-9]|2[0-2])|X|23)$' | \
        sed 's/^chr//' | \
        sed 's/^23$/X/' | \
        tr '\n' ' '
}

get_chr_range() {
    if [[ -n "$CHR_ONLY" ]]; then
        echo "$CHR_ONLY"
    else
        # Dynamic detection based on input file
        # If MYINPUT is defined and exists, detect; otherwise default to all
        if [[ -f "$MYINPUT" ]]; then
            detect_chromosomes "$MYINPUT"
        else
            echo "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"
        fi
    fi
}

################################################################################
# Step 0: Check VCF Build
################################################################################

if [[ $START_FROM -le 0 ]] && [[ $STOP_AFTER -ge 0 ]]; then
    print_step "Step 0: Checking VCF Genome Build"

    STEP0_OUT="${OUTROOT}/0_check_vcf_build"
    mkdir -p "$STEP0_OUT"

    CHECK_VCF_BUILD="${SCRIPT_DIR}/required_tools/check_vcf_build/check_vcf_build.R"

    if [[ ! -f "$CHECK_VCF_BUILD" ]]; then
        print_error "check_vcf_build.R not found at $CHECK_VCF_BUILD"
        exit 1
    fi

    print_info "Running genome build check..."
    Rscript "$CHECK_VCF_BUILD" "$MYINPUT" > "${STEP0_OUT}/${PREFIX}.BuildChecked"

    print_info "Build check completed: ${STEP0_OUT}/${PREFIX}.BuildChecked"

    BUILD_CHECK="${STEP0_OUT}/${PREFIX}.BuildChecked"
else
    BUILD_CHECK="${OUTROOT}/0_check_vcf_build/${PREFIX}.BuildChecked"
    if [[ ! -f "$BUILD_CHECK" ]]; then
        print_error "Build check file not found: $BUILD_CHECK"
        print_error "You may need to run step 0 first"
        exit 1
    fi
fi

################################################################################
# Step 1: Lift to GRCh37
################################################################################

if [[ $START_FROM -le 1 ]] && [[ $STOP_AFTER -ge 1 ]]; then
    print_step "Step 1: Lifting to GRCh37"

    STEP1_OUT="${OUTROOT}/1_lift"
    mkdir -p "$STEP1_OUT"
    mkdir -p "${STEP1_OUT}/temp"

    cd "$STEP1_OUT"

    TEMP_DIR="${TEMPDIR:-${STEP1_OUT}/temp}"
    export TEMP="$TEMP_DIR"

    print_info "Working directory: $STEP1_OUT"
    print_info "Temporary directory: $TEMP"

    # Determine lifted code from build check
    cfilename=$(grep "Use chain file" "$BUILD_CHECK" | tr -d ' ' | tr ':' '\t' | tr -d '"' | cut -f 2 | sed -e 's/->/ /g')
    checknone=$(grep "Use chain file" "$BUILD_CHECK" | grep "none" | wc -l)

    if [[ $checknone -gt 0 ]]; then
        print_warning "Input is already in GRCh37, processing by chromosome"
        LIFTED_CODE="GRCh37"

        # Still need to split by chromosome and convert to BED format
        CHR_RANGE=$(get_chr_range)
        for chr in $(eval echo $CHR_RANGE); do
            print_info "Processing chromosome $chr..."

            name="${PREFIX}.chr${chr}"

            # Extract chromosome
            bcftools view "$MYINPUT" -r "$chr" -Oz -o "${name}.sorted.vcf.gz"
            tabix -p vcf "${name}.sorted.vcf.gz"

            # Remove multi-allelic variants
            bcftools view "${name}.sorted.vcf.gz" -M 2 -m 2 | bcftools norm /dev/stdin -d both -Oz -o "${name}.sorted.bi.vcf.gz"
            tabix -p vcf "${name}.sorted.bi.vcf.gz"

            # Convert to BED format
            bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "${name}.sorted.bi.vcf.gz" > "${name}.sorted.bi.pos"
            $PLINK --vcf "${name}.sorted.bi.vcf.gz" --make-bed --a1-allele "${name}.sorted.bi.pos" 5 3 \
                   --biallelic-only strict --set-missing-var-ids @:#:\$1:\$2 --vcf-half-call missing \
                   --double-id --id-delim '_' --out "${PREFIX}.${LIFTED_CODE}.chr${chr}"
        done
    else
        print_info "Chain file to use: $cfilename"

        # Liftover process for each chromosome
        CHR_RANGE=$(get_chr_range)
        for chr in $(eval echo $CHR_RANGE); do
            print_info "Processing chromosome $chr..."

            name="${PREFIX}.chr${chr}"

            # Extract chromosome if needed
            if [[ ! -f "${name}.sorted.vcf.gz" ]]; then
                bcftools view "$MYINPUT" -r "$chr" -Oz -o "${name}.sorted.vcf.gz"
                tabix -p vcf "${name}.sorted.vcf.gz"
            fi

            # Remove multi-allelic variants
            bcftools view "${name}.sorted.vcf.gz" -M 2 -m 2 | bcftools norm /dev/stdin -d both -Oz -o "${name}.sorted.bi.vcf.gz"
            tabix -p vcf "${name}.sorted.bi.vcf.gz"

            # Convert to ped/map
            bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "${name}.sorted.bi.vcf.gz" > "${name}.sorted.bi.pos"
            $PLINK --vcf "${name}.sorted.bi.vcf.gz" --make-bed --a1-allele "${name}.sorted.bi.pos" 5 3 \
                   --biallelic-only strict --set-missing-var-ids @:#:\$1:\$2 --vcf-half-call missing \
                   --double-id --recode ped --id-delim '_' --out "$name"

            # Liftover
            liftname_base=$(echo "$cfilename" | sed 's/\.chain.*//g')

            nchains=$(echo "$cfilename" | tr ' ' '\n' | wc -l | awk '{print $1}')

            if [ "$nchains" -eq 1 ]; then
                liftname=$(echo $cfilename |  sed -e 's/\.chain.*//g')
                python "$LIFT" -p "${name}.ped" -m "${name}.map" -c "$CPATH/$cfilename" -o "${name}.lifted_${liftname}"
                LIFTED_CODE="lifted_${liftname}"
            elif [ "$nchains" -eq 2 ]; then
                first=$(echo "$cfilename" | tr ' ' '\n' | head -n 1)
                second=$(echo "$cfilename" | tr ' ' '\n' | tail -n 1)
                liftname1=$(echo "$first" | sed -e 's/\.chain.*//g')
                python "$LIFT" -m "${name}.map" -c "$CPATH/$first" -o "${name}.lifted_${liftname1}"
                liftname2=$(echo "$second" | sed -e 's/\.chain.*//g')
                python "$LIFT" -m "${name}.lifted_${liftname1}.map" -c "$CPATH/$second" -o "${name}.lifted_${liftname2}"
                LIFTED_CODE="lifted_${liftname2}"
            fi

            # Convert back to bed
            $PLINK --file "${name}.${LIFTED_CODE}" --make-bed --out "${name}.${LIFTED_CODE}.sorted.with_dup"

            # Remove duplicate markers
            cut -f 2 "${name}.${LIFTED_CODE}.sorted.with_dup.bim" | sort | uniq -d > "${name}.list_multi_a_markers.txt"
            ndup=$(wc -l "${name}.list_multi_a_markers.txt" | awk '{print $1}')

            dupflag=""
            if [ "$ndup" -ge 1 ]; then
                print_info "Found $ndup duplicate markers, removing them."
                dupflag="--exclude ${name}.list_multi_a_markers.txt"
            fi
            $PLINK --bfile "${name}.${LIFTED_CODE}.sorted.with_dup" $dupflag --make-bed --out "${PREFIX}.${LIFTED_CODE}.chr${chr}"
        done
    fi

    print_info "Liftover completed."
fi

# Get lifted code for subsequent steps
if [[ -z "${LIFTED_CODE:-}" ]]; then
    LIFTED_CODE=$(ls -1 "${OUTROOT}/1_lift/" | grep "${PREFIX}" | grep -E 'lifted|GRCh37' | head -1 | sed -e "s/${PREFIX}\.//" -e "s/\.chr.*//" || echo "GRCh37")
fi

################################################################################
# Rename chromosome 23 to X in PLINK files (check both chr23 and chrX files)
################################################################################

cd "${OUTROOT}/1_lift"

# Check both chr23 and chrX files for chromosome 23 in BIM
for chrname in chr23 chrX; do
    BED_FILE="${PREFIX}.${LIFTED_CODE}.${chrname}.bed"
    BIM_FILE="${PREFIX}.${LIFTED_CODE}.${chrname}.bim"

    if [[ -f "$BED_FILE" ]] && [[ -f "$BIM_FILE" ]]; then
        # Check if BIM file contains chromosome 23 (first column)
        if grep -q "^23[[:space:]]" "$BIM_FILE"; then
            print_info "Detected chromosome 23 in ${chrname} BIM file, renaming to X..."

            # Use awk to replace chromosome 23 with X in the BIM file
            awk 'BEGIN {OFS="\t"} {if ($1 == "23") $1 = "X"; print}' "$BIM_FILE" > "${PREFIX}.${LIFTED_CODE}.chrX.temp.bim"

            # Copy other files
            cp "${PREFIX}.${LIFTED_CODE}.${chrname}.bed" "${PREFIX}.${LIFTED_CODE}.chrX.temp.bed"
            cp "${PREFIX}.${LIFTED_CODE}.${chrname}.fam" "${PREFIX}.${LIFTED_CODE}.chrX.temp.fam"

            if [[ $? -eq 0 ]]; then
                print_success "Chromosome 23 successfully renamed to X in ${chrname} files"

                # Replace original files with corrected ones
                mv "${PREFIX}.${LIFTED_CODE}.chrX.temp.bed" "${PREFIX}.${LIFTED_CODE}.chrX.bed"
                mv "${PREFIX}.${LIFTED_CODE}.chrX.temp.bim" "${PREFIX}.${LIFTED_CODE}.chrX.bim"
                mv "${PREFIX}.${LIFTED_CODE}.chrX.temp.fam" "${PREFIX}.${LIFTED_CODE}.chrX.fam"

                # Remove old chr23 files if they exist and are different from chrX
                if [[ "$chrname" == "chr23" ]]; then
                    rm -f "${PREFIX}.${LIFTED_CODE}.chr23.bed" \
                          "${PREFIX}.${LIFTED_CODE}.chr23.bim" \
                          "${PREFIX}.${LIFTED_CODE}.chr23.fam" \
                          "${PREFIX}.${LIFTED_CODE}.chr23.log"
                fi

                rm -f chr_change.map "${PREFIX}.${LIFTED_CODE}.chrX.temp.log"
            else
                print_error "Failed to rename chromosome 23 to X in ${chrname} files"
                exit 1
            fi

            break  # Only process once
        fi
    fi
done

cd - > /dev/null


################################################################################
# Step 2: Genotype Harmonizer QC1
################################################################################

if [[ $START_FROM -le 2 ]] && [[ $STOP_AFTER -ge 2 ]]; then
    print_step "Step 2: Genotype Harmonizer and QC1"

    STEP2_OUT="${OUTROOT}/2_GH"
    mkdir -p "$STEP2_OUT"
    cd "$STEP2_OUT"

    MYINPUT_STEP2="${OUTROOT}/1_lift/${PREFIX}.${LIFTED_CODE}"

    # Reference path - use custom path if provided, otherwise environment variable, otherwise default
    if [[ -n "$REF_PATH_CUSTOM" ]]; then
        REF_PATH="$REF_PATH_CUSTOM"
        print_info "Using custom reference path: $REF_PATH"
    else
        print_error "Custom reference path is required for Genotype Harmonizer step (--ref-path)"
    fi

    FASTA_REF="/reference_prepared/fasta/human_g1k_v37.fasta.gz"

    if [ ! -f "$FASTA_REF" ]; then
        print_warning "Reference FASTA file not found at $FASTA_REF. Skipping fixref step."
        FASTA_REF=""
    fi


    print_info "Running Genotype Harmonizer on chromosomes $(get_chr_range | tr ' ' ',')..."

    # Function to run GH for one chromosome
    run_gh() {
        local chr=$1
        print_info "Processing chromosome ${chr}..."

        if [[ "$chr" == "X" || "$chr" -eq 23 ]]; then
            java -Xmx16g -jar "$GH" --keep --input "${MYINPUT_STEP2}.chr${chr}" \
                --ref "${REF_PATH}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz" \
                --inputType PLINK_BED --callRateFilter 0.90 --output "./${PREFIX}.${LIFTED_CODE}.GH.chr${chr}" \
                --debug
        else
            java -Xmx16g -jar "$GH" --keep --input "${MYINPUT_STEP2}.chr${chr}" \
                --ref "${REF_PATH}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" \
                --inputType PLINK_BED --callRateFilter 0.90 --output "./${PREFIX}.${LIFTED_CODE}.GH.chr${chr}"
        fi
    }

    # Process chromosomes in parallel
    export -f run_gh
    export -f print_info
    export GH PREFIX LIFTED_CODE MYINPUT_STEP2 REF_PATH

    CHR_RANGE=$(get_chr_range)
    # if command -v parallel &> /dev/null; then
    #     print_info "Running Genotype Harmonizer in parallel using $THREADS threads..."
    #     parallel -j "$THREADS" run_gh ::: $(eval echo $CHR_RANGE)
    # else
    print_info "Running Genotype Harmonizer sequentially..."
    for chr in $(eval echo $CHR_RANGE); do
        run_gh "$chr"
    done
    # fi

    print_info "Converting harmonized BED files to VCF and fixing reference alleles..."

    # Convert BED to VCF for each chromosome and run fixref
    fixref_and_convert() {
        local chr=$1
        local outname="${PREFIX}.${LIFTED_CODE}.GH.chr${chr}"
        local VCF_OUT_PREFIX="${outname}" # Define VCF_OUT_PREFIX based on outname
        if [[ -f "${outname}.bed" ]]; then
            print_info "Converting and fixing chromosome $chr..."
            $PLINK2 --bfile "${outname}" \
                --export vcf-4.2 bgz \
                --set-missing-var-ids @:#:\$1:\$2 \
                --out "${outname}.0"
            tabix -p vcf "${outname}.0.vcf.gz"

            if [ -n "$FASTA_REF" ]; then
                print_info "Running fixref for chr $chr"
                # Fix ref allele. The -- -f part is crucial for passing the fasta file to the plugin.
                # Added -m ref-alt to ensure the file is modified, not just stats reported.
                bcftools +fixref "$VCF_OUT_PREFIX.0.vcf.gz" -O z -o "$VCF_OUT_PREFIX.1.vcf.gz" -- -f "$FASTA_REF" -m ref-alt
                bcftools index -t "$VCF_OUT_PREFIX.1.vcf.gz"
                rm "$VCF_OUT_PREFIX.0.vcf.gz"*

                # Now rename the fixed file to the final name
                mv "$VCF_OUT_PREFIX.1.vcf.gz" "$VCF_OUT_PREFIX.vcf.gz"
                mv "$VCF_OUT_PREFIX.1.vcf.gz.tbi" "$VCF_OUT_PREFIX.vcf.gz.tbi"
            else
                # If no reference fasta is provided, just rename the file
                mv "${outname}.0.vcf.gz" "${outname}.vcf.gz"
                mv "${outname}.0.vcf.gz.tbi" "${outname}.vcf.gz.tbi"
            fi
        fi
    }

    export -f fixref_and_convert
    export -f print_info
    export PLINK2 PREFIX LIFTED_CODE THREADS FASTA_REF

    CHR_RANGE=$(get_chr_range)
    # if command -v parallel &> /dev/null; then
    #     parallel -j "$THREADS" fixref_and_convert ::: $(eval echo $CHR_RANGE)
    # else
    for chr in $(eval echo $CHR_RANGE); do
        fixref_and_convert "$chr"
    done
    # fi

    print_info "Genotype Harmonizer and fixref completed"
fi

################################################################################
# Step 3: QC2 (Single Sample Mode - No Ancestry Splitting)
################################################################################

if [[ $START_FROM -le 3 ]] && [[ $STOP_AFTER -ge 3 ]]; then
    print_step "Step 3: Quality Control 2 (Missingness and HWE)"

    STEP3_OUT="${OUTROOT}/3_QC2"
    mkdir -p "$STEP3_OUT"

    cd "$STEP3_OUT"

    # Input is directly from Step 2 (no ancestry splitting)
    MYINPUT_STEP3="${OUTROOT}/2_GH/${PREFIX}.${LIFTED_CODE}.GH"

    # Set filters (HWE disabled for single samples)
    HWEFLAG=""
    if [[ "$HWE_MODE" == "yes" ]]; then
        print_warning "HWE filtering disabled for single sample processing"
    fi

    GENOFLAG="--geno 0.1"
    MINDFLAG="--mind 0.05"

    print_info "Applying filters: $GENOFLAG $MINDFLAG"

    # Process each chromosome
    qc2_chr() {
        local chr=$1

        # FIX: Always regenerate BED from the VCF to ensure we capture the 'fixref' corrections from Step 2
        # Previous logic checked 'if [[ ! -f ... ]]' which used the stale/intermediate BED file.

        print_info "Converting corrected VCF to BED for chromosome $chr..."
        $PLINK2 --vcf "${MYINPUT_STEP3}.chr${chr}.vcf.gz" \
            --make-bed --id-delim '_' --out "${MYINPUT_STEP3}.chr${chr}"

        # Apply QC filters
        $PLINK2 --bfile "${MYINPUT_STEP3}.chr${chr}" \
            --make-bed $GENOFLAG $MINDFLAG \
            --out "${PREFIX}.${LIFTED_CODE}.GH.QC2.chr${chr}"
    }

    export -f qc2_chr
    export -f print_info
    export PLINK2 PREFIX LIFTED_CODE GENOFLAG MINDFLAG MYINPUT_STEP3

    CHR_RANGE=$(get_chr_range)
    if command -v parallel &> /dev/null; then
        parallel -j "$THREADS" qc2_chr ::: $(eval echo $CHR_RANGE)
    else
        for chr in $(eval echo $CHR_RANGE); do
            qc2_chr "$chr"
        done
    fi

    print_info "QC2 completed"
fi

################################################################################
# Step 4: Phasing (Single Sample Mode)
################################################################################

if [[ $START_FROM -le 4 ]] && [[ $STOP_AFTER -ge 4 ]]; then
    print_step "Step 4: Phasing with Eagle"

    STEP4_OUT="${OUTROOT}/4_phase"
    mkdir -p "$STEP4_OUT"
    cd "$STEP4_OUT"

    MYMAP="${SCRIPT_DIR}/required_tools/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz"

    # Set reference - use custom BCF path if provided
    if [[ -n "$REF_BCF_CUSTOM" ]]; then
        REF_BASE="$REF_BCF_CUSTOM"
        print_info "Using custom BCF reference path: $REF_BASE"
    else
        print_error "Custom BCF reference path is required for phasing step (--ref-bcf)"
        exit 1
    fi

    # Phase each chromosome
    phase_chr() {
        local chr=$1

        MYINPUT_STEP4="${OUTROOT}/3_QC2/${PREFIX}.${LIFTED_CODE}.GH.QC2.chr${chr}"

        if [[ ! -f "${MYINPUT_STEP4}.bed" ]]; then
            print_warning "No BED file found for chromosome $chr, skipping..."
            return
        fi

        INPREFIX="${PREFIX}.${LIFTED_CODE}.GH.QC2.chr${chr}"

        print_info "Phasing chromosome $chr..."

        # Convert to VCF
        $PLINK2 --bfile "$MYINPUT_STEP4" --export vcf-4.2 bgz \
            --set-missing-var-ids @:#:\$1:\$2 --out "$INPREFIX"

        tabix -f -p vcf "${INPREFIX}.vcf.gz"

        # Set reference file (handle X chromosome)
        if [[ "$chr" == "X" ]]; then
            MYREF="${REF_BASE}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.bcf"
        else
            MYREF="${REF_BASE}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf"
        fi

        # Run Eagle
        $EAGLE --vcfTarget="${INPREFIX}.vcf.gz" \
            --vcfRef="$MYREF" \
            --noImpMissing \
            --geneticMapFile="$MYMAP" \
            --Kpbwt=100000 --numThreads="$THREADS" \
            --chrom="$chr" --allowRefAltSwap \
            --outPrefix="${INPREFIX}.phased"

        # Index the phased output for minimac4
        tabix -p vcf "${INPREFIX}.phased.vcf.gz"
    }

    export -f phase_chr
    export -f print_info
    export -f print_warning
    export PLINK2 EAGLE PREFIX LIFTED_CODE OUTROOT REF_BASE MYMAP SCRIPT_DIR THREADS

    CHR_RANGE=$(get_chr_range)
    if command -v parallel &> /dev/null; then
        # Ensure at least 1 job is run (integer division 1/4 = 0)
        N_JOBS=$(($THREADS / 4))
        if [[ $N_JOBS -lt 1 ]]; then N_JOBS=1; fi
        parallel -j $N_JOBS phase_chr ::: $(eval echo $CHR_RANGE)
    else
        for chr in $(eval echo $CHR_RANGE); do
            phase_chr "$chr"
        done
    fi

    print_info "Phasing completed"
fi

################################################################################
# Step 5: Imputation (Single Sample Mode)
################################################################################

if [[ $START_FROM -le 5 ]] && [[ $STOP_AFTER -ge 5 ]]; then
    print_step "Step 5: Imputation with Minimac4"

    STEP5_OUT="${OUTROOT}/5_impute"
    mkdir -p "$STEP5_OUT"
    cd "$STEP5_OUT"

    # Set reference - use custom MSAV path if provided
    if [[ -n "$REF_MSAV_CUSTOM" ]]; then
        REF_BASE="$REF_MSAV_CUSTOM"
        print_info "Using custom MSAV reference path: $REF_BASE"
    else
        print_error "Custom MSAV reference path is required for imputation step (--ref-msav)"
        exit 1
    fi

    # Impute each chromosome
    impute_chr() {
        local chr=$1

        MYINPUT_STEP5="${OUTROOT}/4_phase/${PREFIX}.${LIFTED_CODE}.GH.QC2.chr${chr}.phased.vcf.gz"

        if [[ ! -f "$MYINPUT_STEP5" ]]; then
            print_warning "No phased file found for chromosome $chr, skipping..."
            return
        fi

        print_info "Imputing chromosome $chr..."

        # Set reference file (handle X chromosome)
        if [[ "$chr" == "X" ]]; then
            for j in PAR1 PAR2 nonPAR; do
                MYREF="${REF_BASE}/ALL.chrX_${j}.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.msav"
                # Run Minimac4
                minimac4 "$MYREF" "$MYINPUT_STEP5" \
                    --output "${PREFIX}.imputed.chr${chr}_${j}.dose.vcf.gz" \
                    --output-format vcf.gz \
                    --min-ratio 0.001 --chunk 30000000 --overlap 3000000 --threads "$THREADS"
            done
        else
            MYREF="${REF_BASE}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.msav"
            # Run Minimac4
            minimac4 "$MYREF" "$MYINPUT_STEP5" \
                --output "${PREFIX}.imputed.chr${chr}.dose.vcf.gz" \
                --output-format vcf.gz \
                --min-ratio 0.001 --chunk 30000000 --overlap 3000000 --threads "$THREADS"
        fi
    }

    export -f impute_chr
    export -f print_info
    export -f print_warning
    export PREFIX LIFTED_CODE OUTROOT REF_BASE THREADS

    CHR_RANGE=$(get_chr_range)
    for chr in $(eval echo $CHR_RANGE); do
        impute_chr "$chr"
    done

    print_info "Imputation completed"
fi

################################################################################
# Step 6: Convert VCF to GEN format (Single Sample Mode)
################################################################################

if [[ $START_FROM -le 6 ]] && [[ $STOP_AFTER -ge 6 ]]; then
    print_step "Step 6: Converting VCF to GEN format with genotype probabilities"

    STEP6_GEN_OUT="${OUTROOT}/6_gen_format"
    mkdir -p "$STEP6_GEN_OUT"
    cd "$STEP6_GEN_OUT"

    STEP5_OUT="${OUTROOT}/5_impute"

    # Convert each chromosome
    convert_to_gen() {
        local chr=$1

        # Handle X chromosome with PAR regions
        if [[ "$chr" == "X" ]]; then
            print_info "Converting chromosome X (merging PAR1, PAR2, and nonPAR) to GEN format..."

            # Convert each region to GEN format
            for j in PAR1 PAR2 nonPAR; do
                VCF_INPUT="${STEP5_OUT}/${PREFIX}.imputed.chr${chr}_${j}.dose.vcf.gz"

                if [[ ! -f "$VCF_INPUT" ]]; then
                    print_warning "Missing ${chr}_${j} VCF file, skipping..."
                    continue
                fi

                GEN_OUTPUT_TEMP="${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}_${j}.temp"

                # Convert to GEN format
                bcftools convert "$VCF_INPUT" \
                    --gensample "${GEN_OUTPUT_TEMP}" \
                    --tag GP
            done

            # Merge the three regions into one file (in genomic order: PAR1, nonPAR, PAR2)
            GEN_OUTPUT="${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}"

            # Check which files exist
            PAR1_EXISTS=false
            NONPAR_EXISTS=false
            PAR2_EXISTS=false

            [[ -f "${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}_PAR1.temp.gen" ]] && PAR1_EXISTS=true
            [[ -f "${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}_nonPAR.temp.gen" ]] && NONPAR_EXISTS=true
            [[ -f "${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}_PAR2.temp.gen" ]] && PAR2_EXISTS=true

            # Concatenate GEN files in order (PAR1 -> nonPAR -> PAR2)
            > "${GEN_OUTPUT}.gen"  # Create empty file

            if $PAR1_EXISTS; then
                cat "${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}_PAR1.temp.gen" >> "${GEN_OUTPUT}.gen"
                rm "${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}_PAR1.temp.gen"
            fi

            if $NONPAR_EXISTS; then
                cat "${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}_nonPAR.temp.gen" >> "${GEN_OUTPUT}.gen"
                rm "${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}_nonPAR.temp.gen"
            fi

            if $PAR2_EXISTS; then
                cat "${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}_PAR2.temp.gen" >> "${GEN_OUTPUT}.gen"
                rm "${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}_PAR2.temp.gen"
            fi

            # Copy the sample file (they should all be identical, so just use the first one found)
            for j in PAR1 nonPAR PAR2; do
                if [[ -f "${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}_${j}.temp.sample" ]]; then
                    cp "${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}_${j}.temp.sample" "${GEN_OUTPUT}.sample"
                    rm "${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}_${j}.temp.sample"
                    break
                fi
            done

            # Compress the merged GEN file
            gzip -f "${GEN_OUTPUT}.gen"

            print_success "Converted chromosome X (merged PAR regions) to GEN format"
        else
            VCF_INPUT="${STEP5_OUT}/${PREFIX}.imputed.chr${chr}.dose.vcf.gz"

            if [[ ! -f "$VCF_INPUT" ]]; then
                print_warning "No imputed VCF for chromosome $chr, skipping..."
                return
            fi

            print_info "Converting chromosome $chr to GEN format..."

            GEN_OUTPUT="${STEP6_GEN_OUT}/${PREFIX}.imputed.chr${chr}"

            # Use bcftools to convert VCF to GEN format with genotype probabilities
            bcftools convert "$VCF_INPUT" \
                --gensample "${GEN_OUTPUT}" \
                --tag GP

            # Compress the output files
            gzip -f "${GEN_OUTPUT}.gen"

            print_success "Converted chromosome $chr to GEN format"
        fi
    }

    export -f convert_to_gen
    export -f print_info
    export -f print_success
    export -f print_warning
    export PREFIX STEP5_OUT STEP6_GEN_OUT

    CHR_RANGE=$(get_chr_range)
    for chr in $(eval echo $CHR_RANGE); do
        convert_to_gen "$chr"
    done

    print_info "GEN format conversion completed"
fi

################################################################################
# Step 7: Create ZIP archives (Single Sample Mode)
################################################################################

if [[ $START_FROM -le 7 ]] && [[ $STOP_AFTER -ge 7 ]]; then
    print_step "Step 7: Creating ZIP archives of results"

    STEP7_OUT="${OUTROOT}/7_archives"
    mkdir -p "$STEP7_OUT"

    STEP5_OUT="${OUTROOT}/5_impute"
    STEP6_GEN_OUT="${OUTROOT}/6_gen_format"

    # Create VCF archive
    if [[ -d "${STEP5_OUT}" ]]; then
        VCF_FILES=$(find "${STEP5_OUT}" -name "*.vcf.gz" 2>/dev/null)

        if [[ -n "$VCF_FILES" ]]; then
            print_info "Archiving VCF files..."
            cd "${STEP5_OUT}"
            zip -r "${STEP7_OUT}/${PREFIX}_imputed_vcf.zip" *.vcf.gz *.vcf.gz.tbi 2>/dev/null || \
                zip -r "${STEP7_OUT}/${PREFIX}_imputed_vcf.zip" *.vcf.gz
            cd - > /dev/null
            print_success "VCF archive created: ${PREFIX}_imputed_vcf.zip"
        else
            print_warning "No VCF files found"
        fi
    fi

    # Create GEN archive
    if [[ -d "${STEP6_GEN_OUT}" ]]; then
        GEN_FILES=$(find "${STEP6_GEN_OUT}" -name "*.gen.gz" 2>/dev/null)

        if [[ -n "$GEN_FILES" ]]; then
            print_info "Archiving GEN files..."
            cd "${STEP6_GEN_OUT}"
            zip -r "${STEP7_OUT}/${PREFIX}_imputed_gen.zip" *.gen.gz *.sample 2>/dev/null || \
                zip -r "${STEP7_OUT}/${PREFIX}_imputed_gen.zip" *.gen.gz
            cd - > /dev/null
            print_success "GEN archive created: ${PREFIX}_imputed_gen.zip"
        else
            print_warning "No GEN files found"
        fi
    fi

    print_info "Archive creation completed"
    print_info "Archives saved to: $STEP7_OUT"
fi

################################################################################
# Complete
################################################################################

print_step "Pipeline Complete!"

print_info "All steps completed successfully"
print_info "Results are in: $OUTROOT"
print_info ""
print_info "Output structure (single sample mode - no ancestry splitting):"
print_info "  0_check_vcf_build/ - Genome build detection"
print_info "  1_lift/            - Lifted to GRCh37"
print_info "  2_GH/              - Harmonized genotypes"
print_info "  3_QC2/             - QC filtered data"
print_info "  4_phase/           - Phased haplotypes"
print_info "  5_impute/          - Final imputed genotypes (VCF)"
print_info "  6_gen_format/      - GEN format files"
print_info "  7_archives/        - ZIP archives"

echo