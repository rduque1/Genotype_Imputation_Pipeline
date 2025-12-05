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
  --end, -e NUM           Last step to run (default: 6)
  --chr, -c NUM           Process single chromosome only (1-22 or X, default: all)
  --threads, -t NUM       Number of threads to use (default: 8)
  --temp PATH             Custom temporary directory
  --ref-path PATH         Custom reference panel path (VCF files)
  --ref-bcf PATH          Custom reference BCF path (for phasing)
  --ref-m3vcf PATH        Custom reference M3VCF path (for imputation)
  --wgs, -w               Enable WGS mode (variant down-sampling in step 3)
  --hwe, -h               Enable HWE filtering in step 4 (for non-mixed populations)
  --help                  Show this help message

Pipeline Steps:
  0. Check VCF genome build version
  1. Lift to GRCh37 (hg19)
  2. Quality control 1: LD-based fixes, strand corrections
  3. Ancestry analysis and sample splitting
  4. Quality control 2: missingness and HWE filtering
  5. Phasing with Eagle
  6. Imputation with Minimac4

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
STOP_AFTER=6
CHR_ONLY=""
THREADS=8
TEMPDIR=""
REF_PATH_CUSTOM=""
REF_BCF_CUSTOM=""
REF_M3VCF_CUSTOM=""
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
        --ref-m3vcf)
            REF_M3VCF_CUSTOM="$2"
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
        # Default to autosomes + X chromosome (output as space-separated list for iteration)
        # echo "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"
        echo "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
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
                --inputType PLINK_BED --callRateFilter 0.90 --output "./${PREFIX}.${LIFTED_CODE}.GH.chr${chr}"
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
# Step 3: Ancestry Analysis
################################################################################

if [[ $START_FROM -le 3 ]] && [[ $STOP_AFTER -ge 3 ]]; then
    print_step "Step 3: Ancestry Analysis and Sample Splitting"

    STEP3_OUT="${OUTROOT}/3_ancestry"
    mkdir -p "$STEP3_OUT"

    OUTSUBDIR="${PREFIX}"
    mkdir -p "${STEP3_OUT}/${OUTSUBDIR}"
    cd "${STEP3_OUT}/${OUTSUBDIR}"

    MYINPUT_STEP3="${OUTROOT}/2_GH/${PREFIX}.${LIFTED_CODE}.GH"
    REF_POP="${SCRIPT_DIR}/required_tools/1000G_P3_super_pop.pop"

    # Check number of samples
    NUM_SAMPLES=$(bcftools query -l "${MYINPUT_STEP3}.chr$(eval echo $(get_chr_range) | awk '{print $1}').vcf.gz" | wc -l)

    if [[ $NUM_SAMPLES -lt 50 ]]; then
        print_warning "Only $NUM_SAMPLES sample(s) detected. Skipping LD pruning (requires 50+ samples)."
        print_info "Copying VCF files without pruning..."

        CHR_RANGE=$(get_chr_range)
        for chr in $(eval echo $CHR_RANGE); do
            print_info "Processing chromosome $chr..."
            cp "${MYINPUT_STEP3}.chr${chr}.vcf.gz" "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.pruned.vcf.gz"
            tabix -p vcf "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.pruned.vcf.gz" 2>/dev/null || true
        done
    else
        print_info "Pruning markers with LD threshold 0.05..."

        # Prune function
        prune_chr() {
            local chr=$1
            print_info "Pruning chromosome $chr..."

            if [[ "$WGS_MODE" == "yes" ]]; then
                bcftools view -R "${SCRIPT_DIR}/required_tools/chr_pos_23andme.txt" \
                    "${MYINPUT_STEP3}.chr${chr}.vcf.gz" -Ov -o "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.23andMe_pos.vcf"
                $PLINK2 --vcf "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.23andMe_pos.vcf" \
                    --indep-pairwise 100 10 0.05 --out "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}"
            else
                $PLINK2 --vcf "${MYINPUT_STEP3}.chr${chr}.vcf.gz" \
                    --indep-pairwise 100 10 0.05 --out "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}"
            fi

            vcftools --gzvcf "${MYINPUT_STEP3}.chr${chr}.vcf.gz" \
                --snps "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.prune.in" \
                --recode --recode-INFO-all --out "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.pruned"

            bgzip -c "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.pruned.recode.vcf" > "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.pruned.vcf.gz"
            tabix -p vcf "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.pruned.vcf.gz"
            rm "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.pruned.recode.vcf"
        }

        export -f prune_chr
        export -f print_info
        export PREFIX LIFTED_CODE MYINPUT_STEP3 PLINK2 SCRIPT_DIR WGS_MODE

        CHR_RANGE=$(get_chr_range)
        if command -v parallel &> /dev/null; then
            parallel -j "$THREADS" prune_chr ::: $(eval echo $CHR_RANGE)
        else
            for chr in $(eval echo $CHR_RANGE); do
                prune_chr "$chr"
            done
        fi
    fi

    print_info "Intersecting with 1000 Genomes reference..."

    # Get reference path from custom or default
    if [[ -n "$REF_PATH_CUSTOM" ]]; then
        REF_VCF_PATH="$REF_PATH_CUSTOM"
    else
        print_error "Custom reference path is required for ancestry analysis step (--ref-path)"
        exit 1
    fi

    # Intersect with reference and MERGE them together (like archived pipeline)
    CHR_RANGE=$(get_chr_range)
    FIRST_CHR=$(eval echo $CHR_RANGE | awk '{print $1}')

    for chr in $(eval echo $CHR_RANGE); do
        print_info "Intersecting chromosome $chr with 1KG reference..."
        bcftools isec \
            "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.pruned.vcf.gz" \
            "${REF_VCF_PATH}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" \
            -p "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}_tmp" -n =2 -w 1,2 -Oz

        # Save reference sample IDs from the first chromosome before cleanup
        if [[ "$chr" == "$FIRST_CHR" ]]; then
            print_info "Saving reference sample IDs from chromosome $chr..."
            bcftools query -l "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}_tmp/0001.vcf.gz" > 0001.txt
        fi

        # Merge user samples with 1KG reference samples
        bcftools merge "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}_tmp/0000.vcf.gz" \
            "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}_tmp/0001.vcf.gz" \
            -Oz -o "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.pruned.intersect1KG.vcf.gz"

        tabix -f -p vcf "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.pruned.intersect1KG.vcf.gz"
        rm -rf "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}_tmp"
    done

    print_info "Concatenating all chromosomes..."
    ls -v "${PREFIX}.${LIFTED_CODE}.GH.chr"*.pruned.intersect1KG.vcf.gz > "${PREFIX}.pruned.intersect1KG.vcflist"
    bcftools concat -f "${PREFIX}.pruned.intersect1KG.vcflist" -Ov -o "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.vcf"

    # Rename 1000G sample IDs to have FID_IID format (required for PLINK id-delim)
    print_info "Renaming reference sample IDs..."
    # Reference sample IDs were already saved in 0001.txt during the intersection loop
    if [[ ! -f "0001.txt" ]]; then
        print_error "Reference sample IDs file (0001.txt) not found. This should have been created during intersection."
        exit 1
    fi

    while read line; do printf "${line}\t${line}_${line}\n"; done < 0001.txt > 0001.rename
    bcftools reheader "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.vcf" -s 0001.rename > "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.id-delim.vcf"
    rm "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.vcf"

    print_info "Converting merged intersection to BED format..."
    $PLINK2 --vcf "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.id-delim.vcf" \
        --make-bed --id-delim '_' --out "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG"
    rm "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.id-delim.vcf"

    # Count user samples (samples before the 1KG reference samples)
    NUM_USER_SAMPLES=$(bcftools query -l "${PREFIX}.${LIFTED_CODE}.GH.chr$(eval echo $(get_chr_range) | awk '{print $1}').pruned.vcf.gz" | wc -l | tr -d ' ')

    print_info "Creating population file for ADMIXTURE..."
    (
        for i in $(seq 1 $NUM_USER_SAMPLES); do echo "-"; done
        cat "$REF_POP"
    ) > "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.pop"

    print_info "Running ADMIXTURE with supervised learning..."
    admixture --supervised "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.bed" 5 -j${THREADS}

    print_info "Extracting sample IDs and ancestry results..."
    # Get subject IDs with ancestry results (both user and reference)
    cat "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.fam" | awk '{print $1 "_" $2}' > "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.subjectIDs"
    paste "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.subjectIDs" "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.5.Q" | tr '\t' ' ' > "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.5.Q.IDs"

    # Get user sample IDs (FID_IID format) from original pruned VCFs
    bcftools query -l "${PREFIX}.${LIFTED_CODE}.GH.chr$(eval echo $(get_chr_range) | awk '{print $1}').pruned.vcf.gz" | \
        tr '_' '\t' | awk '{print$1"_"$2"\t"$1"\t"$2}' > "${PREFIX}.${LIFTED_CODE}.GH.pruned.subjectIDs"

    # Create readable version with headers
    (
        echo -e "SampleID\tAFR\tAMR\tEAS\tEUR\tSAS"
        head -n $NUM_USER_SAMPLES "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.5.Q.IDs"
    ) > "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.5.Q.readable"

    print_info "Ancestry proportions saved to: ${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.5.Q.readable"
    print_info "Splitting samples by ancestry (threshold 0.95)..."
    Rscript "$SPLIT_BY_ANCESTRY" \
        "${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.5.Q.IDs" \
        "${PREFIX}.${LIFTED_CODE}.GH.pruned.subjectIDs" \
        0.95

    # Now split the original harmonized genotype files by ancestry
    print_info "Creating ancestry-specific BED files from original genotypes..."

    # First convert original VCFs to BED format per chromosome
    for chr in $(eval echo $(get_chr_range)); do
        if [[ -f "${MYINPUT_STEP3}.chr${chr}.vcf.gz" ]]; then
            $PLINK2 --vcf "${MYINPUT_STEP3}.chr${chr}.vcf.gz" \
                --make-bed --id-delim '_' --out "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}"
        fi
    done

    # Now split by ancestry using the .ids files
    CHR_RANGE=$(get_chr_range)
    for anc in 1 2 3 4 5 mixed; do
        IDS_FILE="${PREFIX}.${LIFTED_CODE}.GH.pruned.intersect1KG.5.Q.IDs.${anc}.ids"

        if [[ -f "$IDS_FILE" ]]; then
            print_info "Splitting ancestry group $anc..."

            for chr in $(eval echo $CHR_RANGE); do
                if [[ -f "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}.bed" ]]; then
                    $PLINK2 --bfile "${PREFIX}.${LIFTED_CODE}.GH.chr${chr}" \
                        --keep "$IDS_FILE" \
                        --make-bed \
                        --out "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}"
                fi
            done
        else
            print_info "No samples found for ancestry group $anc"
        fi
    done

    print_info "Ancestry analysis completed"
fi

################################################################################
# Step 4: Split and QC2
################################################################################

if [[ $START_FROM -le 4 ]] && [[ $STOP_AFTER -ge 4 ]]; then
    print_step "Step 4: Quality Control 2 (Missingness and HWE)"

    STEP4_OUT="${OUTROOT}/4_split_QC2"
    mkdir -p "$STEP4_OUT"

    OUTSUBDIR="${PREFIX}"
    mkdir -p "${STEP4_OUT}/${OUTSUBDIR}"

    # Process each ancestry
    for anc in 1 2 3 4 5 mixed; do
        print_info "Processing ancestry group: $anc"

        MYINPUT_STEP4="${OUTROOT}/3_ancestry/${OUTSUBDIR}/${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}"

        if [[ ! -f "${MYINPUT_STEP4}.chr$(eval echo $(get_chr_range) | awk '{print $1}').bed" ]]; then
            print_warning "Ancestry group $anc has no files, skipping..."
            continue
        fi

        cd "${STEP4_OUT}/${OUTSUBDIR}"

        # Set filters
        HWEFLAG=""
        if [[ "$HWE_MODE" == "yes" ]] && [[ "$anc" != "mixed" ]]; then
            HWEFLAG="--hwe 1e-10"
        fi

        GENOFLAG="--geno 0.1"
        MINDFLAG="--mind 0.05"

        print_info "Applying filters: $HWEFLAG $GENOFLAG $MINDFLAG"

        # Process each chromosome
        qc2_chr() {
            local chr=$1
            local anc=$2
            $PLINK2 --bfile "${MYINPUT_STEP4}.chr${chr}" \
                --make-bed $HWEFLAG $GENOFLAG $MINDFLAG \
                --out "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}.tmp"
        }

        export -f qc2_chr
        export PLINK2 PREFIX LIFTED_CODE HWEFLAG GENOFLAG MINDFLAG MYINPUT_STEP4

        CHR_RANGE=$(get_chr_range)
        if command -v parallel &> /dev/null; then
            parallel -j "$THREADS" qc2_chr ::: $(eval echo $CHR_RANGE) ::: "$anc"
        else
            for chr in $(eval echo $CHR_RANGE); do
                qc2_chr "$chr" "$anc"
            done
        fi

        # Remove samples with missingness issues across all chromosomes
        cat ${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr*.tmp.mindrem > "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.tmp.mindrem.id" 2>/dev/null || true

        if [[ -s "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.tmp.mindrem.id" ]]; then
            print_info "Removing samples that failed missingness checks..."
            sort -u "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.tmp.mindrem.id" > "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.tmp.unique.mindrem.id"

            CHR_RANGE=$(get_chr_range)
            for chr in $(eval echo $CHR_RANGE); do
                if [[ -f "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}.tmp.bed" ]]; then
                    $PLINK2 --bfile "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}.tmp" \
                        --remove "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.tmp.unique.mindrem.id" \
                        --make-bed --out "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}"
                fi
            done
            rm -f ${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr*.tmp.*
        else
            print_info "No samples failed missingness checks."
            CHR_RANGE=$(get_chr_range)
            for chr in $(eval echo $CHR_RANGE); do
                if [[ -f "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}.tmp.bed" ]]; then
                    mv "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}.tmp.bed" "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}.bed"
                    mv "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}.tmp.bim" "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}.bim"
                    mv "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}.tmp.fam" "${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}.fam"
                fi
            done
        fi
        # Clean up intermediate files
        rm -f ${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.tmp.*.id


    done

    print_info "QC2 completed"
fi

################################################################################
# Step 5: Phasing
################################################################################

if [[ $START_FROM -le 5 ]] && [[ $STOP_AFTER -ge 5 ]]; then
    print_step "Step 5: Phasing with Eagle"

    STEP5_OUT="${OUTROOT}/5_phase"
    mkdir -p "$STEP5_OUT"

    OUTSUBDIR="${PREFIX}"
    mkdir -p "${STEP5_OUT}/${OUTSUBDIR}"
    cd "${STEP5_OUT}/${OUTSUBDIR}"

    MYMAP="${SCRIPT_DIR}/required_tools/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz"

    # Set reference - use custom BCF path if provided
    if [[ -n "$REF_BCF_CUSTOM" ]]; then
        REF_BASE="$REF_BCF_CUSTOM"
        print_info "Using custom BCF reference path: $REF_BASE"
    else
        print_error "Custom BCF reference path is required for phasing step (--ref-bcf-path)"
        exit 1
    fi

    # Phase each ancestry and chromosome
    phase_chr() {
        local chr=$1
        local anc=$2

        MYINPUT_STEP5="${OUTROOT}/4_split_QC2/${PREFIX}/${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}"

        if [[ ! -f "${MYINPUT_STEP5}.bed" ]]; then
            return
        fi

        INPREFIX=$(basename "$MYINPUT_STEP5")

        print_info "Phasing ancestry $anc, chromosome $chr..."

        # Convert to VCF
        $PLINK2 --bfile "$MYINPUT_STEP5" --export vcf-4.2 bgz \
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
    export PLINK2 EAGLE PREFIX LIFTED_CODE OUTROOT REF_MODE REF_BASE MYMAP SCRIPT_DIR THREADS

    # Process all ancestry groups and chromosomes
    for anc in 1 2 3 4 5 mixed; do
        print_info "Phasing ancestry group: $anc"

        CHR_RANGE=$(get_chr_range)
        if command -v parallel &> /dev/null; then
            parallel -j $(($THREADS / 4)) phase_chr ::: $(eval echo $CHR_RANGE) ::: "$anc"
        else
            for chr in $(eval echo $CHR_RANGE); do
                phase_chr "$chr" "$anc"
            done
        fi
    done

    print_info "Phasing completed"
fi

################################################################################
# Step 6: Imputation
################################################################################

if [[ $START_FROM -le 6 ]] && [[ $STOP_AFTER -ge 6 ]]; then
    print_step "Step 6: Imputation with Minimac4"

    STEP6_OUT="${OUTROOT}/6_impute_${REF_MODE}"
    mkdir -p "$STEP6_OUT"

    # Set reference - use custom M3VCF path if provided
    if [[ -n "$REF_M3VCF_CUSTOM" ]]; then
        REF_BASE="$REF_M3VCF_CUSTOM"
        print_info "Using custom M3VCF reference path: $REF_BASE"
    else
        print_error "Custom M3VCF reference path is required for imputation step (--ref-m3vcf-path)"
        exit 1
    fi

    # Impute each ancestry and chromosome
    impute_chr() {
        local chr=$1
        local anc=$2

        MYINPUT_STEP6="${OUTROOT}/5_phase/${PREFIX}/${PREFIX}.${LIFTED_CODE}.GH.ancestry-${anc}.chr${chr}.phased.vcf.gz"

        if [[ ! -f "$MYINPUT_STEP6" ]]; then # No phased file for this ancestry/chromosome
            return
        fi

        OUTSUBDIR="${PREFIX}/${anc}"
        mkdir -p "${STEP6_OUT}/${OUTSUBDIR}"
        cd "${STEP6_OUT}/${OUTSUBDIR}"

        print_info "Imputing ancestry $anc, chromosome $chr..."

        # Set reference file (handle X chromosome)
        if [[ "$chr" == "X" ]]; then
            for j in {PAR1,PAR2,nonPAR}; do
                MYREF="${REF_BASE}/ALL.chrX_${j}.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.msav"
                # Run Minimac4
                minimac4 "$MYREF" "$MYINPUT_STEP6" \
                    --output "imputed_${chr}_${j}.dose.vcf.gz" \
                    --output-format vcf.gz \
                    --min-ratio 0.001 --chunk 30000000 --overlap 3000000 --threads "$THREADS"
            done
        else
            MYREF="${REF_BASE}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.msav"
            # Run Minimac4
            minimac4 "$MYREF" "$MYINPUT_STEP6" \
                --output "imputed_${chr}.dose.vcf.gz" \
                --output-format vcf.gz \
                --min-ratio 0.001 --chunk 30000000 --overlap 3000000 --threads "$THREADS"
        fi

        # Sort the file before indexing, as minimac4 output may not be sorted.
        # bcftools sort "imputed_${chr}.dose.vcf.gz" -Oz -o "imputed_${chr}.dose.sorted.vcf.gz"
        # mv "imputed_${chr}.dose.sorted.vcf.gz" "imputed_${chr}.dose.vcf.gz"
        # tabix -p vcf "imputed_${chr}.dose.vcf.gz"
    }

    export -f impute_chr
    export -f print_info
    export PREFIX LIFTED_CODE OUTROOT REF_MODE REF_BASE STEP6_OUT THREADS

    # Process all ancestry groups and chromosomes
    for anc in 1 2 3 4 5 mixed; do
        print_info "Imputing ancestry group: $anc"

        CHR_RANGE=$(get_chr_range)
        # if command -v parallel &> /dev/null; then
        #     parallel -j $(($THREADS / 4)) impute_chr ::: $(eval echo $CHR_RANGE) ::: "$anc"
        # else
        for chr in $(eval echo $CHR_RANGE); do
            impute_chr "$chr" "$anc"
        done
        # fi
    done

    print_info "Imputation completed"
fi

################################################################################
# Complete
################################################################################

print_step "Pipeline Complete!"

print_info "All steps completed successfully"
print_info "Results are in: $OUTROOT"
print_info ""
print_info "Output structure:"
print_info "  0_check_vcf_build/ - Genome build detection"
print_info "  1_lift/            - Lifted to GRCh37"
print_info "  2_GH/              - Harmonized genotypes"
print_info "  3_ancestry/        - Ancestry analysis results"
print_info "  4_split_QC2/       - QC filtered data"
print_info "  5_phase/           - Phased haplotypes"
print_info "  6_impute_${REF_MODE}/   - Final imputed genotypes"

echo
