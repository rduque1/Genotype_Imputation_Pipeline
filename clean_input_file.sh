#!/usr/bin/env bash

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/utils.sh"

usage() {
    echo "Usage: $0 --vcf <input_file> --out <output_dir>"
    exit 1
}

# Args
MYINPUT=""
OUTROOT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf|-v) MYINPUT="$2"; shift 2 ;;
        --out|-o) OUTROOT="$2"; shift 2 ;;
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

# Setup
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Fallback for PLINK2
if [[ ! -f "plink2" ]]; then
    if command -v plink2 &> /dev/null; then
        PLINK2="plink2"
    else
        print_warning "plink2 not found at plink2 and not in PATH."
    fi
fi

mkdir -p "$OUTROOT"
INFILE=$(basename "$MYINPUT")
# Define prefix
PREFIX=$(echo "$INFILE" | sed -E 's/\.vcf(\.gz|\.zip)?$//; s/\.txt(\.gz|\.zip)?$//; s/\.tsv(\.gz|\.zip)?$//')

# Logic
inputFile="$MYINPUT"
TEMP_FILE="${OUTROOT}/tempFile.txt"
SORTED_DATA="${OUTROOT}/sorted_rawData.txt"
X_DATA="${OUTROOT}/X_rawData.txt"

if [[ -s "${SORTED_DATA}.dedup.txt" ]]; then
    print_info "Found existing cleaned file, skipping clean_input_file: ${SORTED_DATA}.dedup.txt"
    exit 0
fi

# Convert VCF to 23andMe format if requested
if [[ "$MYINPUT" =~ \.vcf(\.gz|\.zip)?$ ]]; then
    print_step "Converting VCF to 23andMe (txt) format"

    PLINK_DIR="${OUTROOT}/00_plink"
    mkdir -p "$PLINK_DIR"

    VCF_TO_CONVERT="$MYINPUT"

    # Handle zip files
    if [[ "$MYINPUT" =~ \.zip$ ]]; then
        print_info "Unzipping compressed VCF..."
        unzip -p "$MYINPUT" > "${PLINK_DIR}/${PREFIX}.vcf"
        VCF_TO_CONVERT="${PLINK_DIR}/${PREFIX}.vcf"
    fi

    # Run PLINK export to 23andme
    plink2 --max-alleles 2 --vcf "$VCF_TO_CONVERT" --make-bed --out temp_data
    plink --bfile temp_data --recode 23 --snps-only just-acgt --out "${PLINK_DIR}/${PREFIX}"

    if [[ -f "${PLINK_DIR}/${PREFIX}.txt" ]]; then
        print_success "Converted to 23andMe format: ${PLINK_DIR}/${PREFIX}.txt"
        inputFile="${PLINK_DIR}/${PREFIX}.txt"
    else
        print_warning "PLINK conversion may have failed, check log."
    fi
else
    # Handle compressed input files (zip or gzip)
    if [[ "$inputFile" =~ \.gz$ ]]; then
        print_step "Decompressing gzipped input file"
        DECOMPRESSED_FILE="${OUTROOT}/${PREFIX}.txt"
        gunzip -c "$inputFile" > "$DECOMPRESSED_FILE"
        inputFile="$DECOMPRESSED_FILE"
        print_success "Decompressed to: $inputFile"
    elif [[ "$inputFile" =~ \.zip$ ]]; then
        print_step "Extracting zipped input file"
        DECOMPRESSED_FILE="${OUTROOT}/${PREFIX}.txt"
        unzip -p "$inputFile" > "$DECOMPRESSED_FILE"
        inputFile="$DECOMPRESSED_FILE"
        print_success "Extracted to: $inputFile"
    fi
    # Count delimiters, and pick the most frequent one
    SEPARATOR=$(awk '
        {
            t += gsub(/\t/, "");
            c += gsub(/,/, "");
            s += gsub(/;/, "")
        }
        END {
            if (t > c && t > s) print "\t";
            else if (c > t && c > s) print ",";
            else if (s > t && s > c) print ";";
            else print "\t"; # Default to tab
        }' $inputFile)

    echo "Detected Separator: [$SEPARATOR]"
    # If separator is not tab, convert to tab
    if [[ "$SEPARATOR" != $'\t' ]]; then
        print_step "Converting input file to tab-delimited format"
        CONVERTED_FILE="${OUTROOT}/${PREFIX}_tab.txt"
        sed "s/${SEPARATOR}/\t/g" "$inputFile" > "$CONVERTED_FILE"
        inputFile="$CONVERTED_FILE"
        print_success "Converted to tab-delimited: $inputFile"
    fi
    print_info "$(wc -l < "$inputFile" | tr -d ' ') lines in input file."
    cp "$inputFile" "${OUTROOT}/original_input.txt"
    print_info "Removing quotes from input file if any"
    awk '{ gsub(/"/, ""); print}' $inputFile > "${inputFile}.noquotes"
    print_info "$(wc -l < "${inputFile}.noquotes" | tr -d ' ') lines in input file no quotes."
    mv "${inputFile}.noquotes" "$inputFile"
    print_info "$(wc -l < "$inputFile" | tr -d ' ') lines in input file."
fi

print_step "Cleaning and Formatting Input File"
print_info "Processing: $inputFile"
print_info "Output: $SORTED_DATA"

# Remove #, rsid
# Replace comma with tabs
# Replace quotes
# Replace chromosome 23 with X, 24 with Y, 25 with XY, 26 with MT
# Concat column 4 and 5 in case of Ancestry.dna
grep -v "#" "$inputFile" | \
# Replace commas with tabs
sed "s/,/\t/g" | \
# Remove quotes
sed "s/\"//g" | \
# Remove header line
grep -Piv "rsid\t" | \
# Replace 23 with X
sed -E "s/([a-zA-Z0-9\.]+)\t23\t/\1\tX\t/g" | \
# Replace 24 with Y
sed -E "s/([a-zA-Z0-9\.]+)\t24\t/\1\tY\t/g" | \
# Replace 25 with XY
sed -E "s/([a-zA-Z0-9\.]+)\t25\t/\1\tXY\t/g" | \
# Replace 26 with MT
sed -E "s/([a-zA-Z0-9\.]+)\t26\t/\1\tMT\t/g" | \
awk '{ print $1 "\t" $2 "\t" $3 "\t" $4 $5}' > "$TEMP_FILE"

# Remove lines with single allele for chr1-22
sed -Ei "/^[a-zA-Z0-9\.\-\;]+\t([1-9]|[0-9][0-9])+\t[0-9]+\t[A-Z][^A-Z]/d" "$TEMP_FILE"

# Sort content by chromosome for PLINK
set +e
grep "^#" "$inputFile" > "$SORTED_DATA"
set -e

# Numeric chromosomes 1-22
for i in {1..22}; do
    grep -P "^[a-zA-Z0-9\.]+\t$i\t" "$TEMP_FILE" >> "$SORTED_DATA" || true;
done

set +e
# Handle X chromosome
grep -P "^[a-zA-Z0-9\.]+\tX\t" "$TEMP_FILE" > "$X_DATA"

# Append X data or dummy data
if [[ -s "$X_DATA" ]]; then
    cat "$X_DATA" >> "$SORTED_DATA"
else
    # Dummy X data if needed
    DUMMY_X="${SCRIPT_DIR}/required_tools/X_dummy_data.txt"
    [[ ! -f "$DUMMY_X" ]] && DUMMY_X="/scripts/data/X_dummy_data.txt"

    if [[ -f "$DUMMY_X" ]]; then
        print_warning "No X chromosome data found, using dummy data."
        cat "$DUMMY_X" >> "$SORTED_DATA"
    else
        print_warning "No X chromosome data found and no dummy data available at $DUMMY_X."
    fi
fi

grep -P "^[a-zA-Z0-9\.]+\tY\t" "$TEMP_FILE" >> "$SORTED_DATA" || true  # Append Y data
grep -P "^[a-zA-Z0-9\.]+\tXY\t" "$TEMP_FILE" >> "$SORTED_DATA" || true # Append XY data
grep -P "^[a-zA-Z0-9\.]+\tMT\t" "$TEMP_FILE" >> "$SORTED_DATA" || true # Append MT data

rm -f "$TEMP_FILE" "$X_DATA"

# Remove duplicate entries while keeping headers
awk '!/^#/ && !seen[$1]++ || /^#/' "$SORTED_DATA" > "${SORTED_DATA}.dedup.txt"

print_success "Cleaning complete. Output: ${SORTED_DATA}.dedup.txt"

print_step "Converting to PLINK binary format (.bed, .bim, .fam) removing multiallelics"
plink --23file "${SORTED_DATA}.dedup.txt" FAM001 ID001 --snps-only --output-chr M --make-bed --out temp_binary
plink2 --bfile temp_binary --max-alleles 2 --make-bed --out "${OUTROOT}/sorted_cleaned_data"

print_success "PLINK binary files created: ${OUTROOT}/sorted_cleaned_data.bed, ${OUTROOT}/sorted_cleaned_data.bim, ${OUTROOT}/sorted_cleaned_data.fam"