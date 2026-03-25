#!/usr/bin/env bash

####################################
##                                ##
##    Imputation / QC Pipeline    ##
##          Sequential Version    ##
##                                ##
##  Last modified: Sequential Version     ##
##                                ##
####################################

set -e  # Exit on error
set -u  # Exit on undefined variable

# Get script directory (where this script is located)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Default values
myinput=""
outroot=""
ref_mode="HRC"
start_from=0
stop_after=6
tempdir=""
wgs_mode=false
hwe_mode=false
confirm_run=false
ref_path=""
fasta_ref=""
ref_vcf_dir=""
ref_bcf_dir=""
ref_minimac_dir=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf|-v)
            myinput="$2"
            shift 2
            ;;
        --out|-o)
            outroot="$2"
            shift 2
            ;;
        --ref|-r)
            ref_mode="$2"
            shift 2
            ;;
        --start|-s)
            start_from="$2"
            shift 2
            ;;
        --end|-e)
            stop_after="$2"
            shift 2
            ;;
        --temp|-t)
            tempdir="$2"
            shift 2
            ;;
        --wgs|-w)
            wgs_mode=true
            shift
            ;;
        --hwe|-h)
            hwe_mode=true
            shift
            ;;
        --confirm|-c)
            confirm_run=true
            shift
            ;;
        --ref-path|--ref_path)
            ref_path="$2"
            shift 2
            ;;
        --fasta-ref|--fasta_ref)
            fasta_ref="$2"
            shift 2
            ;;
        --ref-vcf-dir|--ref_vcf_dir)
            ref_vcf_dir="$2"
            shift 2
            ;;
        --ref-bcf-dir|--ref_bcf_dir)
            ref_bcf_dir="$2"
            shift 2
            ;;
        --ref-minimac-dir|--ref_minimac_dir)
            ref_minimac_dir="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate arguments
if [ -z "$myinput" ]; then
    echo "ERROR: Please provide --vcf input file path"
    exit 1
fi

if [ -z "$outroot" ]; then
    outroot="$PWD/IMP_QC"
    echo "INFO: No OUT_ROOT provided, using default: $outroot"
fi

if [ "$ref_mode" != "HRC" ] && [ "$ref_mode" != "1KG" ]; then
    echo "ERROR: Invalid reference panel. Use HRC or 1KG"
    exit 1
fi

if ! [[ "$start_from" =~ ^[0-9]+$ ]] || ! [[ "$stop_after" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Start and end steps must be numbers"
    exit 1
fi

if [ $start_from -gt $stop_after ]; then
    echo "ERROR: Start step ($start_from) must be <= end step ($stop_after)"
    exit 1
fi

# Set variables
wgs=${wgs_mode:+yes}
wgs=${wgs:-no}
hwe=${hwe_mode:+yes}
hwe=${hwe:-no}

echo
echo "=========================================="
echo "Sequential Pipeline Configuration"
echo "=========================================="
echo "Input VCF: $myinput"
echo "Output Root: $outroot"
echo "Reference: $ref_mode"
echo "Steps: $start_from to $stop_after"
echo "WGS Mode: $wgs"
echo "HWE Filtering: $hwe"
echo "Reference Path: $ref_path"
echo "FASTA Reference: $fasta_ref"
echo "VCF Reference Directory: $ref_vcf_dir"
echo "BCF Reference Directory: $ref_bcf_dir"
echo "MINIMAC Reference Directory: $ref_minimac_dir"
echo "=========================================="
echo

if [ "$confirm_run" != "true" ]; then
    echo "WARNING: Running in preview mode. Use --confirm to actually run the pipeline."
    echo
fi

# Parse input file info
indir=$(dirname "$myinput")
infile=$(basename "$myinput")
prefix=$(echo "$infile" | sed -e 's/\.vcf.gz$//g' | sed -e 's/\.txt$//g' | sed -e 's/\.vcf$//g')

# Export tools paths (assuming required_tools is in script directory)
export plink="$SCRIPT_DIR/required_tools/plink"
export plink2="$SCRIPT_DIR/required_tools/plink2"
export GH="$SCRIPT_DIR/required_tools/GenotypeHarmonizer/GenotypeHarmonizer.jar"
export lift="$SCRIPT_DIR/required_tools/lift/LiftMap.py"
export cpath="$SCRIPT_DIR/required_tools/chainfiles"
export eagle="$SCRIPT_DIR/required_tools/Eagle_v2.4.1/eagle"
export check_vcf_build="$SCRIPT_DIR/required_tools/check_vcf_build/check_vcf_build.R"
export split_by_ancestry="$SCRIPT_DIR/required_tools/split_by_ancestry/split_by_ancestry.R"
export ref_pop="$SCRIPT_DIR/required_tools/1000G_P3_super_pop.pop"
export mymap="$SCRIPT_DIR/required_tools/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz"

# Set reference paths (use provided arguments or fall back to defaults)
if [ -z "$ref_path" ]; then
    ref_path="/mnt/stsi/stsi3/External/1000G/ref_panel/hg19"
fi
if [ -z "$fasta_ref" ]; then
    fasta_ref="/mnt/stsi/stsi3/External/reference_fasta/hg19/human_g1k_v37.fasta"
fi
if [ -z "$ref_vcf_dir" ]; then
    ref_vcf_dir="/mnt/stsi/stsi3/Internal/1000G/ref_panel/hg19/vcf"
fi
if [ -z "$ref_bcf_dir" ]; then
    if [ "$ref_mode" == "1KG" ]; then
        ref_bcf_dir="/mnt/stsi/stsi3/Internal/1000G/ref_panel/hg19/bcf"
    else
        ref_bcf_dir="/mnt/stsi/stsi3/Internal/HRC/ref_panel/hg19/bcf"
    fi
fi
if [ -z "$ref_minimac_dir" ]; then
    if [ "$ref_mode" == "1KG" ]; then
        ref_minimac_dir="/mnt/stsi/stsi3/Internal/1000G/ref_panel/hg19/m3vcf_erate_rec"
    else
        ref_minimac_dir="/mnt/stsi/stsi3/Internal/HRC/ref_panel/hg19/m3vcf_erate_rec"
    fi
fi

export ref_path
export fasta_ref

# Create output directories
mkdir -p "$outroot/0_check_vcf_build"
mkdir -p "$outroot/1_lift"
mkdir -p "$outroot/2_GH"
mkdir -p "$outroot/3_ancestry"
mkdir -p "$outroot/4_split_QC2"
mkdir -p "$outroot/5_phase"
mkdir -p "$outroot/6_impute_${ref_mode}"

# Set temp directory
if [ -z "$tempdir" ]; then
    if [ -n "${TMPDIR:-}" ]; then
        TEMP="$TMPDIR"
    elif [ -n "${TMP:-}" ]; then
        TEMP="$TMP"
    else
        TEMP="/tmp"
    fi
else
    TEMP="$tempdir"
fi
export TEMP

# Function to run a step if it's in range
run_step() {
    local step=$1
    if [ $step -ge $start_from ] && [ $step -le $stop_after ]; then
        return 0
    else
        echo "# Skipping step $step (not in range $start_from-$stop_after)"
        return 1
    fi
}

# Function to check if file exists
check_file() {
    if [ ! -f "$1" ]; then
        echo "ERROR: Required file not found: $1"
        exit 1
    fi
}

# Function to check if tool exists
check_tool() {
    if [ ! -f "$1" ] && [ ! -x "$1" ]; then
        echo "WARNING: Tool not found or not executable: $1"
        echo "Please ensure all required tools are in the required_tools directory"
    fi
}

# Check required tools
check_tool "$plink"
check_tool "$plink2"
check_tool "$eagle"
[ -f "$GH" ] || echo "WARNING: GenotypeHarmonizer.jar not found at $GH"
[ -f "$check_vcf_build" ] || echo "WARNING: check_vcf_build.R not found at $check_vcf_build"

echo "Starting pipeline execution..."
echo

###############################################################################
# STEP 0: Check VCF build
###############################################################################
if run_step 0; then
    echo "=========================================="
    echo "STEP 0: Checking VCF build version"
    echo "=========================================="
    
    myoutput="$outroot/0_check_vcf_build/${prefix}.BuildChecked"
    
    if [[ "$myinput" == *.txt ]]; then
        echo "Input autosome list detected"
        my_chr1=$(cat "$myinput" | awk '$1==1 {print$2}')
        mydirname=$(dirname "$my_chr1")
        filename=$(basename "$my_chr1")
    elif [[ "$myinput" == *.vcf.gz ]]; then
        echo "Merged vcf.gz file detected"
        mydirname=$(dirname "$myinput")
        filename=$(basename "$myinput")
    else
        echo "ERROR: Invalid file format"
        exit 1
    fi
    
    if [ "$confirm_run" == "true" ]; then
        Rscript "$check_vcf_build" "$mydirname/$filename" > "$myoutput"
        echo "Build check complete: $myoutput"
    else
        echo "Would run: Rscript $check_vcf_build $mydirname/$filename > $myoutput"
    fi
    echo
fi

###############################################################################
# STEP 1: LiftOver to GRCh37
###############################################################################
if run_step 1; then
    echo "=========================================="
    echo "STEP 1: Lifting VCFs to GRCh37"
    echo "=========================================="
    
    buildcheck="$outroot/0_check_vcf_build/${prefix}.BuildChecked"
    myoutdir="$outroot/1_lift"
    
    if [ "$confirm_run" != "true" ]; then
        echo "Would run Step 1 liftOver"
        echo "  Input: $myinput"
        echo "  Build check: $buildcheck"
        echo "  Output dir: $myoutdir"
        echo
        # We need the lifted_code for next steps, so we'll need to run this
        # For now, skip to show structure
    else
        check_file "$buildcheck"
        
        export filename=$(basename "$buildcheck")
        export inprefix=${filename/.BuildChecked/}
        
        cd "$myoutdir"
        
        # Create temp directory
        TEMP_DIR="$TEMP/$$"
        mkdir -p "$TEMP_DIR"
        export TEMP="$TEMP_DIR"
        
        # Source the liftOver functions from step 1 script
        # We'll inline the key functions here
        
        # Check if input is split or merged
        if [[ "$myinput" == *.txt ]]; then
            echo "Input autosome list detected, processing in parallel"
            # Process each chromosome from the list
            while read -r line; do
                chrom=$(echo "$line" | awk '{print$1}')
                fp=$(echo "$line" | awk '{print$2}')
                
                # Check if tabixed
                if [ ! -e "${fp}.tbi" ]; then
                    echo "Tabixing chromosome $chrom..."
                    tabix -p vcf "$fp" || bcftools sort "$fp" -T "$TEMP" -Oz -o "./${inprefix}.chr${chrom}.sorted.vcf.gz"
                    tabix -p vcf "./${inprefix}.chr${chrom}.sorted.vcf.gz"
                else
                    cp "$fp" "./${inprefix}.chr${chrom}.sorted.vcf.gz"
                    cp "${fp}.tbi" "./${inprefix}.chr${chrom}.sorted.vcf.gz.tbi"
                fi
            done < "$myinput"
        elif [[ "$myinput" == *.vcf.gz ]]; then
            echo "Input file is merged; splitting by chromosome"
            # Sort and tabix if needed
            if [ ! -e "${myinput}.tbi" ]; then
                bcftools sort "$myinput" -T "$TEMP" -Oz -o "./${inprefix}.chrALL.sorted.vcf.gz"
                tabix -p vcf "./${inprefix}.chrALL.sorted.vcf.gz"
            else
                cp "$myinput" "./${inprefix}.chrALL.sorted.vcf.gz"
                cp "${myinput}.tbi" "./${inprefix}.chrALL.sorted.vcf.gz.tbi"
            fi
            
            # Split by chromosome
            for chrom in {1..22}; do
                bcftools view "./${inprefix}.chrALL.sorted.vcf.gz" -r "$chrom" -Oz -o "./${inprefix}.chr${chrom}.sorted.vcf.gz"
                tabix -p vcf "./${inprefix}.chr${chrom}.sorted.vcf.gz"
            done
        fi
        
        # Determine chain file from buildcheck
        cfilename=$(grep "Use chain file" "$buildcheck" | tr -d ' ' | tr ':' '\t' | tr -d '"' | cut -f 2 | sed -e 's/->/ /g')
        checknone=$(grep "Use chain file" "$buildcheck" | grep "none" | wc -l)
        
        # LiftOver each chromosome
        for chrom in {1..22}; do
            name="${inprefix}.chr${chrom}"
            
            echo "Processing chromosome $chrom..."
            
            # Remove multi-allelic variants
            bcftools view "$name.sorted.vcf.gz" -M 2 -m 2 | bcftools norm /dev/stdin -d both -Oz -o "$name.sorted.bi.vcf.gz"
            tabix -p vcf "$name.sorted.bi.vcf.gz"
            
            # Convert to ped/map
            bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "$name.sorted.bi.vcf.gz" > "$name.sorted.bi.pos"
            
            "$plink" --vcf "$name.sorted.bi.vcf.gz" --make-bed --a1-allele "$name.sorted.bi.pos" 5 3 --biallelic-only strict \
                --set-missing-var-ids @:#:\$1:\$2 --vcf-half-call missing --double-id --recode ped --id-delim '_' --out "$name"
            
            if [ $checknone -eq 1 ]; then
                echo "Dataset already in GRCh37, converting format"
                "$plink2" --vcf "$name.sorted.bi.vcf.gz" --make-bed --id-delim '_' \
                    --set-missing-var-ids @:#:\$1:\$2 --allow-extra-chr --out "$name.lifted_already_GRCh37.sorted.with_dup"
                
                cut -f 2 "$name.lifted_already_GRCh37.sorted.with_dup.bim" | sort | uniq -d > "$name.list_multi_a_markers.txt"
                ndup=$(wc -l < "$name.list_multi_a_markers.txt")
                
                if [ "$ndup" -ge 1 ]; then
                    "$plink2" --bfile "$name.lifted_already_GRCh37.sorted.with_dup" \
                        --exclude "$name.list_multi_a_markers.txt" --make-bed --out "$name.lifted_already_GRCh37"
                else
                    "$plink2" --bfile "$name.lifted_already_GRCh37.sorted.with_dup" --make-bed --out "$name.lifted_already_GRCh37"
                    for bfile in bed bim fam; do
                        mv "$name.lifted_already_GRCh37.${bfile}" "${inprefix}.lifted_already_GRCh37.chr${chrom}.${bfile}"
                    done
                fi
            else
                # Perform liftOver
                nchains=$(echo "$cfilename" | tr ' ' '\n' | wc -l | awk '{print $1}')
                
                if [ "$nchains" -eq 1 ]; then
                    liftname=$(echo "$cfilename" | sed -e 's/\.chain.*//g')
                    "$lift" -p "$name.ped" -m "$name.map" -c "$cpath/$cfilename" -o "$name.lifted_$liftname"
                elif [ "$nchains" -eq 2 ]; then
                    first=$(echo "$cfilename" | tr ' ' '\n' | head -n 1)
                    second=$(echo "$cfilename" | tr ' ' '\n' | tail -n 1)
                    liftname=$(echo "$first" | sed -e 's/\.chain.*//g')
                    "$lift" -m "$name.map" -c "$cpath/$first" -o "$name.lifted_$liftname"
                    liftname=$(echo "$second" | sed -e 's/\.chain.*//g')
                    "$lift" -m "$name.lifted_$liftname.map" -c "$cpath/$second" -o "$name.lifted_$liftname"
                fi
                
                # Convert lifted output to bed/bim/fam
                "$plink" --file "$name.lifted_$liftname" --a1-allele "$name.sorted.bi.pos" 5 3 --double-id \
                    --set-missing-var-ids @:#:\$1:\$2 --allow-extra-chr --make-bed --out "$name.lifted_$liftname.sorted.with_dup"
                
                cut -f 2 "$name.lifted_$liftname.sorted.with_dup.bim" | sort | uniq -d > "$name.list_multi_a_markers.txt"
                ndup=$(wc -l < "$name.list_multi_a_markers.txt")
                
                if [ "$ndup" -ge 1 ]; then
                    "$plink2" --bfile "$name.lifted_$liftname.sorted.with_dup" --exclude "$name.list_multi_a_markers.txt" \
                        --make-bed --id-delim '_' --allow-extra-chr --out "$name.lifted_$liftname"
                else
                    "$plink2" --bfile "$name.lifted_$liftname.sorted.with_dup" --make-bed --id-delim '_' \
                        --allow-extra-chr --out "$name.lifted_$liftname"
                    for bfile in bed bim fam; do
                        mv "$name.lifted_$liftname.sorted.with_dup.${bfile}" "${inprefix}.lifted_${liftname}.chr${chrom}.${bfile}"
                    done
                fi
            fi
            
            # Cleanup
            rm -f "$name.map" "$name.ped" *with_dup* "$name.list_multi_a_markers.txt" "$name.sorted.bi."*
        done
        
        # Determine lifted_code for next steps from actual output files
        first_output=$(ls "$myoutdir/${prefix}".lifted_*.chr1.bed 2>/dev/null | head -1)
        if [ -n "$first_output" ]; then
            lifted_code=$(basename "$first_output" | sed "s/^${prefix}\.//" | sed "s/\.chr1\.bed$//")
            echo "Lifted code determined: $lifted_code"
        else
            echo "WARNING: Could not find output files to determine lifted_code"
        fi
    fi
    echo
fi

# Get lifted_code (needed for subsequent steps)
if [ -z "${lifted_code:-}" ] && [ -d "$outroot/1_lift" ]; then
    first_bed=$(ls "$outroot/1_lift/${prefix}".lifted_*.chr1.bed 2>/dev/null | head -1)
    if [ -n "$first_bed" ]; then
        lifted_code=$(basename "$first_bed" | sed "s/^${prefix}\.//" | sed "s/\.chr1\.bed$//")
        echo "Detected lifted_code from existing files: $lifted_code"
    fi
fi

# If still no lifted_code and we need it, try to infer or exit
if [ -z "${lifted_code:-}" ] && [ $start_from -le 1 ] && [ $stop_after -ge 2 ]; then
    echo "WARNING: Could not determine lifted_code automatically."
    echo "Please ensure step 1 has completed successfully, or manually set the lifted_code pattern."
    echo "Looking for files matching: ${prefix}.lifted_*.chr*.bed"
    
    # Try to find any lifted file
    any_lifted=$(find "$outroot/1_lift" -name "${prefix}.lifted_*.chr*.bed" 2>/dev/null | head -1)
    if [ -n "$any_lifted" ]; then
        lifted_code=$(basename "$any_lifted" | sed "s/^${prefix}\.//" | sed "s/\.chr.*\.bed$//")
        echo "Found lifted_code from file: $lifted_code"
    else
        echo "ERROR: Cannot proceed without lifted_code. Please run steps 0-1 first."
        exit 1
    fi
fi

###############################################################################
# STEP 2: Genotype Harmonizer QC1
###############################################################################
if run_step 2; then
    echo "=========================================="
    echo "STEP 2: Genotype Harmonizer QC1"
    echo "=========================================="
    
    if [ -z "$lifted_code" ]; then
        echo "ERROR: Cannot proceed without lifted_code. Please run steps 0-1 first."
        exit 1
    fi
    
    # Find the first bed file - use lifted_code if available, otherwise search
    if [ -n "${lifted_code:-}" ]; then
        first_bed=$(ls "$outroot/1_lift/${prefix}.${lifted_code}.chr1.bed" 2>/dev/null | head -1)
    fi
    
    if [ -z "${first_bed:-}" ]; then
        first_bed=$(ls "$outroot/1_lift/${prefix}".lifted_*.chr1.bed 2>/dev/null | head -1)
    fi
    
    if [ -z "${first_bed:-}" ]; then
        echo "ERROR: Could not find lifted bed files in $outroot/1_lift/"
        echo "Please ensure step 1 completed successfully."
        exit 1
    fi
    
    # Update lifted_code from actual file if not set
    if [ -z "${lifted_code:-}" ]; then
        lifted_code=$(basename "$first_bed" | sed "s/^${prefix}\.//" | sed "s/\.chr1\.bed$//")
    fi
    
    myinput_step2=$(echo "$first_bed" | sed 's/\.chr1\.bed$//')
    myoutdir="$outroot/2_GH"
    
    if [ "$confirm_run" != "true" ]; then
        echo "Would run Step 2 Genotype Harmonizer"
        echo "  Input prefix: $myinput_step2"
        echo "  Output dir: $myoutdir"
        echo
    else
        export inprefix=$(basename "$myinput_step2")
        cd "$myoutdir"
        
        export outname="${inprefix}.GH"
        
        # Run GenotypeHarmonizer for each chromosome
        for chrom in {1..22}; do
            echo "Processing chromosome $chrom with GenotypeHarmonizer..."
            
            if [ "$chrom" -eq 23 ]; then
                java -Xmx16g -jar "$GH" --keep --input "${myinput_step2}.chr${chrom}" \
                    --ref "$ref_path/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz" \
                    --inputType PLINK_BED --callRateFilter 0.90 --output "./${outname}.chr${chrom}"
            else
                java -Xmx16g -jar "$GH" --keep --input "${myinput_step2}.chr${chrom}" \
                    --ref "$ref_path/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" \
                    --inputType PLINK_BED --callRateFilter 0.90 --output "./${outname}.chr${chrom}"
            fi
        done
        
        # Fix reference alleles and normalize
        for chrom in {1..22}; do
            echo "Fixing reference alleles for chromosome $chrom..."
            
            "$plink2" --bfile "${outname}.chr${chrom}" --max-alleles 2 \
                --set-missing-var-ids @:#:\$1:\$2 --export vcf-4.2 bgz --out "${outname}.chr${chrom}.0"
            
            tabix -p vcf "${outname}.chr${chrom}.0.vcf.gz"
            
            bcftools +fixref "${outname}.chr${chrom}.0.vcf.gz" --threads 16 -Oz -o "${outname}.chr${chrom}.1.vcf.gz" \
                -- -f "$fasta_ref" -m flip -d
            
            tabix -p vcf "${outname}.chr${chrom}.1.vcf.gz"
            
            bcftools view -Ou -c 2 "${outname}.chr${chrom}.1.vcf.gz" | bcftools norm -m -any | \
                bcftools norm --threads 16 -Oz -o "${outname}.chr${chrom}.vcf.gz" -d both -f "$fasta_ref"
            
            tabix -p vcf "${outname}.chr${chrom}.vcf.gz"
            
            # Cleanup intermediate files
            rm -f "${outname}.chr${chrom}.bed" "${outname}.chr${chrom}.bim" "${outname}.chr${chrom}.fam" \
                  "${outname}.chr${chrom}.0.vcf.gz" "${outname}.chr${chrom}.1.vcf.gz"
        done
    fi
    echo
fi

###############################################################################
# STEP 3: Ancestry Analysis
###############################################################################
if run_step 3; then
    echo "=========================================="
    echo "STEP 3: Ancestry Analysis"
    echo "=========================================="
    
    if [ -z "${lifted_code:-}" ]; then
        echo "ERROR: Cannot proceed without lifted_code. Please run steps 0-2 first."
        exit 1
    fi
    
    myinput_step3="$outroot/2_GH/${prefix}.${lifted_code}.GH"
    myoutdir="$outroot/3_ancestry"
    
    if [ "$confirm_run" != "true" ]; then
        echo "Would run Step 3 Ancestry Analysis"
        echo "  Input prefix: $myinput_step3"
        echo "  Output dir: $myoutdir"
        echo
    else
        export filename=$(basename "$myinput_step3")
        export inprefix=$(basename "$myinput_step3" | sed -e 's/\.vcf.gz$//g')
        export outsubdir=$(basename "$myinput_step3" | sed -e 's~\.lifted.*~~g')
        export indir=$(dirname "$myinput_step3")
        
        mkdir -p "$myoutdir/$outsubdir"
        cd "$myoutdir/$outsubdir"
        
        # Prune markers
        for chrom in {1..22}; do
            echo "Pruning chromosome $chrom..."
            
            if [ "$wgs" == "yes" ]; then
                bcftools view -R "$SCRIPT_DIR/required_tools/chr_pos_23andme.txt" "${myinput_step3}.chr${chrom}.vcf.gz" \
                    -Ov -o "${inprefix}.chr${chrom}.23andMe_pos.vcf"
                "$plink2" --vcf "${inprefix}.chr${chrom}.23andMe_pos.vcf" --indep-pairwise 100 10 0.05 \
                    --out "$inprefix.chr${chrom}"
            else
                "$plink2" --vcf "${myinput_step3}.chr${chrom}.vcf.gz" --indep-pairwise 100 10 0.05 \
                    --out "$inprefix.chr${chrom}"
            fi
            
            vcftools --gzvcf "${myinput_step3}.chr${chrom}.vcf.gz" --snps "$inprefix.chr${chrom}.prune.in" \
                --recode --recode-INFO-all --out "$inprefix.chr${chrom}.pruned"
            
            bgzip -c "$inprefix.chr${chrom}.pruned.recode.vcf" > "$inprefix.chr${chrom}.pruned.vcf.gz"
            tabix -p vcf "$inprefix.chr${chrom}.pruned.vcf.gz"
            rm "$inprefix.chr${chrom}.pruned.recode.vcf"
        done
        
        # Intersect with 1000G reference
        for chrom in {1..22}; do
            echo "Intersecting chromosome $chrom with 1000G..."
            # Find VCF file in reference directory (supports both ALL.chr*.vcf.gz and HRC.chr*.vcf.gz patterns)
            ref_panel_file=$(find "$ref_vcf_dir" -maxdepth 1 -name "*chr${chrom}*.vcf.gz" | head -1)
            if [ -z "$ref_panel_file" ] || [ ! -f "$ref_panel_file" ]; then
                echo "ERROR: Could not find reference VCF file for chromosome $chrom in $ref_vcf_dir"
                echo "Looking for pattern: *chr${chrom}*.vcf.gz"
                exit 1
            fi
            bcftools isec "${inprefix}.chr${chrom}.pruned.vcf.gz" \
                "$ref_panel_file" \
                -p "${inprefix}.chr${chrom}_tmp" -n =2 -w 1,2 -Oz
            
            bcftools merge "${inprefix}.chr${chrom}_tmp/0000.vcf.gz" "${inprefix}.chr${chrom}_tmp/0001.vcf.gz" \
                -Oz -o "$inprefix.chr${chrom}.pruned.intersect1KG.vcf.gz"
            tabix -f -p vcf "$inprefix.chr${chrom}.pruned.intersect1KG.vcf.gz"
        done
        
        # Concatenate
        ls -v "$inprefix.chr"*.pruned.intersect1KG.vcf.gz > "$inprefix.pruned.intersect1KG.vcflist"
        bcftools concat -f "$inprefix.pruned.intersect1KG.vcflist" -Ov -o "$inprefix.pruned.intersect1KG.vcf"
        
        # Rename reference panel IDs
        bcftools query -l "${inprefix}.chr1_tmp/0001.vcf.gz" > 0001.txt
        while read -r line; do printf "${line}\t${line}_${line}\n"; done < 0001.txt > 0001.rename
        bcftools reheader "$inprefix.pruned.intersect1KG.vcf" -s 0001.rename > "$inprefix.pruned.intersect1KG.id-delim.vcf"
        rm "$inprefix.pruned.intersect1KG.vcf"
        
        # Convert to bed
        "$plink2" --vcf "$inprefix.pruned.intersect1KG.id-delim.vcf" --make-bed --id-delim '_' \
            --out "$inprefix.pruned.intersect1KG"
        rm "$inprefix.pruned.intersect1KG.id-delim.vcf"
        
        # Generate population file
        inputN=$(bcftools query -l "$inprefix.chr22.pruned.vcf.gz" | wc -l)
        > "$inprefix.pruned.intersect1KG.pop"
        for i in $(seq 1 1 $inputN); do
            echo "" >> "$inprefix.pruned.intersect1KG.pop"
        done
        cat "$ref_pop" >> "$inprefix.pruned.intersect1KG.pop"
        
        # Run ADMIXTURE
        echo "Running ADMIXTURE..."
        admixture --supervised "$inprefix.pruned.intersect1KG.bed" 5 -j16
        
        # Process results
        cat "$inprefix.pruned.intersect1KG.fam" | awk '{print $1 "_" $2}' > "$inprefix.pruned.intersect1KG.subjectIDs"
        paste "$inprefix.pruned.intersect1KG.subjectIDs" "$inprefix.pruned.intersect1KG.5.Q" | tr '\t' ' ' > "$inprefix.pruned.intersect1KG.5.Q.IDs"
        
        bcftools query -l "$inprefix.chr1.pruned.vcf.gz" | tr '_' '\t' | awk '{print$1"_"$2"\t"$1"\t"$2}' > "$inprefix.pruned.subjectIDs"
        
        # Split by ancestry
        Rscript "$split_by_ancestry" "$inprefix.pruned.intersect1KG.5.Q.IDs" "$inprefix.pruned.subjectIDs" 0.95
        
        # Split input by ancestry
        for chrom in {1..22}; do
            "$plink2" --vcf "${myinput_step3}.chr${chrom}.vcf.gz" --make-bed --id-delim '_' --out "$inprefix.chr${chrom}"
            
            for i in {1..5} mixed; do
                if [ -f "$inprefix.pruned.intersect1KG.5.Q.IDs.$i.ids" ]; then
                    "$plink2" --bfile "$inprefix.chr${chrom}" --keep "$inprefix.pruned.intersect1KG.5.Q.IDs.$i.ids" \
                        --make-bed --out "$inprefix.ancestry-$i.chr${chrom}"
                    echo "ancestry-${i}.chr${chrom} done"
                fi
            done
        done
    fi
    echo
fi

###############################################################################
# STEP 4: QC2 - Split and filter
###############################################################################
if run_step 4; then
    echo "=========================================="
    echo "STEP 4: Quality Control 2 (Split & Filter)"
    echo "=========================================="
    
    # Find ancestry-split files from step 3
    ancestry_files=$(find "$outroot/3_ancestry" -name "*.ancestry-*.chr1.bed" 2>/dev/null || echo "")
    
    if [ -z "$ancestry_files" ]; then
        echo "WARNING: No ancestry-split files found. You may need to run step 3 first."
        echo "Attempting to use files from step 2..."
        # Fallback: use step 2 output
        myinput_step4="$outroot/2_GH/${prefix}.${lifted_code}.GH.chr1.bed"
        myinput_step4=$(echo "$myinput_step4" | sed 's/\.chr1\.bed$//')
    else
        # Process each ancestry group
        for anc_file in $ancestry_files; do
            myinput_step4=$(echo "$anc_file" | sed 's/\.chr1\.bed$//')
            break  # Process first one for now
        done
    fi
    
    myoutdir="$outroot/4_split_QC2"
    
    # Set filtering flags
    if [ -z "${hwe:-}" ] || [ "$hwe" == "no" ]; then
        hweflag=""
    else
        hweflag="--hwe 1e-10"
    fi
    
    if [ -z "${mind:-}" ]; then
        mindflag=""
    else
        mindflag="--mind ${mind:-0.05}"
    fi
    
    if [ -z "${geno:-}" ]; then
        genoflag=""
    else
        genoflag="--geno ${geno:-0.1}"
    fi
    
    if [ "$confirm_run" != "true" ]; then
        echo "Would run Step 4 QC2"
        echo "  Input prefix: $myinput_step4"
        echo "  Output dir: $myoutdir"
        echo "  Filters: $hweflag $mindflag $genoflag"
        echo
    else
        export inprefix=$(basename "$myinput_step4")
        export indir=$(dirname "$myinput_step4")
        outsubdir=$(basename "$myinput_step4" | sed -e 's~\.lifted.*~~g')
        
        mkdir -p "$myoutdir/$outsubdir"
        cd "$myoutdir/$outsubdir"
        
        # Filter each chromosome
        for chrom in {1..22}; do
            echo "Filtering chromosome $chrom..."
            "$plink2" --bfile "${indir}/${inprefix}.chr${chrom}" --make-bed $hweflag $mindflag $genoflag \
                --out "$inprefix.chr${chrom}.tmp"
        done
        
        # Remove samples with incomplete autosomes
        cat "$inprefix.chr"*.tmp.mindrem.id > "$inprefix.tmp.mindrem.id" 2>/dev/null || true
        if [ -f "$inprefix.tmp.mindrem.id" ]; then
            cat "$inprefix.tmp.mindrem.id" | sort -r | uniq > "$inprefix.tmp.unique.mindrem.id"
            for chrom in {1..22}; do
                if [ -f "$inprefix.chr${chrom}.tmp.bed" ]; then
                    "$plink2" --bfile "$inprefix.chr${chrom}.tmp" --remove "$inprefix.tmp.unique.mindrem.id" \
                        --make-bed --out "$inprefix.chr${chrom}"
                fi
            done
        else
            # No samples to remove, just rename
            for chrom in {1..22}; do
                if [ -f "$inprefix.chr${chrom}.tmp.bed" ]; then
                    mv "$inprefix.chr${chrom}.tmp.bed" "$inprefix.chr${chrom}.bed"
                    mv "$inprefix.chr${chrom}.tmp.bim" "$inprefix.chr${chrom}.bim"
                    mv "$inprefix.chr${chrom}.tmp.fam" "$inprefix.chr${chrom}.fam"
                fi
            done
        fi
    fi
    echo
fi

###############################################################################
# STEP 5: Phasing
###############################################################################
if run_step 5; then
    echo "=========================================="
    echo "STEP 5: Phasing"
    echo "=========================================="
    
    # Find QC2 output files
    qc2_files=$(find "$outroot/4_split_QC2" -name "*.chr1.bed" 2>/dev/null | head -1)
    
    if [ -z "$qc2_files" ]; then
        echo "ERROR: No QC2 output files found. Please run step 4 first."
        exit 1
    fi
    
    myoutdir="$outroot/5_phase"
    
    if [ "$confirm_run" != "true" ]; then
        echo "Would run Step 5 Phasing"
        echo "  Input files from: $outroot/4_split_QC2"
        echo "  Output dir: $myoutdir"
        echo
    else
        # Process each ancestry group and chromosome
        for bed_file in $(find "$outroot/4_split_QC2" -name "*.chr*.bed"); do
            inprefix=$(basename "$bed_file" | sed -e 's/\.bed$//g')
            indir=$(dirname "$bed_file")
            mychr=$(echo "$inprefix" | sed -e 's/.*\.chr//g')
            outsubdir=$(basename "$bed_file" | sed -e 's/\.lifted.*//g')
            
            mkdir -p "$myoutdir/$outsubdir"
            cd "$myoutdir/$outsubdir"
            
            echo "Phasing chromosome $mychr..."
            
            # Convert to VCF
            "$plink2" --bfile "$indir/$inprefix" --export vcf-4.2 bgz \
                --set-missing-var-ids @:#:\$1:\$2 --out "$inprefix"
            
            tabix -f -p vcf "$inprefix.vcf.gz"
            
            # Find BCF reference file in directory
            myref=$(find "$ref_bcf_dir" -maxdepth 1 -name "*chr${mychr}*.bcf" | head -1)
            if [ -z "$myref" ] || [ ! -f "$myref" ]; then
                echo "ERROR: Could not find reference BCF file for chromosome $mychr in $ref_bcf_dir"
                echo "Looking for pattern: *chr${mychr}*.bcf"
                exit 1
            fi
            
            # Run Eagle
            "$eagle" --vcfTarget="$inprefix.vcf.gz" \
                --vcfRef="$myref" \
                --noImpMissing \
                --geneticMapFile="$mymap" \
                --Kpbwt=100000 --numThreads=16 \
                --chrom=$mychr --allowRefAltSwap \
                --outPrefix="$inprefix.phased"
        done
    fi
    echo
fi

###############################################################################
# STEP 6: Imputation
###############################################################################
if run_step 6; then
    echo "=========================================="
    echo "STEP 6: Imputation"
    echo "=========================================="
    
    myoutdir="$outroot/6_impute_${ref_mode}"
    
    if [ "$confirm_run" != "true" ]; then
        echo "Would run Step 6 Imputation"
        echo "  Input files from: $outroot/5_phase"
        echo "  Output dir: $myoutdir"
        echo "  Reference: $ref_mode"
        echo
    else
        # Find phased VCF files
        phased_files=$(find "$outroot/5_phase" -name "*.phased.vcf.gz" 2>/dev/null)
        
        if [ -z "$phased_files" ]; then
            echo "ERROR: No phased VCF files found. Please run step 5 first."
            exit 1
        fi
        
        for phased_file in $phased_files; do
            inprefix=$(basename "$phased_file" | sed -e 's/\.phased\.vcf\.gz//g')
            mychr=$(echo "$inprefix" | sed -e 's/.*\.chr//g')
            subdir=$(basename "$phased_file" | sed -e 's/\.lifted.*//g')
            subanc=$(basename "$phased_file" | tr '.' '\n' | grep "ancestry" || echo "")
            
            if [ -n "$subanc" ]; then
                outsubdir="$subdir/$subanc"
            else
                outsubdir="$subdir"
            fi
            
            mkdir -p "$myoutdir/$outsubdir"
            cd "$myoutdir/$outsubdir"
            
            echo "Imputing chromosome $mychr..."
            
            # Find m3vcf reference file in directory
            myref=$(find "$ref_minimac_dir" -maxdepth 1 -name "*chr${mychr}*.m3vcf.gz" | head -1)
            if [ -z "$myref" ] || [ ! -f "$myref" ]; then
                echo "ERROR: Could not find reference m3vcf file for chromosome $mychr in $ref_minimac_dir"
                echo "Looking for pattern: *chr${mychr}*.m3vcf.gz"
                exit 1
            fi
            
            # Run minimac4
            minimac4 --refHaps "$myref" \
                --haps "$phased_file" \
                --prefix "imputed_${mychr}" \
                --ignoreDuplicates \
                --minRatio 0.01 --ChunkLengthMb 30.00 --ChunkOverlapMb 3.00 --cpus 16
            
            # Index output
            tabix -p vcf "imputed_${mychr}.dose.vcf.gz"
        done
    fi
    echo
fi

echo "=========================================="
echo "Pipeline execution complete!"
echo "=========================================="
echo "Output directory: $outroot"
echo "Steps executed: $start_from to $stop_after"
