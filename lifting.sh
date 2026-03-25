#!/usr/bin/env bash

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/utils.sh"


# Get hg18 SNP positions directly
# wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg18/database/snp130.txt.gz"
# wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg17/database/snp125.txt.gz"
# wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg16/database/snp.txt.gz"

# awk 'BEGIN{OFS="\t"} $2 !~ /_/ {print $5, $2, $3+1}' snp130.txt > hg18.tsv
# awk 'BEGIN{OFS="\t"} $2 !~ /_/ {print $5, $2, $3+1}' snp125.txt > hg17.tsv
# awk 'BEGIN{OFS="\t"} $2 !~ /_/ {print $5, $2, $3+1}' snp.txt > hg16.tsv

# awk '{print $1}' build-37-2022.txt > snp_ids_37.txt

# for i in hg18.tsv hg17.tsv hg16.tsv; do
#   # sed -i 's/^chr//g' $i
#   awk 'BEGIN{OFS="\t"} {gsub(/^chr/,"",$2); print}' $i > ${i%.tsv}_nchr.tsv
#   mv ${i%.tsv}_nchr.tsv $i
#   awk 'BEGIN {
#     while ((getline < "snp_ids_37.txt") > 0) keep[$1]=1
#   }
#   $1 in keep {
#     print $0
#   }' OFS="\t" $i > ${i%.tsv}_filtered.tsv
# done

#!/usr/bin/env bash
# Detect genome build and lift to GRCh37 if needed
# Usage: source lifting.sh (expects $rawdata and $CPATH variables set)

# Get rawdata from command line argument
rawdata=${1:-$rawdata}

if [[ -z "$rawdata" ]]; then
    echo "Error: rawdata file not specified. Usage: $0 <rawdata_file>"
    exit 1
fi

if [[ ! -f "$rawdata" ]]; then
    echo "Error: File '$rawdata' not found"
    exit 1
fi

if [[ -f "result.bed" && -f "result.bim" && -f "result.fam" ]]; then
  print_info "Found existing lifting outputs (result.bed/.bim/.fam), skipping lifting step"
  exit 0
fi

# Check if file is compressed and decompress if needed
if [[ "$rawdata" == *.zip ]]; then
  echo "Detected ZIP file, extracting..."
  unzip -p "$rawdata" > "${rawdata%.zip}.txt"
  rawdata="${rawdata%.zip}.txt"
elif [[ "$rawdata" == *.gz ]]; then
  echo "Detected GZIP file, decompressing..."
  gunzip -c "$rawdata" > "${rawdata%.gz}"
  rawdata="${rawdata%.gz}"
fi


# Count matches for each build
count37=$(awk '/^#/ {next} NR==FNR{a[$1,$3];next} ($1,$3) in a' /${SCRIPT_DIR}/required_tools/chainfiles/build-37-2022.txt ${rawdata} 2>/dev/null | wc -l) || count37=0
count38=$(awk '/^#/ {next} NR==FNR{a[$1,$3];next} ($1,$3) in a' /${SCRIPT_DIR}/required_tools/chainfiles/build-38.txt ${rawdata} 2>/dev/null | wc -l)  || count38=0
count36=$(awk '/^#/ {next} NR==FNR{a[$1,$3];next} ($1,$3) in a' /${SCRIPT_DIR}/required_tools/chainfiles/build-36.txt ${rawdata} 2>/dev/null | wc -l) || count36=0
count35=$(awk '/^#/ {next} NR==FNR{a[$1,$3];next} ($1,$3) in a' /${SCRIPT_DIR}/required_tools/chainfiles/build-35.txt ${rawdata} 2>/dev/null | wc -l) || count35=0
count34=$(awk '/^#/ {next} NR==FNR{a[$1,$3];next} ($1,$3) in a' /${SCRIPT_DIR}/required_tools/chainfiles/build-34.txt ${rawdata} 2>/dev/null | wc -l) || count34=0

echo "Build detection: 37=${count37} 38=${count38} 36=${count36} 35=${count35} 34=${count34}"

# Find best match
best_count=$count37
best_build="37"
chain_file=""

if [[ $count38 -gt $best_count ]]; then best_count=$count38; best_build="38"; chain_file="GRCh38_to_GRCh37.chain"; fi
if [[ $count36 -gt $best_count ]]; then best_count=$count36; best_build="36"; chain_file="NCBI36_to_GRCh37.chain"; fi
if [[ $count35 -gt $best_count ]]; then best_count=$count35; best_build="35"; chain_file="NCBI35_to_GRCh37.chain"; fi
if [[ $count34 -gt $best_count ]]; then best_count=$count34; best_build="34"; chain_file="NCBI34_to_GRCh37.chain"; fi

echo "Detected build: ${best_build}"

if [[ -z "$chain_file" ]]; then
  echo "Already GRCh37, no lifting needed"
  cat ${rawdata} > result.txt
  plink --23file ${rawdata} --make-bed --biallelic-only strict --snps-only --output-chr M --double-id --out result
else
  echo "Importing build ${best_build} data for lifting..."
  # Import the raw data into PLINK format first
  plink --23file "${rawdata}" --make-bed --output-chr M --double-id --out genome_to_lift

  echo "LIFTING from build ${best_build} to GRCh37 using ${chain_file}"
  # Recode to ped/map for the LiftMap.py script
  plink --bfile genome_to_lift --recode --out genome_to_lift_ped

  # Run your liftover script
  python "${SCRIPT_DIR}/required_tools/lift/LiftMap.py" \
         -p "genome_to_lift_ped.ped" \
         -m "genome_to_lift_ped.map" \
         -c "$chain_file" \
         -o "result_lifted"

  plink --file result_lifted --make-bed --biallelic-only strict --snps-only --output-chr M --double-id --out result
fi
