#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
ENV_PATH="$REPO_ROOT/.env"
CONFIG_PATH="$REPO_ROOT/pipeline.local.toml"

if [[ -f "$ENV_PATH" ]]; then
  set -a
  # shellcheck disable=SC1090
  source "$ENV_PATH"
  set +a
fi

if [[ -f "$CONFIG_PATH" ]]; then
  eval "$(python3 "$SCRIPT_DIR/load_pipeline_defaults.py" ns3 "$CONFIG_PATH" "$REPO_ROOT")"
fi

PYTHON_BIN="${PYTHON_BIN:-$REPO_ROOT/.venv/bin/python}"

usage() {
  cat <<'EOF'
Usage:
  EXCEL_FILE=/path/to/HCV_BlastHits.xlsx SHEET_NAME=NS3_PtGT0_Check FASTA_POOL=/path/to/FASTA [GENBANK_DIR=/path/to/genbank_seq_files] hcv-ns3-build-workflow/scripts/run_ns3_pipeline.sh

Optional environment variables:
  SHEET_NAME
  OUTPUT_DIR
  REFERENCE_FASTA
  SUBTYPE_JSON
  GT_AA_JSON
  PYTHON_BIN
  TEMP_ROOT
  ACCESSIONS_METADATA_CSV

Defaults can also be provided in .env and pipeline.local.toml at the repository root.
EOF
}

EXCEL_FILE="${EXCEL_FILE:-}"
FASTA_POOL="${FASTA_POOL:-}"
GENBANK_DIR="${GENBANK_DIR:-}"
SHEET_NAME="${SHEET_NAME:-}"
OUTPUT_DIR="${OUTPUT_DIR:-$REPO_ROOT/outputs/local_alignment}"
REFERENCE_FASTA="${REFERENCE_FASTA:-$REPO_ROOT/HCV_GT_RefSeqs.fasta}"
SUBTYPE_JSON="${SUBTYPE_JSON:-$REPO_ROOT/HCV_Subtype_Refs_By_Genome_NA.json}"
GT_AA_JSON="${GT_AA_JSON:-$REPO_ROOT/HCV_GT_Refs_By_Gene_AA.json}"
ACCESSIONS_METADATA_CSV="${ACCESSIONS_METADATA_CSV:-$REPO_ROOT/HCVData/Accessions_metadata.csv}"
SKILL_NAME="hcv-ns3-build-workflow"
TEMP_ROOT="${TEMP_ROOT:-$REPO_ROOT/temp/$SKILL_NAME/$(basename "$0" .sh)}"

if [[ -z "$EXCEL_FILE" || -z "$FASTA_POOL" || -z "$SHEET_NAME" ]]; then
  usage
  exit 1
fi

MATCHED_TXT="$OUTPUT_DIR/NS3_matched_fasta_files.txt"
INCLUDED_FASTA_DIR="$TEMP_ROOT/included_refid_fastas"
DISCOVERY_TMP="$TEMP_ROOT/find_refid_fastas"
DISCOVERY_JSON="$DISCOVERY_TMP/discovery_ns3.json"
SKILL_TEMP_ROOT="$REPO_ROOT/temp/$SKILL_NAME"
GT_ALLSTUDIES_JSON="$SKILL_TEMP_ROOT/build_ns3_gt_allstudies/last_run_summary.json"
SOURCEFEATURES_JSON="$SKILL_TEMP_ROOT/build_ns3_sourcefeatures_csv/last_run_summary.json"
SOURCEFEATURES_GROUPED_JSON="$SKILL_TEMP_ROOT/build_ns3_sourcefeatures_grouped_csv/last_run_summary.json"
SUBTYPE_ALLSTUDIES_JSON="$SKILL_TEMP_ROOT/build_ns3_subtype_allstudies_wseqs/last_run_summary.json"
SUBTYPE_WITH_GT_AA_JSON="$SKILL_TEMP_ROOT/build_ns3_subtype_with_gt_aa/last_run_summary.json"
COMPLETEPROFILES_JSON="$SKILL_TEMP_ROOT/build_ns3_completeprofiles_tabspergt/last_run_summary.json"
GT_RAS_JSON="$SKILL_TEMP_ROOT/build_ns3_gt_ras_profiles/last_run_summary.json"
SUBTYPE_RAS_JSON="$SKILL_TEMP_ROOT/build_ns3_subtype_ras_profiles/last_run_summary.json"
COMBINED_RAS_JSON="$SKILL_TEMP_ROOT/build_ns3_combined_ras_profiles/last_run_summary.json"
GT_AA_DISTANCE_SUMMARY="$SKILL_TEMP_ROOT/build_ns3_gt_aa_distance_matrix/last_run_summary.txt"
SUBTYPE_AA_DISTANCE_SUMMARY="$SKILL_TEMP_ROOT/build_ns3_subtype_aa_distance_matrices/last_run_summary.txt"
RAS_ENTROPY_SUMMARY="$SKILL_TEMP_ROOT/build_ns3_ras_entropy/last_run_summary.txt"
REFID_METADATA_DIR="$TEMP_ROOT/refid_metadata"
mkdir -p "$TEMP_ROOT"
mkdir -p "$(dirname "$GT_ALLSTUDIES_JSON")" "$(dirname "$SOURCEFEATURES_JSON")" "$(dirname "$SOURCEFEATURES_GROUPED_JSON")"
mkdir -p "$(dirname "$SUBTYPE_ALLSTUDIES_JSON")" "$(dirname "$SUBTYPE_WITH_GT_AA_JSON")"
mkdir -p "$(dirname "$COMPLETEPROFILES_JSON")"
mkdir -p "$(dirname "$GT_RAS_JSON")" "$(dirname "$SUBTYPE_RAS_JSON")"
mkdir -p "$(dirname "$COMBINED_RAS_JSON")"
mkdir -p "$(dirname "$GT_AA_DISTANCE_SUMMARY")" "$(dirname "$SUBTYPE_AA_DISTANCE_SUMMARY")" "$(dirname "$RAS_ENTROPY_SUMMARY")"

mkdir -p "$OUTPUT_DIR"

cleanup() {
  if [[ -n "${AA_TMP_WORKBOOK:-}" && -f "${AA_TMP_WORKBOOK:-}" ]]; then
    rm -f "$AA_TMP_WORKBOOK"
  fi
}
trap cleanup EXIT

announce_step() {
  local number="$1"
  local name="$2"
  local input="$3"
  local output="$4"
  input="${input//$REPO_ROOT\//}"
  output="${output//$REPO_ROOT\//}"
  printf '\n=== Step %s: %s ===\nInput: %s\nOutput: %s\n' "$number" "$name" "$input" "$output"
}

announce_step 1 "Prepare working directories" \
  "configuration: $CONFIG_PATH" \
  "temporary files: $TEMP_ROOT; final outputs: $OUTPUT_DIR"
rm -rf "$INCLUDED_FASTA_DIR"
mkdir -p "$INCLUDED_FASTA_DIR"
rm -rf "$DISCOVERY_TMP"
mkdir -p "$DISCOVERY_TMP"
rm -f "$SKILL_TEMP_ROOT/build_ns3_sourcefeatures_csv/NS3_SourceFeatures.csv"
rm -f "$SKILL_TEMP_ROOT/build_ns3_sourcefeatures_grouped_csv/NS3_SourceFeatures_Grouped.csv"

announce_step 2 "Find RefID FASTA files" \
  "workbook: $EXCEL_FILE (sheet: $SHEET_NAME); FASTA pool: $FASTA_POOL" \
  "discovery report: $DISCOVERY_TMP; matched file list: $MATCHED_TXT"
"$PYTHON_BIN" "$SCRIPT_DIR/find_refid_fastas.py" \
  --excel-file "$EXCEL_FILE" \
  --sheet "$SHEET_NAME" \
  --fasta-dir "$FASTA_POOL" \
  --output-dir "$DISCOVERY_TMP" \
  --numpatients-column 'Num Pts' \
  > "$DISCOVERY_JSON"

DISCOVERY_DIR="$(find "$DISCOVERY_TMP" -maxdepth 1 -type d -name 'refid_fasta_*' | head -n 1)"
cp "$DISCOVERY_DIR/matched_fasta_files.txt" "$MATCHED_TXT"

announce_step 3 "Stage matched RefID FASTA files" \
  "matched file list: $MATCHED_TXT" \
  "staged FASTA directory: $INCLUDED_FASTA_DIR"
while IFS= read -r src; do
  [[ -n "$src" ]] || continue
  cp "$src" "$INCLUDED_FASTA_DIR/"
done < "$MATCHED_TXT"

SOURCE_ACCESSION_COUNT="$(while IFS= read -r src; do
  [[ -n "$src" ]] || continue
  awk '/^>/ {sub(/^>/, ""); split($0, fields, /[[:space:]]+/); print fields[1]}' "$src"
done < "$MATCHED_TXT" | LC_ALL=C sort -u | wc -l | tr -d ' ')"
STAGED_ACCESSION_COUNT="$(awk '/^>/ {sub(/^>/, ""); split($0, fields, /[[:space:]]+/); print fields[1]}' "$INCLUDED_FASTA_DIR"/*.fasta | LC_ALL=C sort -u | wc -l | tr -d ' ')"
echo "staged_fasta_accessions_input=$SOURCE_ACCESSION_COUNT"
echo "staged_fasta_accessions_output=$STAGED_ACCESSION_COUNT"

announce_step 4 "Filter master accession metadata" \
  "staged FASTAs: $INCLUDED_FASTA_DIR; master metadata: $ACCESSIONS_METADATA_CSV" \
  "filtered metadata: $TEMP_ROOT/included_accessions_metadata.csv; genotype/subtype metadata: $TEMP_ROOT/included_accessions_genotype_subtype.csv"
"$PYTHON_BIN" "$SCRIPT_DIR/filter_accessions_metadata_by_fasta.py" \
  --fasta-dir "$INCLUDED_FASTA_DIR" \
  --metadata-csv "$ACCESSIONS_METADATA_CSV" \
  --output-dir "$TEMP_ROOT"

rm -rf "$REFID_METADATA_DIR"
announce_step 5 "Apply RefID-specific metadata rules" \
  "filtered metadata: $TEMP_ROOT/included_accessions_metadata.csv" \
  "per-RefID metadata: $REFID_METADATA_DIR"
"$PYTHON_BIN" "$SCRIPT_DIR/split_refid_metadata_csv.py" \
  --input-csv "$TEMP_ROOT/included_accessions_metadata.csv" \
  --output-dir "$REFID_METADATA_DIR" \
  > /dev/null

announce_step 6 "Filter staged FASTA files by RefID metadata" \
  "staged FASTAs: $INCLUDED_FASTA_DIR; per-RefID metadata: $REFID_METADATA_DIR" \
  "filtered staged FASTAs: $INCLUDED_FASTA_DIR"
"$PYTHON_BIN" "$SCRIPT_DIR/filter_refid_fastas_by_metadata.py" \
  --metadata-dir "$REFID_METADATA_DIR" \
  --fasta-dir "$INCLUDED_FASTA_DIR"

announce_step 7 "Build genotype study workbook" \
  "workbook: $EXCEL_FILE; filtered FASTAs: $INCLUDED_FASTA_DIR" \
  "genotype workbook: $OUTPUT_DIR/NS3_GT_AllStudies.xlsx"
echo "Action: align each filtered study sequence to the HCV genotype reference sequences."
echo "Result: group sequences by RefID/study, assign the best genotype, and record alignment results in the genotype workbook."
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_gt_allstudies.py" \
  --excel-file "$EXCEL_FILE" \
  --sheet "$SHEET_NAME" \
  --fasta-dir "$INCLUDED_FASTA_DIR" \
  --reference-fasta "$REFERENCE_FASTA" \
  --genotype-subtype-csv "$TEMP_ROOT/included_accessions_genotype_subtype.csv" \
  --output-dir "$OUTPUT_DIR" \
  --refid-column RefID \
  --refname-column RefName \
  --numpatients-column 'Num Pts' \
  > "$GT_ALLSTUDIES_JSON"

echo "Skipping disabled NS3 source-feature extraction and grouped-summary steps"
# if [[ -n "$GENBANK_DIR" ]]; then
#   "$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_sourcefeatures_csv.py" \
#     --matched-fasta-report "$MATCHED_TXT" \
#     --genbank-dir "$GENBANK_DIR" \
#     > "$SOURCEFEATURES_JSON"
#
#   "$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_sourcefeatures_grouped_csv.py" \
#     --gt-workbook "$OUTPUT_DIR/NS3_GT_AllStudies.xlsx" \
#     --summary-xlsx "$OUTPUT_DIR/NS3_NumSeqs_Naive_1PP_CoversRAS_ByStudy.xlsx" \
#     > "$SOURCEFEATURES_GROUPED_JSON"
# else
#   echo "GENBANK_DIR not provided; skipping NS3 source-feature extraction and grouped summary steps"
# fi

announce_step 8 "Build subtype study workbook" \
  "genotype workbook: $OUTPUT_DIR/NS3_GT_AllStudies.xlsx; filtered FASTAs: $INCLUDED_FASTA_DIR" \
  "subtype workbook: $OUTPUT_DIR/NS3_Subtype_AllStudies_WSeqs.xlsx"
echo "Action: use the genotype assignments and filtered sequences to compare each study with subtype reference sequences."
echo "Result: assign the best subtype within each genotype and write per-sequence amino-acid data to the subtype workbook."
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_subtype_allstudies_wseqs.py" \
  --combined-workbook "$OUTPUT_DIR/NS3_GT_AllStudies.xlsx" \
  --fasta-dir "$INCLUDED_FASTA_DIR" \
  --subtype-json "$SUBTYPE_JSON" \
  --genotype-subtype-csv "$TEMP_ROOT/included_accessions_genotype_subtype.csv" \
  --output-dir "$OUTPUT_DIR" \
  > "$SUBTYPE_ALLSTUDIES_JSON"

announce_step 9 "Add genotype amino-acid references" \
  "subtype workbook: $OUTPUT_DIR/NS3_Subtype_AllStudies_WSeqs.xlsx" \
  "intermediate amino-acid workbook: recorded in $SUBTYPE_WITH_GT_AA_JSON"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_subtype_with_gt_aa.py" \
  --subtype-workbook "$OUTPUT_DIR/NS3_Subtype_AllStudies_WSeqs.xlsx" \
  --fasta-dir "$INCLUDED_FASTA_DIR" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$SUBTYPE_WITH_GT_AA_JSON"

AA_TMP_WORKBOOK="$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["output_workbook"])' "$SUBTYPE_WITH_GT_AA_JSON")"

announce_step 10 "Build complete profile workbooks" \
  "amino-acid workbook: $AA_TMP_WORKBOOK" \
  "profile workbooks: $OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx; $OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx; profile accessions: $OUTPUT_DIR/NS3_Profile_Input_Accessions.csv"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_completeprofiles_tabspergt.py" \
  --input-workbook "$AA_TMP_WORKBOOK" \
  --output-dir "$OUTPUT_DIR" \
  --profile-accessions-csv "$OUTPUT_DIR/NS3_Profile_Input_Accessions.csv" \
  > "$COMPLETEPROFILES_JSON"

announce_step 11 "Export consensus FASTA files" \
  "profile workbooks: $OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx; $OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  "consensus FASTAs: $OUTPUT_DIR/NS3_GT_Consensus.fasta; $OUTPUT_DIR/NS3_Subtype_Consensus.fasta"
"$PYTHON_BIN" "$SCRIPT_DIR/export_ns3_consensus_fasta.py" \
  --gt-profile-workbook "$OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx" \
  --subtype-profile-workbook "$OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  --output-dir "$OUTPUT_DIR" \
  > /dev/null

announce_step 12 "Build genotype RAS profile" \
  "genotype profile workbook: $OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx" \
  "RAS profile: $OUTPUT_DIR/NS3_GT_RAS_Profiles.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_gt_ras_profiles.py" \
  --gt-profile-workbook "$OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$GT_RAS_JSON"

announce_step 13 "Build subtype RAS profile" \
  "subtype profile workbook: $OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  "RAS profile: $OUTPUT_DIR/NS3_Subtype_RAS_Profiles.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_subtype_ras_profiles.py" \
  --subtype-profile-workbook "$OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$SUBTYPE_RAS_JSON"

announce_step 14 "Combine genotype and subtype RAS profiles" \
  "GT RAS profile: $OUTPUT_DIR/NS3_GT_RAS_Profiles.xlsx; subtype RAS profile: $OUTPUT_DIR/NS3_Subtype_RAS_Profiles.xlsx" \
  "combined RAS profile: $OUTPUT_DIR/NS3_Combined_RAS_Profiles.xlsx"
echo "Action: place each genotype RAS consensus row before its subtype rows."
echo "Result: subtype cells retain amino-acid variants strictly above 10%; each genotype block ends with a blank row."
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_combined_ras_profiles.py" \
  --gt-ras-profile-workbook "$OUTPUT_DIR/NS3_GT_RAS_Profiles.xlsx" \
  --subtype-ras-profile-workbook "$OUTPUT_DIR/NS3_Subtype_RAS_Profiles.xlsx" \
  --output-xlsx "$OUTPUT_DIR/NS3_Combined_RAS_Profiles.xlsx" \
  > "$COMBINED_RAS_JSON"

announce_step 15 "Build genotype amino-acid distance matrix" \
  "consensus FASTA: $OUTPUT_DIR/NS3_GT_Consensus.fasta" \
  "distance workbook: $OUTPUT_DIR/NS3_GT_AA_Distance_Pos36_175.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_gt_aa_distance_matrix.py" \
  --input-fasta "$OUTPUT_DIR/NS3_GT_Consensus.fasta" \
  --aligned-fasta "$OUTPUT_DIR/NS3_GT_Consensus_aligned.fasta" \
  --output-xlsx "$OUTPUT_DIR/NS3_GT_AA_Distance_Pos36_175.xlsx" \
  --details-xlsx "$SKILL_TEMP_ROOT/build_ns3_gt_aa_distance_matrix/NS3_GT_AA_Distance_Pos36_175_details.xlsx" \
  --start 36 \
  --end 175 \
  > "$GT_AA_DISTANCE_SUMMARY"

announce_step 16 "Build subtype amino-acid distance matrices" \
  "consensus FASTA: $OUTPUT_DIR/NS3_Subtype_Consensus.fasta" \
  "distance workbook: $OUTPUT_DIR/NS3_Subtype_AA_Distance_Pos36_175.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_subtype_aa_distance_matrices.py" \
  --input-fasta "$OUTPUT_DIR/NS3_Subtype_Consensus.fasta" \
  --output-xlsx "$OUTPUT_DIR/NS3_Subtype_AA_Distance_Pos36_175.xlsx" \
  --temp-dir "$SKILL_TEMP_ROOT/build_ns3_subtype_aa_distance_matrices" \
  --start 36 \
  --end 175 \
  > "$SUBTYPE_AA_DISTANCE_SUMMARY"

announce_step 17 "Build RAS entropy reports" \
  "profile workbooks: $OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx; $OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  "entropy workbooks: $OUTPUT_DIR/NS3_GT_RAS_Entropy.xlsx; $OUTPUT_DIR/NS3_Subtype_RAS_Entropy.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_ras_entropy.py" \
  --gt-profile-workbook "$OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx" \
  --subtype-profile-workbook "$OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  --gt-output-xlsx "$OUTPUT_DIR/NS3_GT_RAS_Entropy.xlsx" \
  --subtype-output-xlsx "$OUTPUT_DIR/NS3_Subtype_RAS_Entropy.xlsx" \
  > "$RAS_ENTROPY_SUMMARY"

echo "NS3 pipeline complete"
echo "matched_fasta_report=${MATCHED_TXT#$REPO_ROOT/}"
echo "included_fasta_dir=${INCLUDED_FASTA_DIR#$REPO_ROOT/}"
echo "output_dir=${OUTPUT_DIR#$REPO_ROOT/}"
echo "temp_root=${TEMP_ROOT#$REPO_ROOT/}"
