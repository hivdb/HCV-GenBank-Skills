#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
ENV_PATH="$REPO_ROOT/.env"
CONFIG_PATH="$REPO_ROOT/pipeline.local.toml"

if [[ -f "$ENV_PATH" ]]; then
  set -a
  # shellcheck disable=SC1090
  source "$ENV_PATH"
  set +a
fi

if [[ -f "$CONFIG_PATH" ]]; then
  eval "$(python3 "$SCRIPT_DIR/load_pipeline_defaults.py" ns3_comet "$CONFIG_PATH" "$REPO_ROOT")"
fi

PYTHON_BIN="${PYTHON_BIN:-$REPO_ROOT/.venv/bin/python}"

usage() {
  cat <<'EOF'
Usage:
  EXCEL_FILE=/path/to/HCV_BlastHits.xlsx SHEET_NAME=NS3_PtGT0_Check FASTA_POOL=/path/to/FASTA [GENBANK_DIR=/path/to/genbank_seq_files] hcv-workflow/hcv-ns3-comet-build-workflow/scripts/run_ns3_pipeline.sh

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
OUTPUT_DIR="${OUTPUT_DIR:-$REPO_ROOT/outputs/comet}"
REFERENCE_FASTA="${REFERENCE_FASTA:-$REPO_ROOT/HCVData/HCV_GT_RefSeqs.fasta}"
SUBTYPE_JSON="${SUBTYPE_JSON:-$REPO_ROOT/HCVData/HCV_Subtype_Refs_By_Genome_NA.json}"
GT_AA_JSON="${GT_AA_JSON:-$REPO_ROOT/HCVData/HCV_GT_Refs_By_Gene_AA.json}"
ACCESSIONS_METADATA_CSV="${ACCESSIONS_METADATA_CSV:-$REPO_ROOT/HCVData/Accessions_metadata.csv}"
COMET_SUBTYPING_CSV="${COMET_SUBTYPING_CSV:-$REPO_ROOT/HCVData/Comet Subtyping/NS3.csv}"
NONCOMET_SUBTYPE_WORKBOOK="${NONCOMET_SUBTYPE_WORKBOOK:-$REPO_ROOT/outputs/local_alignment/NS3_Subtype_AllStudies_WSeqs.xlsx}"
SKILL_NAME="hcv-ns3-comet-build-workflow"
TEMP_ROOT="${TEMP_ROOT:-$REPO_ROOT/outputs/temp/$SKILL_NAME/$(basename "$0" .sh)}"

if [[ -z "$EXCEL_FILE" || -z "$FASTA_POOL" || -z "$SHEET_NAME" ]]; then
  usage
  exit 1
fi
if [[ ! -f "$NONCOMET_SUBTYPE_WORKBOOK" ]]; then
  echo "Missing non-COMET subtype workbook required for 1d and genotype 8 subtype overrides: $NONCOMET_SUBTYPE_WORKBOOK" >&2
  exit 1
fi

MATCHED_TXT="$OUTPUT_DIR/NS3_matched_fasta_files.txt"
INCLUDED_FASTA_DIR="$TEMP_ROOT/included_refid_fastas"
DISCOVERY_TMP="$TEMP_ROOT/find_refid_fastas"
DISCOVERY_JSON="$DISCOVERY_TMP/discovery_ns3.json"
SKILL_TEMP_ROOT="$REPO_ROOT/outputs/temp/$SKILL_NAME"
GT_ALLSTUDIES_JSON="$SKILL_TEMP_ROOT/build_ns3_gt_allstudies/last_run_summary.json"
SOURCEFEATURES_JSON="$SKILL_TEMP_ROOT/build_ns3_sourcefeatures_csv/last_run_summary.json"
SOURCEFEATURES_GROUPED_JSON="$SKILL_TEMP_ROOT/build_ns3_sourcefeatures_grouped_csv/last_run_summary.json"
SUBTYPE_ALLSTUDIES_JSON="$SKILL_TEMP_ROOT/build_ns3_subtype_allstudies_wseqs/last_run_summary.json"
SUBTYPE_WITH_GT_AA_JSON="$SKILL_TEMP_ROOT/build_ns3_subtype_with_gt_aa/last_run_summary.json"
COMPLETEPROFILES_JSON="$SKILL_TEMP_ROOT/build_ns3_completeprofiles_tabspergt/last_run_summary.json"
GT_RAS_JSON="$SKILL_TEMP_ROOT/build_ns3_gt_ras_profiles/last_run_summary.json"
SUBTYPE_RAS_JSON="$SKILL_TEMP_ROOT/build_ns3_subtype_ras_profiles/last_run_summary.json"
COMBINED_RAS_JSON="$SKILL_TEMP_ROOT/build_ns3_combined_ras_profiles/last_run_summary.json"
RAS_CONSENSUS_DIFFERENCE_JSON="$SKILL_TEMP_ROOT/build_ns3_subtype_ras_consensus_difference_summary/last_run_summary.json"
GT_AA_DISTANCE_SUMMARY="$SKILL_TEMP_ROOT/build_ns3_gt_aa_distance_matrix/last_run_summary.txt"
SUBTYPE_AA_DISTANCE_SUMMARY="$SKILL_TEMP_ROOT/build_ns3_subtype_aa_distance_matrices/last_run_summary.txt"
RAS_ENTROPY_SUMMARY="$SKILL_TEMP_ROOT/build_ns3_ras_entropy/last_run_summary.txt"
REFID_METADATA_DIR="$TEMP_ROOT/refid_metadata"
COMET_GENOTYPE_CSV="$TEMP_ROOT/comet_ns3_genotype_assignments.csv"
COMET_SUBTYPE_CSV="$TEMP_ROOT/comet_ns3_subtype_assignments.csv"
COMET_NOT_FOUND_CSV="$TEMP_ROOT/comet_ns3_not_found_or_unassigned.csv"
COMET_NOT_FOUND_FASTA="$TEMP_ROOT/comet_ns3_not_found_or_unassigned.fasta"
mkdir -p "$TEMP_ROOT"
mkdir -p "$(dirname "$GT_ALLSTUDIES_JSON")" "$(dirname "$SOURCEFEATURES_JSON")" "$(dirname "$SOURCEFEATURES_GROUPED_JSON")"
mkdir -p "$(dirname "$SUBTYPE_ALLSTUDIES_JSON")" "$(dirname "$SUBTYPE_WITH_GT_AA_JSON")"
mkdir -p "$(dirname "$COMPLETEPROFILES_JSON")"
mkdir -p "$(dirname "$GT_RAS_JSON")" "$(dirname "$SUBTYPE_RAS_JSON")"
mkdir -p "$(dirname "$COMBINED_RAS_JSON")"
mkdir -p "$(dirname "$RAS_CONSENSUS_DIFFERENCE_JSON")"
mkdir -p "$(dirname "$GT_AA_DISTANCE_SUMMARY")" "$(dirname "$SUBTYPE_AA_DISTANCE_SUMMARY")" "$(dirname "$RAS_ENTROPY_SUMMARY")"

mkdir -p "$OUTPUT_DIR"

cleanup() {
  :
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

announce_step 7 "Report and filter Comet status after RefID filtering" \
  "RefID-filtered FASTAs: $INCLUDED_FASTA_DIR; Comet calls: $COMET_SUBTYPING_CSV" \
  "Comet assignments: $COMET_GENOTYPE_CSV and $COMET_SUBTYPE_CSV; missing/unassigned raw FASTA: $COMET_NOT_FOUND_FASTA; retained FASTAs: $INCLUDED_FASTA_DIR"
echo "Action 1: print assigned, missing, and unassigned counts for the RefID-filtered FASTA accessions."
echo "Action 2: remove only accessions missing from Comet or marked unassigned."
"$PYTHON_BIN" "$SCRIPT_DIR/prepare_comet_ns3_assignments.py" \
  --comet-csv "$COMET_SUBTYPING_CSV" \
  --fasta-dir "$INCLUDED_FASTA_DIR" \
  --genotype-output-csv "$COMET_GENOTYPE_CSV" \
  --subtype-output-csv "$COMET_SUBTYPE_CSV" \
  --not-found-output-csv "$COMET_NOT_FOUND_CSV" \
  --not-found-fasta-output "$COMET_NOT_FOUND_FASTA" \
  --remove-unassigned

announce_step 8 "Build genotype study workbook from Comet" \
  "Comet genotype assignments: $COMET_GENOTYPE_CSV; filtered FASTAs: $INCLUDED_FASTA_DIR" \
  "genotype workbook: $OUTPUT_DIR/NS3_GT_AllStudies.xlsx"
echo "Action: write the Comet genotype for each retained accession; paired distance calculations run later."
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_comet_gt_allstudies.py" \
  --fasta-dir "$INCLUDED_FASTA_DIR" \
  --comet-genotype-csv "$COMET_GENOTYPE_CSV" \
  --output-dir "$OUTPUT_DIR" \
  > "$GT_ALLSTUDIES_JSON"

announce_step 8a "Add genotype-count summary" \
  "genotype workbook: $OUTPUT_DIR/NS3_GT_AllStudies.xlsx" \
  "updated workbook: $OUTPUT_DIR/NS3_GT_AllStudies.xlsx (Genotype_Counts sheet)"
"$PYTHON_BIN" "$SCRIPT_DIR/add_gt_counts_sheet.py" --workbook "$OUTPUT_DIR/NS3_GT_AllStudies.xlsx"

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

announce_step 9 "Build subtype study workbook from Comet" \
  "genotype workbook: $OUTPUT_DIR/NS3_GT_AllStudies.xlsx; Comet subtype assignments: $COMET_SUBTYPE_CSV; non-COMET 1d and genotype 7/8 subtype assignments: $NONCOMET_SUBTYPE_WORKBOOK" \
  "subtype workbook: $OUTPUT_DIR/NS3_Subtype_AllStudies_WSeqs.xlsx"
echo "Action: write the Comet subtype for each retained accession; no subtype alignment is run."
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_comet_subtype_allstudies.py" \
  --genotype-workbook "$OUTPUT_DIR/NS3_GT_AllStudies.xlsx" \
  --comet-subtype-csv "$COMET_SUBTYPE_CSV" \
  --noncomet-subtype-workbook "$NONCOMET_SUBTYPE_WORKBOOK" \
  --output-dir "$OUTPUT_DIR" \
  > "$SUBTYPE_ALLSTUDIES_JSON"

announce_step 10 "Extract genotype-position amino-acid sequences" \
  "subtype workbook: $OUTPUT_DIR/NS3_Subtype_AllStudies_WSeqs.xlsx; FASTA pool: $FASTA_POOL" \
  "amino-acid workbook: $OUTPUT_DIR/NS3_Profile_Input_Source.xlsx"
echo "Action: use blastx once against the Comet genotype reference to map and translate each raw sequence."
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_subtype_with_gt_aa.py" \
  --subtype-workbook "$OUTPUT_DIR/NS3_Subtype_AllStudies_WSeqs.xlsx" \
  --fasta-dir "$FASTA_POOL" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  --output-workbook "$OUTPUT_DIR/NS3_Profile_Input_Source.xlsx" \
  > "$SUBTYPE_WITH_GT_AA_JSON"

AA_TMP_WORKBOOK="$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["output_workbook"])' "$SUBTYPE_WITH_GT_AA_JSON")"
announce_step 10b "Validate NS3 profile alignment coordinates" \
  "amino-acid extraction: $AA_TMP_WORKBOOK" \
  "marked profile input: $OUTPUT_DIR/NS3_Profile_Input_Alignment_QC.xlsx; flagged accessions: $OUTPUT_DIR/NS3_Profile_Alignment_QC_Flagged_Accessions.csv"
"$PYTHON_BIN" "$SCRIPT_DIR/validate_ns3_profile_alignment.py" \
  --input-workbook "$AA_TMP_WORKBOOK" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-workbook "$OUTPUT_DIR/NS3_Profile_Input_Alignment_QC.xlsx" \
  --flagged-accessions-csv "$OUTPUT_DIR/NS3_Profile_Alignment_QC_Flagged_Accessions.csv" \
  > "$TEMP_ROOT/validate_ns3_profile_alignment.json"
AA_TMP_WORKBOOK="$OUTPUT_DIR/NS3_Profile_Input_Alignment_QC.xlsx"
announce_step 10c "Summarize QC-passed NS3 genotype mutation burden" \
  "QC profile input: $AA_TMP_WORKBOOK" \
  "genotype mutation-burden summary: $OUTPUT_DIR/NS3_QC_Passed_Genotype_Mutation_Burden_Summary.csv"
"$PYTHON_BIN" "$REPO_ROOT/build_qc_passed_genotype_mutation_burden_summary/build_qc_passed_genotype_mutation_burden_summary.py" \
  --input-workbook "$AA_TMP_WORKBOOK" \
  --output-csv "$OUTPUT_DIR/NS3_QC_Passed_Genotype_Mutation_Burden_Summary.csv"
PROFILE_INPUT_COUNTS="$("$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_completeprofiles_tabspergt.py" --input-workbook "$AA_TMP_WORKBOOK" --report-only)"
PROFILE_INPUT_INCLUDED_ACCESSIONS="$("$PYTHON_BIN" -c 'import json,sys; print(json.loads(sys.argv[1])["included_accession_count"])' "$PROFILE_INPUT_COUNTS")"
PROFILE_INPUT_WITH_SUBTYPE="$("$PYTHON_BIN" -c 'import json,sys; print(json.loads(sys.argv[1])["accessions_with_comet_subtype_count"])' "$PROFILE_INPUT_COUNTS")"
PROFILE_INPUT_WITHOUT_SUBTYPE="$("$PYTHON_BIN" -c 'import json,sys; print(json.loads(sys.argv[1])["accessions_without_comet_subtype_count"])' "$PROFILE_INPUT_COUNTS")"
PROFILE_INPUT_UNASSIGNED_GT="$("$PYTHON_BIN" -c 'import json,sys; print(json.loads(sys.argv[1])["ignored_unassigned_genotype_accession_count"])' "$PROFILE_INPUT_COUNTS")"
PROFILE_INPUT_UNASSIGNED_SUBTYPE="$("$PYTHON_BIN" -c 'import json,sys; print(json.loads(sys.argv[1])["ignored_unassigned_subtype_accession_count"])' "$PROFILE_INPUT_COUNTS")"
PROFILE_INPUT_QC_COORDINATE="$("$PYTHON_BIN" -c 'import json,sys; print(json.loads(sys.argv[1]).get("excluded_qc_flagged_accession_count", 0))' "$PROFILE_INPUT_COUNTS")"
PROFILE_INPUT_QC_HIGH_DIVERGENCE="$("$PYTHON_BIN" -c 'import json,sys; print(json.loads(sys.argv[1]).get("excluded_qc_high_divergence_accession_count", 0))' "$PROFILE_INPUT_COUNTS")"
PROFILE_INPUT_QC_NO_AA="$("$PYTHON_BIN" -c 'import json,sys; print(json.loads(sys.argv[1]).get("excluded_qc_no_aa_sequence_accession_count", 0))' "$PROFILE_INPUT_COUNTS")"

announce_step 11 "Build complete profile workbooks" \
  "amino-acid workbook: $AA_TMP_WORKBOOK" \
  "profile workbooks: $OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx; $OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx"
echo "profile_input_included_accession_count=$PROFILE_INPUT_INCLUDED_ACCESSIONS"
echo "profile_input_accessions_with_comet_subtype_count=$PROFILE_INPUT_WITH_SUBTYPE"
echo "profile_input_accessions_without_comet_subtype_count=$PROFILE_INPUT_WITHOUT_SUBTYPE"
echo "profile_input_ignored_unassigned_genotype_accession_count=$PROFILE_INPUT_UNASSIGNED_GT"
echo "profile_input_ignored_unassigned_subtype_accession_count=$PROFILE_INPUT_UNASSIGNED_SUBTYPE"
echo "profile_input_excluded_qc_coordinate_span_accession_count=$PROFILE_INPUT_QC_COORDINATE"
echo "profile_input_excluded_qc_high_divergence_accession_count=$PROFILE_INPUT_QC_HIGH_DIVERGENCE"
echo "profile_input_excluded_qc_no_aa_sequence_accession_count=$PROFILE_INPUT_QC_NO_AA"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_completeprofiles_tabspergt.py" \
  --input-workbook "$AA_TMP_WORKBOOK" \
  --output-dir "$OUTPUT_DIR" \
  --profile-accessions-csv "$OUTPUT_DIR/NS3_Profile_Accessions.csv" \
  > "$COMPLETEPROFILES_JSON"
PROFILE_INCLUDED_ACCESSIONS="$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["included_accession_count"])' "$COMPLETEPROFILES_JSON")"
PROFILE_WITH_SUBTYPE_ACCESSIONS="$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["accessions_with_comet_subtype_count"])' "$COMPLETEPROFILES_JSON")"
echo "complete_profile_included_accession_count=$PROFILE_INCLUDED_ACCESSIONS"
echo "complete_profile_accessions_with_subtype_count=$PROFILE_WITH_SUBTYPE_ACCESSIONS"

"$PYTHON_BIN" "$REPO_ROOT/export_noncomet_priority_profile_accessions/export_noncomet_priority_profile_accessions.py" \
  --profile-accessions-csv "$OUTPUT_DIR/NS3_Profile_Accessions.csv" \
  --comet-subtype-csv "$COMET_SUBTYPING_CSV" \
  --noncomet-subtype-workbook "$NONCOMET_SUBTYPE_WORKBOOK" \
  --output-csv "$OUTPUT_DIR/NS3_NonComet_Priority_Profile_Accessions.csv"

announce_step 11a "Build nucleotide distance matrices" \
  "profile accession set: $OUTPUT_DIR/NS3_Profile_Accessions.csv; classified NS3 nucleotide sequences: $AA_TMP_WORKBOOK" \
  "RAS-only distance workbooks: $OUTPUT_DIR/NS3_GT_NA_Distance_RAS.xlsx; $OUTPUT_DIR/NS3_Subtype_NA_Distance_RAS.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_na_distance_matrices.py" \
  --input-workbook "$AA_TMP_WORKBOOK" \
  --profile-accessions-csv "$OUTPUT_DIR/NS3_Profile_Accessions.csv" \
  --gt-output-xlsx "$OUTPUT_DIR/NS3_GT_NA_Distance_RAS.xlsx" \
  --subtype-output-xlsx "$OUTPUT_DIR/NS3_Subtype_NA_Distance_RAS.xlsx" \
  --min-subtype-sequences 10

announce_step 11aa "Build nucleotide position-range distance matrices" \
  "profile accession set: $OUTPUT_DIR/NS3_Profile_Accessions.csv; classified nucleotide sequences: $AA_TMP_WORKBOOK" \
  "distance workbooks: $OUTPUT_DIR/NS3_GT_NA_Distance_Pos36_175.xlsx; $OUTPUT_DIR/NS3_Subtype_NA_Distance_Pos36_175.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_na_distance_matrices.py" \
  --input-workbook "$AA_TMP_WORKBOOK" \
  --profile-accessions-csv "$OUTPUT_DIR/NS3_Profile_Accessions.csv" \
  --gt-output-xlsx "$OUTPUT_DIR/NS3_GT_NA_Distance_Pos36_175.xlsx" \
  --subtype-output-xlsx "$OUTPUT_DIR/NS3_Subtype_NA_Distance_Pos36_175.xlsx" \
  --min-subtype-sequences 10 --start 36 --end 175

announce_step 11b "Build amino-acid RAS distance matrices" \
  "profile accession set: $OUTPUT_DIR/NS3_Profile_Accessions.csv; classified NS3 amino-acid sequences: $AA_TMP_WORKBOOK" \
  "RAS-only distance workbooks: $OUTPUT_DIR/NS3_GT_AA_Distance_RAS.xlsx; $OUTPUT_DIR/NS3_Subtype_AA_Distance_RAS.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_aa_distance_matrices.py" \
  --input-workbook "$AA_TMP_WORKBOOK" \
  --profile-accessions-csv "$OUTPUT_DIR/NS3_Profile_Accessions.csv" \
  --gt-output-xlsx "$OUTPUT_DIR/NS3_GT_AA_Distance_RAS.xlsx" \
  --subtype-output-xlsx "$OUTPUT_DIR/NS3_Subtype_AA_Distance_RAS.xlsx" \
  --min-subtype-sequences 10

announce_step 12 "Export consensus FASTA files" \
  "profile workbooks: $OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx; $OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  "consensus FASTAs: $OUTPUT_DIR/NS3_GT_Consensus.fasta; $OUTPUT_DIR/NS3_Subtype_Consensus.fasta"
"$PYTHON_BIN" "$SCRIPT_DIR/export_ns3_consensus_fasta.py" \
  --gt-profile-workbook "$OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx" \
  --subtype-profile-workbook "$OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  --output-dir "$OUTPUT_DIR" \
  > /dev/null

announce_step 12a "Compare subtype consensus to subtype references" \
  "subtype profile: $OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx; references: $SUBTYPE_JSON" \
  "RAS-only distance report: $OUTPUT_DIR/NS3_Subtype_Consensus_Reference_AA_Distance_RAS.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_subtype_consensus_reference_distance.py" \
  --subtype-profile-workbook "$OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  --subtype-json "$SUBTYPE_JSON" \
  --output-xlsx "$OUTPUT_DIR/NS3_Subtype_Consensus_Reference_AA_Distance_RAS.xlsx"

announce_step 13 "Build genotype RAS profile" \
  "genotype profile workbook: $OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx" \
  "RAS profile: $OUTPUT_DIR/NS3_GT_RAS_Profiles.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_gt_ras_profiles.py" \
  --gt-profile-workbook "$OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$GT_RAS_JSON"

announce_step 14 "Build subtype RAS profile" \
  "subtype profile workbook: $OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  "RAS profile: $OUTPUT_DIR/NS3_Subtype_RAS_Profiles.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_subtype_ras_profiles.py" \
  --subtype-profile-workbook "$OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$SUBTYPE_RAS_JSON"

announce_step 15 "Combine genotype and subtype RAS profiles" \
  "GT RAS profile: $OUTPUT_DIR/NS3_GT_RAS_Profiles.xlsx; subtype RAS profile: $OUTPUT_DIR/NS3_Subtype_RAS_Profiles.xlsx" \
  "combined RAS profile: $OUTPUT_DIR/NS3_Combined_RAS_Profiles.xlsx"
echo "Action: place each genotype RAS consensus row before its subtype rows."
echo "Result: subtype cells retain amino-acid variants strictly above 10%; each genotype block ends with a blank row."
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_combined_ras_profiles.py" \
  --gt-ras-profile-workbook "$OUTPUT_DIR/NS3_GT_RAS_Profiles.xlsx" \
  --subtype-ras-profile-workbook "$OUTPUT_DIR/NS3_Subtype_RAS_Profiles.xlsx" \
  --output-xlsx "$OUTPUT_DIR/NS3_Combined_RAS_Profiles.xlsx" \
  > "$COMBINED_RAS_JSON"
"$PYTHON_BIN" "$REPO_ROOT/replace_comet_profile_coverage_range_with_mean_diff/replace_comet_profile_coverage_range_with_mean_diff.py" \
  --combined-profile-workbook "$OUTPUT_DIR/NS3_Combined_RAS_Profiles.xlsx" \
  --profile-input-workbook "$AA_TMP_WORKBOOK" \
  --profile-accessions-csv "$OUTPUT_DIR/NS3_Profile_Accessions.csv"
"$PYTHON_BIN" "$SCRIPT_DIR/build_comet_subtype_ras_coverage_report.py" \
  --gene NS3 --combined-profile-workbook "$OUTPUT_DIR/NS3_Combined_RAS_Profiles.xlsx" \
  --profile-input-workbook "$AA_TMP_WORKBOOK" --profile-accessions-csv "$OUTPUT_DIR/NS3_Profile_Accessions.csv" \
  --output-xlsx "$OUTPUT_DIR/NS3_Subtype_RAS_Coverage_Report.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_comet_workflow_sequence_audit.py" \
  --gene NS3 --selection-workbook "$EXCEL_FILE" --selection-sheet "$SHEET_NAME" --fasta-dir "$FASTA_POOL" \
  --metadata-csv "$ACCESSIONS_METADATA_CSV" --comet-csv "$COMET_SUBTYPING_CSV" --qc-workbook "$AA_TMP_WORKBOOK" \
  --profile-accessions-csv "$OUTPUT_DIR/NS3_Profile_Accessions.csv" --combined-profile-workbook "$OUTPUT_DIR/NS3_Combined_RAS_Profiles.xlsx" \
  --output-xlsx "$OUTPUT_DIR/NS3_Workflow_Sequence_Audit.xlsx"

announce_step 15a "Summarize subtype RAS differences from genotype consensus" \
  "combined RAS profile: $OUTPUT_DIR/NS3_Combined_RAS_Profiles.xlsx; profile sequences: $AA_TMP_WORKBOOK" \
  "per-subtype mean and median difference workbook: $OUTPUT_DIR/NS3_Subtype_RAS_Consensus_Difference_Summary.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_subtype_ras_consensus_difference_summary.py" \
  --combined-profile-workbook "$OUTPUT_DIR/NS3_Combined_RAS_Profiles.xlsx" \
  --profile-input-workbook "$AA_TMP_WORKBOOK" \
  --profile-accessions-csv "$OUTPUT_DIR/NS3_Profile_Accessions.csv" \
  --genotype-consensus-fasta "$OUTPUT_DIR/NS3_GT_Consensus.fasta" \
  --output-xlsx "$OUTPUT_DIR/NS3_Subtype_RAS_Consensus_Difference_Summary.xlsx" \
  --output-png "$OUTPUT_DIR/NS3_Subtype_RAS_Consensus_Difference_Trend.png" \
  > "$RAS_CONSENSUS_DIFFERENCE_JSON"

announce_step 15b "Build subtype 5a profile coverage report" \
  "profile accession set: $OUTPUT_DIR/NS3_Profile_Accessions.csv; QC-approved amino-acid sequences: $AA_TMP_WORKBOOK" \
  "accession-level, per-position coverage, and chart: $OUTPUT_DIR/NS3_GT5_5a_Profile_Coverage.xlsx; $OUTPUT_DIR/NS3_GT5_5a_Profile_Position_Coverage.csv; $OUTPUT_DIR/NS3_GT5_5a_Profile_Position_Coverage.png"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_subtype_profile_coverage_report.py" \
  --profile-input-workbook "$AA_TMP_WORKBOOK" \
  --profile-accessions-csv "$OUTPUT_DIR/NS3_Profile_Accessions.csv" \
  --genotype 5 \
  --subtype 5a \
  --range-start 1 \
  --range-end 631 \
  --output-xlsx "$OUTPUT_DIR/NS3_GT5_5a_Profile_Coverage.xlsx" \
  --position-coverage-csv "$OUTPUT_DIR/NS3_GT5_5a_Profile_Position_Coverage.csv" \
  --position-coverage-png "$OUTPUT_DIR/NS3_GT5_5a_Profile_Position_Coverage.png"

announce_step 16 "Build paired amino-acid distance matrices" \
  "profile accession set: $OUTPUT_DIR/NS3_Profile_Accessions.csv; classified amino-acid sequences: $AA_TMP_WORKBOOK" \
  "distance workbooks: $OUTPUT_DIR/NS3_GT_AA_Distance_Pos36_175.xlsx; $OUTPUT_DIR/NS3_Subtype_AA_Distance_Pos36_175.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_aa_distance_matrices.py" \
  --input-workbook "$AA_TMP_WORKBOOK" \
  --profile-accessions-csv "$OUTPUT_DIR/NS3_Profile_Accessions.csv" \
  --gt-output-xlsx "$OUTPUT_DIR/NS3_GT_AA_Distance_Pos36_175.xlsx" \
  --subtype-output-xlsx "$OUTPUT_DIR/NS3_Subtype_AA_Distance_Pos36_175.xlsx" \
  --min-subtype-sequences 10 --start 36 --end 175 \
  > "$GT_AA_DISTANCE_SUMMARY"

announce_step 18 "Build RAS entropy reports" \
  "profile workbooks: $OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx; $OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  "entropy workbooks: $OUTPUT_DIR/NS3_GT_RAS_Entropy.xlsx; $OUTPUT_DIR/NS3_Subtype_RAS_Entropy.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_ras_entropy.py" \
  --gt-profile-workbook "$OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx" \
  --subtype-profile-workbook "$OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  --gt-output-xlsx "$OUTPUT_DIR/NS3_GT_RAS_Entropy.xlsx" \
  --subtype-output-xlsx "$OUTPUT_DIR/NS3_Subtype_RAS_Entropy.xlsx" \
  > "$RAS_ENTROPY_SUMMARY"

announce_step 19 "Add combined-profile genotype-consensus difference row" \
  "combined RAS profile: $OUTPUT_DIR/NS3_Combined_RAS_Profiles.xlsx; QC-passed profile sequences: $AA_TMP_WORKBOOK" \
  "updated combined RAS profile: $OUTPUT_DIR/NS3_Combined_RAS_Profiles.xlsx"
"$PYTHON_BIN" "$SCRIPT_DIR/add_combined_profile_nonconsensus_row.py" \
  --combined-profile-workbook "$OUTPUT_DIR/NS3_Combined_RAS_Profiles.xlsx" --profile-input-workbook "$AA_TMP_WORKBOOK" \
  --profile-accessions-csv "$OUTPUT_DIR/NS3_Profile_Accessions.csv" --genotype-consensus-fasta "$OUTPUT_DIR/NS3_GT_Consensus.fasta"

announce_step 20 "Publish shared ICTV reference/consensus comparison report" \
  "all COMET subtype profiles and reference comparison workbooks" \
  "shared report: outputs/shared_report/ICTV_ref_local_cons_compare"
"$PYTHON_BIN" "$SCRIPT_DIR/add_subtype_consensus_mutation_summaries.py"

echo "NS3 pipeline complete"
echo "matched_fasta_report=${MATCHED_TXT#$REPO_ROOT/}"
echo "included_fasta_dir=${INCLUDED_FASTA_DIR#$REPO_ROOT/}"
echo "output_dir=${OUTPUT_DIR#$REPO_ROOT/}"
echo "temp_root=${TEMP_ROOT#$REPO_ROOT/}"
