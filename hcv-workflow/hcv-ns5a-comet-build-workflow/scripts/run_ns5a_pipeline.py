#!/usr/bin/env python3
"""Run named stages of the HCV NS5A COMET build workflow.

Examples:
  scripts/run_ns5a_pipeline.py --list-steps
  scripts/run_ns5a_pipeline.py
  scripts/run_ns5a_pipeline.py --step discover-refid-fastas
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tomllib
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
SKILL_NAME = "hcv-ns5a-comet-build-workflow"
RAS_POSITIONS = "24,26,28,29,30,31,32,38,58,62,92,93"
RANGE_POSITIONS = ",".join(str(position) for position in range(24, 94))
OPTIONAL_STEPS = {"build-source-features", "build-grouped-source-features"}


@dataclass(frozen=True)
class Step:
    name: str
    description: str
    action: Callable[[], None]


def parse_dotenv(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[7:].lstrip()
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
            value = value[1:-1]
        if key:
            values[key] = value
    return values


def load_config(config_path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if config_path.is_file():
        with config_path.open("rb") as handle:
            config = tomllib.load(handle)
        for section_name in ("common", "ns5a_comet"):
            section = config.get(section_name, {})
            if isinstance(section, dict):
                values.update({key.upper(): str(value) for key, value in section.items()})
    return values


def path_value(value: str, *, default: Path | None = None) -> Path:
    path = Path(value).expanduser() if value else default
    if path is None:
        raise RuntimeError("A required path value is missing")
    return path if path.is_absolute() else REPO_ROOT / path


def executable_value(value: str) -> str:
    candidate = Path(value).expanduser()
    if candidate.is_absolute() or "/" in value:
        return str(candidate if candidate.is_absolute() else REPO_ROOT / candidate)
    return value


def repository_relative(value: str | Path) -> str:
    path = Path(value).expanduser()
    if not path.is_absolute():
        return str(path)
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list-steps", action="store_true", help="Print available step names and exit")
    parser.add_argument("--step", action="append", help="Run only this named step; repeat to run multiple steps")
    parser.add_argument("--include-source-features", action="store_true", help="Include disabled source-feature steps in a full run")
    parser.add_argument("--config", type=Path, default=REPO_ROOT / "pipeline.local.toml")
    parser.add_argument("--excel-file")
    parser.add_argument("--sheet-name")
    parser.add_argument("--fasta-pool")
    parser.add_argument("--genbank-dir")
    parser.add_argument("--output-dir")
    parser.add_argument("--reference-fasta")
    parser.add_argument("--subtype-json")
    parser.add_argument("--gt-aa-json")
    parser.add_argument("--accessions-metadata-csv")
    parser.add_argument("--comet-subtyping-csv")
    parser.add_argument("--noncomet-subtype-workbook")
    parser.add_argument("--temp-root")
    parser.add_argument("--python-bin")
    return parser.parse_args()


class Pipeline:
    def __init__(self, args: argparse.Namespace) -> None:
        dotenv = parse_dotenv(REPO_ROOT / ".env")
        config = load_config(args.config.expanduser())
        values = {**config, **dotenv, **os.environ}

        def value(argument: str, environment: str, default: str = "") -> str:
            supplied = getattr(args, argument)
            return supplied if supplied is not None else values.get(environment, default)

        self.python_bin = executable_value(value("python_bin", "PYTHON_BIN", str(REPO_ROOT / ".venv/bin/python")))
        self.excel_file = path_value(value("excel_file", "EXCEL_FILE"))
        self.sheet_name = value("sheet_name", "SHEET_NAME")
        self.fasta_pool = path_value(value("fasta_pool", "FASTA_POOL"))
        self.genbank_dir = path_value(value("genbank_dir", "GENBANK_DIR"))
        self.output_dir = path_value(value("output_dir", "OUTPUT_DIR", "outputs/comet"))
        self.reference_fasta = path_value(value("reference_fasta", "REFERENCE_FASTA", "HCVData/HCV_GT_RefSeqs.fasta"))
        self.subtype_json = path_value(value("subtype_json", "SUBTYPE_JSON", "HCVData/HCV_Subtype_Refs_By_Genome_NA.json"))
        self.gt_aa_json = path_value(value("gt_aa_json", "GT_AA_JSON", "HCVData/HCV_GT_Refs_By_Gene_AA.json"))
        self.metadata_csv = path_value(value("accessions_metadata_csv", "ACCESSIONS_METADATA_CSV", "HCVData/Accessions_metadata.csv"))
        self.comet_csv = path_value(value("comet_subtyping_csv", "COMET_SUBTYPING_CSV", "HCVData/Comet Subtyping/NS5A.csv"))
        self.noncomet_workbook = path_value(value("noncomet_subtype_workbook", "NONCOMET_SUBTYPE_WORKBOOK", "outputs/local_alignment/NS5A_Subtype_AllStudies_WSeqs.xlsx"))
        self.temp_root = path_value(value("temp_root", "TEMP_ROOT", f"outputs/temp/{SKILL_NAME}/run_ns5a_pipeline"))
        self.included_fasta_dir = self.temp_root / "included_refid_fastas"
        self.discovery_tmp = self.temp_root / "find_refid_fastas"
        self.skill_temp_root = REPO_ROOT / "outputs/temp" / SKILL_NAME
        self.refid_metadata_dir = self.temp_root / "refid_metadata"
        self.matched_txt = self.output_dir / "NS5A_matched_fasta_files.txt"
        self.aa_workbook = self.output_dir / "NS5A_Profile_Input_Alignment_QC.xlsx"
        self.profile_accessions_csv = self.output_dir / "NS5A_Profile_Accessions.csv"
    def ensure_summary_directories(self) -> None:
        self.temp_root.mkdir(parents=True, exist_ok=True)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        for name in (
            "build_ns5a_gt_allstudies", "build_ns5a_sourcefeatures_csv", "build_ns5a_sourcefeatures_grouped_csv",
            "build_ns5a_subtype_allstudies_wseqs", "build_ns5a_subtype_with_gt_aa", "build_ns5a_completeprofiles_tabspergt",
            "build_ns5a_gt_ras_profiles", "build_ns5a_subtype_ras_profiles", "build_ns5a_combined_ras_profiles",
            "build_ns5a_subtype_ras_consensus_difference_summary", "build_ns5a_gt_aa_distance_matrix",
            "build_ns5a_subtype_aa_distance_matrices",
        ):
            (self.skill_temp_root / name).mkdir(parents=True, exist_ok=True)

    def run(self, script: str, *arguments: str, stdout_path: Path | None = None) -> str:
        command = [self.python_bin, repository_relative(SCRIPT_DIR / script), *(repository_relative(argument) for argument in arguments)]
        if stdout_path is None:
            subprocess.run(command, check=True, cwd=REPO_ROOT)
            return ""
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        result = subprocess.run(command, check=True, text=True, capture_output=True, cwd=REPO_ROOT)
        stdout_path.write_text(result.stdout, encoding="utf-8")
        if result.stderr:
            print(result.stderr, end="", file=sys.stderr)
        return result.stdout

    def announce(self, name: str, description: str) -> None:
        print(f"\n=== {name}: {description} ===")

    def steps(self) -> list[Step]:
        summary = lambda directory, filename="last_run_summary.json": self.skill_temp_root / directory / filename
        gt_workbook = self.output_dir / "NS5A_GT_AllStudies.xlsx"
        subtype_workbook = self.output_dir / "NS5A_Subtype_AllStudies_WSeqs.xlsx"
        gt_profile = self.output_dir / "NS5A_GT_CompleteProfiles_TabsPerGT.xlsx"
        subtype_profile = self.output_dir / "NS5A_Subtype_CompleteProfiles_TabsPerGT.xlsx"
        gt_consensus = self.output_dir / "NS5A_GT_Consensus.fasta"
        subtype_consensus = self.output_dir / "NS5A_Subtype_Consensus.fasta"
        gt_ras = self.output_dir / "NS5A_GT_RAS_Profiles.xlsx"
        subtype_ras = self.output_dir / "NS5A_Subtype_RAS_Profiles.xlsx"
        combined_ras = self.output_dir / "NS5A_Combined_RAS_Profiles.xlsx"
        comet_gt_csv = self.temp_root / "comet_ns5a_genotype_assignments.csv"
        comet_subtype_csv = self.temp_root / "comet_ns5a_subtype_assignments.csv"

        def prepare() -> None:
            self.run("prepare_ns5a_pipeline_workdirs.py", "--clean-dir", self.included_fasta_dir, "--clean-dir", self.discovery_tmp, "--remove-file", summary("build_ns5a_sourcefeatures_csv", "NS5A_SourceFeatures.csv"), "--remove-file", summary("build_ns5a_sourcefeatures_grouped_csv", "NS5A_SourceFeatures_Grouped.csv"))

        def discover() -> None:
            output = self.run("find_refid_fastas.py", "--excel-file", self.excel_file, "--sheet", self.sheet_name, "--fasta-dir", self.fasta_pool, "--output-dir", self.discovery_tmp, stdout_path=self.discovery_tmp / "discovery_ns5a.json")
            del output
            matches = sorted(self.discovery_tmp.glob("refid_fasta_*/matched_fasta_files.txt"))
            if len(matches) != 1:
                raise RuntimeError(f"Expected one matched FASTA file list under {self.discovery_tmp}, found {len(matches)}")
            self.matched_txt.write_text(matches[0].read_text(encoding="utf-8"), encoding="utf-8")

        def prepare_comet() -> None:
            self.run("prepare_comet_ns5a_assignments.py", "--comet-csv", self.comet_csv, "--fasta-dir", self.included_fasta_dir, "--genotype-output-csv", comet_gt_csv, "--subtype-output-csv", comet_subtype_csv, "--not-found-output-csv", self.temp_root / "comet_ns5a_not_found_or_unassigned.csv", "--not-found-fasta-output", self.temp_root / "comet_ns5a_not_found_or_unassigned.fasta", "--remove-unassigned")

        def extract_aa() -> None:
            self.run("build_ns5a_subtype_with_gt_aa.py", "--subtype-workbook", subtype_workbook, "--fasta-dir", self.fasta_pool, "--gt-aa-json", self.gt_aa_json, "--output-dir", self.output_dir, "--output-workbook", self.output_dir / "NS5A_Profile_Input_Source.xlsx", stdout_path=summary("build_ns5a_subtype_with_gt_aa"))

        def report_profile_counts() -> None:
            output = self.run("build_ns5a_completeprofiles_tabspergt.py", "--input-workbook", self.aa_workbook, "--report-only", stdout_path=self.temp_root / "profile_input_counts.json")
            print(output, end="")

        def paired_distances() -> None:
            for sequence_type in ("aa", "na"):
                upper = sequence_type.upper()
                script = f"build_ns5a_{sequence_type}_distance_matrices.py"
                for positions, suffix in ((RAS_POSITIONS, "RAS"), (RANGE_POSITIONS, "Pos24_93")):
                    self.run(script, "--input-workbook", self.aa_workbook, "--profile-accessions-csv", self.profile_accessions_csv, "--sequence-type", sequence_type, "--positions", positions, "--gt-output-xlsx", self.output_dir / f"NS5A_GT_{upper}_Distance_{suffix}.xlsx", "--subtype-output-xlsx", self.output_dir / f"NS5A_Subtype_{upper}_Distance_{suffix}.xlsx", "--min-subtype-sequences", "10")

        return [
            Step("prepare-workdirs", "recreate run directories and clear stale source-feature files", prepare),
            Step("discover-refid-fastas", "discover selected RefID FASTA files", discover),
            Step("stage-refid-fastas", "copy discovered FASTA files into the staging directory", lambda: self.run("stage_matched_refid_fastas.py", "--matched-files", self.matched_txt, "--output-dir", self.included_fasta_dir)),
            Step("filter-accession-metadata", "filter master metadata to staged FASTA accessions", lambda: self.run("filter_accessions_metadata_by_fasta.py", "--fasta-dir", self.included_fasta_dir, "--metadata-csv", self.metadata_csv, "--output-dir", self.temp_root)),
            Step("split-refid-metadata", "create per-RefID metadata filters", lambda: (self.run("prepare_ns5a_pipeline_workdirs.py", "--clean-dir", self.refid_metadata_dir), self.run("split_refid_metadata_csv.py", "--input-csv", self.temp_root / "included_accessions_metadata.csv", "--output-dir", self.refid_metadata_dir))),
            Step("filter-refid-fastas", "filter staged FASTA records by RefID metadata", lambda: self.run("filter_refid_fastas_by_metadata.py", "--metadata-dir", self.refid_metadata_dir, "--fasta-dir", self.included_fasta_dir)),
            Step("prepare-comet-assignments", "create COMET calls and remove unassigned records", prepare_comet),
            Step("build-genotype-workbook", "build the COMET genotype workbook", lambda: self.run("build_ns5a_comet_gt_allstudies.py", "--fasta-dir", self.included_fasta_dir, "--comet-genotype-csv", comet_gt_csv, "--output-dir", self.output_dir, stdout_path=summary("build_ns5a_gt_allstudies"))),
            Step("add-genotype-counts", "add the genotype-count worksheet", lambda: self.run("add_gt_counts_sheet.py", "--workbook", gt_workbook)),
            Step("build-source-features", "build optional source-feature CSV output", lambda: self.run("build_ns5a_sourcefeatures_csv.py", "--matched-fasta-report", self.matched_txt, "--genbank-dir", self.genbank_dir, stdout_path=summary("build_ns5a_sourcefeatures_csv"))),
            Step("build-grouped-source-features", "build optional grouped source-feature output", lambda: self.run("build_ns5a_sourcefeatures_grouped_csv.py", "--gt-workbook", gt_workbook, "--summary-xlsx", self.output_dir / "NS5A_NumSeqs_Naive_1PP_CoversRAS_ByStudy.xlsx", stdout_path=summary("build_ns5a_sourcefeatures_grouped_csv"))),
            Step("build-subtype-workbook", "build COMET subtype workbook with non-COMET priorities", lambda: self.run("build_ns5a_comet_subtype_allstudies.py", "--genotype-workbook", gt_workbook, "--comet-subtype-csv", comet_subtype_csv, "--noncomet-subtype-workbook", self.noncomet_workbook, "--output-dir", self.output_dir, stdout_path=summary("build_ns5a_subtype_allstudies_wseqs"))),
            Step("extract-profile-aa", "extract genotype-position amino-acid sequences", extract_aa),
            Step("validate-profile-alignment", "write alignment QC columns and flagged-accession report", lambda: self.run("validate_ns5a_profile_alignment.py", "--input-workbook", self.output_dir / "NS5A_Profile_Input_Source.xlsx", "--gt-aa-json", self.gt_aa_json, "--output-workbook", self.aa_workbook, "--flagged-accessions-csv", self.output_dir / "NS5A_Profile_Alignment_QC_Flagged_Accessions.csv", stdout_path=self.skill_temp_root / "validate_ns5a_profile_alignment.json")),
            Step("summarize-qc-mutation-burden", "summarize QC-passed genotype mutation burden", lambda: self.run("build_qc_passed_genotype_mutation_burden_summary.py", "--input-workbook", self.aa_workbook, "--output-csv", self.output_dir / "NS5A_QC_Passed_Genotype_Mutation_Burden_Summary.csv")),
            Step("report-profile-input-counts", "report profile-input inclusion counts", report_profile_counts),
            Step("build-complete-profiles", "build genotype and subtype complete profiles", lambda: self.run("build_ns5a_completeprofiles_tabspergt.py", "--input-workbook", self.aa_workbook, "--output-dir", self.output_dir, "--profile-accessions-csv", self.profile_accessions_csv, stdout_path=summary("build_ns5a_completeprofiles_tabspergt"))),
            Step("export-noncomet-priority-accessions", "export non-COMET priority profile accessions", lambda: self.run("export_noncomet_priority_profile_accessions.py", "--profile-accessions-csv", self.profile_accessions_csv, "--comet-subtype-csv", self.comet_csv, "--noncomet-subtype-workbook", self.noncomet_workbook, "--output-csv", self.output_dir / "NS5A_NonComet_Priority_Profile_Accessions.csv")),
            Step("export-consensus-fastas", "export genotype and subtype consensus FASTA files", lambda: self.run("export_ns5a_consensus_fasta.py", "--gt-profile-workbook", gt_profile, "--subtype-profile-workbook", subtype_profile, "--output-dir", self.output_dir)),
            Step("compare-reference-consensus", "compare COMET consensus sequences to references", lambda: (self.run("export_gt_reference_consensus_differences.py", "--gene", "NS5A_NTD", "--reference-fasta", REPO_ROOT / "HCVData/HCV_GT_Refs_NS3_NS5A_NTD_NS5B_AA.fasta", "--subtype-reference-dir", REPO_ROOT / "HCVData/Reference_seqs"), self.run("build_ns5a_subtype_consensus_reference_distance.py", "--subtype-profile-workbook", subtype_profile, "--subtype-json", self.subtype_json, "--output-xlsx", self.output_dir / "NS5A_Subtype_Consensus_Reference_AA_Distance_RAS.xlsx"))),
            Step("build-genotype-ras-profile", "build genotype RAS profile", lambda: self.run("build_ns5a_gt_ras_profiles.py", "--gt-profile-workbook", gt_profile, "--gt-aa-json", self.gt_aa_json, "--output-dir", self.output_dir, stdout_path=summary("build_ns5a_gt_ras_profiles"))),
            Step("build-subtype-ras-profile", "build subtype RAS profile", lambda: self.run("build_ns5a_subtype_ras_profiles.py", "--subtype-profile-workbook", subtype_profile, "--gt-aa-json", self.gt_aa_json, "--output-dir", self.output_dir, stdout_path=summary("build_ns5a_subtype_ras_profiles"))),
            Step("build-combined-ras-reports", "build combined RAS, coverage, and sequence-audit reports", lambda: (self.run("build_ns5a_combined_ras_profiles.py", "--gt-ras-profile-workbook", gt_ras, "--subtype-ras-profile-workbook", subtype_ras, "--output-xlsx", combined_ras, stdout_path=summary("build_ns5a_combined_ras_profiles")), self.run("build_comet_subtype_ras_coverage_report.py", "--gene", "NS5A", "--combined-profile-workbook", combined_ras, "--profile-input-workbook", self.aa_workbook, "--profile-accessions-csv", self.profile_accessions_csv, "--output-xlsx", self.output_dir / "NS5A_Subtype_RAS_Coverage_Report.xlsx"), self.run("build_comet_workflow_sequence_audit.py", "--gene", "NS5A", "--selection-workbook", self.excel_file, "--selection-sheet", self.sheet_name, "--fasta-dir", self.fasta_pool, "--metadata-csv", self.metadata_csv, "--comet-csv", self.comet_csv, "--qc-workbook", self.aa_workbook, "--profile-accessions-csv", self.profile_accessions_csv, "--combined-profile-workbook", combined_ras, "--combined-subtype-cutoff", "5", "--combined-subtype-cutoff-exclusive", "--output-xlsx", self.output_dir / "NS5A_Workflow_Sequence_Audit.xlsx"))),
            Step("summarize-subtype-ras-differences", "summarize subtype RAS differences from genotype consensus", lambda: self.run("build_ns5a_subtype_ras_consensus_difference_summary.py", "--combined-profile-workbook", combined_ras, "--profile-input-workbook", self.aa_workbook, "--profile-accessions-csv", self.profile_accessions_csv, "--genotype-consensus-fasta", gt_consensus, "--output-xlsx", self.output_dir / "NS5A_Subtype_RAS_Consensus_Difference_Summary.xlsx", "--output-png", self.output_dir / "NS5A_Subtype_RAS_Consensus_Difference_Trend.png", stdout_path=summary("build_ns5a_subtype_ras_consensus_difference_summary"))),
            Step("build-genotype-aa-distance", "build genotype amino-acid distance matrix", lambda: self.run("build_ns5a_gt_aa_distance_matrix.py", "--input-fasta", gt_consensus, "--aligned-fasta", self.output_dir / "NS5A_GT_Consensus_aligned.fasta", "--output-xlsx", self.output_dir / "NS5A_GT_AA_Distance_Pos24_93.xlsx", "--details-xlsx", summary("build_ns5a_gt_aa_distance_matrix", "NS5A_GT_AA_Distance_Pos24_93_details.xlsx"), "--start", "24", "--end", "93", stdout_path=summary("build_ns5a_gt_aa_distance_matrix", "last_run_summary.txt"))),
            Step("build-subtype-aa-distance", "build subtype amino-acid distance matrices", lambda: self.run("build_ns5a_subtype_aa_distance_matrices.py", "--input-fasta", subtype_consensus, "--subtype-profile-workbook", subtype_profile, "--min-subtype-sequences", "10", "--output-xlsx", self.output_dir / "NS5A_Subtype_AA_Distance_Pos24_93.xlsx", "--temp-dir", summary("build_ns5a_subtype_aa_distance_matrices").parent, "--start", "24", "--end", "93", stdout_path=summary("build_ns5a_subtype_aa_distance_matrices", "last_run_summary.txt"))),
            Step("build-paired-distance-matrices", "build paired AA and NA RAS/range distance matrices", paired_distances),
            Step("build-ras-entropy", "build genotype and subtype RAS entropy reports", lambda: self.run("build_ns5a_ras_entropy.py", "--gt-profile-workbook", gt_profile, "--subtype-profile-workbook", subtype_profile, "--gt-output-xlsx", self.output_dir / "NS5A_GT_RAS_Entropy.xlsx", "--subtype-output-xlsx", self.output_dir / "NS5A_Subtype_RAS_Entropy.xlsx", stdout_path=summary("build_ns5a_ras_entropy", "last_run_summary.txt"))),
            Step("add-nonconsensus-row", "add combined-profile genotype-consensus difference rows", lambda: self.run("add_combined_profile_nonconsensus_row.py", "--combined-profile-workbook", combined_ras, "--profile-input-workbook", self.aa_workbook, "--profile-accessions-csv", self.profile_accessions_csv, "--genotype-consensus-fasta", gt_consensus)),
            Step("update-coverage-labels", "replace combined-profile coverage labels", lambda: self.run("replace_comet_profile_coverage_range_with_mean_diff.py", "--combined-profile-workbook", combined_ras, "--profile-input-workbook", self.aa_workbook, "--profile-accessions-csv", self.profile_accessions_csv)),
            Step("publish-ictv-report", "publish the shared ICTV comparison report", lambda: self.run("add_subtype_consensus_mutation_summaries.py")),
        ]


def main() -> int:
    args = parse_args()
    pipeline = Pipeline(args)
    steps = pipeline.steps()
    step_by_name = {step.name: step for step in steps}
    if args.list_steps:
        for step in steps:
            print(step.name)
        return 0

    requested = args.step or [step.name for step in steps if args.include_source_features or step.name not in OPTIONAL_STEPS]
    unknown = [name for name in requested if name not in step_by_name]
    if unknown:
        raise SystemExit(f"Unknown step name(s): {', '.join(unknown)}. Use --list-steps.")
    if not args.excel_file and not pipeline.excel_file.is_file() and "discover-refid-fastas" in requested:
        raise SystemExit(f"Excel file was not found: {pipeline.excel_file}")
    if not pipeline.sheet_name and "discover-refid-fastas" in requested:
        raise SystemExit("SHEET_NAME or --sheet-name is required for discover-refid-fastas")
    pipeline.ensure_summary_directories()
    for name in requested:
        step = step_by_name[name]
        pipeline.announce(step.name, step.description)
        step.action()
    print("NS5A pipeline complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
