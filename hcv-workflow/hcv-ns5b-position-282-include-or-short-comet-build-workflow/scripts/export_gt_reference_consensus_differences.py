#!/usr/bin/env python3
"""Export aligned amino acids between GT references and Comet GT consensuses."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Font


GENES = ("NS3", "NS5A_NTD", "NS5B")
VALID_AAS = set("ACDEFGHIKLMNPQRSTVWY*")
RAS_POSITIONS = {
    "NS3": (36, 41, 43, 54, 55, 56, 80, 122, 155, 156, 158, 166, 168, 170, 175),
    "NS5A_NTD": (24, 26, 28, 29, 30, 31, 32, 38, 58, 62, 92, 93),
    "NS5B": (150, 159, 206, 282, 316, 320, 321),
}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks).upper()
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line)
    if header is not None:
        records[header] = "".join(chunks).upper()
    return records


def load_references(path: Path) -> dict[tuple[str, str], str]:
    references: dict[tuple[str, str], str] = {}
    for header, sequence in read_fasta(path).items():
        match = re.fullmatch(r"genotype ([1-8]) \| (NS3|NS5A_NTD|NS5B)", header)
        if match:
            references[(match.group(1), match.group(2))] = sequence
    missing = [
        f"GT{gt} {gene}"
        for gt in "12345678"
        for gene in GENES
        if (gt, gene) not in references
    ]
    if missing:
        raise ValueError(f"Missing reference sequences: {', '.join(missing)}")
    return references


def covered_alignment_pairs(
    reference: str, consensus: str
) -> dict[int, tuple[str, str]]:
    """Compare fixed gene coordinates; inputs are already position-aligned."""
    return {
        position: (reference_aa, consensus_aa)
        for position, (reference_aa, consensus_aa) in enumerate(
            zip(reference, consensus), start=1
        )
        if reference_aa in VALID_AAS and consensus_aa in VALID_AAS
    }


def header_fields(header: str) -> dict[str, str]:
    return dict(field.split("=", 1) for field in header.split("|") if "=" in field)


def write_excel(output_path: Path, rows: list[list[str | int]]) -> None:
    """Write a RAS-only comparison workbook."""
    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "RAS comparison"
    for row in rows:
        worksheet.append(row)

    for row in worksheet.iter_rows():
        if row[0].value == "Gene":
            for cell in row:
                cell.font = Font(bold=True)
    worksheet.freeze_panes = "A2"
    for column_cells in worksheet.columns:
        width = max((len(str(cell.value or "")) for cell in column_cells), default=0)
        worksheet.column_dimensions[column_cells[0].column_letter].width = min(
            width + 2, 30
        )
    workbook.save(output_path)


def write_subtype_alignment(
    gene: str, reference_path: Path, consensus_path: Path, output_path: Path
) -> tuple[int, list[tuple[str, str, str, str]], list[tuple[str, str, int, int]]]:
    consensuses = read_fasta(consensus_path)
    rows: list[tuple[str, str, str, dict[int, tuple[str, str]]]] = []
    missing_consensuses: list[tuple[str, str, str, str]] = []
    matched_lengths: list[tuple[str, str, int, int]] = []
    seen_subtypes: set[tuple[str, str]] = set()
    for header, reference in read_fasta(reference_path).items():
        fields = header_fields(header)
        genotype = fields.get("genotype", "")
        subtype = fields.get("subtype", "")
        accession = fields.get("accession", "")
        subtype_key = (genotype, subtype)
        if genotype and subtype and subtype_key in seen_subtypes:
            continue
        if genotype and subtype:
            seen_subtypes.add(subtype_key)
        consensus_name = f"GT{genotype}_{subtype}"
        consensus = consensuses.get(consensus_name)
        if not genotype or not subtype or consensus is None:
            missing_consensuses.append((gene, f"GT{genotype}", subtype, accession))
            continue
        rows.append(
            (
                genotype,
                subtype,
                accession,
                covered_alignment_pairs(reference, consensus),
            )
        )
        matched_lengths.append((header, consensus_name, len(reference), len(consensus)))
    rows.sort(key=lambda row: (int(row[0]), row[1], row[2]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    excel_rows: list[list[str | int]] = []
    for genotype, subtype, _accession, pairs in rows:
        positions = RAS_POSITIONS[gene]
        excel_rows.append([])
        excel_rows.append(["Gene", "Genotype", "Subtype", "Sequence"] + list(positions))
        excel_rows.append(
            [gene, f"GT{genotype}", subtype, "Reference"]
            + [
                pairs[position][0] if position in pairs else ""
                for position in positions
            ]
        )
        excel_rows.append(
            [gene, f"GT{genotype}", subtype, "Consensus"]
            + [
                pairs[position][1]
                if position in pairs and pairs[position][1] != pairs[position][0]
                else ""
                for position in positions
            ]
        )
    write_excel(output_path, excel_rows)
    return len(rows), missing_consensuses, matched_lengths


def write_subtype_length_csv(path: Path, rows: list[tuple[str, str, int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ref_name", "cons_name", "ref_aa_length", "cons_aa_length"])
        writer.writerows(rows)


def write_missing_subtype_consensus_report(
    output_path: Path, rows: list[tuple[str, str, str, str]]
) -> None:
    """Write all subtype references that lack a profile-derived Comet consensus."""
    sorted_rows = sorted(
        rows,
        key=lambda row: (row[0], int(row[1].removeprefix("GT") or 0), row[2], row[3]),
    )
    write_excel(
        output_path,
        [["Gene", "Genotype", "Subtype", "ReferenceAccession"]]
        + [list(row) for row in sorted_rows],
    )


def write_differences(
    gene: str, refs: dict[tuple[str, str], str], consensus_path: Path, output_path: Path
) -> int:
    consensuses = read_fasta(consensus_path)
    sequences: dict[str, dict[int, tuple[str, str]]] = {}
    for genotype in "12345678":
        reference = refs[(genotype, gene)]
        consensus = consensuses.get(f"GT{genotype}")
        if consensus is None:
            raise ValueError(f"{consensus_path} has no GT{genotype} record")
        pairs = covered_alignment_pairs(reference, consensus)
        sequences[genotype] = pairs
    output_path.parent.mkdir(parents=True, exist_ok=True)
    excel_rows: list[list[str | int]] = []
    for genotype in "12345678":
        pairs = sequences[genotype]
        positions = RAS_POSITIONS[gene]
        excel_rows.append([])
        excel_rows.append(["Gene", "Genotype", "Sequence"] + list(positions))
        values = [
            pairs[position][0] if position in pairs else "" for position in positions
        ]
        excel_rows.append([gene, f"GT{genotype}", "Reference"] + values)
        values = [
            pairs[position][1]
            if position in pairs and pairs[position][1] != pairs[position][0]
            else ""
            for position in positions
        ]
        excel_rows.append([gene, f"GT{genotype}", "Consensus"] + values)
    write_excel(output_path, excel_rows)
    return sum(len(pairs) for pairs in sequences.values())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument(
        "--gene",
        choices=GENES,
        action="append",
        help="Gene to process. Repeat to process several; default processes all genes.",
    )
    parser.add_argument(
        "--consensus-dir",
        type=Path,
        default=Path("outputs/comet"),
        help="Directory containing the current Comet consensus FASTAs (default: outputs/comet).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/comet-NS5B-position-282-include-or-short"),
        help="Directory for generated comparison workbooks (default: outputs/comet-NS5B-position-282).",
    )
    parser.add_argument(
        "--subtype-reference-dir",
        type=Path,
        help="Optional directory containing HCV_Subtype_Refs_<gene>_AA.fasta files.",
    )
    parser.add_argument(
        "--subtype-reference-fasta",
        type=Path,
        help="Explicit subtype reference FASTA for a single selected gene.",
    )
    parser.add_argument(
        "--subtype-consensus-fasta",
        type=Path,
        help="Explicit subtype consensus FASTA for a single selected gene.",
    )
    parser.add_argument(
        "--subtype-length-csv",
        type=Path,
        help="Write matched subtype reference/consensus names and amino-acid lengths to this CSV.",
    )
    args = parser.parse_args()

    refs = load_references(args.reference_fasta)
    genes = args.gene or GENES
    if (
        args.subtype_reference_fasta
        or args.subtype_consensus_fasta
        or args.subtype_length_csv
    ) and len(genes) != 1:
        raise ValueError(
            "Explicit subtype FASTA and length-CSV options require exactly one --gene value"
        )
    for gene in genes:
        consensus_gene = "NS5A" if gene == "NS5A_NTD" else gene
        output = (
            args.output_dir
            / f"HCV_GT_Ref_vs_Comet_GT_Consensus_Differences_{gene}.xlsx"
        )
        count = write_differences(
            gene,
            refs,
            args.consensus_dir / f"{consensus_gene}_GT_Consensus.fasta",
            output,
        )
        print(f"{output}: {count} aligned reference-consensus amino-acid pairs")
    if args.subtype_reference_dir or args.subtype_reference_fasta:
        missing_consensuses: list[tuple[str, str, str, str]] = []
        matched_lengths: list[tuple[str, str, int, int]] = []
        for gene in genes:
            consensus_gene = "NS5A" if gene == "NS5A_NTD" else gene
            reference_path = (
                args.subtype_reference_fasta
                or args.subtype_reference_dir / f"HCV_Subtype_Refs_{gene}_AA.fasta"
            )
            output = (
                args.output_dir
                / f"HCV_Subtype_Ref_vs_Comet_Subtype_Consensus_Aligned_{gene}.xlsx"
            )
            consensus_path = (
                args.subtype_consensus_fasta
                or args.consensus_dir / f"{consensus_gene}_Subtype_Consensus.fasta"
            )
            matched, missing, lengths = write_subtype_alignment(
                gene,
                reference_path,
                consensus_path,
                output,
            )
            missing_consensuses.extend(missing)
            matched_lengths.extend(lengths)
            print(
                f"{output}: {matched} matched subtype references; {len(missing)} without a Comet consensus"
            )
        suffix = "" if len(genes) == len(GENES) else f"_{genes[0]}"
        missing_output = (
            args.output_dir / f"HCV_Subtype_Refs_Without_Comet_Consensus{suffix}.xlsx"
        )
        write_missing_subtype_consensus_report(missing_output, missing_consensuses)
        print(
            f"{missing_output}: {len(missing_consensuses)} subtype references without a Comet consensus"
        )
        if args.subtype_length_csv:
            write_subtype_length_csv(args.subtype_length_csv, matched_lengths)
            print(
                f"{args.subtype_length_csv}: {len(matched_lengths)} matched subtype reference/consensus lengths"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
