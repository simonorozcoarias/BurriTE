#!/usr/bin/env python3
"""Command-line adapters used by the BurriTE v4 Nextflow workflow.

The biological algorithms remain in :mod:`burrite_core`, which is the last
validated BurriTE v3 Python implementation.  This module exposes small,
deterministic commands so Nextflow can schedule the independent work per
assembly and per assembly/chromosome.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd
from Bio import SeqIO

import burrite_core as core


def _mkdir(path: str | Path) -> Path:
    result = Path(path)
    result.mkdir(parents=True, exist_ok=True)
    return result


def _copy_files(files: Iterable[str], destination: Path) -> list[Path]:
    destination.mkdir(parents=True, exist_ok=True)
    copied = []
    for source in files:
        src = Path(source)
        dst = destination / src.name
        shutil.copy2(src, dst)
        copied.append(dst)
    return copied


def _escape_anchor_gff_id(value: object) -> str:
    """Escape delimiters before storing an original anchor ID in its map.

    TE names reconstructed by One Code To Find Them All can contain commas.
    Preserve valid percent escapes, escape literal percent signs, and encode
    commas without changing the embedded ``_ID=`` tokens used by the validated
    BurriTE v3 downstream parser.  Liftoff itself receives only short opaque
    IDs; these escaped originals are restored after mapping.
    """
    text = str(value)
    text = re.sub(r"%(?![0-9A-Fa-f]{2})", "%25", text)
    return text.replace(",", "%2C")


def cmd_normalize_samplesheet(args: argparse.Namespace) -> None:
    input_path = Path(args.input).resolve()
    base_dir = Path(args.base_dir).resolve()
    rows = []
    with input_path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if header is None or len(header) < 2:
            raise ValueError("The assemblies CSV must contain a header and at least two columns")
        for line_number, row in enumerate(reader, 2):
            if not row or all(not value.strip() for value in row):
                continue
            if len(row) < 2:
                raise ValueError(f"Invalid assemblies CSV row {line_number}: expected path,sample")
            raw_path, sample = row[0].strip(), row[1].strip()
            if not re.fullmatch(r"[A-Za-z0-9_-]+", sample) or "." in sample:
                raise ValueError(
                    f"Invalid sample name at row {line_number}: {sample!r}. "
                    "Use only letters, numbers, '_' and '-', without dots."
                )
            candidate = Path(raw_path).expanduser()
            if not candidate.is_absolute():
                from_launch = base_dir / candidate
                from_sheet = input_path.parent / candidate
                candidate = from_launch if from_launch.exists() else from_sheet
            candidate = candidate.resolve()
            if not candidate.is_file():
                raise FileNotFoundError(f"Assembly FASTA not found at row {line_number}: {candidate}")
            rows.append((str(candidate), sample))

    samples = [sample for _, sample in rows]
    if len(samples) != len(set(samples)):
        duplicates = sorted({sample for sample in samples if samples.count(sample) > 1})
        raise ValueError("Duplicated assembly IDs: " + ", ".join(duplicates))

    with Path(args.output).open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([header[0], header[1]])
        writer.writerows(rows)


def cmd_graffite_vcf_to_rm(args: argparse.Namespace) -> None:
    samples = [x for x in args.samples.split(",") if x]
    allowed = set(samples) | {args.reference_name}
    records: dict[str, list[str]] = {sample: [] for sample in allowed}

    with Path(args.vcf).open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            id_tokens = fields[2].split(".")
            if len(id_tokens) < 2:
                continue
            sv_type = id_tokens[-2]
            if sv_type == "DEL":
                sample = args.reference_name
            elif sv_type == "INS":
                sample = id_tokens[0]
            else:
                continue
            if sample not in allowed:
                continue

            contig = fields[0]
            start = int(fields[1])
            info = fields[7]
            length_match = re.search(r"(?:^|;)SVLEN=(-?\d+)", info)
            if length_match is None:
                continue
            length = abs(int(length_match.group(1)))
            end = start + length
            repeat_match = re.search(r"(?:^|;)repeat_ids=([^;]+)", info)
            class_match = re.search(r"(?:^|;)matching_classes=([^;]+)", info)
            strand_match = re.search(r"(?:^|;)STRANDS=([+-])", info)
            repeat_name = repeat_match.group(1) if repeat_match else "Unknown"
            repeat_class = class_match.group(1) if class_match else "Unknown"
            strand = strand_match.group(1) if strand_match else "+"
            records[sample].append(
                f"{sv_type}\t0\t0\t0\t{contig}\t{start}\t{end}\t{length}\t{strand}\t"
                f"{repeat_name}\t{repeat_class}\t0\t0\t0\t0\t0\tNo_ref_available\n"
            )

    output_dir = _mkdir(args.output_dir)
    for sample in [*samples, args.reference_name]:
        with (output_dir / f"{sample}_GraffiTE.out").open("w") as handle:
            handle.writelines(records.get(sample, []))


RM_COLUMNS = [
    "SW_score", "perc_div", "perc_del", "perc_ins", "query_seq", "query_start",
    "query_end", "length", "sense", "element", "family", "Pos_Repeat_Beg",
    "Pos_Repeat_End", "Pos_Repeat_Left", "ID", "Num_Assembled", "%_of_Ref",
]


def _rm_table_to_gff(input_path: str, output_path: str) -> None:
    path = Path(input_path)
    if path.stat().st_size == 0:
        Path(output_path).write_text("")
        return
    try:
        frame = pd.read_csv(path, sep="\t", header=None, comment="#")
    except pd.errors.EmptyDataError:
        Path(output_path).write_text("")
        return
    if frame.empty:
        Path(output_path).write_text("")
        return
    if frame.shape[1] < len(RM_COLUMNS):
        raise ValueError(f"{path} has {frame.shape[1]} columns; expected at least 17")
    frame = frame.iloc[:, : len(RM_COLUMNS)]
    frame.columns = RM_COLUMNS
    frame["start"] = frame[["query_start", "query_end"]].min(axis=1)
    frame["end"] = frame[["query_start", "query_end"]].max(axis=1)
    frame["sense"] = frame["sense"].replace({"C": "-"}).fillna("+")
    frame["attributes"] = (
        "ID=TE" + frame.index.astype(str)
        + ";Name=" + frame["element"].astype(str)
        + ";Class=" + frame["family"].astype(str)
    )
    gff = pd.DataFrame(
        {
            "seqid": frame["query_seq"],
            # Retain the v3 convention. GraffiTE annotations transferred with
            # Liftoff acquire source=Liftoff in the next stage.
            "source": "RepeatMasker",
            "type": "repeat_region",
            "start": frame["start"],
            "end": frame["end"],
            "score": ".",
            "strand": frame["sense"],
            "phase": ".",
            "attributes": frame["attributes"],
        }
    )
    gff.to_csv(output_path, sep="\t", header=False, index=False)


def cmd_annotation_gffs(args: argparse.Namespace) -> None:
    _rm_table_to_gff(args.repeatmasker, f"{args.sample}_RepeatMasker.gff")
    _rm_table_to_gff(args.graffite, f"{args.sample}_GraffiTE.gff")


def _read_gff(path: str) -> pd.DataFrame:
    columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    if Path(path).stat().st_size == 0:
        return pd.DataFrame(columns=columns)
    try:
        frame = pd.read_csv(path, sep="\t", comment="#", header=None, names=columns)
    except pd.errors.EmptyDataError:
        frame = pd.DataFrame(columns=columns)
    return frame


def cmd_join_annotations(args: argparse.Namespace) -> None:
    rm = _read_gff(args.repeatmasker_gff)
    graffite = _read_gff(args.graffite_gff)
    joined = pd.concat([rm, graffite], ignore_index=True)
    if not joined.empty:
        joined = joined.sort_values(by=["seqid", "start", "score"], kind="stable")
    joined.to_csv(f"{args.sample}_final_annot.gff", sep="\t", header=False, index=False)


def cmd_filter_bed(args: argparse.Namespace) -> None:
    frame = _read_gff(args.annotation)
    if frame.empty:
        Path(f"{args.sample}_final_annot.bed").write_text("")
        return
    lengths = frame["end"] - frame["start"]
    keep = lengths >= args.min_length
    if args.max_length is not None:
        keep &= lengths <= args.max_length
    frame = frame.loc[keep].reset_index(drop=True)
    ids = frame.apply(
        lambda row: (
            f"{row.name + 1}_{row['attributes']}_{row['start']}_{row['end']}"
            .replace("Name=", "")
            .replace("Class=", "-")
            .replace(";", "-")
        ),
        axis=1,
    )
    output = pd.DataFrame({0: frame["seqid"], 1: frame["start"], 2: frame["end"], 3: ids})
    output.to_csv(f"{args.sample}_final_annot.bed", sep="\t", header=False, index=False)


def cmd_make_anchors(args: argparse.Namespace) -> None:
    bed_path = Path(args.bed)
    output = Path(f"{args.sample}_anchors{args.flank_size}.gff")
    mapping_output = Path(f"{args.sample}_anchor_id_map.tsv")
    if bed_path.stat().st_size == 0:
        output.write_text("##gff-version 3\n")
        mapping_output.write_text("short_id\toriginal_id\n")
        return
    bed = pd.read_csv(bed_path, sep="\t", header=None, names=["chr", "start", "end", "ID"])
    bed = bed.reset_index(drop=True)
    bed["ID"] = bed["ID"].astype(str)
    bed["short_base"] = [f"BTE{index:08d}" for index in range(1, len(bed) + 1)]
    bed["F1_start"] = (bed["start"] - args.flank_size).clip(lower=1)
    bed["F2_end"] = bed["end"] + args.flank_size
    f1 = pd.DataFrame(
        {
            "chr": bed["chr"], "start": bed["F1_start"], "end": bed["start"],
            "original_id": (bed["ID"] + "_F1").map(_escape_anchor_gff_id),
            "short_id": bed["short_base"] + "_F1",
        }
    )
    f2 = pd.DataFrame(
        {
            "chr": bed["chr"], "start": bed["end"], "end": bed["F2_end"],
            "original_id": (bed["ID"] + "_F2").map(_escape_anchor_gff_id),
            "short_id": bed["short_base"] + "_F2",
        }
    )
    combined = pd.concat([f1, f2], ignore_index=True).sort_values("short_id")
    if combined["short_id"].str.len().max() > 254:
        raise ValueError("Internal BurriTE anchor IDs exceed the SAM QNAME limit")
    gff = pd.DataFrame(
        {
            "seqid": combined["chr"], "source": "Mmd", "type": "repeat_region",
            "start": combined["start"], "end": combined["end"], "score": ".",
            "strand": "+", "phase": ".", "attributes": "ID=" + combined["short_id"],
        }
    )
    with output.open("w") as handle:
        handle.write("##gff-version 3\n#!gff-spec-version 1.21\n#!processor BurriTE v4\n")
        gff.to_csv(handle, sep="\t", index=False, header=False)
    combined[["short_id", "original_id"]].to_csv(
        mapping_output, sep="\t", index=False
    )


def cmd_restore_anchor_ids(args: argparse.Namespace) -> None:
    """Restore descriptive IDs after Liftoff has finished using short SAM names."""
    with Path(args.mapping).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["short_id", "original_id"]:
            raise ValueError(f"Invalid anchor ID map header in {args.mapping}")
        mapping: dict[str, str] = {}
        for row_number, row in enumerate(reader, 2):
            short_id = row["short_id"]
            original_id = row["original_id"]
            if not re.fullmatch(r"BTE[0-9]{8}_F[12]", short_id):
                raise ValueError(
                    f"Invalid short anchor ID {short_id!r} at {args.mapping}:{row_number}"
                )
            if short_id in mapping:
                raise ValueError(f"Duplicated short anchor ID in {args.mapping}: {short_id}")
            mapping[short_id] = original_id

    short_with_copy = re.compile(r"^(BTE[0-9]{8}_F[12])(_[0-9]+)?$")
    restored_attributes = 0
    data_records = 0
    unresolved: set[str] = set()

    def restore_value(value: str) -> str:
        nonlocal restored_attributes
        match = short_with_copy.fullmatch(value)
        if not match:
            return value
        short_id, copy_suffix = match.groups()
        original_id = mapping.get(short_id)
        if original_id is None:
            unresolved.add(short_id)
            return value
        restored_attributes += 1
        return original_id + (copy_suffix or "")

    with Path(args.gff).open() as source, Path(args.output).open("w") as destination:
        for line_number, line in enumerate(source, 1):
            if not line.strip() or line.startswith("#"):
                destination.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Invalid Liftoff GFF record at {args.gff}:{line_number}: "
                    f"expected 9 columns, found {len(fields)}"
                )
            data_records += 1
            attributes = []
            for item in fields[8].split(";"):
                if "=" not in item:
                    attributes.append(item)
                    continue
                key, value = item.split("=", 1)
                if key in {"ID", "copy_num_ID"}:
                    value = restore_value(value)
                attributes.append(f"{key}={value}")
            fields[8] = ";".join(attributes)
            destination.write("\t".join(fields) + "\n")

    if unresolved:
        preview = ", ".join(sorted(unresolved)[:10])
        raise ValueError(f"Unresolved short anchor IDs in {args.gff}: {preview}")
    if data_records and restored_attributes == 0:
        raise ValueError(f"No short anchor IDs were restored in non-empty GFF {args.gff}")
    print(
        f"Restored {restored_attributes} anchor ID attributes across "
        f"{data_records} Liftoff records"
    )


def cmd_rename_by_paf(args: argparse.Namespace) -> None:
    scaffold_map: dict[str, str] = {}
    paf_path = Path(args.paf)
    if paf_path.stat().st_size:
        names = ["query", "qlen", "qstart", "qend", "strand", "target", "tlen", "tstart", "tend", "nmatch", "alen", "mapq"]
        frame = pd.read_csv(paf_path, sep="\t", header=None, usecols=range(12), names=names, comment="#")
        frame["aligned_len"] = frame["qend"] - frame["qstart"]
        frame = frame[frame["aligned_len"] > 0]
        coverage = frame.groupby(["target", "query"])["aligned_len"].sum().reset_index()
        scaffold_map = (
            coverage.sort_values("aligned_len", ascending=False)
            .drop_duplicates("target")
            .set_index("query")["target"]
            .to_dict()
        )

    renamed_records = []
    for record in SeqIO.parse(args.fasta, "fasta"):
        new_id = scaffold_map.get(record.id, record.id)
        record.id = new_id
        record.name = new_id
        record.description = ""
        renamed_records.append(record)
    SeqIO.write(renamed_records, f"{args.sample}_chrRenamed.fasta", "fasta")

    columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    if Path(args.anchors).stat().st_size:
        try:
            anchors = pd.read_csv(args.anchors, sep="\t", comment="#", header=None, names=columns)
        except pd.errors.EmptyDataError:
            anchors = pd.DataFrame(columns=columns)
        anchors["seqid"] = anchors["seqid"].map(scaffold_map).fillna(anchors["seqid"])
    else:
        anchors = pd.DataFrame(columns=columns)
    anchors.to_csv(f"{args.sample}_anchors_chrRenamed.gff", sep="\t", header=False, index=False)


def _extract_one_chromosome(fasta: str, chromosome: str, output: str) -> None:
    found = [record for record in SeqIO.parse(fasta, "fasta") if record.id == chromosome]
    if len(found) != 1:
        raise ValueError(f"Expected exactly one sequence named {chromosome!r} in {fasta}; found {len(found)}")
    SeqIO.write(found, output, "fasta")


def cmd_extract_reference_chrom(args: argparse.Namespace) -> None:
    _extract_one_chromosome(args.reference, args.chromosome, f"reference_{args.chromosome}.fasta")


def cmd_prepare_alt_chrom(args: argparse.Namespace) -> None:
    _extract_one_chromosome(args.fasta, args.chromosome, f"{args.sample}_{args.chromosome}.fasta")
    columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    if Path(args.anchors).stat().st_size:
        try:
            anchors = pd.read_csv(args.anchors, sep="\t", comment="#", header=None, names=columns)
        except pd.errors.EmptyDataError:
            anchors = pd.DataFrame(columns=columns)
    else:
        anchors = pd.DataFrame(columns=columns)
    anchors = anchors[anchors["seqid"] == args.chromosome]
    anchors.to_csv(f"{args.sample}_{args.chromosome}_anchors.gff", sep="\t", header=False, index=False)
    Path(f"{args.chromosome}_repeat_region.features").write_text("repeat_region\n")


def cmd_sample_bed(args: argparse.Namespace) -> None:
    work = _mkdir(f"transfer_{args.sample}")
    combined = work / f"{args.sample}_combined_liffout.gff_polished"
    with combined.open("w") as output:
        for path in sorted(args.gff):
            with Path(path).open() as source:
                shutil.copyfileobj(source, output)

    with combined.open() as handle:
        has_records = any(line.strip() and not line.startswith("#") for line in handle)
    species_tag = args.sample.replace("_", "")
    if has_records:
        produced = Path(core.te_transfer_per_species(str(combined), species_tag, str(work)))
        shutil.copy2(produced, f"{args.sample}.bed")
    else:
        (work / f"{species_tag}_F1F2_coordinates.txt").write_text("")
        (work / f"{species_tag}_TEinsertions.txt").write_text("")
        (work / f"Table_{species_tag}_500filtered.txt").write_text(
            "Chrm\tF1_StartR\tF1_EndR\tF2_StartR\tF2_EndR\tTE_Family\tLengthT\tLengthR\tF1_TE_F2\tState\n"
        )
        Path(f"{args.sample}.bed").write_text("")


def cmd_reference_bed(args: argparse.Namespace) -> None:
    slots = [None] * 12
    slots[0] = args.reference
    slots[7] = args.annotation
    file_dict = {args.sample: slots}
    out_dir = _mkdir("reference_bed")
    produced = core.build_reference_bed(
        file_dict, args.sample, args.min_length, args.max_length, str(out_dir)
    )
    shutil.copy2(produced, f"{args.sample}.bed")


def _expected_samples(raw: str) -> list[str]:
    samples = [sample.strip() for sample in raw.split(",") if sample.strip()]
    if not samples:
        raise ValueError("The expected sample list is empty")
    if len(samples) != len(set(samples)):
        raise ValueError("The expected sample list contains duplicates")
    return samples


def _data_line_count(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def cmd_validate_sample_beds(args: argparse.Namespace) -> None:
    """Require exactly one BED for every genome before locus collapsing.

    Nextflow channel mistakes can otherwise remove a genome without causing a
    process failure.  This explicit barrier also rejects header-only/empty BEDs
    by default, while ``--allow-empty`` is available for intentional edge cases.
    """
    expected = _expected_samples(args.expected)
    expected_set = set(expected)
    by_sample: dict[str, list[Path]] = defaultdict(list)
    for raw_path in args.bed:
        path = Path(raw_path)
        by_sample[path.stem].append(path)

    errors: list[str] = []
    unexpected = sorted(set(by_sample) - expected_set)
    if unexpected:
        errors.append("unexpected BED sample IDs: " + ", ".join(unexpected))

    output_dir = _mkdir(args.output_dir)
    report_rows: list[dict[str, object]] = []
    for sample in expected:
        matches = by_sample.get(sample, [])
        if not matches:
            report_rows.append(
                {"sample": sample, "bed": "", "rows": 0, "status": "MISSING"}
            )
            errors.append(f"missing BED for {sample}")
            continue
        if len(matches) > 1:
            report_rows.append(
                {
                    "sample": sample,
                    "bed": ",".join(str(path) for path in matches),
                    "rows": "",
                    "status": "DUPLICATE",
                }
            )
            errors.append(f"multiple BED files for {sample}")
            continue

        path = matches[0]
        rows = _data_line_count(path)
        if rows == 0 and not args.allow_empty:
            status = "EMPTY"
            errors.append(f"BED for {sample} contains no TE loci")
        elif rows == 0:
            status = "EMPTY_ALLOWED"
        else:
            status = "OK"
        shutil.copy2(path, output_dir / f"{sample}.bed")
        report_rows.append(
            {"sample": sample, "bed": path.name, "rows": rows, "status": status}
        )

    with Path(args.report).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sample", "bed", "rows", "status"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(report_rows)

    if errors:
        raise ValueError(
            "Per-genome BED validation failed: " + "; ".join(errors)
            + f". See {args.report}."
        )


def cmd_collapse(args: argparse.Namespace) -> None:
    if not args.bed:
        raise ValueError("At least one per-sample BED is required")
    pairs = sorted(((Path(path).stem, Path(path)) for path in args.bed), key=lambda x: x[0])
    samples = [sample for sample, _ in pairs]
    beds = [str(path) for _, path in pairs]
    out_dir = _mkdir("07_merge")

    if args.method == "incremental":
        collapsed = core.incremental_collapse(
            beds, samples, str(out_dir), "annots_incremental",
            args.lenthr, args.minlen, args.bedmin, merge_by=args.merge_by,
        )
    else:
        bed_dir = _mkdir(out_dir / "bed_files")
        _copy_files(beds, bed_dir)
        merged = out_dir / "annots_merged.bed"
        grouped = out_dir / "annots_groupby.bed"
        collapsed = out_dir / "annots_groupby_collapsed.bed"
        core.merge_beds_to_manaus_nhv(str(bed_dir), str(merged))
        core.run_bedtools_groupby(str(merged), str(grouped))
        core.collapse_and_filter_grouped_bed(str(grouped), str(collapsed))
    shutil.copy2(collapsed, out_dir / "collapsed.bed")


def cmd_consensus_flanks(args: argparse.Namespace) -> None:
    chromosomes = [line.strip() for line in Path(args.chromosome_list).read_text().splitlines() if line.strip()]
    core.generate_flanking_gffs(args.collapsed_bed, ".", args.flank_size, chromosomes)


def cmd_evaluate_sample(args: argparse.Namespace) -> None:
    combined = Path(f"{args.sample}_liffout_polished.gff")
    with combined.open("w") as output:
        for path in sorted(args.gff):
            with Path(path).open() as source:
                shutil.copyfileobj(source, output)
    core.evaluate_te_transfer(args.collapsed_bed, str(combined), f"{args.sample}_Transfers_TE.txt")


def cmd_evaluate_reference(args: argparse.Namespace) -> None:
    core.evaluate_reference_transfer(
        args.collapsed_bed, f"{args.sample}_Transfers_TE.txt", args.flank_size
    )


def cmd_validate_transfers(args: argparse.Namespace) -> None:
    """Check that every genome has one complete locus-genotype table."""
    expected = _expected_samples(args.expected)
    expected_set = set(expected)
    suffix = "_Transfers_TE"
    by_sample: dict[str, list[Path]] = defaultdict(list)
    for raw_path in args.transfer:
        path = Path(raw_path)
        stem = path.stem
        sample = stem[: -len(suffix)] if stem.endswith(suffix) else stem
        by_sample[sample].append(path)

    expected_loci = _data_line_count(Path(args.collapsed_bed))
    errors: list[str] = []
    unexpected = sorted(set(by_sample) - expected_set)
    if unexpected:
        errors.append("unexpected transfer sample IDs: " + ", ".join(unexpected))

    output_dir = _mkdir(args.output_dir)
    report_fields = [
        "sample", "transfer", "rows", "expected_loci", "Present", "Absent",
        "Ambiguous", "NA", "Other", "status",
    ]
    report_rows: list[dict[str, object]] = []
    for sample in expected:
        matches = by_sample.get(sample, [])
        if not matches:
            report_rows.append(
                {
                    "sample": sample, "transfer": "", "rows": 0,
                    "expected_loci": expected_loci, "Present": 0, "Absent": 0,
                    "Ambiguous": 0, "NA": 0, "Other": 0, "status": "MISSING",
                }
            )
            errors.append(f"missing transfer table for {sample}")
            continue
        if len(matches) > 1:
            report_rows.append(
                {
                    "sample": sample,
                    "transfer": ",".join(str(path) for path in matches),
                    "rows": "", "expected_loci": expected_loci, "Present": "",
                    "Absent": "", "Ambiguous": "", "NA": "", "Other": "",
                    "status": "DUPLICATE",
                }
            )
            errors.append(f"multiple transfer tables for {sample}")
            continue

        path = matches[0]
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or "State" not in reader.fieldnames:
                states: Counter[str] = Counter()
                rows = 0
                status = "INVALID_HEADER"
                errors.append(f"transfer table for {sample} has no State column")
            else:
                states = Counter((row.get("State") or "").strip() for row in reader)
                rows = sum(states.values())
                if rows != expected_loci:
                    status = "ROW_MISMATCH"
                    errors.append(
                        f"transfer table for {sample} has {rows} rows; expected {expected_loci}"
                    )
                else:
                    status = "OK"

        known = sum(states.get(state, 0) for state in ("Present", "Absent", "Ambiguous", "NA"))
        shutil.copy2(path, output_dir / f"{sample}_Transfers_TE.txt")
        report_rows.append(
            {
                "sample": sample, "transfer": path.name, "rows": rows,
                "expected_loci": expected_loci, "Present": states.get("Present", 0),
                "Absent": states.get("Absent", 0),
                "Ambiguous": states.get("Ambiguous", 0), "NA": states.get("NA", 0),
                "Other": rows - known, "status": status,
            }
        )

    with Path(args.report).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=report_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(report_rows)

    if errors:
        raise ValueError(
            "Per-genome transfer validation failed: " + "; ".join(errors)
            + f". See {args.report}."
        )


def cmd_matrix(args: argparse.Namespace) -> None:
    transfer_dir = _mkdir("transfers")
    _copy_files(args.transfer, transfer_dir)
    samples = [sample for sample in args.samples.split(",") if sample]
    contrib = core.build_contributor_map(args.method, args.merge_dir, args.collapsed_bed)
    loc_meta = core.build_locus_meta(args.merge_dir) if args.method == "incremental" else None
    core.build_presence_absence_matrix(
        args.collapsed_bed, samples, str(transfer_dir), "Presence_Absence_matrix.txt",
        contrib=contrib, loc_meta=loc_meta,
    )


def cmd_final_gff(args: argparse.Namespace) -> None:
    transfer_dir = _mkdir(f"transfers_{args.sample}")
    _copy_files(args.transfer, transfer_dir)
    slots = [None] * 12
    slots[0] = args.original_fasta
    slots[7] = args.annotation
    slots[8] = args.filtered_bed
    slots[11] = args.renamed_fasta
    file_dict = {args.sample: slots}
    output_dir = _mkdir(f"final_{args.sample}")

    if args.mode == "all":
        core.write_per_genome_all_gffs(
            args.matrix, str(transfer_dir), str(output_dir), file_dict,
            args.min_length, args.max_length, args.merge_dir,
            overlap_frac=args.overlap,
        )
        source = output_dir / f"{args.sample}_all_TEs.gff"
    else:
        core.write_per_genome_gffs(
            args.matrix, str(transfer_dir), str(output_dir), file_dict, args.merge_dir
        )
        source = output_dir / f"{args.sample}_TEs.gff"
    if not source.is_file():
        raise FileNotFoundError(f"The final GFF was not generated for {args.sample}: {source}")
    shutil.copy2(source, source.name)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="burrite_stage.py")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("normalize-samplesheet")
    p.add_argument("--input", required=True); p.add_argument("--base-dir", required=True)
    p.add_argument("--output", required=True); p.set_defaults(func=cmd_normalize_samplesheet)

    p = sub.add_parser("graffite-vcf-to-rm")
    p.add_argument("--vcf", required=True); p.add_argument("--samples", required=True)
    p.add_argument("--reference-name", required=True); p.add_argument("--output-dir", required=True)
    p.set_defaults(func=cmd_graffite_vcf_to_rm)

    p = sub.add_parser("annotation-gffs")
    p.add_argument("--sample", required=True); p.add_argument("--repeatmasker", required=True)
    p.add_argument("--graffite", required=True); p.set_defaults(func=cmd_annotation_gffs)

    p = sub.add_parser("join-annotations")
    p.add_argument("--sample", required=True); p.add_argument("--repeatmasker-gff", required=True)
    p.add_argument("--graffite-gff", required=True); p.set_defaults(func=cmd_join_annotations)

    p = sub.add_parser("filter-bed")
    p.add_argument("--sample", required=True); p.add_argument("--annotation", required=True)
    p.add_argument("--min-length", type=int, default=0); p.add_argument("--max-length", type=int)
    p.set_defaults(func=cmd_filter_bed)

    p = sub.add_parser("make-anchors")
    p.add_argument("--sample", required=True); p.add_argument("--bed", required=True)
    p.add_argument("--flank-size", type=int, required=True)
    p.add_argument("--gff-id-encoding", choices=["short-map-v2"], default="short-map-v2")
    p.set_defaults(func=cmd_make_anchors)

    p = sub.add_parser("restore-anchor-ids")
    p.add_argument("--gff", required=True); p.add_argument("--mapping", required=True)
    p.add_argument("--output", required=True); p.set_defaults(func=cmd_restore_anchor_ids)

    p = sub.add_parser("rename-by-paf")
    p.add_argument("--sample", required=True); p.add_argument("--fasta", required=True)
    p.add_argument("--anchors", required=True); p.add_argument("--paf", required=True)
    p.set_defaults(func=cmd_rename_by_paf)

    p = sub.add_parser("extract-reference-chrom")
    p.add_argument("--reference", required=True); p.add_argument("--chromosome", required=True)
    p.set_defaults(func=cmd_extract_reference_chrom)

    p = sub.add_parser("prepare-alt-chrom")
    p.add_argument("--sample", required=True); p.add_argument("--chromosome", required=True)
    p.add_argument("--fasta", required=True); p.add_argument("--anchors", required=True)
    p.set_defaults(func=cmd_prepare_alt_chrom)

    p = sub.add_parser("sample-bed")
    p.add_argument("--sample", required=True); p.add_argument("--gff", action="append", required=True)
    p.set_defaults(func=cmd_sample_bed)

    p = sub.add_parser("reference-bed")
    p.add_argument("--sample", required=True); p.add_argument("--reference", required=True)
    p.add_argument("--annotation", required=True); p.add_argument("--min-length", type=int, default=0)
    p.add_argument("--max-length", type=int); p.set_defaults(func=cmd_reference_bed)

    p = sub.add_parser("validate-sample-beds")
    p.add_argument("--bed", action="append", required=True); p.add_argument("--expected", required=True)
    p.add_argument("--output-dir", required=True); p.add_argument("--report", required=True)
    p.add_argument("--allow-empty", action="store_true"); p.set_defaults(func=cmd_validate_sample_beds)

    p = sub.add_parser("collapse")
    p.add_argument("--bed", action="append", required=True)
    p.add_argument("--method", choices=["incremental", "groupby"], required=True)
    p.add_argument("--merge-by", choices=["family", "superfamily"], default="superfamily")
    p.add_argument("--lenthr", type=int, default=100); p.add_argument("--minlen", type=int, default=0)
    p.add_argument("--bedmin", type=int, default=100); p.set_defaults(func=cmd_collapse)

    p = sub.add_parser("consensus-flanks")
    p.add_argument("--collapsed-bed", required=True); p.add_argument("--chromosome-list", required=True)
    p.add_argument("--flank-size", type=int, required=True); p.set_defaults(func=cmd_consensus_flanks)

    p = sub.add_parser("evaluate-sample")
    p.add_argument("--sample", required=True); p.add_argument("--collapsed-bed", required=True)
    p.add_argument("--gff", action="append", required=True); p.set_defaults(func=cmd_evaluate_sample)

    p = sub.add_parser("evaluate-reference")
    p.add_argument("--sample", required=True); p.add_argument("--collapsed-bed", required=True)
    p.add_argument("--flank-size", type=int, required=True); p.set_defaults(func=cmd_evaluate_reference)

    p = sub.add_parser("validate-transfers")
    p.add_argument("--transfer", action="append", required=True); p.add_argument("--expected", required=True)
    p.add_argument("--collapsed-bed", required=True); p.add_argument("--output-dir", required=True)
    p.add_argument("--report", required=True); p.set_defaults(func=cmd_validate_transfers)

    p = sub.add_parser("matrix")
    p.add_argument("--collapsed-bed", required=True); p.add_argument("--merge-dir", required=True)
    p.add_argument("--transfer", action="append", required=True); p.add_argument("--samples", required=True)
    p.add_argument("--method", choices=["incremental", "groupby"], required=True)
    p.set_defaults(func=cmd_matrix)

    p = sub.add_parser("final-gff")
    p.add_argument("--sample", required=True); p.add_argument("--matrix", required=True)
    p.add_argument("--merge-dir", required=True); p.add_argument("--transfer", action="append", required=True)
    p.add_argument("--original-fasta", required=True); p.add_argument("--renamed-fasta", required=True)
    p.add_argument("--annotation", required=True); p.add_argument("--filtered-bed", required=True)
    p.add_argument("--mode", choices=["burrite", "all"], default="burrite")
    p.add_argument("--min-length", type=int, default=0); p.add_argument("--max-length", type=int)
    p.add_argument("--overlap", type=float, default=0.5); p.set_defaults(func=cmd_final_gff)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
