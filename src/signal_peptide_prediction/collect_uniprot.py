"""Collect and audit the preliminary UniProtKB datasets."""

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
from textwrap import fill

if __package__:
    from . import selection, uniprot_api
else:
    import selection
    import uniprot_api


# UniProt taxonomy ID 2759 is Eukaryota.
QUERY = "reviewed:true AND taxonomy_id:2759"


def local_entries(directory):
    paths = sorted(directory.glob("*.json"))
    if not paths:
        raise ValueError(f"No UniProt entry JSON files in {directory}")
    for path in paths:
        yield json.loads(path.read_text(encoding="utf-8"))


def write_fasta(file, accession, sequence):
    file.write(f">{accession}\n{fill(sequence, width=60)}\n")


def collect(entries, output):
    counts = Counter()
    accessions = set()
    columns = ("accession", "organism", "kingdom", "protein_length")
    with (
        (output / "positives.tsv").open("w", encoding="utf-8", newline="") as positive_tsv,
        (output / "negatives.tsv").open("w", encoding="utf-8", newline="") as negative_tsv,
        (output / "positives.fasta").open("w", encoding="utf-8") as positive_fasta,
        (output / "negatives.fasta").open("w", encoding="utf-8") as negative_fasta,
    ):
        positives = csv.writer(positive_tsv, delimiter="\t", lineterminator="\n")
        negatives = csv.writer(negative_tsv, delimiter="\t", lineterminator="\n")
        positives.writerow((*columns, "cleavage_after"))
        negatives.writerow((*columns, "early_tm_90"))

        for entry in entries:
            accession = entry["primaryAccession"]
            if accession in accessions:
                raise ValueError(f"Duplicate accession: {accession}")
            accessions.add(accession)

            row = selection.shared_fields(entry, counts)
            if row is None:
                continue
            _, _, _, protein_length = row
            sequence = entry["sequence"]["value"]
            if len(sequence) != protein_length:
                raise ValueError(f"Sequence length mismatch for {accession}")

            features = entry.get("features", [])
            signals = [feature for feature in features if feature["type"] == "Signal"]
            # A predicted Signal annotation also rules out a negative label.
            if signals:
                cut = selection.positive_cleavage(signals, protein_length, counts)
                if cut is not None:
                    positives.writerow((*row, cut))
                    write_fasta(positive_fasta, accession, sequence)
            elif selection.keep_negative(entry, counts):
                early_tm = selection.has_early_tm(features)
                counts["negatives_with_early_tm"] += early_tm
                negatives.writerow((*row, str(early_tm).lower()))
                write_fasta(negative_fasta, accession, sequence)

    return counts


SHARED_STEPS = (
    "records_read", "reviewed_eukaryote", "non_fragment",
    "evidence_at_protein_level", "protein_length_at_least_40",
)
POSITIVE_STEPS = (
    "has_signal_annotation", "has_experimental_signal",
    "has_exact_cleavage_site", "signal_length_at_least_14",
    "one_distinct_eligible_cleavage_site",
)
NEGATIVE_STEPS = (
    "no_signal_annotation", "has_allowed_location",
    "has_experimental_evidence_for_allowed_location",
)


def filter_report(counts):
    yield "stage", "filter", "input", "kept", "dropped"
    # Both class branches start with the proteins that passed all shared rules.
    shared_total = counts[SHARED_STEPS[-1]]
    for stage, steps, previous in (
        ("shared", SHARED_STEPS, counts["records_read"]),
        ("positive", POSITIVE_STEPS, shared_total),
        ("negative", NEGATIVE_STEPS, shared_total),
    ):
        for step in steps:
            kept = counts[step]
            yield stage, step, previous, kept, previous - kept
            previous = kept


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path, help="new output directory")
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--max-pages", type=int, help=f"preview the first N pages of {uniprot_api.PAGE_SIZE} entries"
    )
    source.add_argument("--input-json-dir", type=Path, help="read individual entry JSON files")
    args = parser.parse_args()
    if args.max_pages is not None and args.max_pages < 1:
        parser.error("--max-pages must be at least 1")

    destination = args.output_dir
    if destination.exists():
        parser.error(f"Output directory already exists: {destination}")

    if args.input_json_dir:
        metadata = {"source": str(args.input_json_dir), "complete": False}
        entries = local_entries(args.input_json_dir)
    else:
        metadata = {}
        entries = uniprot_api.uniprot_entries(QUERY, args.max_pages, metadata)

    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(prefix=".uniprot-", dir=destination.parent) as staging:
            output = Path(staging)
            counts = collect(entries, output)
            rows = list(filter_report(counts))
            with (output / "filter_counts.tsv").open("w", encoding="utf-8", newline="") as file:
                csv.writer(file, delimiter="\t", lineterminator="\n").writerows(rows)
            metadata.update(
                records_read=counts["records_read"],
                positives=counts[POSITIVE_STEPS[-1]],
                negatives=counts[NEGATIVE_STEPS[-1]],
                negatives_with_early_tm=counts["negatives_with_early_tm"],
            )
            (output / "run.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
            output.rename(destination)
    except (OSError, ValueError) as error:
        parser.exit(1, f"Collection failed: {error}\n")

    csv.writer(sys.stdout, delimiter="\t", lineterminator="\n").writerows(rows)
    print(f"Negatives with an early transmembrane helix: {metadata['negatives_with_early_tm']:,}")
    print(f"Saved {'full collection' if metadata['complete'] else 'preview'} to {destination}")


if __name__ == "__main__":
    main()
