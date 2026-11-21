"""Decide which UniProt entries qualify for the positive and negative datasets."""

REVIEWED = "UniProtKB reviewed (Swiss-Prot)"
EXPERIMENTAL_EVIDENCE = "ECO:0000269"
MIN_PROTEIN_LENGTH = 40
MIN_SIGNAL_LENGTH = 14
EARLY_TM_END = 90
NEGATIVE_COMPARTMENTS = (
    "cytosol", "cytoplasm", "nucleus", "mitochondrion",
    "plastid", "chloroplast", "peroxisome", "cell membrane",
)
KINGDOMS = {"Metazoa": "Metazoa", "Fungi": "Fungi", "Viridiplantae": "Plants"}

def shared_fields(entry, counts):
    """Return TSV fields only for records that pass the shared criteria."""
    counts["records_read"] += 1
    organism = entry["organism"]
    lineage = organism["lineage"]
    if entry.get("entryType") != REVIEWED or "Eukaryota" not in lineage:
        return None
    counts["reviewed_eukaryote"] += 1

    flag = entry.get("proteinDescription", {}).get("flag") or ""
    if "fragment" in flag.casefold():
        return None
    counts["non_fragment"] += 1

    if entry.get("proteinExistence") != "1: Evidence at protein level":
        return None
    counts["evidence_at_protein_level"] += 1

    length = entry["sequence"]["length"]
    if length < MIN_PROTEIN_LENGTH:
        return None
    counts["protein_length_at_least_40"] += 1

    kingdom = next((KINGDOMS[taxon] for taxon in lineage if taxon in KINGDOMS), "Other")
    return (entry["primaryAccession"], organism["scientificName"], kingdom, length)


def has_experimental_evidence(annotation):
    return any(
        evidence["evidenceCode"] == EXPERIMENTAL_EVIDENCE
        for evidence in annotation.get("evidences", [])
    )


def known_cleavage(signal, protein_length):
    start = signal["location"].get("start") or {}
    end = signal["location"].get("end") or {}
    # UniProt gives the last signal peptide residue. The cut is after it.
    position = end.get("value")
    if (
        start.get("value") == 1
        and start.get("modifier") == end.get("modifier") == "EXACT"
        and isinstance(position, int)
        and 1 <= position < protein_length
        and "not cleaved" not in (signal.get("description") or "").casefold()
    ):
        return position
    return None


def positive_cleavage(signals, protein_length, counts):
    counts["has_signal_annotation"] += 1
    experimental_signals = [signal for signal in signals if has_experimental_evidence(signal)]
    if not experimental_signals:
        return None
    counts["has_experimental_signal"] += 1

    cuts = {known_cleavage(signal, protein_length) for signal in experimental_signals}
    cuts.discard(None)
    if not cuts:
        return None
    counts["has_exact_cleavage_site"] += 1

    eligible_cuts = {cut for cut in cuts if cut >= MIN_SIGNAL_LENGTH}
    if not eligible_cuts:
        return None
    counts["signal_length_at_least_14"] += 1

    if len(eligible_cuts) != 1:
        return None
    counts["one_distinct_eligible_cleavage_site"] += 1
    return eligible_cuts.pop()


def is_allowed_compartment(name):
    name = name.casefold()
    return any(
        name == compartment or name.startswith((f"{compartment},", f"{compartment} "))
        for compartment in NEGATIVE_COMPARTMENTS
    )


def keep_negative(entry, counts):
    counts["no_signal_annotation"] += 1
    locations = []
    for comment in entry.get("comments", []):
        # A named isoform may differ from the reference sequence we export.
        if comment["commentType"] != "SUBCELLULAR LOCATION" or "molecule" in comment:
            continue
        for item in comment.get("subcellularLocations", []):
            if "location" in item and is_allowed_compartment(item["location"]["value"]):
                locations.append(item["location"])

    if not locations:
        return False
    counts["has_allowed_location"] += 1

    # The experimental evidence must support one of the allowed locations.
    if not any(has_experimental_evidence(location) for location in locations):
        return False
    counts["has_experimental_evidence_for_allowed_location"] += 1
    return True


def has_early_tm(features):
    for feature in features:
        if feature["type"] != "Transmembrane":
            continue
        start = feature["location"].get("start") or {}
        position = start.get("value")
        if (
            start.get("modifier") == "EXACT"
            and isinstance(position, int)
            and 1 <= position <= EARLY_TM_END
        ):
            return True
    return False
