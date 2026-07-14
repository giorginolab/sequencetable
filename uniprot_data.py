import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from urllib.error import URLError
from urllib.parse import quote
from urllib.request import urlopen

import pandas as pd

UNIPROT_NAMESPACE = {"up": "http://uniprot.org/uniprot"}
UNIPROT_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
FIXTURE_DIRECTORY = Path(__file__).resolve().parent / "test"

DATAFRAME_COLUMNS = [
    "Residue Number",
    "Residue code",
    "Secondary structure",
    "Domain",
    "Pfam entries (protein-level)",
    "Disorder",
    "Disulfide bridges",
    "Glycosylation sites",
    "Phosphorylation sites",
    "active sites",
    "Binding sites",
    "modified",
]


def _exact_position(element):
    """Return a positive position, or None for unknown/invalid locations."""
    if element is None:
        return None

    value = element.get("position")
    if value is None:
        return None

    try:
        position = int(value)
    except (TypeError, ValueError):
        return None
    return position if position > 0 else None


def _feature_description(feature):
    description = feature.get("description", "").strip()
    if description:
        return description

    ligand = feature.find("up:ligand", UNIPROT_NAMESPACE)
    if ligand is None:
        return ""

    name = ligand.findtext("up:name", default="", namespaces=UNIPROT_NAMESPACE).strip()
    note = ligand.findtext("up:note", default="", namespaces=UNIPROT_NAMESPACE).strip()
    if name and note:
        return f"{name} ({note})"
    return name or note


def _normalized_feature_type(feature_type, description):
    feature_type = feature_type.lower()
    description = description.lower()

    if feature_type == "modified residue" and "phospho" in description:
        return "phosphorylation site"
    if feature_type in {"region", "region of interest"} and "disorder" in description:
        return "disorder"
    if feature_type == "region" and "pfam" in description:
        return "pfam"
    return feature_type


def _parse_uniprot_xml(response):
    try:
        root = ET.fromstring(response)
    except (ET.ParseError, TypeError, ValueError):
        return None, {}, "UniProt returned data that could not be parsed."

    sequence_elem = root.find("./up:entry/up:sequence", UNIPROT_NAMESPACE)
    if sequence_elem is None or not sequence_elem.text:
        return None, {}, "Could not find a protein sequence in the UniProt response."

    protein_sequence = "".join(sequence_elem.text.split())
    if not protein_sequence:
        return None, {}, "Could not find a protein sequence in the UniProt response."

    annotations = defaultdict(list)
    for feature in root.findall(".//up:feature", UNIPROT_NAMESPACE):
        feature_type = feature.get("type")
        if not feature_type:
            continue

        location = feature.find("up:location", UNIPROT_NAMESPACE)
        if location is None or location.get("sequence") is not None:
            # Locations on non-canonical sequences do not index this table.
            continue

        position = _exact_position(location.find("up:position", UNIPROT_NAMESPACE))
        begin = _exact_position(location.find("up:begin", UNIPROT_NAMESPACE))
        end = _exact_position(location.find("up:end", UNIPROT_NAMESPACE))

        if position is not None:
            coordinates = {"position": position}
        elif begin is not None and end is not None:
            coordinates = {"begin": begin, "end": end}
        else:
            # UniProt permits uncertain and unknown locations without a position.
            continue

        description = _feature_description(feature)
        normalized_type = _normalized_feature_type(feature_type, description)
        annotations[normalized_type].append(
            {**coordinates, "description": description}
        )

    # UniProt XML exposes Pfam as protein-level cross-references, without residue
    # coordinates. Label that scope explicitly instead of presenting blank or
    # misleading residue-level domain ranges.
    for reference in root.findall(
        "./up:entry/up:dbReference[@type='Pfam']", UNIPROT_NAMESPACE
    ):
        reference_id = reference.get("id", "").strip()
        entry_name = ""
        for prop in reference.findall("up:property", UNIPROT_NAMESPACE):
            if prop.get("type") == "entry name":
                entry_name = prop.get("value", "").strip()
                break

        description = reference_id
        if reference_id and entry_name:
            description = f"{reference_id} ({entry_name})"
        elif entry_name:
            description = entry_name

        if description:
            annotations["pfam"].append(
                {
                    "begin": 1,
                    "end": len(protein_sequence),
                    "description": description,
                }
            )

    return protein_sequence, dict(annotations), None


def get_uniprot_data(uniprot_id):
    """Fetch a UniProt sequence and its annotations.

    Returns a consistent ``(sequence, annotations, error_message)`` tuple.
    ``error_message`` is None on success.
    """
    if not isinstance(uniprot_id, str):
        return None, {}, "Please enter a valid UniProt ID."

    uniprot_id = uniprot_id.strip()
    if not uniprot_id or not UNIPROT_ID_PATTERN.fullmatch(uniprot_id):
        return None, {}, "Please enter a valid UniProt ID."

    local_file_path = FIXTURE_DIRECTORY / f"{uniprot_id}.xml"
    try:
        if local_file_path.is_file():
            response = local_file_path.read_bytes()
        else:
            encoded_id = quote(uniprot_id, safe="")
            url = f"https://rest.uniprot.org/uniprotkb/{encoded_id}?format=xml"
            with urlopen(url, timeout=30) as remote_response:
                response = remote_response.read()
    except (OSError, TimeoutError, URLError):
        return None, {}, f"Could not retrieve UniProt entry {uniprot_id}."

    return _parse_uniprot_xml(response)


def _annotation_range(item, sequence_length):
    start = item.get("begin", item.get("position"))
    end = item.get("end", item.get("position"))
    try:
        start = int(start)
        end = int(end)
    except (TypeError, ValueError):
        return None

    if start < 1 or end < start or end > sequence_length:
        return None
    return start - 1, end


def _append_value(df, row, column, value):
    if not value:
        return

    current = df.at[row, column]
    if not current:
        df.at[row, column] = value
    elif current != value:
        df.at[row, column] = f"{current}; {value}"


def create_dataframe(protein_sequence, annotations):
    """Create a residue-level dataframe from a sequence and annotations."""
    rows = []
    for residue_number, residue in enumerate(protein_sequence, start=1):
        rows.append(
            {
                "Residue Number": residue_number,
                "Residue code": residue,
                **{column: "" for column in DATAFRAME_COLUMNS[2:]},
            }
        )

    df = pd.DataFrame(rows, columns=DATAFRAME_COLUMNS)
    sequence_length = len(protein_sequence)

    feature_mapping = {
        "strand": "Secondary structure",
        "helix": "Secondary structure",
        "turn": "Secondary structure",
        "domain": "Domain",
        "pfam": "Pfam entries (protein-level)",
        "disorder": "Disorder",
        "glycosylation site": "Glycosylation sites",
        "phosphorylation site": "Phosphorylation sites",
        "modified residue": "modified",
        "active site": "active sites",
        "binding site": "Binding sites",
    }

    for feature_type, values in annotations.items():
        feature_type = feature_type.lower()

        if feature_type == "disulfide bond":
            for item in values:
                annotation_range = _annotation_range(item, sequence_length)
                if annotation_range is None:
                    continue
                start, end = annotation_range
                _append_value(df, start, "Disulfide bridges", f"Cys-{end + 1}")
                _append_value(df, end - 1, "Disulfide bridges", f"Cys-{start + 1}")
            continue

        column = feature_mapping.get(feature_type)
        if column is None:
            continue

        for item in values:
            annotation_range = _annotation_range(item, sequence_length)
            if annotation_range is None:
                continue

            start, end = annotation_range
            value = item.get("description", "")
            if column == "Secondary structure":
                value = feature_type.upper()

            for row in range(start, end):
                _append_value(df, row, column, value)

    return df
