from pathlib import Path
from urllib.error import HTTPError

import pandas as pd

import app
from uniprot_data import _parse_uniprot_xml, create_dataframe, get_uniprot_data


def test_bundled_entry_has_consistent_success_result():
    sequence, annotations, error = get_uniprot_data("P01308")

    assert error is None
    assert len(sequence) == 110
    assert annotations["disulfide bond"]


def test_binding_ligands_are_included():
    sequence, annotations, error = get_uniprot_data("P07550")
    frame = create_dataframe(sequence, annotations)

    assert error is None
    assert frame.loc[112, "Binding sites"] == (
        "(S)-carazolol (inverse agonist); (S)-timolol (inverse agonist)"
    )


def test_annotation_types_are_mapped_to_their_columns():
    sequence, annotations, error = get_uniprot_data("P07550")
    frame = create_dataframe(sequence, annotations)

    assert error is None
    assert frame.loc[140, "Phosphorylation sites"] == "Phosphotyrosine"
    assert frame.loc[381, "modified"] == "4-hydroxyproline"
    assert (frame.loc[391:412, "Disorder"] == "Disordered").all()
    assert (frame["Pfam entries (protein-level)"] == "PF00001 (7tm_1)").all()


def test_sequence_without_annotations_is_still_processed(monkeypatch, tmp_path):
    monkeypatch.setattr(app, "get_uniprot_data", lambda _: ("AC", {}, None))
    monkeypatch.setattr(app, "get_upload_folder", lambda: str(tmp_path))

    frame, output_path = app.process_uniprot_id("TEST")

    assert isinstance(frame, pd.DataFrame)
    assert frame["Residue code"].tolist() == ["A", "C"]
    assert Path(output_path).parent == tmp_path
    assert Path(output_path).is_file()


def test_missing_sequence_returns_a_displayable_error(monkeypatch):
    result = _parse_uniprot_xml(
        '<uniprot xmlns="http://uniprot.org/uniprot"><entry /></uniprot>'
    )
    monkeypatch.setattr(app, "get_uniprot_data", lambda _: result)

    message, download = app.process_uniprot_id("TEST")

    assert result == (
        None,
        {},
        "Could not find a protein sequence in the UniProt response.",
    )
    assert "Could not find" in message
    assert download is None


def test_network_errors_return_a_displayable_error(monkeypatch):
    error = HTTPError("https://example.test", 404, "Not Found", {}, None)

    def raise_error(*args, **kwargs):
        raise error

    monkeypatch.setattr("uniprot_data.urlopen", raise_error)

    sequence, annotations, message = get_uniprot_data("DOES_NOT_EXIST")

    assert sequence is None
    assert annotations == {}
    assert message == "Could not retrieve UniProt entry DOES_NOT_EXIST."


def test_unknown_and_out_of_range_positions_do_not_corrupt_the_table():
    sequence, annotations, error = _parse_uniprot_xml(
        """\
        <uniprot xmlns="http://uniprot.org/uniprot">
          <entry>
            <feature type="active site" description="unknown">
              <location><position status="unknown" /></location>
            </feature>
            <sequence>AC</sequence>
          </entry>
        </uniprot>
        """
    )
    annotations["glycosylation site"] = [
        {"position": 4, "description": "outside sequence"}
    ]

    frame = create_dataframe(sequence, annotations)

    assert error is None
    assert len(frame) == 2
    assert not frame["active sites"].astype(bool).any()
    assert not frame["Glycosylation sites"].astype(bool).any()
