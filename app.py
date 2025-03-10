import gradio as gr
import pandas as pd
from io import StringIO
from urllib.request import urlopen
import re
import xml.etree.ElementTree as ET

def get_uniprot_data(uniprot_id):
    """
    Fetches protein sequence and annotation data from UniProt in XML format.

    Args:
        uniprot_id: The UniProt ID of the protein.

    Returns:
        A tuple containing:
        - protein_sequence: The protein sequence as a string.
        - annotations: A dictionary containing annotations.
        - error_message: An error message if something goes wrong, otherwise None
    """
    try:
        # Fetch XML data
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
        response = urlopen(url).read().decode('utf-8')
        
        # Parse XML
        root = ET.fromstring(response)
        
        # Get sequence
        sequence_elem = root.find(".//{http://uniprot.org/uniprot}sequence")
        if sequence_elem is None:
            return None, None, "Could not find sequence in UniProt response"
        protein_sequence = sequence_elem.text.strip()
        
        # Get feature annotations
        annotations = {}
        for feature in root.findall(".//{http://uniprot.org/uniprot}feature"):
            feature_type = feature.get('type')
            description = feature.get('description', '')
            
            # Get position information
            location = feature.find("{http://uniprot.org/uniprot}location")
            if location is None:
                continue
                
            position = location.find("{http://uniprot.org/uniprot}position")
            if position is not None:
                start = end = int(position.get('position'))
            else:
                begin = location.find("{http://uniprot.org/uniprot}begin")
                end_elem = location.find("{http://uniprot.org/uniprot}end")
                if begin is None or end_elem is None:
                    continue
                start = int(begin.get('position'))
                end = int(end_elem.get('position'))
            
            if feature_type not in annotations:
                annotations[feature_type] = []
            annotations[feature_type].append((start, end, description))
            
        return protein_sequence, annotations, None
        
    except Exception as e:
        return None, None, f"Error fetching or processing data from UniProt: {e}"

def create_dataframe(protein_sequence, annotations):
    """
    Creates a Pandas DataFrame from protein sequence and annotations.
    """
    data = []
    for i, residue in enumerate(protein_sequence):
        residue_number = i + 1
        row = {
            "Residue Number": residue_number,
            "Residue code": residue,
            "Secondary structure": "",
            "Domain": "",
            "Pfam domain": "",
            "Disorder": "",
            "Disulfide bridges": "",
            "Glycosylation sites": "",
            "Phosphorylation sites": "",
            "active sites": "",
            "metal binding sites": "",
            "DNA binding sites": "",
            "RNA binding sites": "",
            "ligand binding sites": "",
            "modified": ""
        }
        data.append(row)

    df = pd.DataFrame(data)

    # Map UniProt feature types to our column names
    feature_mapping = {
        'strand': 'Secondary structure',
        'helix': 'Secondary structure', 
        'turn': 'Secondary structure',
        'domain': 'Domain',
        'region': ['Pfam domain', 'Disorder'],
        'disulfide bond': 'Disulfide bridges',
        'glycosylation site': 'Glycosylation sites',
        'modified residue': 'modified',
        'active site': 'active sites',
        'binding site': ['metal binding sites', 'DNA binding sites', 'RNA binding sites', 'ligand binding sites'],
        'site': 'Phosphorylation sites'
    }

    for feature_type, values in annotations.items():
        feature_type = feature_type.lower()  # Convert to lowercase for matching
        for start, end, description in values:
            # Get the corresponding column(s)
            column = feature_mapping.get(feature_type)
            if not column:
                continue
                
            # Handle cases where one feature type maps to multiple possible columns
            if isinstance(column, list):
                if feature_type == 'region':
                    if 'Pfam' in description:
                        column = 'Pfam domain'
                    elif 'disorder' in description.lower():
                        column = 'Disorder'
                    else:
                        continue
                elif feature_type == 'binding site':
                    if 'metal' in description.lower():
                        column = 'metal binding sites'
                    elif 'DNA' in description:
                        column = 'DNA binding sites'
                    elif 'RNA' in description:
                        column = 'RNA binding sites'
                    else:
                        column = 'ligand binding sites'
                        
            # Fill in the annotation
            for i in range(start - 1, end):
                if i < len(df):
                    if column == 'Secondary structure':
                        df.loc[i, column] = feature_type.upper()  # Use uppercase for secondary structure
                    else:
                        current_value = df.loc[i, column]
                        if current_value:
                            df.loc[i, column] = f"{current_value}; {description}"
                        else:
                            df.loc[i, column] = description

    return df

def process_uniprot_id(uniprot_id):
    """
    Main function to process a UniProt ID.

    Args:
        uniprot_id: The UniProt ID.

    Returns:
        A Pandas DataFrame or an error message.
    """
    protein_sequence, annotations, error_message = get_uniprot_data(uniprot_id)

    if error_message:
        return error_message

    if protein_sequence and annotations:
        df = create_dataframe(protein_sequence, annotations)
        return df
    else:
        return "Could not retrieve or process data for the given Uniprot ID"


# Gradio Interface
iface = gr.Interface(
    fn=process_uniprot_id,
    inputs=gr.Textbox(label="UniProt ID", placeholder="e.g., P04637"),
    outputs=gr.Dataframe(label="Protein Sequence and Annotations"),
    title="UniProt Protein Sequence and Annotation Viewer",
    description="Enter a UniProt ID to view the protein sequence and its annotations in a DataFrame."
)

iface.launch()
