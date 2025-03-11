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
    # Fetch XML data
    local_file_path = os.path.join("test", f"{uniprot_id}.xml")
    if os.path.exists(local_file_path):
        with open(local_file_path, "r", encoding="utf-8") as file:
            response = file.read()
    else:
        # Fetch XML data from UniProt
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
        response = urlopen(url).read().decode('utf-8')
    
    # Parse XML with namespace
    root = ET.fromstring(response)
    ns = {'up': 'http://uniprot.org/uniprot'}
    
    # Get sequence
    sequence_elem = root.find("./up:entry/up:sequence", ns)
    if sequence_elem is None:
        return None, None, "Could not find sequence in UniProt response"
    protein_sequence = sequence_elem.text.strip()
    
    # Get feature annotations
    annotations = {}
    for feature in root.findall(".//up:feature", ns):
        feature_type = feature.get('type')
        description = feature.get('description', '')
        
        # Get position information
        location = feature.find("up:location", ns)
        if location is None:
            continue
        
        # Handle different types of position elements
        position = location.find("up:position", ns)
        begin = location.find("up:begin", ns)
        end_elem = location.find("up:end", ns)
        
        if position is not None:
            pos = int(position.get('position'))
            # For single position features
            if feature_type not in annotations:
                annotations[feature_type] = []
            annotations[feature_type].append({
                'position': pos,
                'description': description
            })
        elif begin is not None and end_elem is not None:
            start = int(begin.get('position'))
            end = int(end_elem.get('position'))
            # For range features and disulfide bonds
            if feature_type not in annotations:
                annotations[feature_type] = []
            annotations[feature_type].append({
                'begin': start,
                'end': end,
                'description': description
            })
        
    return protein_sequence, annotations

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
            "Binding sites": "",  # Combined binding sites column
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
        'disulfide bond': 'Disulfide bridges',
        'glycosylation site': 'Glycosylation sites',
        'modified residue': 'modified',
        'active site': 'active sites',
        'site': 'Phosphorylation sites'
    }

    # Special mappings that need additional processing
    region_mapping = {
        'pfam': 'Pfam domain',
        'disorder': 'Disorder'
    }

    for feature_type, values in annotations.items():
        feature_type = feature_type.lower()
        
        # Handle disulfide bond pairs
        if feature_type == 'disulfide bond':
            for item in values:
                start = item['begin']
                end = item['end']
                desc = f"Cys-{end}"
                df.at[start-1, 'Disulfide bridges'] = desc
                desc = f"Cys-{start}"
                df.at[end-1, 'Disulfide bridges'] = desc
                
        # Handle glycosylation sites
        elif feature_type == 'glycosylation site':
            for item in values:
                pos = item['position'] - 1
                df.at[pos, 'Glycosylation sites'] = item['description']
                
        # Handle region features
        elif feature_type == 'region':
            for item in values:
                start = item.get('begin', item.get('position'))
                end = item.get('end', item.get('position'))
                if not start:
                    continue
                    
                start = int(start)
                end = int(end) if end else start
                desc = item['description'].lower()
                
                # Map to appropriate column based on description
                column = None
                if 'pfam' in desc:
                    column = 'Pfam domain'
                elif 'disorder' in desc:
                    column = 'Disorder'
                    
                if column:
                    for i in range(start - 1, end):
                        if i >= len(df):
                            continue
                        current = df.at[i, column]
                        if isinstance(current, str) and current != "" and desc:
                            df.at[i, column] = f"{current}; {desc}"
                        elif desc:
                            df.at[i, column] = desc
                            
        # Handle binding site features
        elif feature_type == 'binding site':
            for item in values:
                start = item.get('begin', item.get('position'))
                end = item.get('end', item.get('position'))
                if not start:
                    continue
                    
                start = int(start)
                end = int(end) if end else start
                desc = item['description']
                
                for i in range(start - 1, end):
                    if i >= len(df):
                        continue
                    current = df.at[i, 'Binding sites']
                    if isinstance(current, str) and current != "" and desc:
                        df.at[i, 'Binding sites'] = f"{current}; {desc}"
                    elif desc:
                        df.at[i, 'Binding sites'] = desc
        
        # Handle other features
        else:
            column = feature_mapping.get(feature_type)
            if not column:
                continue
                
            for item in values:
                start = item.get('begin', item.get('position'))
                end = item.get('end', item.get('position'))
                if not start:
                    continue
                    
                start = int(start)
                end = int(end) if end else start
                
                for i in range(start - 1, end):
                    if i >= len(df):
                        continue
                    if column == 'Secondary structure':
                        df.at[i, column] = feature_type.upper()
                    else:
                        current = df.at[i, column]
                        desc = item['description']
                        if isinstance(current, str) and current != "" and desc:
                            df.at[i, column] = f"{current}; {desc}"
                        elif desc:
                            df.at[i, column] = desc

    return df

def process_uniprot_id(uniprot_id):
    """
    Main function to process a UniProt ID.

    Args:
        uniprot_id: The UniProt ID.

    Returns:
        A Pandas DataFrame or an error message.
    """
    protein_sequence, annotations = get_uniprot_data(uniprot_id)

    if protein_sequence and annotations:
        df = create_dataframe(protein_sequence, annotations)
        return df
    else:
        return "Could not retrieve or process data for the given Uniprot ID"


# Gradio Interface
with gr.Blocks() as demo:
    with gr.Column():
        gr.Markdown("# Protein Sequence Analysis")
        input_text = gr.Textbox(
            label="UniProt ID", 
            placeholder="Enter UniProt ID (e.g., P53_HUMAN)",
            value="",  # Empty default value
        )
        submit_btn = gr.Button("Submit")
        
        # Add examples
        gr.Examples(
            examples=[
                ["P06280"],  # Alpha-galactosidase A
                ["P07550"],  # beta-2 AR
                ["P01308"],  # Insulin
                ["Q8WZ42"],  # Titin
                ["P0DTC2"],  # SARS-CoV-2 Spike protein
            ],
            example_labels=["Alpha-galactosidase A", "Beta-2 adrenergic receptor", "Insulin", "Titin", "SARS-CoV-2 Spike protein"],
            inputs=input_text,
            label="Example UniProt IDs"
        )

        output_df = gr.Dataframe(interactive=False)
        
        submit_btn.click(
            fn=process_uniprot_id,
            inputs=input_text,
            outputs=output_df
        )

if __name__ == "__main__":
    demo.launch()
