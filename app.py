import io
import gradio as gr
import tempfile

from uniprot_data import create_dataframe, get_uniprot_data


def process_uniprot_id(uniprot_id):
    """
    Main function to process a UniProt ID.

    Args:
        uniprot_id: The UniProt ID.

    Returns:
        Tuple of (DataFrame, downloadable file) or error message
    """
    protein_sequence, annotations = get_uniprot_data(uniprot_id)

    if protein_sequence and annotations:
        df = create_dataframe(protein_sequence, annotations)
        # Create Excel file
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx')
        df.to_excel(temp_file.name, index=False)
        return df, temp_file.name
    return "Could not retrieve or process data for the given Uniprot ID", None


# Gradio Interface
with gr.Blocks() as demo:
    with gr.Column():
        gr.Markdown("# Protein Sequence Analysis")
        gr.Markdown(
            "This app fetches protein sequence and annotation data from UniProt using a UniProt ID and prints a copy-pasteable table for note-taking. **DO NOT TRUST, THIS IS A CODING EXPERIMENT**"
        )
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
            example_labels=[
                "Alpha-galactosidase A",
                "Beta-2 adrenergic receptor",
                "Insulin",
                "Titin",
                "SARS-CoV-2 Spike protein",
            ],
            inputs=input_text,
            label="Example UniProt IDs",
        )

        output_df = gr.Dataframe(interactive=False)
        download_btn = gr.DownloadButton(label="Download Excel")

        outputs = submit_btn.click(
            fn=process_uniprot_id, 
            inputs=input_text, 
            outputs=[output_df, download_btn]
        )

if __name__ == "__main__":
    demo.launch()
