---
title: Protein Sequence Table
colorFrom: pink
colorTo: gray
sdk: gradio
app_file: app.py
pinned: false
---


# Protein Sequence Table

A Gradio-based web application that reformats protein sequences based on UniProt IDs and displays detailed annotations in a structured format.

## Features

The application retrieves protein data from UniProt and presents the following information for each residue:
- Position number in the sequence
- Amino acid (single-letter code)
- Secondary structure annotation
- Protein-level Pfam cross-references
- Disorder prediction
- Participation in disulfide bridges
- Post-translational modifications:
  * Glycosylation sites
  * Phosphorylation sites
- Functional annotations:
  * Active sites
  * Metal binding sites
  * DNA binding regions
  * RNA binding regions
  * Ligand binding sites
  * Other modifications

## Usage

1. Install [uv](https://docs.astral.sh/uv/) if it is not already available.
2. Create the project environment and install the locked dependencies:

   ```bash
   uv sync
   ```

3. Launch the application:

   ```bash
   uv run gradio app.py
   ```

4. Enter a valid UniProt ID (e.g., P53_HUMAN) in the input field.
5. Click "Submit" to generate the analysis.
6. Results will be displayed in an interactive data frame format.

To run the linter:

```bash
uv run ruff check .
```

To run the test suite:

```bash
uv run pytest
```

## Requirements

- Python 3.11+
- Dependencies defined in `pyproject.toml` and locked in `uv.lock`

## Note

The application processes UniProt's XML format to extract annotations. 
