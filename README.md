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
- Associated Pfam domain
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

1. Launch the application
2. Enter a valid UniProt ID (e.g., P53_HUMAN) in the input field
3. Click "Submit" to generate the analysis
4. Results will be displayed in a interactive data frame format

## Requirements

- Python 3.7+
- Gradio
- Pandas
- Requests
- XML parsing libraries

## Note

The application processes UniProt's XML format to extract annotations. 

