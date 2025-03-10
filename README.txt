Create a Gradio app which requests the uniprot ID of a protein and returns the protein sequence as a data frame.
Use the XML file returned. Parse the annotations and returns the following columns in the dataframe:

 - Residue Number
 - Residue code (1 letter) for the wild type residue
 - Secondary structure (taken from the annotations)
 - Pfam domain  
 - Disorder  
 - Disulfilde bridges 
 - Glycosylation sites
 - Phosphorylation sites 
 - active sites
 - metal binding sites
 - DNA binding sites
 - RNA binding sites
 - ligand binding sites
 - modified
 
 