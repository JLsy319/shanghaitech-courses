## Q1

Please use PDB website to genreate a ligand report for human tyrosine-protein kinase ABL1 (Hint: First, find its UniProt; then, use the UniProt to search at PDB website). In the report, please only keep PDB entries with structure resolution higher than 2.5 Angstroms (Note: the jargon 'high resolution' means better resolution, that means the resolution number is smaller than 2.5 Angstroms), and please do not keep the ligands with Molecular Weight < 250 since they are likely to be solvent, please also remove ATP and ATP-like ligands (i.e., only keep drug like ligands).

## Q2

Please (use RDKit or other tools to) generate 2 different binary fingerprints (e.g., MACSS, ECFP4) for the ligands of human tyrosine-protein kinase ABL1 in your report. Please submit an excel file which should have 5 columns (PDB-ID, 3-digit alphanumerial ID of ligand, smiles of ligand, fp1 of ligand, fp2 of ligand).

[Hint: You may start from the 3-digit alphanumerial ID of each ligand (e.g. STI is a ligand of human tyrosine-protein kinase ABL1 from PDB 2HYY, the PDB page for the ligand is https://www.rcsb.org/ligand/STI, you can find the smiles representation/format of STI in this page. Once you have its smiles, you can readily convert it to binary fingerprints.]

## Q3

Use Tanimoto coeff cutoff = 0.5 and 0.3 to group the ligands in your report, respectively. Please submit your group results.  

## Q4

Find the ChEMBL target ID for human tyrosine-protein kinase ABL. At ChEMBL website, go to the webpage of human tyrosine-protein kinase ABL, use the navigation panel to help you figure out 

1. How many associated compounds have been tested against this target.
2. In these compounds, how many of them have IC50 data? Choose those with IC50 data, can you filter out the part that with IC50 <= 0.1 nM (Hint, use pIC50 to filer)? How do you realize it? Download the CSV table for these highly potent compounds as your answer.

## Q5

1IEP is the first crystal structure of mouse tyrosine-protein kinase ABL in complex with the famous drug imatinib (伊马替尼，商品名：格列卫） https://www.rcsb.org/structure/1IEP ). Please follow the UCSF ChimeraX tutorial for protein-ligand binding site (https://www.cgl.ucsf.edu/chimerax/docs/user/tutorials/binding-sites.html), generate a picture for protein-ligand binding for 1IEP. Please submit your final picture, and the command lines (in time line order) you used to generate the picture.