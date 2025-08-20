# Codon Metrics Shiny App

This repository contains a Python [Shiny](https://shiny.posit.co/py/) web application for computing codon usage metrics from coding sequence (CDS) data.  
The app allows you to score sequences for **CAI**, **tAI**, **ENC**, and **Copt** metrics, visualize the results, and export them as tables or plots.

---

## Features

- **Input options**
  - Upload a FASTA file of CDS sequences.
  - Paste sequences directly into a text area (FASTA format or one sequence per line).
  - Upload a CSV table with sequence data (select ID, sequence, and optional species columns).

- **Built-in weights**
  - *E. coli* (CAI, tAI, Copt ratios).
  - *S. cerevisiae* (CAI, tAI, Copt ratios).

- **Metrics**
  - **CAI** (Codon Adaptation Index)
  - **tAI** (tRNA Adaptation Index)
  - **ENC** (Effective Number of Codons)
  - **Copt**
    - Gene-level ratio (geometric mean of codon High/Low ratios)
    - Gene-level log2 mean
    - Binary % optimal codons (≥ 1.0 threshold)

- **Outputs**
  - Interactive metrics table
  - Copy-pasteable TSV export
  - Boxplot visualizations (per metric or combined)
  - CSV download of computed metrics
  - Augmented CSV download (for Table mode: merges metrics into the original input table)

---

## Installation

1. Clone or download this repository.
2. Install dependencies (Python ≥3.9 recommended):

   ```bash
   pip install shiny biopython pandas matplotlib seaborn scikit-learn
Running the App
From the repository directory, run:

bash
Copy
Edit
shiny run --reload app.py
This will start a local server. By default, the app is accessible at:

cpp
Copy
Edit
http://127.0.0.1:8000
Usage
Choose input mode (sidebar):

Upload FASTA: supply .fa, .fasta, or .fna file.

Paste sequences: copy/paste directly into text box.

Upload table (CSV): load a CSV containing CDS sequences.

If using CSV table mode:

Select the ID, Sequence, and optional Species columns.

Or set a fixed species for all rows.

Species choice:

For FASTA/paste modes, select E. coli or S. cerevisiae for CAI/tAI weights.

Copt ratios are built in for both.

Outputs:

Codon Metrics Table: view results interactively.

Copy-paste table (TSV): quickly grab results for spreadsheets or R/Python.

Boxplots: visualize distributions of CAI, tAI, ENC, and Copt.

Downloads:

Metrics table as CSV

Boxplots as PNG

Augmented CSV (Table mode only: original table + computed metrics)

Example
Input CSV:

tx_id	species	CDS
gfp	ecoli	ATGAGTAAAGGAGAAGAACTTT…
mcherry	scer	ATGGTGAGCAAGGGCGAGGAG…

Output CSV (augmented):

tx_id	species	CDS	CAI	tAI	ENC	Copt_ratio (Zhou)	Copt_log2 (Zhou)	Copt (%)
gfp	ecoli	…	0.73	0.68	42.1	1.23	0.30	64.5
mcherry	scer	…	0.81	0.74	39.7	0.95	-0.07	48.9

Notes
Input sequences are automatically normalized (uppercase, U→T, trimmed to full codons).

Single-codon amino acids (Met, Trp) are treated as neutral (ratio = 1.0) in Copt calculations.

Exponential/log transforms require strictly positive values (e.g., when comparing metrics downstream).

Built-in weights currently cover E. coli and S. cerevisiae only, but additional species can be added by extending the weight dictionaries.

License
MIT License – feel free to adapt and extend.

Reference
Sharp et al., The codon Adaptation Index--a measure of directional synonymous codon usage bias, and its potential applications, Nucleic Acids Research (1987).
Wright, The ‘effective number of codons’ used in a gene, Gene (1990).
dos Reis et al., Solving the riddle of codon usage preferences: a test for translational selection, Nucleic Acids Research (2004).
Zhou et al., Measuring codon usage bias in microbial genomes, Nucleic Acids Research (2016).
