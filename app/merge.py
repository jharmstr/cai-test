#!/usr/bin/env python3
import pandas as pd

TE_PATH = "te_table_with_mean.csv"   # must contain 'gene' and 'locus_tag'
CDS_PATH = "te_with_codons - Li.csv"           # produced from your GenBank parsing
OUT_PATH = "te_with_cds.csv"

# 1) Load
te = pd.read_csv(TE_PATH)
cds = pd.read_csv(CDS_PATH)

# 2) Decide which sequence column exists in CDS: prefer 'cds_seq', fallback to 'codon_seq'
seq_col = None
for candidate in ["cds_seq", "codon_seq"]:
    if candidate in cds.columns:
        seq_col = candidate
        break
if seq_col is None:
    raise ValueError(
        f"No sequence column found in {CDS_PATH}. "
        "Expected one of: 'cds_seq' or 'codon_seq'. "
        f"Found columns: {list(cds.columns)}"
    )

# 3) Build a slim CDS table with only needed columns (and only those that actually exist)
keep_cols = ["gene", "locus_tag", seq_col, "cds_len_nt", "protein_id", 'TE']
keep_cols = [c for c in keep_cols if c in cds.columns]
cds_slim = cds[keep_cols].copy()

# 4) Primary merge on locus_tag (dedupe on that key first)
if "locus_tag" in cds_slim.columns:
    cds_by_locus = cds_slim.dropna(subset=["locus_tag"]).drop_duplicates(subset=["locus_tag"])
    merged = te.merge(cds_by_locus, on="locus_tag", how="left", validate="m:1")
else:
    # If locus_tag isn't present in CDS, start with a gene merge directly
    merged = te.copy()

# 5) Fallback merge on gene and coalesce any missing values
if "gene" in te.columns and "gene" in cds_slim.columns:
    cds_by_gene = cds_slim.dropna(subset=["gene"]).drop_duplicates(subset=["gene"]).copy()

    # rename value columns so we can combine_first() cleanly
    value_cols = [c for c in [seq_col, "cds_len_nt", "protein_id"] if c in cds_by_gene.columns]
    rename_map = {c: f"{c}_gene" for c in value_cols}
    cds_by_gene = cds_by_gene.rename(columns=rename_map)

    # merge on gene
    merged_gene = te.merge(cds_by_gene, on="gene", how="left", validate="m:1")

    # coalesce: fill missing from gene-based merge
    for c in value_cols:
        gene_c = f"{c}_gene"
        if c in merged.columns and gene_c in merged_gene.columns:
            merged[c] = merged[c].combine_first(merged_gene[gene_c])
        elif c not in merged.columns and gene_c in merged_gene.columns:
            merged[c] = merged_gene[gene_c]  # if primary merge never created the column
else:
    print("Warning: could not run gene-based fallback; 'gene' not present in both tables.")

# 6) Report unmatched after both passes
n_unmatched = merged[seq_col].isna().sum()
print(f"Final unmatched rows (no sequence): {n_unmatched}")

# 7) Save
merged.to_csv(OUT_PATH, index=False)
print(f"Wrote {OUT_PATH} with sequence column '{seq_col}'")