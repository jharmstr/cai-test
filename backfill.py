import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO

def parse_args():
    p = argparse.ArgumentParser(description="Backfill missing CDS sequences from GenBank by locus_tag (b-number).")
    p.add_argument("--in", dest="te_with_cds_clean.csv", required=True, help="Input CSV with columns: gene, locus_tag, and a sequence column (e.g., cds_seq).")
    p.add_argument("--genbank", dest="GCF_000005845.2_ASM584v2_genomic.gbff", required=True, help="GenBank file (.gb or .gbff) for E. coli K-12 MG1655 (e.g., GCF_000005845.2...gbff).")
    p.add_argument("--out", dest="te_with_cds_backfill.csv", required=True, help="Output CSV path.")
    p.add_argument("--seq-col", dest="seq_col", default="cds_seq", help="Name of the sequence column to fill (default: cds_seq).")
    p.add_argument("--include-stop", dest="include_stop", default="true",
                   help="Whether to keep terminal stop codon if present: true/false (default: true)")
    return p.parse_args()

def maybe_drop_stop(seq: str, include_stop: bool) -> str:
    seq = (seq or "").upper()
    if not include_stop and len(seq) >= 3 and len(seq) % 3 == 0:
        last = seq[-3:]
        if last in ("TAA", "TAG", "TGA"):
            return seq[:-3]
    return seq

def build_locus_to_cds(gb_path: str, include_stop: bool = True):
    """
    Returns dict: locus_tag -> dict(cds_seq, cds_len_nt, protein_id, gene)
    Skips pseudogenes.
    If multiple CDS share a locus_tag, keeps the longest CDS.
    """
    locus_map = {}
    n_records = 0
    n_cds = 0
    n_pseudo = 0
    for rec in SeqIO.parse(gb_path, "genbank"):
        n_records += 1
        for f in rec.features:
            if f.type != "CDS":
                continue
            quals = f.qualifiers
            if "pseudo" in quals:
                n_pseudo += 1
                continue
            locus_tag = (quals.get("locus_tag", [None])[0] or "").strip()
            if not locus_tag:
                continue
            try:
                seq_str = str(f.extract(rec.seq)).upper()
            except Exception:
                continue
            seq_str = maybe_drop_stop(seq_str, include_stop=include_stop)
            entry = {
                "cds_seq": seq_str,
                "cds_len_nt": len(seq_str),
                "protein_id": quals.get("protein_id", [None])[0],
                "gene": quals.get("gene", [None])[0],
            }
            # keep longest if duplicate locus_tag appears
            if (locus_tag not in locus_map) or (entry["cds_len_nt"] > locus_map[locus_tag]["cds_len_nt"]):
                locus_map[locus_tag] = entry
                n_cds += 1
    print(f"[INFO] Parsed {n_records} record(s). Kept {n_cds} CDS entries. Skipped {n_pseudo} pseudogenes.")
    return locus_map

def main():
    args = parse_args()
    include_stop = str(args.include_stop).strip().lower() in {"1","true","t","yes","y"}
    df = pd.read_csv(args.in_csv)

    # Ensure expected columns exist
    for col in ["locus_tag", args.seq_col]:
        if col not in df.columns:
            raise ValueError(f"Input CSV missing required column: {col}")

    # Identify rows needing backfill
    missing_mask = df[args.seq_col].isna() | (df[args.seq_col].astype(str).str.len() == 0)
    n_missing = int(missing_mask.sum())
    print(f"[INFO] Rows missing {args.seq_col}: {n_missing}")

    if n_missing == 0:
        print("[INFO] Nothing to backfill; writing a copy.")
        df.to_csv(args.out_csv, index=False)
        print(f"[OK] Wrote {args.out_csv}")
        return

    # Build locus_tag -> CDS map from GenBank
    locus_map = build_locus_to_cds(args.gb_path, include_stop=include_stop)

    filled = 0
    updated_len = 0
    # Also fill cds_len_nt and protein_id if present (create columns if absent)
    added_cols = []
    for col in ["cds_len_nt", "protein_id", "gene_from_gb"]:
        if col not in df.columns:
            df[col] = np.nan
            added_cols.append(col)

    # Backfill per row where missing
    for idx in df.index[missing_mask]:
        btag = str(df.at[idx, "locus_tag"]) if pd.notna(df.at[idx, "locus_tag"]) else ""
        btag = btag.strip()
        if not btag:
            continue
        if btag in locus_map:
            entry = locus_map[btag]
            df.at[idx, args.seq_col] = entry["cds_seq"]
            df.at[idx, "cds_len_nt"] = entry["cds_len_nt"]
            df.at[idx, "protein_id"] = entry["protein_id"]
            df.at[idx, "gene_from_gb"] = entry["gene"]
            filled += 1

    # Optional: if some rows had a shorter/empty seq, you can choose to overwrite only when missing.
    # (Current behavior only fills previously-missing rows.)

    remaining = int(df[args.seq_col].isna().sum() | (df[args.seq_col].astype(str).str.len() == 0).sum())
    print(f"[INFO] Filled sequences for {filled} row(s). Still missing: {remaining}")

    df.to_csv(args.out_csv, index=False)
    print(f"[OK] Wrote {args.out_csv}")
    if added_cols:
        print(f"[INFO] Added columns: {', '.join(added_cols)}")

if __name__ == "__main__":
    main()