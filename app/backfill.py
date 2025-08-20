import pandas as pd
import numpy as np
from Bio import SeqIO

# -----------------------------
# Hardcoded paths and settings
# -----------------------------
INPUT_CSV    = "TE_with_cds_clean.csv"
GENBANK_GB   = "GCF_000005845.2_ASM584v2_genomic.gbff"
OUTPUT_CSV   = "TE_with_cds_backfill.csv"
SEQ_COL      = "cds_seq"   # column to fill
INCLUDE_STOP = True        # keep terminal stop codon if present

UNRESOLVED_LIST = "unresolved_btags.txt"
STILL_MISSING_CSV = "still_missing_rows.csv"
# -----------------------------

def maybe_drop_stop(seq: str, include_stop: bool) -> str:
    seq = (seq or "").upper()
    if not include_stop and len(seq) >= 3 and len(seq) % 3 == 0:
        if seq[-3:] in ("TAA", "TAG", "TGA"):
            return seq[:-3]
    return seq

def build_locus_to_cds(gb_path: str, include_stop: bool = True):
    """Return dict: locus_tag -> {cds_seq, cds_len_nt, protein_id, gene} (skip pseudogenes).
       If duplicates exist, keep the longest CDS.
    """
    locus_map = {}
    for rec in SeqIO.parse(gb_path, "genbank"):
        for f in rec.features:
            if f.type != "CDS":
                continue
            if "pseudo" in f.qualifiers:
                continue
            locus_tag = (f.qualifiers.get("locus_tag", [None])[0] or "").strip()
            if not locus_tag:
                continue
            try:
                seq_str = str(f.extract(rec.seq)).upper()
            except Exception:
                continue
            seq_str = maybe_drop_stop(seq_str, include_stop)
            entry = {
                "cds_seq": seq_str,
                "cds_len_nt": len(seq_str),
                "protein_id": f.qualifiers.get("protein_id", [None])[0],
                "gene": f.qualifiers.get("gene", [None])[0],
            }
            if (locus_tag not in locus_map) or (entry["cds_len_nt"] > locus_map[locus_tag]["cds_len_nt"]):
                locus_map[locus_tag] = entry
    return locus_map

def main():
    df = pd.read_csv(INPUT_CSV)

    # Ensure required columns exist
    if "locus_tag" not in df.columns:
        raise ValueError("Input CSV must contain a 'locus_tag' column.")
    if SEQ_COL not in df.columns:
        df[SEQ_COL] = np.nan

    # Identify rows missing sequence
    missing_mask = df[SEQ_COL].isna() | (df[SEQ_COL].astype(str).str.len() == 0)
    n_missing_init = int(missing_mask.sum())
    print(f"[INFO] Rows missing {SEQ_COL} before backfill: {n_missing_init}")

    if n_missing_init == 0:
        df.to_csv(OUTPUT_CSV, index=False)
        # Empty reports
        open(UNRESOLVED_LIST, "w").close()
        pd.DataFrame(columns=["gene","locus_tag"]).to_csv(STILL_MISSING_CSV, index=False)
        print(f"[OK] Nothing to fill. Wrote {OUTPUT_CSV}")
        return

    # Build map from GenBank
    locus_map = build_locus_to_cds(GENBANK_GB, include_stop=INCLUDE_STOP)

    # Ensure metadata cols exist
    for col in ["cds_len_nt", "protein_id", "gene_from_gb"]:
        if col not in df.columns:
            df[col] = np.nan

    # Backfill by locus_tag
    filled = 0
    unresolved_btags = set()

    for idx in df.index[missing_mask]:
        btag = df.at[idx, "locus_tag"]
        if pd.isna(btag):
            unresolved_btags.add("")  # blank/missing btag
            continue
        btag = str(btag).strip()
        if not btag or btag not in locus_map:
            unresolved_btags.add(btag)
            continue

        entry = locus_map[btag]
        df.at[idx, SEQ_COL]       = entry["cds_seq"]
        df.at[idx, "cds_len_nt"]  = entry["cds_len_nt"]
        df.at[idx, "protein_id"]  = entry["protein_id"]
        df.at[idx, "gene_from_gb"]= entry["gene"]
        filled += 1

    # Recompute missing after fill (don’t double-count)
    missing_mask_after = df[SEQ_COL].isna() | (df[SEQ_COL].astype(str).str.len() == 0)
    n_missing_final = int(missing_mask_after.sum())
    print(f"[INFO] Filled {filled} row(s). Still missing: {n_missing_final}")

    # Write outputs
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"[OK] Wrote {OUTPUT_CSV}")

    # Unresolved list (only the ones that were actually missing at start)
    # Filter unresolved_btags to those that correspond to rows still missing
    still_missing_tags = (
        df.loc[missing_mask_after, "locus_tag"]
          .astype(str)
          .str.strip()
          .tolist()
    )
    unresolved_filtered = sorted(set([t for t in unresolved_btags if t in still_missing_tags]))
    with open(UNRESOLVED_LIST, "w") as fh:
        for t in unresolved_filtered:
            if t:  # skip empty
                fh.write(f"{t}\n")
    print(f"[OK] Wrote unresolved b-tags to {UNRESOLVED_LIST} ({len(unresolved_filtered)} item(s)).")

    # Also dump a small CSV with the still-missing rows (handy for inspection)
    cols_to_show = [c for c in ["gene","locus_tag", SEQ_COL, "cds_len_nt", "protein_id"] if c in df.columns]
    df.loc[missing_mask_after, cols_to_show].to_csv(STILL_MISSING_CSV, index=False)
    print(f"[OK] Wrote still-missing rows to {STILL_MISSING_CSV}")

if __name__ == "__main__":
    main()