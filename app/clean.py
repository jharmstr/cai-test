import pandas as pd

te = pd.read_csv("TE_with_cds_backfill.csv")

te = te.drop(columns=['gene_from_gb'])   # assign back

te.to_csv("TE_with_cds_backfill_clean.csv", index=False)