import pandas as pd
import numpy as np

# Load
df = pd.read_csv("te_table.csv")

# 1) Identify TE columns
te_cols = [c for c in df.columns if c.startswith("te")]

# 2) Replace bad values with NaN
df[te_cols] = df[te_cols].replace(["#DIV/0!", 0, 0.0], np.nan).astype(float)

# 3) Count non-missing TE values
df["TE_n_nonmissing"] = df[te_cols].notna().sum(axis=1)

# 4) Row-wise mean excluding NaN (will be NaN if all are NaN)
df["TE_mean"] = df[te_cols].mean(axis=1, skipna=True)

# 5) Explicit flag when all three are NaN
df["TE_all_missing"] = df["TE_n_nonmissing"] == 0

# (optional) if you want TE_mean explicitly set to NaN when all missing (it already will be)
df.loc[df["TE_all_missing"], "TE_mean"] = np.nan

# Save
df.to_csv("te_table_with_mean.csv", index=False)
print("Wrote te_table_with_mean.csv")