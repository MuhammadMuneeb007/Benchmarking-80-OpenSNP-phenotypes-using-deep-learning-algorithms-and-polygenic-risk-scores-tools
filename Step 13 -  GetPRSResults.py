"""
Get PRS results for all three tools: Plink, PRSice, Lassosum
Handles unequal fold sizes by padding shorter folds with NaN,
then averaging per hyperparameter position using only valid folds.

Run:
    python GetPRSResults.py AUC Plink
    python GetPRSResults.py AUC PRSice
    python GetPRSResults.py AUC Lassosum
"""

import pandas as pd
import numpy as np
import os
from os.path import exists
import re
import sys

def sorted_nicely(l):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

# ============================================================
# ARGUMENTS
# ============================================================
metric = sys.argv[1]   # e.g. AUC
tool   = sys.argv[2]   # e.g. Plink, PRSice, Lassosum

TOOL_COLUMN = {
    "Plink":    "PlinkAUC",
    "PRSice":   "PRSiceAUC",
    "Lassosum": "LassosumAUC",
}

if tool not in TOOL_COLUMN:
    print(f"ERROR: tool must be one of {list(TOOL_COLUMN.keys())}")
    sys.exit(1)

col_name    = TOOL_COLUMN[tool]
result_file = f"Results_{tool}_{metric}.txt"
MAX_ROWS    = 675  # maximum possible hyperparameter combinations

print(f"Tool        : {tool}")
print(f"Metric      : {metric}")
print(f"Result file : {result_file}")
print(f"Column name : {col_name}")

# ============================================================
# LOAD PHENOTYPE LIST
# ============================================================
phenotypes = pd.read_csv("allphenotypesname2.txt", header=None)[0].values
print(f"\nTotal phenotypes to process: {len(phenotypes)}")

# ============================================================
# AGGREGATE RESULTS
# ============================================================
disease = []
auc     = []
STD     = []
ALGO    = []

for loop in phenotypes:
    # Find which folds exist
    fold_files = []
    for loop2 in range(1, 6):
        path = f"./{loop}/{loop2}/{result_file}"
        if exists(path):
            fold_files.append((loop2, path))

    if len(fold_files) < 5:
        print(f"  SKIP {loop}: only {len(fold_files)}/5 folds found")
        continue

    try:
        # --------------------------------------------------------
        # Load each fold and pad to MAX_ROWS with NaN
        # This handles the case where different folds completed
        # different numbers of hyperparameter combinations
        # --------------------------------------------------------
        padded_folds = []
        fold_row_counts = []

        for loop2, path in fold_files:
            df = pd.read_csv(path, sep=",")

            if col_name not in df.columns:
                print(f"  WARNING: column '{col_name}' not found in {path}")
                print(f"  Available columns: {df.columns.tolist()}")
                # Pad entirely with NaN
                padded = np.full(MAX_ROWS, np.nan)
            else:
                values = df[col_name].values.astype(float)
                n = len(values)
                fold_row_counts.append(n)

                # Pad shorter folds with NaN up to MAX_ROWS
                if n < MAX_ROWS:
                    padded = np.full(MAX_ROWS, np.nan)
                    padded[:n] = values
                else:
                    padded = values[:MAX_ROWS]

            padded_folds.append(padded)

        if len(fold_row_counts) > 0:
            print(f"  {loop}: fold sizes = {fold_row_counts}")

        # Stack into matrix: shape (5, MAX_ROWS)
        matrix = np.vstack(padded_folds)

        # Replace inf with nan
        matrix[np.isinf(matrix)] = np.nan

        # Count valid folds per hyperparameter position
        valid_counts = np.sum(~np.isnan(matrix), axis=0)

        # Average using nanmean — divides by actual valid count per position
        average = np.nanmean(matrix, axis=0)

        # Standard deviation using nanstd ddof=1
        # Only compute where we have at least 2 valid folds
        std_vals = np.where(
            valid_counts >= 2,
            np.nanstd(matrix, axis=0, ddof=1),
            np.nan
        )

        # Only consider hyperparameter positions with at least 3 valid folds
        valid_mask = valid_counts >= 3
        if not np.any(valid_mask):
            print(f"  SKIP {loop}: no hyperparameter had >= 3 valid folds")
            continue

        # Find best hyperparameter among those with >= 3 valid folds
        masked_average = np.where(valid_mask, average, np.nan)
        best_idx = np.nanargmax(masked_average)
        best_auc = average[best_idx]
        best_std = std_vals[best_idx] if not np.isnan(std_vals[best_idx]) else 0.0
        best_count = valid_counts[best_idx]

        disease.append(loop)
        auc.append(best_auc)
        STD.append(best_std)
        ALGO.append(f"{tool}_{best_idx + 1}")

        print(f"  {loop}: best AUC = {best_auc*100:.2f}% "
              f"(hyperparameter {best_idx+1}, valid in {best_count}/5 folds)")

    except Exception as e:
        print(f"  ERROR processing {loop}: {e}")
        continue

# ============================================================
# SAVE RESULTS
# ============================================================
print(f"\nSuccessfully processed: {len(disease)} phenotypes")

hu = pd.DataFrame()
hu['Phenotype']                            = disease
hu[f'Test {metric} 5 Iterations Average'] = np.array(auc) * 100
hu['Standard Deviation']                  = np.array(STD) * 100
hu['Best parameters index']               = ALGO

out_csv  = f"PRS_{tool}AUC_{metric}basedbechmarking.csv"
out_html = f"PRS_{tool}AUC_{metric}basedbechmarking.html"

hu.to_csv(out_csv, index=False, sep=",")
hu.to_html(out_html)

print(f"\nSaved: {out_csv}")
print(f"Saved: {out_html}")
print(f"\nSummary:")
print(f"  Phenotypes processed : {len(disease)}")
print(f"  Mean AUC             : {hu[f'Test {metric} 5 Iterations Average'].mean():.2f}%")
print(f"  Min AUC              : {hu[f'Test {metric} 5 Iterations Average'].min():.2f}%")
print(f"  Max AUC              : {hu[f'Test {metric} 5 Iterations Average'].max():.2f}%")
print("\nFirst 10 results:")
print(hu.head(10).to_string(index=False))
