import argparse
import pandas as pd
import numpy as np
import os
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

# --------------------------
# 1. Preprocessing
# --------------------------
def preprocess_lcms(file, gene_col, intensity_cols, impute=True):
    """Preprocess LC-MS data: filter, subset, and optionally impute missing values."""
    if file.endswith(".xlsx"):
        df = pd.read_excel(file)
    else:
        df = pd.read_csv(file)
    
    # Filter out contaminants
    df = df[~df[gene_col].str.contains("^contaminant", case=False, na=False)]
    
    # Keep only required columns
    cols_to_keep = [gene_col] + intensity_cols
    df = df[cols_to_keep]
    
    # Impute missing values (optional)
    if impute:
        df[intensity_cols] = df[intensity_cols].fillna(df[intensity_cols].median())
    
    return df

# --------------------------
# 2. Differential analysis
# --------------------------
def differential_analysis(expression_df, metadata_df, gene_col, group_col):
    """Compute log2 fold-change between two groups with t-test and FDR."""
    # Expect exactly two groups
    groups = metadata_df[group_col].unique()
    if len(groups) != 2:
        raise ValueError("Currently only supports two-group comparisons")
    
    grp1_samples = metadata_df[metadata_df[group_col] == groups[0]]['Sample'].tolist()
    grp2_samples = metadata_df[metadata_df[group_col] == groups[1]]['Sample'].tolist()
    
    # Compute mean intensities
    expression_df['mean_grp1'] = expression_df[grp1_samples].mean(axis=1)
    expression_df['mean_grp2'] = expression_df[grp2_samples].mean(axis=1)
    
    # Log2 fold-change
    expression_df['logFC'] = np.log2(expression_df['mean_grp2'] + 1) - np.log2(expression_df['mean_grp1'] + 1)
    
    # T-test p-values
    pvals = []
    for idx, row in expression_df.iterrows():
        tstat, pval = ttest_ind(row[grp1_samples], row[grp2_samples], equal_var=False, nan_policy='omit')
        pvals.append(pval)
    expression_df['pvalue'] = pvals
    
    # FDR correction
    _, fdrs, _, _ = multipletests(expression_df['pvalue'], method='fdr_bh')
    expression_df['FDR'] = fdrs
    
    return expression_df[[gene_col, 'logFC', 'pvalue', 'FDR'] + grp1_samples + grp2_samples]

# --------------------------
# 3. Volcano plot
# --------------------------
def generate_volcano(df, gene_col, fc_col, fdr_col, fc_thresh=1.5, fdr_thresh=0.05, out_file="volcano.png"):
    """Generate a volcano plot from differential expression results."""
    df['-log10FDR'] = -np.log10(df[fdr_col])
    
    # Classification
    df['diffexpressed'] = 'NS'
    df.loc[(df[fc_col] >= fc_thresh) & (df[fdr_col] < fdr_thresh), 'diffexpressed'] = 'Up'
    df.loc[(df[fc_col] <= -fc_thresh) & (df[fdr_col] < fdr_thresh), 'diffexpressed'] = 'Down'
    
    # Plot
    plt.figure(figsize=(6,5))
    sns.scatterplot(data=df, x=fc_col, y='-log10FDR', hue='diffexpressed', palette={'Up':'red', 'Down':'blue', 'NS':'grey'}, edgecolor=None, s=30)
    
    # Label top hits
    for _, row in df.iterrows():
        if row['diffexpressed'] != 'NS':
            plt.text(row[fc_col], row['-log10FDR'], row[gene_col], fontsize=6)
    
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(FDR)")
    plt.title("Volcano Plot")
    plt.legend(title='')
    plt.tight_layout()
    
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"Volcano plot saved to {out_file}")

# --------------------------
# 4. Main
# --------------------------
def main():
    parser = argparse.ArgumentParser(description="General LC-MS preprocessing, differential analysis, and volcano plot")
    parser.add_argument("--input", required=True, help="LC-MS data file (.xlsx or .csv)")
    parser.add_argument("--metadata", required=True, help="Metadata file (CSV with Sample, Group columns)")
    parser.add_argument("--gene_col", required=True, help="Column name for gene/protein IDs")
    parser.add_argument("--intensity_cols", required=True, nargs="+", help="Columns with intensity values")
    parser.add_argument("--group_col", required=True, help="Column in metadata indicating groups")
    parser.add_argument("--out_dir", default="output", help="Output directory")
    parser.add_argument("--fc_thresh", type=float, default=1.5, help="Fold-change threshold for volcano plot")
    parser.add_argument("--fdr_thresh", type=float, default=0.05, help="FDR threshold for volcano plot")
    args = parser.parse_args()
    
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Preprocess
    expression_df = preprocess_lcms(args.input, args.gene_col, args.intensity_cols)
    metadata_df = pd.read_csv(args.metadata)
    
    # Differential analysis
    diff_df = differential_analysis(expression_df, metadata_df, args.gene_col, args.group_col)
    
    # Save table
    table_file = os.path.join(args.out_dir, "diff_expression.tsv")
    diff_df.to_csv(table_file, sep="\t", index=False)
    print(f"Differential expression table saved at {table_file}")
    
    # Volcano plot
    volcano_file = os.path.join(args.out_dir, "volcano.png")
    generate_volcano(diff_df, args.gene_col, 'logFC', 'FDR', args.fc_thresh, args.fdr_thresh, volcano_file)

if __name__ == "__main__":
    main()
