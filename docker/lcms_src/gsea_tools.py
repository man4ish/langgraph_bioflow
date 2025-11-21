#!/usr/bin/env python3
"""
gsea_tools.py

Python reimplementation of gsea_tools.R functionality:
- volcano plotting
- TERM2GENE extraction (msigdbr, Reactome)
- preranked GSEA using gseapy
- dotplot saving (png/pdf)
- combined plots and NES heatmap

Usage examples:
python gsea_tools.py --infile data/diff_expression.tsv --gene-col Gene \
    --fc-col logFC --sample-label SampleA --out-prefix sampleA

python gsea_tools.py --mode batch --batch-file samples_list.tsv --out-prefix study
"""

import argparse
import os
import sys
import math
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# try imports that might not be installed
try:
    import msigdbr
except Exception as e:
    msigdbr = None
try:
    import gseapy as gp
except Exception as e:
    gp = None

# adjustText is useful for label repelling (like ggrepel)
try:
    from adjustText import adjust_text
except Exception:
    adjust_text = None

# -------------------------
# Utilities
# -------------------------
def ensure_dir(path):
    d = os.path.dirname(path)
    if d and not os.path.exists(d):
        os.makedirs(d, exist_ok=True)


# -------------------------
# Volcano plot
# -------------------------
def volcano_plot(df,
                 logfc_col,
                 fdr_col,
                 gene_name_col,
                 fc_thresh_high=1.5,
                 fc_thresh_low=1.0,
                 logfdr_thresh=1.3,
                 plot_title="Volcano Plot",
                 only_positive_logfc=False,
                 out_file=None,
                 figsize=(6, 4)):
    """
    Create a volcano plot. Saves to out_file if provided.
    - df: pandas DataFrame
    - logfc_col: name of column containing log2 fold-change
    - fdr_col: name of column containing FDR / adjusted p-value
    - gene_name_col: name of gene identifier column
    """
    d = df.copy()

    # Ensure numeric
    d[logfc_col] = pd.to_numeric(d[logfc_col], errors="coerce")
    d[fdr_col] = pd.to_numeric(d[fdr_col], errors="coerce")
    d[gene_name_col] = d[gene_name_col].astype(str)

    if only_positive_logfc:
        d = d[d[logfc_col] > 0].copy()

    d["logFDR"] = -np.log10(d[fdr_col].replace(0, np.nan))
    d["logFDR"] = d["logFDR"].fillna(0)
    d[logfc_col] = d[logfc_col].fillna(0)

    # Label criteria
    d["diffexpressed"] = "NO"
    mask_label = (d[logfc_col] > fc_thresh_high) & (d["logFDR"] > logfdr_thresh)
    d.loc[mask_label, "diffexpressed"] = "LABEL"
    d["label"] = ""
    d.loc[mask_label, "label"] = d.loc[mask_label, gene_name_col]

    # Plot
    plt.figure(figsize=figsize)
    sns.scatterplot(data=d, x=logfc_col, y="logFDR", hue="diffexpressed",
                    palette={"NO": "grey", "LABEL": "red"}, legend=False, s=20)
    plt.title(plot_title)
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(FDR)")
    plt.grid(False)

    # add labels using adjustText if available
    texts = []
    label_df = d[d["label"] != ""]
    for _, row in label_df.iterrows():
        texts.append(plt.text(row[logfc_col], row["logFDR"], row["label"], fontsize=8))

    if adjust_text and texts:
        try:
            adjust_text(texts, arrowprops=dict(arrowstyle='-', color='0.5', lw=0.5))
        except Exception:
            pass

    if out_file:
        ensure_dir(out_file)
        plt.savefig(out_file, bbox_inches="tight", dpi=300)
        plt.close()
        print(f"[OK] Volcano plot saved to: {out_file}")
    else:
        plt.show()
        plt.close()

    return d


# -------------------------
# Prepare TERM2GENE (Reactome)
# -------------------------
def prepare_term2gene_reactome():
    if msigdbr is None:
        raise ImportError("msigdbr is required for TERM2GENE (pip install msigdbr)")

    # get C2 CP:REACTOME for Homo sapiens
    df = msigdbr.msigdb(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
    # msigdbr returns columns 'gs_name' and 'gene_symbol'
    term2gene = df[["gs_name", "gene_symbol"]].drop_duplicates()
    term2gene.columns = ["term", "gene"]
    return term2gene


# -------------------------
# Run preranked GSEA using gseapy (fgsea-like)
# -------------------------
def run_gsea_prerank(gene_fc_series, term2gene_df,
                     minGSSize=10, maxGSSize=500,
                     out_prefix="gsea", sample_label="sample", mode="all",
                     top_categories=10, out_dir="out"):
    """
    gene_fc_series: pandas Series indexed by gene name, values are numeric scores (logFC)
    term2gene_df: DataFrame with columns ['term', 'gene']
    mode: "all", "pos", "neg" - determines filtering applied prior to ranking
    """
    if gp is None:
        raise ImportError("gseapy is required (pip install gseapy)")

    s = gene_fc_series.copy().dropna()
    if mode == "pos":
        s = s[s > 0]
    elif mode == "neg":
        s = s[s < 0]

    if s.shape[0] == 0:
        print(f"[WARN] No genes found for mode {mode} in {sample_label}. Skipping.")
        return {"results": pd.DataFrame(), "plot": None}

    # keep strongest per gene already assumed (series indexed uniquely)
    # gseapy.prerank requires a two-column file or series; we'll use prerank with rnk DataFrame
    rnk = s.sort_values(ascending=False)
    # write rnk to temporary file (gseapy accepts series via DataFrame too)
    # build gene sets dict from term2gene
    gs = {}
    for term, grp in term2gene_df.groupby("term"):
        genes = grp["gene"].astype(str).unique().tolist()
        # filter gene sets by size: we'll rely on gseapy min/maxGSSize
        gs[term] = genes

    # gseapy prerank: provide the rnk as a file or DataFrame with gene name and score
    # gseapy.prerank accepts 'rnk' as path or DataFrame (2 columns)
    rnk_df = pd.DataFrame({"Gene": rnk.index, "Score": rnk.values})
    # gseapy expects gene names in first column, scores second, sorted descending
    try:
        pre_res = gp.prerank(rnk=rnk_df,
                             gene_sets=gs,
                             min_size=minGSSize,
                             max_size=maxGSSize,
                             permutation_num=1000,  # higher for stable p-values; change if slow
                             processes=1,
                             outdir=None,  # we will not rely on gseapy's outdir, we capture result
                             seed=123,
                             verbose=False)
    except Exception as e:
        # fallback to calling prerank with temporary file if needed
        tmp_rnk = f"{out_prefix}_{sample_label}_{mode}.rnk.tmp"
        rnk_df.to_csv(tmp_rnk, sep="\t", header=False, index=False)
        pre_res = gp.prerank(rnk=tmp_rnk,
                             gene_sets=gs,
                             min_size=minGSSize,
                             max_size=maxGSSize,
                             permutation_num=1000,
                             processes=1,
                             outdir=None,
                             seed=123,
                             verbose=False)
        try:
            os.remove(tmp_rnk)
        except Exception:
            pass

    # pre_res.res2d is a pandas DataFrame with GSEA results
    if hasattr(pre_res, "res2d"):
        gsea_df = pre_res.res2d.reset_index().rename(columns={"index": "ID"})
    else:
        # gseapy returns object with attribute results; try to convert
        gsea_df = pd.DataFrame(pre_res.results).reset_index().rename(columns={"index": "ID"})

    if gsea_df.shape[0] == 0:
        print(f"[INFO] No significant GSEA results for {sample_label} ({mode}).")
        return {"results": pd.DataFrame(), "plot": None}

    # select columns approx similar to clusterProfiler output
    keep_cols = [c for c in ["ID", "nes", "pval", "fdr", "nes_norm", "nes", "es", "pval", "fdr"] if c in gsea_df.columns]
    # gseapy uses 'nes','pval','fdr' names typically; rename to consistent names
    rename_map = {}
    if "nes" in gsea_df.columns:
        rename_map["nes"] = "NES"
    if "fdr" in gsea_df.columns:
        rename_map["fdr"] = "p.adjust"
    if "pval" in gsea_df.columns:
        rename_map["pval"] = "pvalue"
    if "es" in gsea_df.columns:
        rename_map["es"] = "ES"
    if "lead_genes" in gsea_df.columns:
        # keep but not required
        pass

    gsea_df = gsea_df.rename(columns=rename_map)

    # Ensure columns exist
    for c in ["ID", "NES", "pvalue", "p.adjust"]:
        if c not in gsea_df.columns:
            gsea_df[c] = np.nan

    # add sample label
    gsea_df["sample"] = sample_label

    # Save table
    out_tables_dir = os.path.join(out_dir, "tables")
    os.makedirs(out_tables_dir, exist_ok=True)
    suffix = "" if mode == "all" else f"_{mode}FC"
    tsv_file = os.path.join(out_tables_dir, f"{out_prefix}{suffix}_results.tsv")
    gsea_df.to_csv(tsv_file, sep="\t", index=False)
    print(f"[OK] GSEA results table saved at: {tsv_file}")

    # Dotplot: create a dotplot for top categories by absolute NES or by padj
    # take top_categories by absolute NES
    plot_df = gsea_df.copy().dropna(subset=["NES"]).sort_values("NES", ascending=False)
    # pick top by abs NES
    plot_df["absNES"] = plot_df["NES"].abs()
    top_df = plot_df.sort_values("absNES", ascending=False).head(top_categories)

    if top_df.shape[0] == 0:
        return {"results": gsea_df, "plot": None}

    plt.figure(figsize=(6, max(3, top_df.shape[0] * 0.4)))
    # horizontal dotplot: y = term, x = NES, size = -log10(p.adjust), color = NES sign
    top_df = top_df[::-1]  # so highest NES on top
    y = top_df["ID"].astype(str)
    x = top_df["NES"]
    size = -np.log10(top_df["p.adjust"].replace(0, np.nan)).fillna(0) * 20 + 20
    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    sc = plt.scatter(x=x, y=y, s=size, c=x, cmap=cmap, edgecolors="k")
    plt.colorbar(sc, label="NES")
    plt.xlabel("NES")
    plt.title(f"GSEA: {sample_label} ({mode} FC)")
    plt.tight_layout()

    out_figs_dir = os.path.join(out_dir, "figures")
    os.makedirs(out_figs_dir, exist_ok=True)
    png_file = os.path.join(out_figs_dir, f"{out_prefix}{suffix}_results.png")
    pdf_file = os.path.join(out_figs_dir, f"{out_prefix}{suffix}_results.pdf")
    plt.savefig(png_file, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_file, bbox_inches="tight")
    plt.close()
    print(f"[OK] GSEA dotplot saved at: {png_file} and {pdf_file}")

    return {"results": gsea_df, "plot": top_df}


# -------------------------
# Combine and save plots as grid
# -------------------------
def combine_and_save(plot_dict, out_prefix, suffix, out_dir="out"):
    """
    plot_dict: dict name -> {'results': df, 'plot': top_df}
    We'll create a grid of small dotplots using seaborn/matplotlib.
    """
    names = list(plot_dict.keys())
    plots = [plot_dict[n]["plot"] for n in names]

    n = len(plots)
    if n == 0:
        print("[WARN] No plots to combine.")
        return

    ncol = 3
    nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(6 * ncol, 0.75 * nrow * 6))
    axes = axes.flatten()

    for ax_idx, (name, plot_df) in enumerate(zip(names, plots)):
        ax = axes[ax_idx]
        ax.set_title(name)
        if plot_df is None or plot_df.shape[0] == 0:
            ax.text(0.5, 0.5, name, ha='center', va='center')
            ax.set_axis_off()
            continue
        df = plot_df
        y = df["ID"].astype(str)
        x = df["NES"]
        size = -np.log10(df["p.adjust"].replace(0, np.nan)).fillna(0) * 20 + 20
        sc = ax.scatter(x=x, y=y, s=size, c=x, cmap='coolwarm', edgecolors='k')
        ax.set_xlabel("NES")
        # tidy y labels
        ax.set_yticks(range(len(y)))
        ax.set_yticklabels(list(y))
    # hide unused axes
    for i in range(len(plots), len(axes)):
        axes[i].set_axis_off()

    fig.tight_layout()
    out_figs_dir = os.path.join(out_dir, "figures")
    os.makedirs(out_figs_dir, exist_ok=True)
    png_file = os.path.join(out_figs_dir, f"{out_prefix}_{suffix}.png")
    pdf_file = os.path.join(out_figs_dir, f"{out_prefix}_{suffix}.pdf")
    fig.savefig(png_file, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_file, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Combined {suffix} GSEA plot saved to: {png_file} and {pdf_file}")


# -------------------------
# Plot GSEA heatmap (NES)
# -------------------------
def plot_gsea_heatmap(gsea_all_df, col_order, out_file=None, title="GSEA NES Heatmap", cluster_cols=False):
    """
    gsea_all_df: DataFrame containing columns ['ID', 'sample', 'NES']
    col_order: desired sample column order (list of sample labels)
    """
    if gsea_all_df.shape[0] == 0:
        print("[WARN] Empty GSEA results supplied to heatmap function.")
        return None

    nes_matrix = (gsea_all_df
                  .groupby(["ID", "sample"])["NES"]
                  .mean()
                  .reset_index()
                  .pivot(index="ID", columns="sample", values="NES"))

    # ensure all requested cols present
    for c in col_order:
        if c not in nes_matrix.columns:
            nes_matrix[c] = np.nan

    nes_matrix = nes_matrix[col_order]
    nes_matrix = nes_matrix.fillna(0)

    # plot using seaborn clustermap or heatmap
    plt.figure(figsize=(max(6, 0.75 * len(col_order)), max(6, 0.25 * nes_matrix.shape[0])))
    sns.heatmap(nes_matrix, cmap="RdYlBu_r", center=0, linewidths=0.5)
    plt.title(title)
    plt.tight_layout()

    if out_file:
        ensure_dir(out_file)
        plt.savefig(out_file, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"[OK] Heatmap saved to {out_file}")
    else:
        plt.show()
        plt.close()

    return nes_matrix


# -------------------------
# CLI and main runner
# -------------------------
def parse_args():
    p = argparse.ArgumentParser(description="GSEA tools (Python port of gsea_tools.R)")
    p.add_argument("--infile", help="Input differential expression TSV/CSV (Gene, logFC, FDR etc.)")
    p.add_argument("--gene-col", default="Gene", help="Gene column name (default: Gene)")
    p.add_argument("--fc-col", default="logFC", help="Fold-change (numeric) column (default: logFC)")
    p.add_argument("--fdr-col", default="FDR", help="FDR / padj column name (default: FDR)")
    p.add_argument("--sample-label", default="sample", help="Label for this sample (used for output names)")
    p.add_argument("--out-prefix", default="gsea", help="Prefix for output files")
    p.add_argument("--mode", choices=["all", "pos", "neg"], default="all", help="Run mode for prerank")
    p.add_argument("--run-volcano", action="store_true", help="Produce volcano plot")
    p.add_argument("--run-gsea", action="store_true", help="Run preranked GSEA for the provided file")
    p.add_argument("--min-gs-size", type=int, default=10)
    p.add_argument("--max-gs-size", type=int, default=500)
    p.add_argument("--top-cats", type=int, default=10)
    p.add_argument("--out-dir", default="out", help="Base output directory")
    p.add_argument("--batch-file", default=None, help="If provided, a TSV with columns: infile, sample_label, out_prefix for batch processing")
    return p.parse_args()


def main():
    args = parse_args()

    # If batch mode provided, override single-run flags
    if args.batch_file:
        batch_df = pd.read_csv(args.batch_file, sep=None, engine="python")
        plot_dict = {}
        gsea_all_results = []
        for idx, row in batch_df.iterrows():
            infile = row.get("infile") or row.get("file") or row.get("path")
            sample_label = row.get("sample_label") or row.get("sample") or f"sample{idx}"
            out_prefix = row.get("out_prefix") or sample_label
            print(f"[BATCH] Processing {infile} as {sample_label}")
            df = pd.read_csv(infile, sep=None, engine="python")
            # build gene -> strongest logFC
            tmp = df[[args.gene_col, args.fc_col]].copy()
            tmp[args.fc_col] = pd.to_numeric(tmp[args.fc_col], errors="coerce")
            tmp = tmp.dropna(subset=[args.gene_col]).drop_duplicates(args.gene_col)
            gene_ser = tmp.groupby(args.gene_col)[args.fc_col].apply(lambda s: s.iloc[np.argmax(np.abs(s.values))])
            # run gsea (all/pos/neg) modes and save each
            term2gene = prepare_term2gene_reactome()
            res = run_gsea_prerank(gene_ser, term2gene,
                                   minGSSize=args.min_gs_size, maxGSSize=args.max_gs_size,
                                   out_prefix=out_prefix, sample_label=sample_label, mode=args.mode,
                                   top_categories=args.top_cats, out_dir=args.out_dir)
            plot_dict[sample_label] = res
            if res["results"].shape[0] > 0:
                gsea_all_results.append(res["results"][["ID", "NES", "sample"]])
        # combine plots
        combine_and_save(plot_dict, args.out_prefix, "gsea_combined", out_dir=args.out_dir)
        if gsea_all_results:
            all_gsea_df = pd.concat(gsea_all_results, ignore_index=True)
            heat_out = os.path.join(args.out_dir, "figures", f"{args.out_prefix}_gsea_heatmap.png")
            col_order = sorted(all_gsea_df["sample"].unique().tolist())
            plot_gsea_heatmap(all_gsea_df, col_order, out_file=heat_out)
        sys.exit(0)

    # single-run mode
    if not args.infile:
        print("Either --infile or --batch-file is required. See --help.")
        sys.exit(1)

    df = pd.read_csv(args.infile, sep=None, engine="python")
    # drop contaminant rows if first column begins with "contaminant" (mimic earlier behavior)
    first_col = df.columns[0]
    df = df[~df[first_col].astype(str).str.lower().str.startswith("contaminant")]

    # volcano
    if args.run_volcano:
        volcano_out = os.path.join(args.out_dir, "figures", f"{args.out_prefix}_volcano.png")
        volcano_plot(df, args.fc_col, args.fdr_col, args.gene_col,
                     out_file=volcano_out, plot_title=f"Volcano: {args.sample_label}")

    # prepare gene vector for prerank: choose strongest logFC per gene (by absolute value)
    df_temp = df[[args.gene_col, args.fc_col]].copy()
    df_temp[args.fc_col] = pd.to_numeric(df_temp[args.fc_col], errors="coerce")
    df_temp = df_temp.dropna(subset=[args.gene_col])
    # keep strongest per gene
    gene_ser = df_temp.groupby(args.gene_col)[args.fc_col].apply(lambda s: s.iloc[np.argmax(np.abs(s.values))])
    gene_ser.index = gene_ser.index.astype(str)

    if args.run_gsea:
        term2gene = prepare_term2gene_reactome()
        res = run_gsea_prerank(gene_ser, term2gene,
                               minGSSize=args.min_gs_size, maxGSSize=args.max_gs_size,
                               out_prefix=args.out_prefix, sample_label=args.sample_label, mode=args.mode,
                               top_categories=args.top_cats, out_dir=args.out_dir)
        # res["results"] contains the table; res["plot"] contains the top_df used for dotplot

        # If you want to create combined heatmap across multiple runs, you should run multiple sample files and then call plot_gsea_heatmap on the merged results.

    print("Done.")


if __name__ == "__main__":
    main()
