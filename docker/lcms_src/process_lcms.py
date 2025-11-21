import pandas as pd
import os
import argparse


def process_expression_sheet(input_file, outdir, prefix):
    df = pd.read_excel(input_file, sheet_name="Proteins", header=0)

    # Remove rows starting with "contaminant"
    first_col = df.columns[0]
    df_clean = df[~df[first_col].str.contains("^contaminant", case=False, na=False)]

    # Determine columns using first data row values
    first_row = df_clean.iloc[0].astype(str)
    cols_to_keep = [
        i for i, val in enumerate(first_row)
        if val == "Genes" or val.endswith("_normalized_intensity")
    ]

    df_subset = df_clean.iloc[:, cols_to_keep]

    # Remove the row used for selecting columns
    df_subset = df_subset.iloc[1:].reset_index(drop=True)

    # Rename first column
    df_subset = df_subset.rename(columns={df_subset.columns[0]: "Gene Name"})

    # Remove SEV_HC_EDTA if present
    if "SEV_HC_EDTA" in df_subset.columns:
        df_subset = df_subset.drop(columns=["SEV_HC_EDTA"])

    os.makedirs(outdir, exist_ok=True)

    outfile = os.path.join(outdir, f"{prefix}_expression.tsv")
    df_subset.to_csv(outfile, sep="\t", index=False)

    print(f"[OK] Saved expression data to: {outfile}")


def process_diff_expression_sheet(input_file, outdir, prefix):
    df = pd.read_excel(input_file, sheet_name="DE Analysis", header=0)

    # Remove contaminants
    first_col = df.columns[0]
    df_clean = df[~df[first_col].str.contains("^contaminant", case=False, na=False)]

    # Keep Gene Name + fold-change and p.adj columns
    pattern = r"fold change$|_p\.adj$"
    cols_to_keep = ["Gene Name"] + [
        col for col in df_clean.columns
        if pd.Series(col).str.contains(pattern, regex=True).iloc[0]
    ]

    df_subset = df_clean[cols_to_keep]

    os.makedirs(outdir, exist_ok=True)

    outfile = os.path.join(outdir, f"{prefix}_diff_expression.tsv")
    df_subset.to_csv(outfile, sep="\t", index=False)

    print(f"[OK] Saved differential expression data to: {outfile}")


def main():
    parser = argparse.ArgumentParser(
        description="Process LC-MS Proteomics Excel file into clean expression tables."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Path to input XLSX file."
    )

    parser.add_argument(
        "--outdir",
        default="out/tables",
        help="Output directory (default: out/tables)"
    )

    parser.add_argument(
        "--prefix",
        default="output",
        help="Prefix for output files (default: output)"
    )

    args = parser.parse_args()

    print("Processing LC-MS Excel file...")
    process_expression_sheet(args.input, args.outdir, args.prefix)
    process_diff_expression_sheet(args.input, args.outdir, args.prefix)
    print("Done.")


if __name__ == "__main__":
    main()
