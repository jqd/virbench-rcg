#!/usr/bin/env python3
"""
Merge Kaiju taxonomy results with contig abundance data.
"""
import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(description='Merge Kaiju results with abundance')
    parser.add_argument('--kaiju', required=True, help='Kaiju output TSV file')
    parser.add_argument('--abundance', required=True, help='Abundance TSV file')
    parser.add_argument('--output', required=True, help='Output TSV file')
    args = parser.parse_args()

    # Read kaiju results (only first 3 columns)
    kaiju = pd.read_csv(args.kaiju, sep="\t", header=None, 
                       names=["class", "contig", "taxon_id"], usecols=[0, 1, 2])
    kaiju_results = kaiju[["contig", "taxon_id"]]

    # Read abundance data
    abundance = pd.read_csv(args.abundance, sep="\t", header=None, 
                           names=["contig", "abundance"])

    # Merge by contig
    merged = pd.merge(kaiju_results, abundance, on="contig", how="inner")

    # Collapse by taxon_id and sum abundance
    taxon_abundance = merged.groupby("taxon_id")["abundance"].sum().reset_index()
    taxon_abundance.columns = ["Taxon", "Abundance"]

    # Save results
    taxon_abundance.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
