#!/usr/bin/env python3
"""
Merge MMseqs2 taxonomy results with contig abundance data.
"""
import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(description='Merge MMseqs2 results with abundance')
    parser.add_argument('--mmseqs2', required=True, help='MMseqs2 output TSV file')
    parser.add_argument('--abundance', required=True, help='Abundance TSV file')
    parser.add_argument('--output', required=True, help='Output TSV file')
    args = parser.parse_args()

    # Read mmseqs2 results (all contigs with taxon IDs)
    # MMseqs2 output: contig_id, taxon_id (first two columns)
    mmseqs2 = pd.read_csv(args.mmseqs2, sep="\t", header=None, usecols=[0, 1],
                          names=["contig", "taxon_id"])

    # Read abundance data
    abundance = pd.read_csv(args.abundance, sep="\t", header=None, 
                           names=["contig", "abundance"])

    # Merge by contig
    merged = pd.merge(mmseqs2, abundance, on="contig", how="inner")

    # Collapse by taxon_id and sum abundance
    taxon_abundance = merged.groupby("taxon_id")["abundance"].sum().reset_index()
    taxon_abundance.columns = ["Taxon", "Abundance"]

    # Save results
    taxon_abundance.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
