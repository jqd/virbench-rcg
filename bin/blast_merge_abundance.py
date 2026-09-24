#!/usr/bin/env python3
"""
Merge BLAST taxonomy results with contig abundance and CD-HIT cluster data.
"""
import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(description='Merge BLAST results with abundance and cluster data')
    parser.add_argument('--blast', required=True, help='Processed BLAST TSV file')
    parser.add_argument('--abundance', required=True, help='Abundance TSV file')
    parser.add_argument('--cdhit', required=True, help='CD-HIT cluster TSV file')
    parser.add_argument('--output', required=True, help='Output TSV file')
    args = parser.parse_args()

    # Read BLAST results with expected columns
    blast = pd.read_csv(args.blast, sep="\t")
    
    # Read abundance data (contig, abundance)
    abundance = pd.read_csv(args.abundance, sep="\t", header=None, 
                           names=["contig", "abundance"])

    # Read CD-HIT cluster data (only columns 1, 2, 5: id, clstr, clstr_rep)
    cdhit = pd.read_csv(args.cdhit, sep="\t", usecols=["id", "clstr", "clstr_rep"])
    cdhit.rename(columns={"id": "contig"}, inplace=True)

    # Merge BLAST with abundance on contig
    merged = pd.merge(blast, abundance, on="contig", how="left")
    
    # Merge with CD-HIT clusters
    merged = pd.merge(merged, cdhit, on="contig", how="left")
    
    # Reorder columns to put contig, abundance, clstr, clstr_rep first, then all BLAST columns
    first_cols = ["contig", "abundance", "clstr", "clstr_rep"]
    blast_cols = [col for col in blast.columns if col != "contig"]
    output_cols = first_cols + blast_cols
    # Only keep columns that exist
    output_cols = [col for col in output_cols if col in merged.columns]
    merged = merged[output_cols]

    # Save results
    merged.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
