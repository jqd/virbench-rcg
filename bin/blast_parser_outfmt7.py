#!/usr/bin/env python
import sys, os, re, argparse
import pandas as pd
from collections import Counter
import json
import time
import subprocess
import tempfile

def generate_taxon_lineage_with_taxonkit(taxids, ncbi_taxonomy_dir):
    """
    Use TaxonKit to generate full lineage for a list of taxon IDs.
    
    Parameters:
        taxids (pd.Series or list): List/Series of taxon IDs
        ncbi_taxonomy_dir (str): Path to NCBI taxonomy database directory
        
    Returns:
        pd.DataFrame: DataFrame with columns [Taxid, Lineage, Kingdom, Phylum, Class, Order, Family, Genus, Species]
    """
    # Get unique taxon IDs
    unique_taxids = set(taxids.dropna().unique()) if hasattr(taxids, 'dropna') else set(taxids)
    
    # Write taxon IDs to a temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_file:
        tmp_taxid_file = tmp_file.name
        for taxid in unique_taxids:
            tmp_file.write(f"{int(taxid)}\n")
    
    try:
        # Run taxonkit reformat to get lineage
        cmd = f"cat {tmp_taxid_file} | ~/vir_pipelines/redo_2025/taxonkit_results/taxonkit reformat2 -I 1 --data-dir {ncbi_taxonomy_dir} -r 'Unclassified'"
        
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        
        if result.returncode != 0:
            print(f"TaxonKit error: {result.stderr}")
            raise RuntimeError(f"TaxonKit failed with return code {result.returncode}")
        
        # Parse into DataFrame
        lineage_data = []
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    taxid = int(parts[0])
                    lineage = parts[1]
                    lineage_data.append({'Taxid': taxid, 'Lineage': lineage})
        
        taxon_lineage = pd.DataFrame(lineage_data)
        
        # Split lineage into rank columns
        ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        taxon_lineage[ranks] = taxon_lineage.Lineage.str.split(pat=';', expand=True)
        
        # Save taxon_lineage to file
        taxon_lineage.to_csv("taxon_lineage_generated.tsv", sep='\t', index=False)
        print(f"Saved taxon lineage to taxon_lineage_generated.tsv")
        
        print(f"Generated lineage for {len(taxon_lineage)} taxon IDs")
        return taxon_lineage
        
    finally:
        # Clean up temporary file
        if os.path.exists(tmp_taxid_file):
            os.remove(tmp_taxid_file)


def get_species_and_taxid_lineage(taxon_lineage, taxid, missing_taxid_file="missing_taxids.txt"):
    """
    Retrieve the Species name and a concatenated string of Taxid and Lineage from the taxon_lineage DataFrame based on the input Taxid.
    
    Parameters:
        taxon_lineage (pd.DataFrame): DataFrame containing species and taxon data.
        taxid (int): TaxID to lookup.
        missing_taxid_file (str): File path to log missing taxids.
        
    Returns:
        tuple: Species name and concatenated string of Taxid and Lineage in the format (Species, 'Taxid;Lineage').
    """
    # Filter the DataFrame by the given Taxid
    taxon_row = taxon_lineage[taxon_lineage.Taxid == taxid]

    if taxon_row.empty:
        with open(missing_taxid_file, 'a') as f:
            f.write(f"Warning: Taxid {taxid} not found in taxon_lineage\n")
            f.write(f"Empty row: {taxon_row}\n\n")
        return None, None  # Return tuple of None if the Taxid is not found

    # Extract Species, Taxid, and Lineage
    species = taxon_row['Species'].iloc[0]
    taxon_id = taxon_row['Taxid'].iloc[0]
    lineage = taxon_row['Lineage'].iloc[0]
    
    return species, f"{taxon_id};{lineage}"



def get_most_common_taxon(df, taxon_column, include_kingdoms=None, exclude_kingdoms=None):
    """
    Retrieve the most common taxon (Taxid) from the specified kingdoms or exclude certain kingdoms, 
    using the specified taxon column.

    Parameters:
        df (pd.DataFrame): DataFrame containing species and taxon data.
        taxon_column (str): The name of the column to use for taxon values (e.g., 'Staxids_new').
        include_kingdoms (list of str, optional): The kingdoms to include in the calculation. 
                                                 If None, includes all kingdoms.
        exclude_kingdoms (list of str, optional): The kingdoms to exclude from the calculation. 
                                                 If None, no kingdoms are excluded.
        
    Returns:
        str: The most common taxon (Taxid).
    """
    # Filter by include_kingdoms if specified
    if include_kingdoms:
        kingdom_data = df[df.Kingdom.isin(include_kingdoms)]
    elif exclude_kingdoms:
        # Filter out the kingdoms to exclude
        kingdom_data = df[~df.Kingdom.isin(exclude_kingdoms)]
    else:
        # If no include or exclude list is provided, use the entire DataFrame
        kingdom_data = df
    
    # Use mode to get the most common Taxid (from the specified taxon column)
    # If there is a tie, break it by selecting the one with the highest Bitscore, or by the index
    #most_common_taxon = kingdom_data[taxon_column].mode()[0]

    bitscore_column = 'bitscore'
    # Step 1: Count frequencies of each Taxid (or other taxon)
    taxon_counts = kingdom_data[taxon_column].value_counts()

    # Step 2: Sum the Bitscores for each Taxid (or taxon)
    sum_bitscores = kingdom_data.groupby(taxon_column)[bitscore_column].sum()

    # Step 3: Identify the most common Taxid(s) based on frequency
    most_common_taxons = taxon_counts[taxon_counts == taxon_counts.max()].index

    # Step 4: If there is a tie, break it by selecting the Taxid with the highest summed Bitscore
    if len(most_common_taxons) > 1:
        # Sum the Bitscores for each tied Taxid and pick the one with the highest summed Bitscore
        summed_bitscore = sum_bitscores[most_common_taxons]
        best_taxon = summed_bitscore.idxmax()  # Select the Taxid with the highest summed bitscore
        #best_bitscore = summed_bitscore.max()  # Get the highest summed bitscore
        return best_taxon
    else:
        # No tie, just pick the most common Taxid
        best_taxon = most_common_taxons[0]
        return best_taxon


def append_to_tsv(df, output_tsv_filename):
    """
    Append data to a TSV file. If the file doesn't exist, it creates a new one.
    
    Parameters:
    df (pandas.DataFrame): The DataFrame containing data to append.
    output_tsv_filename (str): The path to the output TSV file.
    """
    try:
        # Check if the file exists. If it doesn't, write the header.
        file_exists = os.path.exists(output_tsv_filename)
        
        # Append to the file, writing the header only if the file doesn't already exist
        df.to_csv(output_tsv_filename, sep='\t', index=False, header=not file_exists, mode='a')
        # print(f"Data successfully appended to {output_tsv_filename}")
    
    except Exception as e:
        print(f"Error appending data to {output_tsv_filename}: {e}")

def convert_tsv_to_json(tsv_filename, json_filename):
    """
    Convert a TSV file to a JSON file.

    Parameters:
    tsv_filename (str): The path to the input TSV file.
    json_filename (str): The path to the output JSON file.
    """
    try:
        # Load the TSV data into a pandas DataFrame
        df = pd.read_csv(tsv_filename, sep='\t')
        
        # Convert the DataFrame to JSON
        df.to_json(json_filename, orient='records', lines=False, indent=4)
        # print(f"Data successfully converted to {json_filename}")
    
    except Exception as e:
        print(f"Error while converting TSV to JSON: {e}")

def assign_tax(contig, df, taxon_lineage):
    #to-do list
    #1: uncultured bacterium: no taxnonomy, may use dada2 to assign taxonomy

    #df.to_csv(str(contig) + '.tsv', sep="\t")
    #The first hit among all hits
    species_top_hit = 'NONE'
    species_top_hit_lineage = 'NONE'

    #The first viral hit among all virus hits
    top_virus_hit = ''
    top_virus_species = "NONE"
    top_virus_species_lineage = 'NONE'

    #The species that occurs most frequently across all hits
    most_common_species = 'NONE'
    most_common_species_lineage = 'NONE'

    #The most frequently occurring species among viral hits
    common_species_viruses = "NONE"
    common_species_viruses_lineage = 'NONE'

    #The most frequently occurring species among non-viral hits
    common_species_non_viruses = "NONE"
    common_species_non_viruses_lineage = 'NONE'

    kingdom_tag = ''
    phage_tag = 'Non_phage' #Defalut
    phages_class = ['Caudoviricetes', 'Malgrandaviricetes', 'Leviviricetes'] #Common phage class, also check "Cystoviridae" family


    df_merged = df.merge(taxon_lineage, left_on='Staxids_new', right_on="Taxid", how='left').sort_values(by='original_index')

    #Filter out N/A or missing kingdoms
    df_merged = df_merged[df_merged.Kingdom != 'Unclassified']

    # File path for the CSV
    file_path = 'output_file.tsv'

    # Check if the file exists
    if os.path.exists(file_path):
        # Append to CSV (without writing the header)
        df_merged.to_csv(file_path, mode='a', header=False, index=False, sep = '\t')
        #print(f"DataFrame successfully appended to {file_path}")
    else:
        # If the file doesn't exist, create it and write the header
        df_merged.to_csv(file_path, mode='w', header=True, index=False, sep = '\t')
        #print(f"DataFrame successfully saved to {file_path}")

    ############Get top hit
    top_hit = df_merged.iloc[0]
    top_hit_taxon = top_hit.Staxids_new
    if top_hit_taxon:
        species_top_hit, species_top_hit_lineage = get_species_and_taxid_lineage(taxon_lineage, top_hit_taxon)

    ############Get the most common species (use the taxon)
    most_common_species_taxon = get_most_common_taxon(df_merged, taxon_column="Staxids_new")
    most_common_species, most_common_species_lineage = get_species_and_taxid_lineage(taxon_lineage, most_common_species_taxon)

    ##########Get the most common kingdom
    kingdoms_counter = Counter(df_merged.Kingdom)
    kingdom_str = ";".join([f"{kingdom}:{count}" for kingdom, count in kingdoms_counter.items()])

    if len(kingdoms_counter) == 1:
        kingdom_tag = kingdoms_counter.most_common(1)[0][0]
        if "Viruses" in kingdoms_counter:
            common_species_viruses = most_common_species
            common_species_viruses_lineage = most_common_species_lineage
            top_virus_species = species_top_hit
            top_virus_species_lineage = species_top_hit_lineage
        else:
            kingdom_tag = 'Non_viral|' + kingdom_str
            common_species_non_viruses = most_common_species
            common_species_non_viruses_lineage = most_common_species_lineage
    elif len(kingdoms_counter) > 1:
        if "Viruses" in kingdoms_counter:
            kingdom_tag = 'mix_with_viruses|' + kingdom_str 

            common_species_viruses_taxon = get_most_common_taxon(df_merged, taxon_column="Staxids_new", include_kingdoms=["Viruses"])
            common_species_viruses, common_species_viruses_lineage = get_species_and_taxid_lineage(taxon_lineage, common_species_viruses_taxon)

            common_species_non_viruses_taxon = get_most_common_taxon(df_merged, taxon_column="Staxids_new", exclude_kingdoms=["Viruses"])
            common_species_non_viruses, common_species_non_viruses_lineage = get_species_and_taxid_lineage(taxon_lineage, common_species_non_viruses_taxon)
     
            virus_hits = df_merged[df_merged.Kingdom == 'Viruses'].sort_values(by='original_index')
            top_virus_hit = virus_hits.iloc[0]
            top_virus_taxon = top_virus_hit.Staxids_new
            top_virus_species, top_virus_species_lineage = get_species_and_taxid_lineage(taxon_lineage, top_virus_taxon)
        else:           
            kingdom_tag = 'mix_without_viruses|' + kingdom_str
            common_species_non_viruses = most_common_species
            common_species_non_viruses_lineage = most_common_species_lineage

    ################phage tag
    if df_merged['Class'].isin(phages_class).any() or 'Cystoviridae' in df_merged['Family'].values:
        phage_tag = 'phage'

    #Creating the list to return
    output_tsv = [
        contig,
        kingdom_tag,
        kingdom_str,
        species_top_hit,
        species_top_hit_lineage,
        most_common_species,
        most_common_species_lineage,
        top_virus_species,
        top_virus_species_lineage,
        common_species_viruses,
        common_species_viruses_lineage,
        common_species_non_viruses,
        common_species_non_viruses_lineage,
        phage_tag
    ]
    # Define the column names corresponding to the output list
    columns = [
        "contig",
        "kingdom_tag",
        "kingdom_str",
        "species_top_hit",
        "species_top_hit_lineage",
        "most_common_species",
        "most_common_species_lineage",
        "top_virus_species",
        "top_virus_species_lineage",
        "common_species_viruses",
        "common_species_viruses_lineage",
        "common_species_non_viruses",
        "common_species_non_viruses_lineage",
        "phage_tag"
    ]

    # Convert the output list to a pandas DataFrame
    df_to_append = pd.DataFrame([output_tsv], columns=columns)

    return(df_to_append)

def main():
    parser = argparse.ArgumentParser()

    # Input file (required) - short: -i, long: --input
    parser.add_argument('-i', '--input', required=True,
                        help='blast format in outfmt 7')

    # E-value cutoff (optional) - short: -e, long: --evalue
    parser.add_argument('-e', '--evalue', type=float, required=False, default=1e-5,
                        help='evalue cutoff (default: 1e-5)')

    # Number of distinct bit scores to consider (optional) - short: -n, long: --number_bitscores
    parser.add_argument('-n', '--number_bitscores', type=int, required=False, default=100,
                        help='Number of distinct bit scores to consider. Hits corresponding to any of these bit scores will be kept (default: 100)')

    # Taxonomy data (optional) - short: -t, long: --taxonomy
    parser.add_argument('-t', '--taxonomy', required=False,
                        help='Pre-computed taxonomy lineage file')

    # NCBI taxonomy database directory for taxonkit (optional)
    parser.add_argument('-d', '--ncbi_taxonomy_dir', required=False,
                        help='Path to NCBI taxonomy database directory for taxonkit')

    # Output filename (optional) - short: --output, long: --output_prefix
    parser.add_argument('-o', '--output_prefix', type=str, default='results',
                        help="Prefix for the output TSV and JSON files.")


    # get the start time
    st = time.time()
    # Accessing the parameters
    args=parser.parse_args()
    evalue =  args.evalue  
    n_bitscores = args.number_bitscores
    output_prefix = args.output_prefix

    # Define the output filenames
    output_tsv_filename = f"{output_prefix}.tsv"
    output_json_filename = f"{output_prefix}.json"

    # Open blast output
    headers = ['query',
            'subject',
            'identity',
            'alignment_length',
            'mismatches',
            'gapopens',
            'Qstart',
            'Qend',
            'Sstart',
            'Send',
            'evalue',
            'bitscore',
            'Skingdom',
            'Staxids',
            'Sscinames',
            'Stitle'
    ]
    df = pd.read_csv(args.input, 
                    sep='\t', 
                    comment='#', 
                    names = headers,
                    dtype={'evalue': float,'bitscore': float, 'Staxids': str})
    
    # Add a column to keep track of the original index
    df['original_index'] = df.index

    # remove hits with unknown kingdom, such as synthetic construct(taxid:32630) / vector
    df = df[df.Skingdom.notnull()] 

    # Process 'Staxids': remove leading zeros and keep the first non-zero taxon ID, for example:0;562;83333;1110693
    df['Staxids_new'] = df['Staxids'].str.lstrip('0;').str.split(';').str[0].astype(int)
    
    # Keep only rows with acceptable e-values
    df = df[df['evalue'] <= evalue]

    # Get taxon lineage - either from pre-computed file or generate with taxonkit
    ranks = [
        "Kingdom", 
        "Phylum", 
        "Class", 
        "Order", 
        "Family", 
        "Genus", 
        "Species"
    ]
    
    if args.taxonomy:
        # Use pre-computed taxonomy file
        taxon_lineage = pd.read_csv(args.taxonomy, sep='\t', usecols=[0,1], names=["Taxid", "Lineage"])
        taxon_lineage[ranks] = taxon_lineage.Lineage.str.split(pat=';', expand=True)
    elif args.ncbi_taxonomy_dir:
        # Generate lineage using taxonkit
        taxon_lineage = generate_taxon_lineage_with_taxonkit(df['Staxids_new'], args.ncbi_taxonomy_dir)
        print(f"Sample row from taxon_lineage:\n{taxon_lineage.head(1)}")
    else:
        print("Error: Either --taxonomy or --ncbi_taxonomy_dir must be provided")
        sys.exit(1)


    # Group by the 'query' column for further analysis
    df_grouped = df.groupby('query', sort=False)

	#Iterate over each contig group in the DataFrame
    for contig, group in df_grouped:
        # Sort the group by 'bitscore' in descending order
        sorted_group = group.sort_values(by='bitscore', ascending=False, kind='stable')
        
        # Select the top 'n_bitscores' based on bitscore
        selected_scores = sorted_group['bitscore'].unique()[:n_bitscores]
        selected_rows = sorted_group[sorted_group['bitscore'].isin(selected_scores)]

        num_columns = 14
        # Check if there are selected rows
        if selected_rows.empty:
            # Create a list with empty strings for the annotations, and add the contig as the first column
            anno = [contig] + [""] * (num_columns - 1)
        else:
            # Assign taxonomy based on the selected rows
            anno = assign_tax(contig, selected_rows, taxon_lineage)
        # Append the annotation data to the JSON file
        #append_to_json_file(anno_jason, output_json_filename)
        append_to_tsv(anno, output_tsv_filename)
        convert_tsv_to_json(output_tsv_filename, output_json_filename)


if __name__ == "__main__":
    main()
