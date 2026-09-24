#!/usr/bin/env Rscript
###R 4.4.0
# Script for calculating metrics in Nextflow pipeline
# Usage:
#   Rscript calculate_metrics.R \
#     --blast <file> \
#     --kaiju <file> \
#     --mmseqs2 <file> \
#     --kraken2 <file> \
#     --fve <file> \
#     --viromescan <file> \
#     --ground_truth <file> \
#     --ncbi_taxonomy_dir <dir> \
#     --tool <assembler> \
#     --output <metrics.tsv> \
#     --FVE_default_db <true|false> \
#     --FVE_acc2taxon <file> \
#     --viromescan_default_db <true|false> \
#     --viromescan_acc2taxon <file>

library(tidyverse)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse named arguments
parse_args <- function(args) {
  arg_list <- list()
  i <- 1
  while (i <= length(args)) {
    if (startsWith(args[i], "--")) {
      key <- sub("^--", "", args[i])
      if (i < length(args) && !startsWith(args[i + 1], "--")) {
        arg_list[[key]] <- args[i + 1]
        i <- i + 2
      } else {
        stop(paste("Missing value for argument:", args[i]))
      }
    } else {
      i <- i + 1
    }
  }
  return(arg_list)
}

params <- parse_args(args)

# Validate required arguments
required <- c("blast", "kaiju", "mmseqs2", "kraken2", "fve", "viromescan", 
              "ground_truth", "ncbi_taxonomy_dir", "tool", "output",
              "FVE_default_db", "FVE_acc2taxon",
              "viromescan_default_db", "viromescan_acc2taxon")
missing <- setdiff(required, names(params))
if (length(missing) > 0) {
  stop(paste("Missing required arguments:", paste(missing, collapse = ", ")))
}

cat("Configuration:\n")
cat("  BLAST file:", params$blast, "\n")
cat("  Kaiju file:", params$kaiju, "\n")
cat("  MMseqs2 file:", params$mmseqs2, "\n")
cat("  Kraken2 file:", params$kraken2, "\n")
cat("  FVE file:", params$fve, "\n")
cat("  FVE default db:", params$FVE_default_db, "\n")
cat("  FVE acc2taxon:", params$FVE_acc2taxon, "\n")
cat("  ViromeScan file:", params$viromescan, "\n")
cat("  Ground truth:", params$ground_truth, "\n")
cat("  NCBI taxonomy dir:", params$ncbi_taxonomy_dir, "\n")
cat("  Assembler tool:", params$tool, "\n")
cat("  Output prefix:", params$output, "\n\n")

# Function definitions
L2 <- function(x, y){
  d <- x - y
  sqrt(sum(d^2, na.rm = TRUE))
}

cal_ROC <- function(target, pred, tool){
  target <- target
  pred <- pred
  tool <- tool
  TP <- intersect(target, pred) %>% nrow()
  FN <- setdiff(target, pred) %>% nrow()
  FP <- setdiff(pred, target) %>% nrow()
  TN <- 0
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1score <- 2 * (precision * recall) / (precision + recall)
  return(list("Tool" = tool,
              "TP" = TP,
              "FP" = FP,
              "FN" = FN,
              "Precision" = precision, 
              "Recall" = recall, 
              "F1score" = F1score))
}

######Truth value
data_real <- read_tsv(params$ground_truth, col_names = TRUE, show_col_types = FALSE)
data_real <- data_real %>% select(Taxon, Abundance)

# Define parameters
assembler <- params$tool  # "megahit" or "metaspades"
rank <- 'Species'
target <- "Real"
tools <- c(target,
           "blast_Tophit_species",
           "blast_Common_Species",
           "blast_top_virus_species",
           "Kraken2",
           "Kaiju",
           "MMseqs2 taxonomy",
           "FastViromeExplorer",
           "ViromeScan")

cat("\n=== Processing", assembler, "===\n")

###blast
# Read blast output: first column is abundance, second is semicolon-delimited lineage (taxon id first)
blast_raw <- read_tsv(params$blast, col_names = TRUE, show_col_types = FALSE)

# Helper function to extract taxon abundance from a lineage column
extract_taxon_abun <- function(df, lineage_col) {
  df %>%
    filter(.data[[lineage_col]] != "NONE") %>%
    mutate(Taxon = as.numeric(str_split_fixed(.data[[lineage_col]], ";", 2)[,1])) %>%
    group_by(Taxon) %>%
    summarise(Abundance = sum(abundance), .groups = "drop") %>%
    filter(!is.na(Taxon))
}

blast_Tophit <- extract_taxon_abun(blast_raw, "species_top_hit_lineage")
blast_Common_Species <- extract_taxon_abun(blast_raw, "most_common_species_lineage")
blast_top_virus_species <- extract_taxon_abun(blast_raw, "top_virus_species_lineage")

# Read other tool results
kaiju_data <- read_tsv(params$kaiju, col_names = TRUE, show_col_types = FALSE)
mmseqs2_data <- read_tsv(params$mmseqs2, col_names = TRUE, show_col_types = FALSE)
kraken2_data <- read_tsv(params$kraken2, col_names = TRUE, show_col_types = FALSE)
fve_data <- read_tsv(params$fve, col_names = TRUE, show_col_types = FALSE)
head(fve_data)
# Only convert if using default database
if (tolower(params$FVE_default_db) == "true") {
  acc2tax <- read_tsv(params$FVE_acc2taxon,
                      col_names = c("Accession", "TaxonName"),
                      show_col_types = FALSE)

  fve_data <- fve_data %>%
    mutate(Taxon = as.character(Taxon)) %>%
    left_join(acc2tax, by = c("Taxon" = "Accession")) %>%
    select(TaxonName, Abundance) %>%
    rename(Taxon = TaxonName)
}

viromescan_data <- read_tsv(params$viromescan, col_names = TRUE, show_col_types = FALSE)
# Only convert if using default database (accession → taxon ID)
if (tolower(params$viromescan_default_db) == "true") {
  viromescan_acc2tax <- read_tsv(params$viromescan_acc2taxon,
                                col_names = c("Accession", "TaxonName"),
                                show_col_types = FALSE)

  viromescan_data <- viromescan_data %>%
    mutate(Taxon = as.character(Taxon)) %>%
    left_join(viromescan_acc2tax, by = c("Taxon" = "Accession")) %>%
    select(TaxonName, Abundance) %>%
    rename(Taxon = TaxonName)
}
# Combine all data sources
dfs <- list(
  "Real" = data_real,
  "blast_Tophit_species" = blast_Tophit,
  "blast_Common_Species" = blast_Common_Species,
  "blast_top_virus_species" = blast_top_virus_species,
  "Kaiju" = kaiju_data,
  "MMseqs2 taxonomy" = mmseqs2_data,
  "Kraken2" = kraken2_data,
  "FastViromeExplorer" = fve_data,
  "ViromeScan" = viromescan_data
)

combined <- dplyr::bind_rows(dfs, .id = "Tool")

###Taxon lineage - generate using taxonkit
# Extract unique taxon IDs
unique_taxons <- combined %>% 
  select(Taxon) %>% 
  distinct() %>% 
  filter(!is.na(Taxon))

# Write taxon IDs to temporary file
temp_taxon_file <- paste0("temp_taxons_", assembler, ".txt")
write_tsv(unique_taxons, temp_taxon_file, col_names = FALSE)

# Call taxonkit to get lineage
taxonkit_cmd <- sprintf(
  "cat %s | ~/vir_pipelines/redo_2025/taxonkit_results/taxonkit reformat2 -I 1 --data-dir %s -r 'Unclassified' > combined_taxons_lineage_%s.tsv",
  temp_taxon_file, params$ncbi_taxonomy_dir, assembler
)

cat("Running taxonkit for", assembler, "...\n")
system(taxonkit_cmd)

# Read the generated lineage file
taxon_lineage <- read_tsv(paste0("combined_taxons_lineage_", assembler, ".tsv"), 
                          col_names = c("Taxon", "Lineage"),
                          show_col_types = FALSE)

# Split lineage into rank columns
ranks <- c("Acellular_root", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxon_lineage <- taxon_lineage %>%
  separate(Lineage, into = ranks, sep = ";", fill = "right")

# Clean up temporary file
file.remove(temp_taxon_file)

combined_with_lineage <- combined %>% left_join(taxon_lineage, by = "Taxon")

cat("Done taxon lineage processing.\n")
# Species-level analysis
species_results <- combined_with_lineage %>% 
  select(Tool, Abundance, all_of(rank)) %>%
  group_by(Tool, .data[[rank]]) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  as.data.frame()


write_tsv(species_results, paste0("species_results_", assembler, ".tsv"))

cat("Done species-level aggregation.\n")


##Keep species with abundance >= 1000, convert to relative abundance
species_results_filtered <- species_results %>%
  group_by(Tool) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  filter(Abundance >= 1000) %>%
  select(Tool, all_of(rank), RelativeAbundance) %>%
  rename(Abundance = RelativeAbundance) %>%
  filter(Species != "Unclassified")

cat("Done filtering species and converting to relative abundance.\n")
species_abun <- split(species_results_filtered, species_results_filtered$Tool)

write_tsv(species_results_filtered, paste0("species_results_filtered_", assembler, ".tsv"))
### Calculate Precision, Recall, F1, L2
species_metrics <- data.frame()
cat("Calculating metrics for each tool...\n")
for (tool in tools){
  if(tool != target){
    dt_target <- species_abun[[ target ]] %>% select(all_of(rank))
    dt_pred <- species_abun[[ tool ]] %>% select(all_of(rank))
    roc_value <- cal_ROC(dt_target, dt_pred, tool)
    dt_temp <- species_abun[[ target ]] %>%
      full_join(species_abun[[ tool ]], by = setNames(rank,rank)) %>%
      select(rank, starts_with("Abundance")) %>%
      mutate(across(-c(rank), ~ replace_na(., 0)))
    
    l2_value <- with(dt_temp, L2(Abundance.x, Abundance.y))
    roc_value[["L2"]] <- l2_value
    species_metrics <- rbind(species_metrics, roc_value %>% as.data.frame())
  }
}
cat("Done calculating metrics.\n")
##Add assembler and method info
species_metrics <- species_metrics %>%
  mutate(Method = case_when(
    Tool %in% c("Kraken2", "FastViromeExplorer", "ViromeScan") ~ "Read-based",
    Tool %in% c("blast_Tophit_species", "blast_Common_Species", "blast_top_virus_species", "Kaiju", "MMseqs2 taxonomy") ~ "Contig-based"
  )) %>%
  mutate(Method = factor(Method, levels = c("Read-based", "Contig-based"))) %>%
  mutate(Tool = factor(Tool, levels = c(
    "Kraken2", "FastViromeExplorer", "ViromeScan",
    "blast_Tophit_species", "blast_Common_Species", "blast_top_virus_species",
    "Kaiju", "MMseqs2 taxonomy"
  ))) %>%
  mutate(Assembler = assembler)

# Plot: F1 Score vs L2 Distance
p1 <- ggplot(species_metrics, aes(x = L2, y = F1score, color = Tool, shape = Method)) +
  geom_point(size = 4) +
  geom_text(aes(label = Tool), vjust = -0.5, hjust = 0.5, size = 3) +
  labs(title = paste("F1 Score vs L2 Distance at Species Level -", assembler),
       x = "L2 Distance (lower is better)",
       y = "F1 Score (higher is better)") +
  theme_bw() +
  theme(legend.position = "right")

# Save plot as PDF (no additional packages required)
plot_file <- sub("\\.tsv$", "_plot.pdf", params$output)
ggsave(plot_file, plot = p1, width = 10, height = 7, units = "in")

# Save results
cat("\n=== Saving results ===\n")
write_tsv(species_metrics, params$output)

cat("Saved:\n")
cat("  -", params$output, "\n")
cat("  -", plot_file, "\n")
cat("\nProcessing complete!\n")
