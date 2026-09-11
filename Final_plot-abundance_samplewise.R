# ============================================================
# eggnog_samplewise_comprehensive_publication.R
# Full Functional Profiling: Multi-Domain Plots & Exhaustive TSVs
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(curl)
})

# ============================================================
# 1. USER SETTINGS & DIRECTORIES
# ============================================================
input_dir <- "/home/gsbtmjrf/plots"
args <- commandArgs(trailingOnly = TRUE)
sample_name <- if (length(args) >= 1) args[[1]] else "Chemical_03"

output_dir <- file.path(input_dir, sample_name)
top_n <- 20

# Create structured output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "02_functional_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "03_abundance_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "04_plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "06_annotated_matrices"), recursive = TRUE, showWarnings = FALSE)

anno_file <- file.path(input_dir, paste0(sample_name, ".emapper.annotations"))

cat("\n============================================================\n")
cat("COMPREHENSIVE SINGLE-SAMPLE FUNCTIONAL ANALYSIS\n")
cat("Sample     :", sample_name, "\n")
cat("Input Path :", anno_file, "\n")
cat("============================================================\n\n")

if (!file.exists(anno_file)) {
  stop("ERROR: Annotation file not found at: ", anno_file)
}

# ============================================================
# 2. READ ANNOTATION FILE
# ============================================================
df <- read_tsv(anno_file, comment = "##", na = c("", "-", "NA"), show_col_types = FALSE)
names(df)[1] <- sub("^#", "", names(df)[1])

# ============================================================
# 3. NCBI COG MAPPING & ABUNDANCE CALCULATION
# ============================================================
cog_def_url <- "https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.def.tab"
cog_def_file <- file.path(input_dir, "cog-24.def.tab")

if (!file.exists(cog_def_file)) {
  tryCatch({
    download.file(cog_def_url, cog_def_file, mode = "wb", quiet = TRUE)
  }, error = function(e) {
    message("Notice: Could not download NCBI COG table. Using internal definitions.")
  })
}

cog_dict <- tribble(
  ~Letter, ~Name,                                                       ~Supergroup,
  "J",     "Translation, ribosomal structure and biogenesis",           "Information storage and processing",
  "A",     "RNA processing and modification",                           "Information storage and processing",
  "K",     "Transcription",                                             "Information storage and processing",
  "L",     "Replication, recombination and repair",                     "Information storage and processing",
  "B",     "Chromatin structure and dynamics",                          "Information storage and processing",
  "D",     "Cell cycle control, cell division, chromosome partitioning","Cellular processes and signaling",
  "Y",     "Nuclear structure",                                         "Cellular processes and signaling",
  "V",     "Defense mechanisms",                                        "Cellular processes and signaling",
  "T",     "Signal transduction mechanisms",                            "Cellular processes and signaling",
  "M",     "Cell wall/membrane/envelope biogenesis",                    "Cellular processes and signaling",
  "N",     "Cell motility",                                             "Cellular processes and signaling",
  "Z",     "Cytoskeleton",                                              "Cellular processes and signaling",
  "W",     "Extracellular structures",                                  "Cellular processes and signaling",
  "U",     "Intracellular trafficking, secretion, and vesicular transport", "Cellular processes and signaling",
  "O",     "Posttranslational modification, protein turnover, chaperones","Cellular processes and signaling",
  "C",     "Energy production and conversion",                          "Metabolism",
  "G",     "Carbohydrate transport and metabolism",                     "Metabolism",
  "E",     "Amino acid transport and metabolism",                       "Metabolism",
  "F",     "Nucleotide transport and metabolism",                       "Metabolism",
  "H",     "Coenzyme transport and metabolism",                         "Metabolism",
  "I",     "Lipid transport and metabolism",                            "Metabolism",
  "P",     "Inorganic ion transport and metabolism",                    "Metabolism",
  "Q",     "Secondary metabolites biosynthesis, transport and catabolism","Metabolism",
  "R",     "General function prediction only",                          "Poorly characterized",
  "S",     "Function unknown",                                          "Poorly characterized"
)

if (file.exists(cog_def_file)) {
  ncbi_map <- read_tsv(cog_def_file, col_names = c("COG_ID", "NCBI_Category", "Name", "Gene", "Pathway", "PubMed", "PDB"), show_col_types = FALSE) %>%
    select(COG_ID, NCBI_Category)
  
  df_cog_mapped <- df %>%
    select(query, COG_category) %>%
    filter(!is.na(COG_category)) %>%
    left_join(ncbi_map, by = c("COG_category" = "COG_ID")) %>%
    mutate(Final_Category = coalesce(NCBI_Category, COG_category))
} else {
  df_cog_mapped <- df %>%
    select(query, COG_category) %>%
    filter(!is.na(COG_category)) %>%
    mutate(Final_Category = COG_category)
}

cog_extracted <- df_cog_mapped %>%
  mutate(Final_Category = str_extract_all(toupper(Final_Category), "[A-Z]")) %>%
  unnest(Final_Category) %>%
  filter(Final_Category %in% cog_dict$Letter) %>%
  distinct(query, Final_Category) %>%
  left_join(cog_dict, by = c("Final_Category" = "Letter")) %>%
  rename(COG_category = Final_Category, COG_Name = Name)

cog_abundance <- cog_extracted %>%
  count(COG_category, COG_Name, Supergroup, name = "Gene_Count") %>%
  arrange(desc(Gene_Count)) %>%
  mutate(Percent = 100 * Gene_Count / sum(Gene_Count))

supergroup_abundance <- cog_extracted %>%
  count(Supergroup, name = "Gene_Count") %>%
  mutate(Percent = 100 * Gene_Count / sum(Gene_Count)) %>%
  arrange(desc(Gene_Count))

# ============================================================
# 4. EXTRACTION HELPER FOR OTHER FUNCTIONAL DATABASES
# ============================================================
extract_field <- function(df, col, fix_map_prefix = FALSE) {
  if (!col %in% names(df)) return(tibble(query = character(), !!col := character()))
  out <- df %>%
    select(query, !!sym(col)) %>%
    filter(!is.na(.data[[col]])) %>%
    separate_rows(!!sym(col), sep = ",")
  
  if (fix_map_prefix) {
    out <- out %>%
      mutate(!!col := if_else(str_detect(.data[[col]], "^[0-9]+$"),
                              paste0("map", .data[[col]]), .data[[col]]))
  }
  distinct(out, query, !!sym(col))
}

abundance_table <- function(tbl, col) {
  if (nrow(tbl) == 0) return(tibble())
  tbl %>%
    count(!!sym(col), name = "Gene_Count") %>%
    arrange(desc(Gene_Count)) %>%
    mutate(Percent = 100 * Gene_Count / sum(Gene_Count))
}

ko      <- extract_field(df, "KEGG_ko")
pathway <- extract_field(df, "KEGG_Pathway", fix_map_prefix = TRUE)
module  <- extract_field(df, "KEGG_Module")
cazy    <- extract_field(df, "CAZy")
pfam    <- extract_field(df, "PFAMs")

ko_abundance      <- abundance_table(ko, "KEGG_ko")
pathway_abundance <- abundance_table(pathway, "KEGG_Pathway")
module_abundance  <- abundance_table(module, "KEGG_Module")
cazy_abundance    <- abundance_table(cazy, "CAZy")
pfam_abundance    <- abundance_table(pfam, "PFAMs")

# ============================================================
# 5. EXPORT COMPLETE, EXHAUSTIVE TABLES & ANNOTATED MATRICES
# ============================================================
tsv_dir <- file.path(output_dir, "03_abundance_tables")
mat_dir <- file.path(output_dir, "06_annotated_matrices")
dir.create(mat_dir, recursive = TRUE, showWarnings = FALSE)

# Summary Abundance Tables (Down to 1 count)
write_csv(cog_abundance, file.path(tsv_dir, paste0(sample_name, "_COG_abundance_all.csv")))
write_csv(supergroup_abundance, file.path(tsv_dir, paste0(sample_name, "_COG_Supergroup_abundance_all.csv")))
if(nrow(ko_abundance) > 0) write_csv(ko_abundance, file.path(tsv_dir, paste0(sample_name, "_KEGG_KO_abundance_all.csv")))
if(nrow(pathway_abundance) > 0) write_csv(pathway_abundance, file.path(tsv_dir, paste0(sample_name, "_KEGG_Pathway_abundance_all.csv")))
if(nrow(module_abundance) > 0) write_csv(module_abundance, file.path(tsv_dir, paste0(sample_name, "_KEGG_Module_abundance_all.csv")))
if(nrow(cazy_abundance) > 0) write_csv(cazy_abundance, file.path(tsv_dir, paste0(sample_name, "_CAZy_abundance_all.csv")))
if(nrow(pfam_abundance) > 0) write_csv(pfam_abundance, file.path(tsv_dir, paste0(sample_name, "_PFAM_abundance_all.csv")))

# Detailed Feature-to-Gene Mapping Ledgers (Separated & Sorted)
write_csv(cog_extracted %>% arrange(COG_category), file.path(mat_dir, paste0(sample_name, "_COG_detailed_annotations.csv")))
if(nrow(ko) > 0) write_csv(ko %>% arrange(KEGG_ko), file.path(mat_dir, paste0(sample_name, "_KEGG_KO_detailed_annotations.csv")))
if(nrow(pathway) > 0) write_csv(pathway %>% arrange(KEGG_Pathway), file.path(mat_dir, paste0(sample_name, "_KEGG_Pathway_detailed_annotations.csv")))
if(nrow(module) > 0) write_csv(module %>% arrange(KEGG_Module), file.path(mat_dir, paste0(sample_name, "_KEGG_Module_detailed_annotations.csv")))
if(nrow(cazy) > 0) write_csv(cazy %>% arrange(CAZy), file.path(mat_dir, paste0(sample_name, "_CAZy_detailed_annotations.csv")))
if(nrow(pfam) > 0) write_csv(pfam %>% arrange(PFAMs), file.path(mat_dir, paste0(sample_name, "_PFAM_detailed_annotations.csv")))

# Global Gene-level Master Table
gene_level_export <- df %>%
  select(any_of(c("query", "seed_ortholog", "evalue", "score", "COG_category", "Preferred_name", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "CAZy", "PFAMs"))) %>%
  arrange(desc(score))

write_csv(gene_level_export, file.path(output_dir, "02_functional_tables", paste0(sample_name, "_gene_level_annotations.csv")))

cat("Successfully exported full sorted abundance tables and separated feature matrices down to 1 count.\n")

# ============================================================
# 6. GENERATE MULTI-DOMAIN PUBLICATION PLOTS
# ============================================================
plot_dir <- file.path(output_dir, "04_plots")
supergroup_colors <- c(
  "Information storage and processing" = "#1b9e77",
  "Cellular processes and signaling"   = "#d95f02",
  "Metabolism"                         = "#7570b3",
  "Poorly characterized"               = "#999999"
)

# 1. COG Profile Plot
p_cog <- cog_abundance %>%
  slice_head(n = top_n) %>%
  mutate(Label = paste0("[", COG_category, "] ", COG_Name)) %>%
  ggplot(aes(x = reorder(Label, Percent), y = Percent, fill = Supergroup)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = supergroup_colors) +
  labs(title = paste0("COG Functional Category Profile - ", sample_name),
       x = NULL, y = "Share of COG-assigned genes (%)", fill = "COG supergroup") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4))

ggsave(file.path(plot_dir, paste0(sample_name, "_COG_profile.png")), p_cog, width = 11, height = 7, dpi = 300)

# 2. KEGG Pathways Plot
if(nrow(pathway_abundance) > 0) {
  p_path <- pathway_abundance %>%
    slice_head(n = top_n) %>%
    ggplot(aes(x = reorder(KEGG_Pathway, Percent), y = Percent)) +
    geom_col(fill = "#2b8cbe", width = 0.75) +
    coord_flip() +
    labs(title = paste0("Top ", top_n, " KEGG Pathways - ", sample_name),
         x = NULL, y = "Share of KEGG-pathway-assigned genes (%)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4))
  
  ggsave(file.path(plot_dir, paste0(sample_name, "_KEGG_Pathways.png")), p_path, width = 10, height = 7, dpi = 300)
}

# 3. CAZy Carbohydrate-Active Enzymes Plot
if(nrow(cazy_abundance) > 0) {
  p_cazy <- cazy_abundance %>%
    slice_head(n = top_n) %>%
    ggplot(aes(x = reorder(CAZy, Percent), y = Percent)) +
    geom_col(fill = "#31a354", width = 0.75) +
    coord_flip() +
    labs(title = paste0("Top ", top_n, " CAZy Families - ", sample_name),
         x = NULL, y = "Share of CAZy-assigned genes (%)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4))
  
  ggsave(file.path(plot_dir, paste0(sample_name, "_CAZy_profile.png")), p_cazy, width = 10, height = 7, dpi = 300)
}

# 4. PFAM Protein Domains Plot
if(nrow(pfam_abundance) > 0) {
  p_pfam <- pfam_abundance %>%
    slice_head(n = top_n) %>%
    ggplot(aes(x = reorder(PFAMs, Percent), y = Percent)) +
    geom_col(fill = "#de2d26", width = 0.75) +
    coord_flip() +
    labs(title = paste0("Top ", top_n, " PFAM Domains - ", sample_name),
         x = NULL, y = "Share of PFAM-assigned genes (%)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4))
  
  ggsave(file.path(plot_dir, paste0(sample_name, "_PFAM_profile.png")), p_pfam, width = 10, height = 7, dpi = 300)
}

cat("\n============================================================\n")
cat("ALL ANALYSES COMPLETED SUCCESSFULLY!\n")
cat("Plots saved to :", plot_dir, "\n")
cat("Full TSVs saved:", tsv_dir, "\n")
cat("Detailed Matrices:", mat_dir, "\n")
cat("============================================================\n")
