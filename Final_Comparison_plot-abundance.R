# ============================================================
# compare_treatments_eggnog.R
# Chemical vs Natural comparative functional analysis
# Builds cross-sample abundance matrices from raw eggNOG
# .emapper.annotations files, then runs diversity + differential
# abundance statistics and comparison plots.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)     # diversity(), vegdist(), adonis2(), betadisper()
  library(ggpubr)    # stat_compare_means(), themes
  library(pheatmap)  # heatmaps
})

# ============================================================
# 1. USER SETTINGS
# ============================================================
input_dir  <- "/home/gsbtmjrf/plots"         # same folder as your per-sample .emapper.annotations files
output_dir <- file.path(input_dir, "COMPARATIVE_Chemical_vs_Natural")
top_n_plot <- 20                             # top N features shown in bar/heatmap plots
alpha_sig  <- 0.05

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "matrices"),   showWarnings = FALSE)
dir.create(file.path(output_dir, "diversity"),  showWarnings = FALSE)
dir.create(file.path(output_dir, "diffabund"),  showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"),      showWarnings = FALSE)

# COG letter -> name/supergroup dictionary
cog_dict <- tribble(
  ~Letter, ~Name,                                                        ~Supergroup,
  "J",     "Translation, ribosomal structure and biogenesis",            "Information storage and processing",
  "A",     "RNA processing and modification",                            "Information storage and processing",
  "K",     "Transcription",                                              "Information storage and processing",
  "L",     "Replication, recombination and repair",                      "Information storage and processing",
  "B",     "Chromatin structure and dynamics",                           "Information storage and processing",
  "D",     "Cell cycle control, cell division, chromosome partitioning","Cellular processes and signaling",
  "Y",     "Nuclear structure",                                          "Cellular processes and signaling",
  "V",     "Defense mechanisms",                                         "Cellular processes and signaling",
  "T",     "Signal transduction mechanisms",                             "Cellular processes and signaling",
  "M",     "Cell wall/membrane/envelope biogenesis",                     "Cellular processes and signaling",
  "N",     "Cell motility",                                              "Cellular processes and signaling",
  "Z",     "Cytoskeleton",                                               "Cellular processes and signaling",
  "W",     "Extracellular structures",                                   "Cellular processes and signaling",
  "U",     "Intracellular trafficking, secretion, and vesicular transport", "Cellular processes and signaling",
  "O",     "Posttranslational modification, protein turnover, chaperones","Cellular processes and signaling",
  "C",     "Energy production and conversion",                           "Metabolism",
  "G",     "Carbohydrate transport and metabolism",                      "Metabolism",
  "E",     "Amino acid transport and metabolism",                        "Metabolism",
  "F",     "Nucleotide transport and metabolism",                        "Metabolism",
  "H",     "Coenzyme transport and metabolism",                          "Metabolism",
  "I",     "Lipid transport and metabolism",                             "Metabolism",
  "P",     "Inorganic ion transport and metabolism",                     "Metabolism",
  "Q",     "Secondary metabolites biosynthesis, transport and catabolism","Metabolism",
  "R",     "General function prediction only",                           "Poorly characterized",
  "S",     "Function unknown",                                           "Poorly characterized"
)

# ============================================================
# 2. DISCOVER SAMPLES & READ ALL ANNOTATION FILES
# ============================================================
anno_files <- list.files(input_dir, pattern = "\\.emapper\\.annotations$", full.names = TRUE)
if (length(anno_files) == 0) stop("No .emapper.annotations files found in ", input_dir)

sample_names <- str_remove(basename(anno_files), "\\.emapper\\.annotations$")
groups <- str_extract(sample_names, "^[A-Za-z]+")   # "Chemical" / "Natural"

metadata <- tibble(Sample = sample_names, Group = groups)
cat("Detected samples:\n"); print(metadata)

group_sizes <- metadata %>% count(Group, name = "n")
cat("\nSamples per group:\n"); print(group_sizes)
min_group_n <- min(group_sizes$n)
STATS_OK <- min_group_n >= 2   

if (!STATS_OK) {
  cat("\nNOTE: At least one group has < 2 replicates (n=1 detected).\n",
      "Statistical tests will be skipped, but matrices, diversity values, and a readable summary report will be generated.\n\n")
}

read_one <- function(f, sname) {
  read_tsv(f, comment = "##", na = c("", "-", "NA"), show_col_types = FALSE) %>%
    rename_with(~ sub("^#", "", .x), 1) %>%
    mutate(Sample = sname)
}

all_df <- map2_dfr(anno_files, sample_names, read_one)
total_genes_per_sample <- all_df %>% count(Sample, name = "Total_Genes")

# ============================================================
# 3. BUILD FEATURE x SAMPLE COUNT MATRICES
# ============================================================
extract_field <- function(df, col, fix_map_prefix = FALSE) {
  if (!col %in% names(df)) return(tibble(query = character(), Sample = character(), !!col := character()))
  out <- df %>%
    select(query, Sample, !!sym(col)) %>%
    filter(!is.na(.data[[col]])) %>%
    separate_rows(!!sym(col), sep = ",")
  if (fix_map_prefix) {
    out <- out %>%
      mutate(!!col := if_else(str_detect(.data[[col]], "^[0-9]+$"),
                            paste0("map", .data[[col]]), .data[[col]]))
  }
  distinct(out, query, Sample, !!sym(col))
}

make_matrix <- function(long_tbl, feature_col) {
  if (nrow(long_tbl) == 0) return(NULL)
  long_tbl %>%
    count(!!sym(feature_col), Sample, name = "n") %>%
    pivot_wider(names_from = Sample, values_from = n, values_fill = 0)
}

cog_long <- all_df %>%
  select(query, Sample, COG_category) %>%
  filter(!is.na(COG_category)) %>%
  mutate(COG_category = str_extract_all(toupper(COG_category), "[A-Z]")) %>%
  unnest(COG_category) %>%
  filter(COG_category %in% cog_dict$Letter) %>%
  distinct(query, Sample, COG_category)

cog_matrix     <- make_matrix(cog_long, "COG_category")
ko_matrix      <- make_matrix(extract_field(all_df, "KEGG_ko"), "KEGG_ko")
pathway_matrix <- make_matrix(extract_field(all_df, "KEGG_Pathway", fix_map_prefix = TRUE), "KEGG_Pathway")
module_matrix  <- make_matrix(extract_field(all_df, "KEGG_Module"), "KEGG_Module")
cazy_matrix    <- make_matrix(extract_field(all_df, "CAZy"), "CAZy")

matrices <- list(COG = cog_matrix, KEGG_KO = ko_matrix, KEGG_Pathway = pathway_matrix,
                 KEGG_Module = module_matrix, CAZy = cazy_matrix)
matrices <- matrices[!map_lgl(matrices, is.null)]

for (nm in names(matrices)) {
  write_csv(matrices[[nm]], file.path(output_dir, "matrices", paste0(nm, "_count_matrix.csv")))
}

to_num_matrix <- function(wide_tbl) {
  feat_col <- names(wide_tbl)[1]
  m <- as.matrix(wide_tbl[, -1, drop = FALSE])
  rownames(m) <- wide_tbl[[feat_col]]
  storage.mode(m) <- "numeric"
  m[, metadata$Sample, drop = FALSE]
}

# ============================================================
# 5. ALPHA DIVERSITY (Shannon) PER FUNCTIONAL LEVEL
# ============================================================
alpha_results <- list()
for (nm in names(matrices)) {
  m <- to_num_matrix(matrices[[nm]])
  shannon <- vegan::diversity(t(m), index = "shannon")
  simpson <- vegan::diversity(t(m), index = "simpson")
  richness <- colSums(m > 0)
  
  div_tbl <- tibble(Sample = names(shannon), Level = nm,
                    Shannon = shannon, Simpson = simpson, Richness = richness) %>%
    left_join(metadata, by = "Sample")
  alpha_results[[nm]] <- div_tbl

  p_alpha <- ggplot(div_tbl, aes(Group, Shannon, fill = Group)) +
    geom_bar(stat = "identity", alpha = 0.7, width = 0.5) +
    geom_text(aes(label = round(Shannon, 3)), vjust = -0.5, size = 4) +
    labs(title = paste0(nm, " Shannon diversity: Chemical vs Natural"), y = "Shannon index", x = NULL) +
    theme_classic(base_size = 12) + theme(legend.position = "none") +
    expand_limits(y = max(div_tbl$Shannon) * 1.15)
  
  if (STATS_OK && length(unique(metadata$Group)) > 1 && min(group_sizes$n) >= 2) {
    p_alpha <- p_alpha + stat_compare_means(method = "wilcox.test", label = "p.format")
  }
  
  ggsave(file.path(output_dir, "diversity", paste0(nm, "_Shannon_boxplot.png")),
         p_alpha, width = 5, height = 5, dpi = 300)
}
alpha_all <- bind_rows(alpha_results)
write_csv(alpha_all, file.path(output_dir, "diversity", "alpha_diversity_all_levels.csv"))

# ============================================================
# 6. BETA DIVERSITY & HUMAN-READABLE REPORT
# ============================================================
report_lines <- c(
  "============================================================",
  "       COMPARATIVE FUNCTIONAL ANALYSIS: CHEMICAL VS NATURAL  ",
  "============================================================\n",
  paste0("Total Samples Analyzed: ", nrow(metadata)),
  paste0("Groups: ", paste(unique(metadata$Group), collapse = ", ")),
  "\n--- 1. ALPHA DIVERSITY (SHANNON INDEX) ---"
)

for (nm in names(alpha_results)) {
  df_sub <- alpha_results[[nm]]
  report_lines <- c(report_lines, paste0("\nLevel: ", nm))
  for (i in seq_len(nrow(df_sub))) {
    report_lines <- c(report_lines, sprintf("  Sample: %s (%s) | Shannon: %.4f | Richness: %d", 
                                            df_sub$Sample[i], df_sub$Group[i], df_sub$Shannon[i], df_sub$Richness[i]))
  }
}

writeLines(report_lines, file.path(output_dir, "comparison_summary_report.txt"))

cat("\n============================================================\n")
cat("COMPARATIVE ANALYSIS COMPLETE\n")
cat("Output directory:", output_dir, "\n")
cat("  - Human-readable report saved to: COMPARATIVE_Chemical_vs_Natural/comparison_summary_report.txt\n")
cat("============================================================\n")
