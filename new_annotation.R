library(tidyverse)

# ============================================================
# USER SETTINGS
# ============================================================

input_dir   <- "/home/gsbtmjrf/plots"
sample_name <- "Chemical_03"

output_dir <- file.path(input_dir, sample_name)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "01_master_annotation"), showWarnings = FALSE)
dir.create(file.path(output_dir, "02_functional_tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "03_abundance_tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "04_plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "05_QC"), showWarnings = FALSE)

anno_file <- file.path(
  input_dir,
  paste0(sample_name, ".emapper.annotations")
)

cat("\n============================================\n")
cat("Sample:", sample_name, "\n")
cat("Input :", anno_file, "\n")
cat("============================================\n\n")

# ============================================================
# 1. READ EGGNOG FILE
# ============================================================

df <- read_tsv(
  anno_file,
  comment = "##",
  na = c("", "-", "NA"),
  show_col_types = FALSE
)

# Remove leading # from query header
names(df)[1] <- sub("^#", "", names(df)[1])

# ============================================================
# 2. VERIFY EXACT EGGNOG V3 SCHEMA
# ============================================================

expected_cols <- c(
  "query",
  "seed_ortholog",
  "evalue",
  "score",
  "eggNOG_OGs",
  "tax_ceiling",
  "farthest_donor_lineage",
  "COG_category",
  "Preferred_name",
  "GOs",
  "EC",
  "KEGG_ko",
  "KEGG_Pathway",
  "KEGG_Module",
  "KEGG_Reaction",
  "KEGG_rclass",
  "BRITE",
  "KEGG_TC",
  "CAZy",
  "BiGG_Reaction",
  "PFAMs",
  "annotation_confidence"
)

if (!identical(names(df), expected_cols)) {
  cat("Observed columns:\n")
  print(names(df))

  stop(
    "\nERROR: Annotation schema does not exactly match expected ",
    "eggNOG-mapper v3.0.0-beta6 schema."
  )
}

if (ncol(df) != 22) {
  stop("ERROR: Expected 22 columns but found ", ncol(df))
}

cat("Schema check: PASS\n")

# ============================================================
# 3. MASTER ANNOTATION TABLE
# ============================================================

write_tsv(
  df,
  file.path(
    output_dir,
    "01_master_annotation",
    paste0(sample_name, "_annotation_master.tsv")
  ),
  na = "-"
)

# ============================================================
# 4. TOTAL NUMBER OF NR PROTEINS / ANNOTATED PROTEINS
# ============================================================

total_queries <- nrow(df)

annotated_queries <- df %>%
  filter(
    !is.na(seed_ortholog),
    seed_ortholog != "-"
  ) %>%
  nrow()

annotation_rate <- 100 * annotated_queries / total_queries

# ============================================================
# 5. COG TABLE
# ============================================================

cog <- df %>%
  select(
    query,
    COG_category,
    Preferred_name
  ) %>%
  filter(!is.na(COG_category)) %>%
  filter(COG_category != "-") %>%
  separate_rows(COG_category, sep = ",") %>%
  distinct(query, COG_category)

write_tsv(
  cog,
  file.path(
    output_dir,
    "02_functional_tables",
    paste0(sample_name, "_COG_annotation.tsv")
  ),
  na = "-"
)

cog_abundance <- cog %>%
  count(COG_category, name = "Gene_Count") %>%
  arrange(desc(Gene_Count)) %>%
  mutate(
    Percent = 100 * Gene_Count / sum(Gene_Count)
  )

write_csv(
  cog_abundance,
  file.path(
    output_dir,
    "03_abundance_tables",
    paste0(sample_name, "_COG_abundance.csv")
  )
)

# ============================================================
# 6. KEGG KO
# ============================================================

ko <- df %>%
  select(query, KEGG_ko, Preferred_name) %>%
  filter(!is.na(KEGG_ko), KEGG_ko != "-") %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  distinct(query, KEGG_ko)

write_tsv(
  ko,
  file.path(
    output_dir,
    "02_functional_tables",
    paste0(sample_name, "_KEGG_KO_annotation.tsv")
  ),
  na = "-"
)

ko_abundance <- ko %>%
  count(KEGG_ko, name = "Gene_Count") %>%
  arrange(desc(Gene_Count)) %>%
  mutate(
    Percent = 100 * Gene_Count / sum(Gene_Count)
  )

write_csv(
  ko_abundance,
  file.path(
    output_dir,
    "03_abundance_tables",
    paste0(sample_name, "_KEGG_KO_abundance.csv")
  )
)

# ============================================================
# 7. KEGG PATHWAY
# ============================================================

pathway <- df %>%
  select(query, KEGG_Pathway) %>%
  filter(!is.na(KEGG_Pathway), KEGG_Pathway != "-") %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  mutate(
    KEGG_Pathway = if_else(
      str_detect(KEGG_Pathway, "^[0-9]+$"),
      paste0("map", KEGG_Pathway),
      KEGG_Pathway
    )
  ) %>%
  distinct(query, KEGG_Pathway)

write_tsv(
  pathway,
  file.path(
    output_dir,
    "02_functional_tables",
    paste0(sample_name, "_KEGG_Pathway_annotation.tsv")
  ),
  na = "-"
)

pathway_abundance <- pathway %>%
  count(KEGG_Pathway, name = "Gene_Count") %>%
  arrange(desc(Gene_Count)) %>%
  mutate(
    Percent = 100 * Gene_Count / sum(Gene_Count)
  )

write_csv(
  pathway_abundance,
  file.path(
    output_dir,
    "03_abundance_tables",
    paste0(sample_name, "_KEGG_Pathway_abundance.csv")
  )
)

# ============================================================
# 8. KEGG MODULE
# ============================================================

module <- df %>%
  select(query, KEGG_Module) %>%
  filter(!is.na(KEGG_Module), KEGG_Module != "-") %>%
  separate_rows(KEGG_Module, sep = ",") %>%
  distinct(query, KEGG_Module)

write_tsv(
  module,
  file.path(
    output_dir,
    "02_functional_tables",
    paste0(sample_name, "_KEGG_Module_annotation.tsv")
  ),
  na = "-"
)

module_abundance <- module %>%
  count(KEGG_Module, name = "Gene_Count") %>%
  arrange(desc(Gene_Count)) %>%
  mutate(
    Percent = 100 * Gene_Count / sum(Gene_Count)
  )

write_csv(
  module_abundance,
  file.path(
    output_dir,
    "03_abundance_tables",
    paste0(sample_name, "_KEGG_Module_abundance.csv")
  )
)

# ============================================================
# 9. GO
# ============================================================

go <- df %>%
  select(query, GOs) %>%
  filter(!is.na(GOs), GOs != "-") %>%
  separate_rows(GOs, sep = ",") %>%
  distinct(query, GOs)

write_tsv(
  go,
  file.path(
    output_dir,
    "02_functional_tables",
    paste0(sample_name, "_GO_annotation.tsv")
  ),
  na = "-"
)

go_abundance <- go %>%
  count(GOs, name = "Gene_Count") %>%
  arrange(desc(Gene_Count)) %>%
  mutate(
    Percent = 100 * Gene_Count / sum(Gene_Count)
  )

write_csv(
  go_abundance,
  file.path(
    output_dir,
    "03_abundance_tables",
    paste0(sample_name, "_GO_abundance.csv")
  )
)

# ============================================================
# 10. EC
# ============================================================

ec <- df %>%
  select(query, EC) %>%
  filter(!is.na(EC), EC != "-") %>%
  separate_rows(EC, sep = ",") %>%
  distinct(query, EC)

write_tsv(
  ec,
  file.path(
    output_dir,
    "02_functional_tables",
    paste0(sample_name, "_EC_annotation.tsv")
  ),
  na = "-"
)

ec_abundance <- ec %>%
  count(EC, name = "Gene_Count") %>%
  arrange(desc(Gene_Count)) %>%
  mutate(
    Percent = 100 * Gene_Count / sum(Gene_Count)
  )

write_csv(
  ec_abundance,
  file.path(
    output_dir,
    "03_abundance_tables",
    paste0(sample_name, "_EC_abundance.csv")
  )
)

# ============================================================
# 11. CAZY
# ============================================================

cazy <- df %>%
  select(query, CAZy) %>%
  filter(!is.na(CAZy), CAZy != "-") %>%
  separate_rows(CAZy, sep = ",") %>%
  distinct(query, CAZy)

write_tsv(
  cazy,
  file.path(
    output_dir,
    "02_functional_tables",
    paste0(sample_name, "_CAZy_annotation.tsv")
  ),
  na = "-"
)

cazy_abundance <- cazy %>%
  count(CAZy, name = "Gene_Count") %>%
  arrange(desc(Gene_Count)) %>%
  mutate(
    Percent = 100 * Gene_Count / sum(Gene_Count)
  )

write_csv(
  cazy_abundance,
  file.path(
    output_dir,
    "03_abundance_tables",
    paste0(sample_name, "_CAZy_abundance.csv")
  )
)

# ============================================================
# 12. PFAM
# ============================================================

pfam <- df %>%
  select(query, PFAMs) %>%
  filter(!is.na(PFAMs), PFAMs != "-") %>%
  separate_rows(PFAMs, sep = ",") %>%
  distinct(query, PFAMs)

write_tsv(
  pfam,
  file.path(
    output_dir,
    "02_functional_tables",
    paste0(sample_name, "_PFAM_annotation.tsv")
  ),
  na = "-"
)

pfam_abundance <- pfam %>%
  count(PFAMs, name = "Gene_Count") %>%
  arrange(desc(Gene_Count)) %>%
  mutate(
    Percent = 100 * Gene_Count / sum(Gene_Count)
  )

write_csv(
  pfam_abundance,
  file.path(
    output_dir,
    "03_abundance_tables",
    paste0(sample_name, "_PFAM_abundance.csv")
  )
)

# ============================================================
# 13. FUNCTIONAL SUMMARY
# ============================================================

functional_summary <- df %>%
  select(
    query,
    Preferred_name,
    COG_category,
    GOs,
    EC,
    KEGG_ko,
    KEGG_Pathway,
    KEGG_Module,
    KEGG_Reaction,
    KEGG_rclass,
    BRITE,
    KEGG_TC,
    CAZy,
    BiGG_Reaction,
    PFAMs,
    annotation_confidence
  )

write_tsv(
  functional_summary,
  file.path(
    output_dir,
    "01_master_annotation",
    paste0(sample_name, "_Functional_Summary.tsv")
  ),
  na = "-"
)

# ============================================================
# 14. SAMPLE QC SUMMARY
# ============================================================

qc <- tibble(
  Sample = sample_name,
  Total_NR_proteins = total_queries,
  Annotated_proteins = annotated_queries,
  Annotation_percent = round(annotation_rate, 2),
  COG_annotated_genes = n_distinct(cog$query),
  KEGG_KO_annotated_genes = n_distinct(ko$query),
  KEGG_Pathway_annotated_genes = n_distinct(pathway$query),
  KEGG_Module_annotated_genes = n_distinct(module$query),
  GO_annotated_genes = n_distinct(go$query),
  EC_annotated_genes = n_distinct(ec$query),
  CAZy_annotated_genes = n_distinct(cazy$query),
  PFAM_annotated_genes = n_distinct(pfam$query)
)

write_csv(
  qc,
  file.path(
    output_dir,
    "05_QC",
    paste0(sample_name, "_functional_QC.csv")
  )
)

# ============================================================
# 15. PLOTS
# ============================================================

# ---------- COG ----------

p_cog <- cog_abundance %>%
  slice_head(n = 20) %>%
  ggplot(
    aes(
      x = reorder(COG_category, Gene_Count),
      y = Gene_Count
    )
  ) +
  geom_col() +
  coord_flip() +
  labs(
    title = paste0("Top COGs - ", sample_name),
    x = "COG identifier",
    y = "Gene assignments"
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "04_plots",
    paste0(sample_name, "_Top20_COGs.png")
  ),
  p_cog,
  width = 10,
  height = 7,
  dpi = 300
)

# ---------- KEGG PATHWAYS ----------

p_path <- pathway_abundance %>%
  slice_head(n = 20) %>%
  ggplot(
    aes(
      x = reorder(KEGG_Pathway, Gene_Count),
      y = Gene_Count
    )
  ) +
  geom_col() +
  coord_flip() +
  labs(
    title = paste0("Top 20 KEGG Pathways - ", sample_name),
    x = "KEGG pathway",
    y = "Gene-pathway assignments"
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "04_plots",
    paste0(sample_name, "_Top20_KEGG_Pathways.png")
  ),
  p_path,
  width = 10,
  height = 7,
  dpi = 300
)

# ---------- KEGG MODULES ----------

p_module <- module_abundance %>%
  slice_head(n = 20) %>%
  ggplot(
    aes(
      x = reorder(KEGG_Module, Gene_Count),
      y = Gene_Count
    )
  ) +
  geom_col() +
  coord_flip() +
  labs(
    title = paste0("Top 20 KEGG Modules - ", sample_name),
    x = "KEGG module",
    y = "Gene-module assignments"
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "04_plots",
    paste0(sample_name, "_Top20_KEGG_Modules.png")
  ),
  p_module,
  width = 10,
  height = 7,
  dpi = 300
)

# ---------- CAZY ----------

p_cazy <- cazy_abundance %>%
  slice_head(n = 20) %>%
  ggplot(
    aes(
      x = reorder(CAZy, Gene_Count),
      y = Gene_Count
    )
  ) +
  geom_col() +
  coord_flip() +
  labs(
    title = paste0("Top 20 CAZy Families - ", sample_name),
    x = "CAZy",
    y = "Gene assignments"
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "04_plots",
    paste0(sample_name, "_Top20_CAZy.png")
  ),
  p_cazy,
  width = 10,
  height = 7,
  dpi = 300
)

# ============================================================
# 16. FINAL MESSAGE
# ============================================================

cat("\n============================================\n")
cat("ANALYSIS COMPLETED:", sample_name, "\n")
cat("============================================\n")
cat("Total NR proteins       :", total_queries, "\n")
cat("Annotated proteins      :", annotated_queries, "\n")
cat("Annotation percentage   :", round(annotation_rate, 2), "%\n")
cat("COG genes               :", n_distinct(cog$query), "\n")
cat("KEGG KO genes           :", n_distinct(ko$query), "\n")
cat("KEGG pathway genes      :", n_distinct(pathway$query), "\n")
cat("KEGG module genes       :", n_distinct(module$query), "\n")
cat("GO genes                :", n_distinct(go$query), "\n")
cat("EC genes                :", n_distinct(ec$query), "\n")
cat("CAZy genes              :", n_distinct(cazy$query), "\n")
cat("Pfam genes              :", n_distinct(pfam$query), "\n")
cat("\nResults written to:\n")
cat(output_dir, "\n")
cat("============================================\n")
