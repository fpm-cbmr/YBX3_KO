library(tidyverse)


# Protein-level analysis --------------------------------------------------

# Load protein quantification data, drop the ProteinGroups column,
# remove duplicate gene entries, remove blank gene names,
# and set gene names as row names
tmp <- vroom::vroom(here::here("data-raw/YBX3_KO_protein_quant.tsv")) |>
  dplyr::select(
    !PG.ProteinGroups
  ) |>
  dplyr::filter(!duplicated(PG.Genes)) |>
  dplyr::filter(!PG.Genes == " ") |>
  tibble::column_to_rownames("PG.Genes")

# Rename columns to meaningful sample identifiers
colnames(tmp) <- c("WT1",
                   "WT2",
                   "WT3",
                   "KOg2_1",
                   "KOg2_2",
                   "KOg2_3",
                   "KOg5_1",
                   "KOg5_2")

# Sample metadata mapping each sample to its experimental group
metadata <- data.frame("sample_id" =  c("WT1",
                                        "WT2",
                                        "WT3",
                                        "KOg2_1",
                                        "KOg2_2",
                                        "KOg2_3",
                                        "KOg5_1",
                                        "KOg5_2"),
                       "grouping" = c("WT", "WT", "WT", "KO_g2", "KO_g2", "KO_g2", "KO_g5", "KO_g5")
)

# Log2-transform intensities, transpose so samples are rows,
# extract Ybx3 column only, and join with metadata
data <- tmp |>
  log2() |>
  t() |>
  as.data.frame() |>
  dplyr::select("Ybx3") |>
  tibble::rownames_to_column("sample_id") |>
  dplyr::inner_join(metadata)

# Boxplot of Ybx3 protein abundance across experimental groups
data |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = grouping,
      y = Ybx3,
      fill = grouping
    )
  ) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_point() +
  ggplot2::theme_classic() +
  ggplot2::ggtitle("Ybx3") +
  ggplot2::ylab("LFQ intensity (log2)") +
  ggplot2::scale_fill_viridis_d(option = "turbo") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.title.x = element_blank(),
    text = element_text(size = 6
  ))

ggplot2::ggsave(here::here("plots/protein_LFQs.pdf"),
                units = "mm",
                height = 60,
                width = 60)


# YBX3 peptide-level analysis ---------------------------------------------

# Load peptide quantification data, keep only Ybx3 peptides,
# drop quality/filter columns, and set precursor ID as row names
tmp <- vroom::vroom(here::here("data-raw/YBX3_KO_peptide_quant.tsv")) |>
  dplyr::select(
    !c(PG.Qvalue, PEP.IsProteotypic, PEP.IsGeneSpecific)
  ) |>
  dplyr::filter(PG.Genes == "Ybx3") |>
  dplyr::select(!PG.Genes) |>
  tibble::column_to_rownames("EG.PrecursorId")

# Rename columns to meaningful sample identifiers
colnames(tmp) <- c("WT1",
                   "WT2",
                   "WT3",
                   "KOg2_1",
                   "KOg2_2",
                   "KOg2_3",
                   "KOg5_1",
                   "KOg5_2")

# Reshape to long format, replace "Filtered" strings with NA,
# convert LFQ values to numeric, log2-transform, and join metadata
tmp2 <- tmp |>
  as.data.frame() |>
  tibble::rownames_to_column("peptide_seq") |>
  tidyr::pivot_longer(cols = !peptide_seq, names_to = "sample_id", values_to = "LFQs") |>
  dplyr::mutate(
    LFQs = case_when(
      LFQs == "Filtered" ~ NA,
      TRUE ~ LFQs
    )
  ) |>
  mutate(
    LFQs = as.numeric(LFQs),
    LFQs = log2(LFQs)
  ) |>
  inner_join(metadata)

# Faceted boxplot showing log2 LFQ intensities per peptide across groups
tmp2 |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = grouping,
      y = LFQs,
      fill = grouping
    )
  ) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_point() +
  ggplot2::theme_classic() +
  ggplot2::ylab("LFQ intensity (log2)") +
  ggplot2::scale_fill_viridis_d(option = "turbo") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.position = "none",
    axis.title.x = element_blank(),
    text = element_text(size = 6)
  ) +
  facet_wrap(~peptide_seq, ncol = 6)

ggplot2::ggsave(here::here("plots/peptide_LFQs.pdf"),
                units = "mm",
                height = 90,
                width = 180)


# LFQs peptide sequence and position --------------------------------------

# Compute per-group mean LFQ for each peptide, then calculate
# log2 fold-change as KO_g2 minus WT (negative = loss of abundance)
LFQ_seq <- tmp2 |>
  group_by(grouping, peptide_seq) |>
  summarise(
    mean_LFQ = mean(LFQs, na.rm = TRUE)
  ) |>
  pivot_wider(names_from = grouping,
              values_from = mean_LFQ) |>
  mutate(
    log2FC = KO_g2 - WT
  )


# log2FC color gradient ---------------------------------------------------

# Bar plot of per-peptide log2FC (KO_g2 vs WT), with bars ordered by
# fold-change magnitude and filled by a continuous gradient.
# Scale runs from 0.1 (near-zero change, white) to -4 (strong loss, darkblue).
# Values outside the range are clamped to the nearest limit.
p_bars <- LFQ_seq |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = reorder(peptide_seq, log2FC),
      y = log2FC,
      fill = log2FC
    )
  ) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_gradient2(
    low  = "darkred",
    mid = "grey",
    high = "darkblue",
    limits = c(-4, 0.1),
    oob = scales::squish,   # clamp out-of-range values instead of turning them grey
    name = "log2FC"
  ) +
  ggplot2::theme_classic() +
  ggplot2::ylab("log2FC (KO_g2 - WT)") +
  ggplot2::xlab("Peptide") +
  ggplot2::theme(
    text                  = ggplot2::element_text(size = 5),
    axis.text.x           = ggplot2::element_text(angle = 45, hjust = 1, size = 4),
    axis.text.y           = ggplot2::element_text(size = 4),
    axis.title            = ggplot2::element_text(size = 5),
    legend.title          = ggplot2::element_text(size = 5),
    legend.text           = ggplot2::element_text(size = 4),
    legend.key.height     = ggplot2::unit(3, "mm"),
    legend.key.width      = ggplot2::unit(2, "mm"),
    plot.margin           = ggplot2::margin(1, 1, 1, 1, "mm")
  )

ggplot2::ggsave(here::here("plots/peptide_LFQs_bars.pdf"),
                plot   = p_bars,
                units  = "mm",
                height = 60,
                width  = 120)

# Extract the legend as a standalone grob and save at a larger size
p_bars_legend <- cowplot::get_legend(
  p_bars +
    ggplot2::theme(
      legend.title     = ggplot2::element_text(size = 8),
      legend.text      = ggplot2::element_text(size = 6),
      legend.key.height = ggplot2::unit(6, "mm"),
      legend.key.width  = ggplot2::unit(3, "mm")
    )
)

ggplot2::ggsave(here::here("plots/peptide_LFQs_bars_legend.pdf"),
                plot   = p_bars_legend,
                units  = "mm",
                height = 60,
                width  = 30)


# Differential abundance --------------------------------------------------

# Log2-transform and quantile-normalize peptide intensities across samples,
# then run a limma linear model to test KO_g2 vs WT
data <- tmp |>
  log2() |>
  limma::normalizeBetweenArrays() |>
  as.data.frame()

# No-intercept design matrix so each group gets its own coefficient
design_matrix <- model.matrix(
  ~ 0 + metadata$grouping,
  data
)

colnames(design_matrix) <- c(
  "KO_g2",
  "KO_g5",
  "WT"
)

# Fit per-peptide linear models
fit <- limma::lmFit(
  data,
  design_matrix
)

# Define the KO_g2 vs WT contrast
contrast_matrix <- limma::makeContrasts(
  "KO_g2-WT" = KO_g2-WT,
  levels = design_matrix
)

# Re-fit on the contrast and apply empirical Bayes shrinkage
tmp <- limma::contrasts.fit(
  fit,
  contrast_matrix
)

tmp <- limma::eBayes(tmp)

# Extract all results sorted by p-value
DE_results <- limma::topTable(tmp,
                              sort.by = "P",
                              n = Inf
)

# Interactive volcano for quick exploration (not saved)
plotly::ggplotly(
  DE_results |>
    dplyr::mutate("significant" = dplyr::case_when(
      adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated",
      adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated",
      TRUE ~ "not significant"
    )) |>
    tibble::rownames_to_column("Genes") |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = logFC,
        y = -log10(P.Value),
        color = significant,
        names = Genes
      )
    ) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_color_manual(values = c(
      "blue",
      "grey",
      "red"
    )) +
    ggplot2::theme_minimal()
)

# Annotate significance, recover gene names as a column,
# and flag specific genes of interest for labelling
DE_results <- DE_results |>
  dplyr::mutate("significant" = dplyr::case_when(
    adj.P.Val <= 0.05 & logFC > 0 ~ "Upregulated",
    adj.P.Val <= 0.05 & logFC < 0 ~ "Downregulated",
    TRUE ~ "not significant"
  )) |>
  tibble::rownames_to_column("Genes") |>
  dplyr::mutate("names" = dplyr::case_when(
    Genes %in% c(
      "Ybx3",
      "Col1a1"
    ) ~ Genes,
    TRUE ~ ""
  ))

# Publication-style volcano plot; only significant hits are coloured,
# and selected genes are labelled with repelled text boxes
DE_results |>
  ggplot2::ggplot(ggplot2::aes(
    x = logFC,
    y = -log10(P.Value),
    color = significant,
    names = Genes
  )) +
  ggplot2::geom_point(
    size = 0.5,
    alpha =0.8
  ) +
  ggplot2::scale_color_manual(values = c("darkblue",
                                         "lightgrey"
                                         )) +
  ggplot2::theme_classic() +
  ggplot2::ggtitle("KO Guide 2 vs Control") +
  ggplot2::xlab("log2FC (KO_g2 - WT)") +
  ggplot2::ylab("-log10(P-value)") +
  ggplot2::theme(
    text = ggplot2::element_text(
      size = 5,
      colour = "black",
      face = "bold"
    ),
    strip.text = ggplot2::element_text(colour = "white"),
    strip.background = ggplot2::element_rect(fill = "black"),
    legend.position = "none",
    plot.title = ggplot2::element_text(hjust = 0.5,
                                       face = "bold")
  ) +
  ggrepel::geom_label_repel(
    data = DE_results |>
      dplyr::filter(!names == ""),
    mapping = ggplot2::aes(
      x = logFC,
      y = -log10(P.Value),
      fill = significant,
      label = names
    ),
    color = "black",
    size = 1.5,
    label.padding = 0.1,
    min.segment.length = 0.1,
    segment.size = 0.2,
    force = 40
  ) +
  ggplot2::scale_fill_manual(values = c(
    "lightblue",
    "#d5d4d4"))

ggplot2::ggsave(here::here("plots/volcano_DE.pdf"),
                units = "mm",
                height = 60,
                width = 60)


