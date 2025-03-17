library(ggplot2)
library(ggvenn)
library(readr)
library(dplyr)
library(glue)
library(patchwork)

source('03-scripts/R/seq_functions.R')

nuclei <- LoadDataset('Dec2024')
activity_colors <- LoadActivityColors()
subclass_colors <- LoadSubclassColors()


# subset object and get genes ----------------------------------------
library(glue)
library(ggvenn)

# A helper function that safely reads a CSV file. If the file doesn't exist,
# returns an empty tibble with the expected columns.
safe_read_csv <- function(filepath) {
  if (file.exists(filepath)) {
    read_csv(filepath, show_col_types = FALSE)
  } else {
    warning(glue("File not found: {filepath}. Using an empty data frame."))
    tibble(gene = character(), classification = character())
  }
}

# Example: load the subclasses to iterate over.
subclasses_raw <- LoadSubclassesToUse(nuclei)
subclasses <- subclasses_raw[-c(2,8,10:14)] |> print() # dont load small subclasses

df_all <- tibble(df_all <- tibble(
  subclass    = character(),
  time        = character(),
  direction   = character(),
  EE_specific = integer(),
  KA_specific = integer(),
  shared      = integer()
))


# loop thru subclasses ------------------------------------------
for (subclass in subclasses) {
  subclass_str <- ShrinkSubclassName(subclass)
  
  file_deg_30m_EE <- glue("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/{subclass_str}__EE30m_vs_SE.csv")
  file_deg_30m_KA <- glue("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/{subclass_str}__KA30m_vs_SE.csv")
  file_deg_6h_EE  <- glue("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/{subclass_str}__EE6h_vs_SE.csv")
  file_deg_6h_KA  <- glue("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/{subclass_str}__KA6h_vs_SE.csv")
  
  deg_30m_EE <- safe_read_csv(file_deg_30m_EE)
  deg_30m_KA <- safe_read_csv(file_deg_30m_KA)
  deg_6h_EE  <- safe_read_csv(file_deg_6h_EE)
  deg_6h_KA  <- safe_read_csv(file_deg_6h_KA)
  
  safe_filter <- function(df, direction) {
    if ("classification" %in% names(df)) {
      df |> filter(classification == direction) |> pull(gene)
    } else {
      character()
    }
  }
  
  # Create a tibble summarizing counts for 30m and 6h for up and down regulated genes.
  df_current_subclass <- tibble(
    subclass = subclass,
    time = c('30m', '30m', '6h', '6h'),
    direction = c('up', 'down', 'up', 'down'),
    EE_specific = c(
      length(setdiff(safe_filter(deg_30m_EE, 'upregulated'),
                     safe_filter(deg_30m_KA, 'upregulated'))),
      length(setdiff(safe_filter(deg_30m_EE, 'downregulated'),
                     safe_filter(deg_30m_KA, 'downregulated'))),
      length(setdiff(safe_filter(deg_6h_EE, 'upregulated'),
                     safe_filter(deg_6h_KA, 'upregulated'))),
      length(setdiff(safe_filter(deg_6h_EE, 'downregulated'),
                     safe_filter(deg_6h_KA, 'downregulated')))
    ),
    KA_specific = c(
      length(setdiff(safe_filter(deg_30m_KA, 'upregulated'),
                     safe_filter(deg_30m_EE, 'upregulated'))),
      length(setdiff(safe_filter(deg_30m_KA, 'downregulated'),
                     safe_filter(deg_30m_EE, 'downregulated'))),
      length(setdiff(safe_filter(deg_6h_KA, 'upregulated'),
                     safe_filter(deg_6h_EE, 'upregulated'))),
      length(setdiff(safe_filter(deg_6h_KA, 'downregulated'),
                     safe_filter(deg_6h_EE, 'downregulated')))
    ),
    shared = c(
      length(intersect(safe_filter(deg_30m_EE, 'upregulated'),
                       safe_filter(deg_30m_KA, 'upregulated'))),
      length(intersect(safe_filter(deg_30m_EE, 'downregulated'),
                       safe_filter(deg_30m_KA, 'downregulated'))),
      length(intersect(safe_filter(deg_6h_EE, 'upregulated'),
                       safe_filter(deg_6h_KA, 'upregulated'))),
      length(intersect(safe_filter(deg_6h_EE, 'downregulated'),
                       safe_filter(deg_6h_KA, 'downregulated')))
    )
  )
  
  # if df_current_subclass is empty, skip it
  if (nrow(df_current_subclass) == 0) {
    next
  }
  
  print(df_current_subclass)
  df_all <- bind_rows(df_all, df_current_subclass)
}


df_all_processed <- df_all |>
  pivot_longer(
    cols = c(EE_specific, shared, KA_specific),
    names_to = "group",
    values_to = "count"
  ) |>
  mutate(
    # for downregulated genes, use negative counts
    count = ifelse(direction == "down", -count, count),
    time = factor(time, levels = c("30m", "6h")),
    subclass = factor(subclass, levels = (subclasses)),
    group = factor(group, levels = c("KA_specific", "shared", "EE_specific"))
  ) |> 
  print()

df_30m <- df_all_processed |> filter(time == "30m")
df_6h  <- df_all_processed |> filter(time == "6h")

gene_category_colors <- c(EE_specific = 'gray20', shared = 'gray50', KA_specific = 'gray80')

# 30m plot
ggplot(df_30m, aes(y = count, x = subclass, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = 'black') +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = gene_category_colors) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_blank(),
  )
si(1000, 500)

# 6h plot
ggplot(df_6h, aes(y = count, x = subclass, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = 'black') +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = gene_category_colors) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_blank(),
  )
si(1000,500)


# make UMAPs for each gigaclass ------------------------------------------------
subclass_colors <- LoadAllenColors('subclass')
supertype_colors <- LoadAllenColors('supertype')
cluster_colors <- LoadAllenColors('cluster')
activity_colors <- LoadActivityColors()

subclass_sets <- LoadSubclassesToUse(nuclei, as_gigaclasses = TRUE)
subclass_sets$excitatory <- subclass_sets$excitatory[-2] # rm CA2
subclass_sets$inhibitory <- subclass_sets$inhibitory[-c(2:6)] # rm low-count interneurons
subclass_sets$glia <- subclass_sets$glia[-1] # rm IMNs

# re-embed nuclei
for (gigatype in names(subclass_sets)) {
  nuclei_gigatype <- nuclei |> 
    subset(subclass_name %in% subclass_sets[[gigatype]]) |> 
    RunUMAP(dims = 1:41)
  p1 <- DimPlot(nuclei_gigatype, group.by = 'subclass_name') +
    scale_color_manual(values = subclass_colors) +
    theme_void() +
    theme(legend.position = 'none',
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_blank())
  p2 <- DimPlot(nuclei_gigatype, group.by = 'activity_condition') +
    scale_color_manual(values = activity_colors) +
    theme_void() +
    theme(legend.position = 'none',
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_blank())
  p3 <- DimPlot(nuclei_gigatype, group.by = 'supertype_name') +
    scale_color_manual(values = supertype_colors) +
    theme_void() +
    theme(legend.position = 'none',
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_blank())
  p4 <- DimPlot(nuclei_gigatype, group.by = 'cluster_name') +
    scale_color_manual(values = cluster_colors) +
    theme_void() +
    theme(legend.position = 'none',
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_blank())
  
  design <- "
  133
  233
  "
  print(p1+p2+p3 + plot_layout(design = design, widths = c(1,0.75,0.75)))
  si(1500,1000,'png')
}

# supertype composition by condition ------------------------------------------------
df_cluster_counts <- nuclei@meta.data |> 
  filter(subclass_name %in% LoadSubclassesToUse(nuclei)) |>
  group_by(activity_condition, subclass_name, supertype_name) |>
  count(cluster_name, name = 'cluster_count')

# excitatory
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$excitatory) |> 
  group_by(activity_condition, subclass_name, supertype_name) |> 
  summarize(supertype_count = sum(cluster_count)) |> 
ggplot() +
  aes(x = activity_condition, y = supertype_count, fill = supertype_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = supertype_colors)
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())

# inhibitory
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$inhibitory) |> 
  group_by(activity_condition, subclass_name, supertype_name) |> 
  summarize(supertype_count = sum(cluster_count)) |> 
ggplot() +
  aes(x = activity_condition, y = supertype_count, fill = supertype_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = supertype_colors) +
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())

# glia
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$glia) |>
  group_by(activity_condition, subclass_name, supertype_name) |> 
  summarize(supertype_count = sum(cluster_count)) |> 
ggplot() +
  aes(x = activity_condition, y = supertype_count, fill = supertype_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = supertype_colors) +
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())


# cluster composition by condition ------------------------------------------------
# excitatory
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$excitatory) |> 
ggplot() +
  aes(x = activity_condition, y = cluster_count, fill = cluster_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = cluster_colors) +
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())

# inhibitory
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$inhibitory) |> 
ggplot() +
  aes(x = activity_condition, y = cluster_count, fill = cluster_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = cluster_colors) +
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())

# glia
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$glia) |>
ggplot() +
  aes(x = activity_condition, y = cluster_count, fill = cluster_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = cluster_colors) +
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())



# 0136 DG Glut_1 0137 DG Glut_2 0138 DG Glut_3 0139 DG Glut_4 
c("#F573FF",      "#009BCC",      "#FF7726",      "#FFBF26")

# 0503 DG Glut_1 0504 DG Glut_1 0505 DG Glut_2 0506 DG Glut_2 0507 DG Glut_2 0508 DG Glut_3 0509 DG Glut_3 0510 DG Glut_4 
c("#5CAECC",      "#FF0700",      "#CC603D",      "#FFD826",      "#99176A",      "#2E9964",      "#CC3DA3",      "#8EFF4D")


002 IT EP-CLA Glut  003 L5/6 IT TPE-ENT Glut        004 L6 IT CTX Glut        005 L5 IT CTX Glut      006 L4/5 IT CTX Glut 
"#1F665D"                 "#C400FF"                 "#C0FF4D"                 "#660F38"                 "#3DCCB7" 
007 L2/3 IT CTX Glut      008 L2/3 IT ENT Glut 009 L2/3 IT PIR-ENTl Glut     010 IT AON-TT-DP Glut     011 L2 IT ENT-po Glut 
"#E9530F"                 "#1F46CC"                 "#940099"                 "#FF004B"                 "#00FFC4" 
012 MEA Slc17a7 Glut      013 COAp Grxcr2 Glut    014 LA-BLA-BMA-PA Glut    015 ENTmv-PA-COAp Glut         016 CA1-ProS Glut 
"#456C99"                 "#61CC5C"                 "#BBFF73"                 "#B6CC7A"                 "#0F3466" 
017 CA3 Glut         018 L2 IT PPP-APr Glut      019 L2/3 IT PPP Glut      020 L2/3 IT RSP Glut       021 L4 RSP-ACA Glut 
"#C973FF"                 "#0F6632"                 "#925C99"                 "#999145"                 "#FFF226" 
022 L5 ET CTX Glut         023 SUB-ProS Glut           024 L5 PPP Glut        025 CA2-FC-IG Glut       028 L6b/CT ENT Glut 
"#99A2FF"                 "#FF9999"                 "#6B995C"                 "#FF7391"                 "#CC0098" 
029 L6b CTX Glut        030 L6 CT CTX Glut           031 CT SUB Glut        032 L5 NP CTX Glut           033 NP SUB Glut 
"#9EFF99"                 "#34661F"                 "#440066"                 "#7B7ACC"                 "#99CAFF" 
034 NP PPP Glut  035 OB Eomes Ms4a15 Glut           036 HPF CR Glut               037 DG Glut         038 DG-PIR Ex IMN 
"#B8FF26"                 "#7F1FCC"                 "#919900"                 "#CC2400"                 "#3D53CC" 
039 OB Meis2 Thsd7b Gaba    045 OB-STR-CTX Inh IMN              046 Vip Gaba             047 Sncg Gaba     048 RHP-COA Ndnf Gaba 
"#FF99F7"                 "#819945"                 "#663D47"                 "#FFE999"                 "#0003FF" 
049 Lamp5 Gaba       050 Lamp5 Lhx6 Gaba 051 Pvalb chandelier Gaba            052 Pvalb Gaba              053 Sst Gaba 
"#FF764D"                 "#26FF3E"                 "#2F992E"                 "#661F3B"                 "#99F2FF" 
054 STR Prox1 Lhx6 Gaba        056 Sst Chodl Gaba   072 LSX Sall3 Lmo1 Gaba        115 MS-SF Bsx Glut     208 SC Lef1 Otx2 Gaba 
"#c60f3a"                 "#EDA32D"                 "#00FFD0"                 "#1F6641"                 "#994F45" 
262 Pineal Crx Glut           318 Astro-NT NN           319 Astro-TE NN          320 Astro-OLF NN     321 Astroependymal NN 
"#0F6654"                 "#AD5CCC"                 "#3DCCB1"                 "#FF73CF"                 "#0833ce" 
322 Tanycyte NN                326 OPC NN              327 Oligo NN               330 VLMC NN               331 Peri NN 
"#FF26A1"                 "#FF26CB"                 "#99FFBC"                 "#653D66"                 "#82992E" 
333 Endo NN          334 Microglia NN                335 BAM NN 
"#994567"                 "#CC1F4E"                 "#66493D" 
