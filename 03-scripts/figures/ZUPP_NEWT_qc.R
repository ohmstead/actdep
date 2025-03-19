library(tidyverse)
library(ggplot2)
library(patchwork)
library(paletteer)
library(Seurat)
library(SeuratDisk)

source("03-scripts/R/seq_functions.R")

# load colors  ----------------------------------------------------------------------
activity_colors <- LoadActivityColors('Dec2024')
subclass_colors <- LoadAllenColors('subclass')
ZT_colors <- LoadZTColors(5)
sex_colors <- LoadSexColors()


# load nuclei ----------------------------------------------------------------------
nuclei <- LoadDataset("Dec2024")


# plot Xist expression ----------------------------------------------------------------------
# let's ascertain the sex of biological samples based on Xist expression
nuclei <- nuclei |> 
  AddMetaData(metadata = nuclei@assays$RNA$counts['Xist',], col = 'Xist')

meta <- nuclei@meta.data

# plot Xist expression in each sample
ggplot(meta) +
  geom_histogram(aes(x = Xist, fill = activity_condition)) +
  labs(title = 'Xist expression') +
  facet_wrap(~sample, nrow = 4) +
  scale_fill_manual(values = activity_colors)


# relevel, filter subclasses ----------------------------------------------------------------------
meta <- nuclei@meta.data

subclasses_to_use <- LoadSubclassesToUse(nuclei)

meta <- meta |> 
  filter(subclass_name %in% subclasses_to_use) |>
  mutate(count = n(), .by = subclass_name) |>
  mutate(subclass_name = fct_reorder(subclass_name, count, .desc = FALSE)) |> 
  select(-count)

# plot ----------------------------------------------------------------------
my_theme <- theme(
  # text = element_text(family = 'Mono'),
  legend.position = 'none',
  plot.margin = margin(0, 0, 0, 0),
  panel.background = element_blank(),
  panel.grid.major.x = element_line(color = 'gray90'),
  panel.grid.minor.x = element_line(color = 'gray93'),
  panel.spacing = unit(0, "lines"),
)

# num cells
pNum <- ggplot(meta) +
  aes(y = subclass_name) +
  geom_bar(aes(fill = subclass_name)) + 
  scale_fill_manual(values = subclass_colors) +
  scale_x_log10(breaks = c(100, 1000, 10000), position = 'bottom') +
  labs(x = '', y = '') +
  my_theme +
  theme(
    axis.text.x.top = element_text(angle = 270, hjust = 1, vjust = 0.5),
    # axis.text.y = element_text(size = 18),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
  )

# condition ratios
pActivity <- meta |> 
  group_by(subclass_name, activity_condition) |>
  summarise(n = n()) |> 
ggplot() +
    aes(y = subclass_name, x = n, fill = activity_condition) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_fill_manual(values = activity_colors) +
    scale_x_continuous(breaks = NULL, position = 'top') +
    labs(y = '', x = '') +
    my_theme +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 270, hjust = 1),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
    )

# sex ratios
pSex <- meta |> 
  group_by(subclass_name, sex) |>
  summarise(n = n()) |> 
ggplot() +
    aes(y = subclass_name, x = n, fill = sex) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_x_continuous(breaks = NULL, position = 'top') +
    scale_fill_manual(values = sex_colors) +
    labs(y = '', x = '') +
    my_theme +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 270, hjust = 1),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
    )

# sublibrary ratio
pSublib <- meta |> 
  group_by(subclass_name, sublibrary) |>
  summarise(n = n()) |> 
ggplot() +
    aes(y = subclass_name, x = n, fill = sublibrary) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_x_continuous(breaks = NULL, position = 'top') +
    scale_fill_paletteer_d("ltc::reading") + 
    labs(y = '', x = '') +
    my_theme +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 270, hjust = 1),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
    )

# ZT ratio
pZT <- meta |> 
  group_by(subclass_name, ZT) |>
  summarise(n = n()) |> 
ggplot() +
    aes(y = subclass_name, x = n, fill = ZT) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_x_continuous(breaks = NULL, position = 'top') +
    scale_fill_manual(values = ZT_colors) +
    labs(y = '', x = 'ZT') +
    my_theme +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 270, hjust = 1),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
    )

# UMI plot
pUMI <- ggplot(meta) +
  aes(y = subclass_name, x = nCount_RNA) +
  geom_violin(fill = 'gray50', width = 2) +
  geom_violin(aes(fill = subclass_name), width = 2) +
  scale_fill_manual(values = subclass_colors) +
  geom_boxplot(fill = 'gray80', width = 0.3, outliers = FALSE) +
  scale_x_log10(position = 'top') +
  labs(y = '', x = 'UMI count') +
  my_theme +
  theme(
    axis.text.x.top = element_text(angle = 270, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  )

# gene plot
pGenes <- ggplot(meta) +
  aes(y = subclass_name, x = nFeature_RNA) +
  geom_violin(aes(fill = subclass_name), width = 2) +
  scale_fill_manual(values = subclass_colors) +
  geom_boxplot(fill = 'gray80', width = 0.3, outliers = FALSE) +
  scale_x_continuous(position = 'top') +
  labs(y = '', x = 'gene count') +
  my_theme +
  theme(
    axis.text.x.top = element_text(angle = 270, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  )

# pct.mt plot
pMt <- ggplot(meta) +
  aes(y = subclass_name, x = percent.mt) +
  geom_violin(aes(fill = subclass_name), width = 2) +
  scale_fill_manual(values = subclass_colors) +
  geom_boxplot(fill = 'gray80', width = 0.3, outliers = FALSE) +
  scale_x_continuous(position = 'top') +
  labs(y = '', x = '% mt reads') +
  my_theme +
  theme(
    axis.text.x.top = element_text(angle = 270, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  )

pNum + pActivity + pSex + pSublib + plot_layout(widths = c(30, 3, 3, 3))
ggsave('quality_by_subclass_short.png', path = "05-results/figure_supp_qc/raw_R_plots", width = 15, height = 10, dpi = 600)
ggsave('quality_by_subclass_short.svg', path = "05-results/figure_supp_qc/raw_R_plots", width = 15, height = 10)

pNum + pActivity + pSex + pSublib + pZT + pUMI + pGenes + pMt + plot_layout(widths = c(5, 3, 3, 3, 3, 25, 25, 25, 10))
ggsave('quality_by_subclass.png', path = "05-results/figure_supp_qc/raw_R_plots", width = 27, height = 9, dpi = 600)


## ---------------------------------------------------------------------------------------------------------------------------------------------
my_theme <- theme(
  legend.position = 'none',
  plot.margin = margin(0, 0, 0, 0),
  panel.background = element_blank(),
  panel.grid.major.x = element_line(color = 'gray90'),
  panel.grid.minor.x = element_line(color = 'gray93'),
  panel.spacing = unit(0, "lines"),
  axis.title.x = element_text(size = 25),
  axis.text.y = element_text(size = 18),
)

# num cells
pNum <- ggplot(meta) +
  aes(y = sample, fill = activity_condition) +
  geom_bar() + 
  scale_fill_manual(values = activity_colors) +
  scale_x_continuous(breaks = c(100, 500, 1000), position = 'top') +
  scale_y_discrete(limits = rev) +
  labs(x = 'cells', y = '') +
  my_theme +
  theme(
    axis.text.x.top = element_text(angle = 270, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 18)
  )

# sex ratios
pSex <- meta |> 
  group_by(sample, sex) |>
  summarise(n = n()) |> 
  ggplot() +
    aes(y = sample, x = n, fill = sex) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_fill_manual(values = sex_colors) +
    scale_x_continuous(breaks = NULL, position = 'top') +
    scale_y_discrete(limits = rev) +
    labs(y = '', x = 'sex') +
    my_theme +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 270, hjust = 1),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
    )

# mread plot
pReads <- ggplot(meta) +
  aes(y = sample, x = mread_count, fill = activity_condition) +
  geom_violin(width = 2) +
  geom_boxplot(width = 0.2, outliers = FALSE) +
  scale_fill_manual(values = activity_colors) +
  scale_x_log10(position = 'top') +
  scale_y_discrete(limits = rev) +
  labs(y = '', x = 'read count') +
  my_theme +
  theme(
    axis.text.x.top = element_text(angle = 270, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  )

# UMI plot
pUMI <- ggplot(meta) +
  aes(y = sample, x = tscp_count, fill = activity_condition) +
  geom_violin(width = 2) +
  geom_boxplot(width = 0.2, outliers = FALSE) +
  scale_fill_manual(values = activity_colors) +
  scale_x_log10(position = 'top') +
  scale_y_discrete(limits = rev) +
  labs(y = '', x = 'UMI count') +
  my_theme +
  theme(
    axis.text.x.top = element_text(angle = 270, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  )

# gene plot
pGenes <- ggplot(meta) +
  aes(y = sample, x = gene_count, fill = activity_condition) +
  geom_violin(width = 2) +
  geom_boxplot(width = 0.2, outliers = FALSE) +
  scale_fill_manual(values = activity_colors) +
  scale_x_log10(position = 'top') +
  scale_y_discrete(limits = rev) +
  labs(y = '', x = 'gene count') +
  my_theme +
  theme(
    axis.text.x.top = element_text(angle = 270, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  )

# pct.mt plot
pMt <- ggplot(meta) +
  aes(y = sample, x = percent.mt, fill = activity_condition) +
  geom_violin(width = 2) +
  geom_boxplot(width = 0.2, outliers = FALSE) +
  scale_fill_manual(values = activity_colors) +
  scale_x_continuous(position = 'top') +
  scale_y_discrete(limits = rev) +
  labs(y = '', x = '% mt reads') +
  my_theme +
  theme(
    axis.text.x.top = element_text(angle = 270, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  )

pNum + pSex + pReads + pUMI + pGenes + pMt + plot_layout(widths = c(5, 3, 25, 25, 25, 10), nrow = 1)

ggsave('quality_by_sample.png', path = "05-results/figure_supp_qc/raw_R_plots", width = 27, height = 14)
ggsave('quality_by_sample.svg', path = "05-results/figure_supp_qc/raw_R_plots", width = 27, height = 14)


# show disto for each class ----------------------------------------------------------------------
class_colors <- LoadAllenColors('class')

meta |> 
  filter(ZT == 'ZT4') |> 
ggplot() +
  aes(y = interaction(class_name, sample), x = gene_count, fill = class_name) +
  # geom_violin(width = 2) +
  geom_boxplot(width = 0.2, outliers = FALSE) +
  scale_x_log10(position = 'top') +
  scale_y_discrete(limits = rev) +
  labs(y = '', x = 'gene count') +
  my_theme +
  theme(
    axis.text.x.top = element_text(angle = 270, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


# distro for reads ----------------------------------------------------------------------
# plot reads by class
pReads <- ggplot(meta) +
  aes(x = sample, y = mread_count, fill = class_name) +
  geom_violin() +
  geom_boxplot(width = 0.2, outliers = FALSE) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_x_discrete() +
  facet_wrap(~class_name, 
             nrow = length(unique(meta$class_name)),
             scales = 'free_y') +
  my_theme + 
  theme(strip.text = element_text(size = 25),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# plot for n cells
pN <- meta |> 
  group_by(sample) |> 
  summarize(n = n(),
            ZT = first(ZT),
            activity_condition = first(activity_condition)) |> 
ggplot() +
  aes(x = sample, y = n) +
  geom_col() +
  scale_y_log10(breaks = c(100, 1000, 2000, 3000),
                labels = scales::label_comma()) +
  my_theme +
  theme(
    axis.text.x.top = element_text(angle = 270, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15, angle = 90, hjust = 1)
  )

pReads / pN + plot_layout(heights = c(10, 1))

ggsave('reads_by_class.png', path = "05-results/figure_supp_qc/raw_R_plots", width = 9, height = 18, dpi = 600)


# distro for UMIs ----------------------------------------------------------------------
# plot reads by class
pUMI <- ggplot(meta) +
  aes(x = sample, y = tscp_count, fill = class_name) +
  geom_violin() +
  geom_boxplot(width = 0.2, outliers = FALSE) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_x_discrete() +
  facet_wrap(~class_name, 
             nrow = length(unique(meta$class_name)),
             scales = 'free_y') +
  my_theme + 
  theme(strip.text = element_text(size = 25),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pUMI / pN + plot_layout(heights = c(10, 1))

ggsave('tscp_by_class.png', path = "05-results/figure_supp_qc/raw_R_plots", width = 9, height = 18, dpi = 600)


# distro for genes ----------------------------------------------------------------------
# plot reads by class
pGenes <- ggplot(meta) +
  aes(x = sample, y = gene_count, fill = class_name) +
  geom_violin() +
  geom_boxplot(width = 0.2, outliers = FALSE) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_x_discrete() +
  facet_wrap(~class_name, 
             nrow = length(unique(meta$class_name)),
             scales = 'free_y') +
  my_theme + 
  theme(strip.text = element_text(size = 25),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pGenes / pN + plot_layout(heights = c(10, 1))

ggsave('genes_by_class.png', path = "05-results/figure_supp_qc/raw_R_plots", width = 9, height = 18, dpi = 600)

