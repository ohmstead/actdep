## ---- Fig1B
library(paletteer)
library(Seurat)

source("03-scripts/R/seq_functions.R")

# only load if not already in environment
if (!exists("nuclei")) {nuclei <- LoadDataset("Dec2024")}

# load colors  ----------------------------------------------------------------------
class_colors <- LoadAllenColors('class')
subclass_colors <- LoadAllenColors('subclass')
activity_colors <- LoadActivityColors()
ZT_colors <- LoadZTColors()
sex_colors <- LoadSexColors()

save_path <- "05-results/Figure1/raw_R_plots"

save_plots_if_requested <- function(plot_args_list) {
  if (!exists("SAVE_PLOTS") || !isTRUE(SAVE_PLOTS)) {
    return(invisible())
  }

  for (plot_args in plot_args_list) {
    do.call(ggsave, plot_args)
  }
}


# relevel, filter subclasses ----------------------------------------------------------------------
meta <- nuclei@meta.data

subclasses_to_use <- LoadSubclassesToUse(nuclei)

meta <- meta |> 
  filter(subclass_name %in% subclasses_to_use) |>
  mutate(count = n(), .by = subclass_name) |>
  mutate(subclass_name = fct_reorder(subclass_name, count, .desc = FALSE)) |> 
  select(-count)

# qc subclass ----------------------------------------------------------------------
my_theme <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle=30, hjust=1),
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
  labs(x = 'num cells') +
  my_theme +
  theme(
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
  )

# condition ratios
pActivity <- meta |> 
  group_by(subclass_name, activity_condition) |>
  summarise(n = n()) |> 
ggplot() +
    aes(y = subclass_name, x = n, fill = activity_condition) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_fill_manual(values = activity_colors) +
    scale_x_continuous(breaks = NULL, position = 'bottom') +
    labs(x = 'activity') +
    my_theme

# sex ratios
pSex <- meta |> 
  group_by(subclass_name, sex) |>
  summarise(n = n()) |> 
ggplot() +
    aes(y = subclass_name, x = n, fill = sex) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_x_continuous(breaks = NULL, position = 'bottom') +
    scale_fill_manual(values = sex_colors) +
    labs(x = 'sex') +
    my_theme

# sublibrary ratio
pSublib <- meta |> 
  group_by(subclass_name, sublibrary) |>
  summarise(n = n()) |> 
ggplot() +
    aes(y = subclass_name, x = n, fill = sublibrary) +
    geom_bar(stat = 'identity', position = 'fill') +
    paletteer::scale_fill_paletteer_d("ltc::reading") + 
    scale_x_continuous(breaks = NULL, position = 'bottom') +
    labs(x = 'sublib') +
    my_theme

# ZT ratio
pZT <- meta |> 
  group_by(subclass_name, ZT) |>
  summarise(n = n()) |> 
ggplot() +
    aes(y = subclass_name, x = n, fill = ZT) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_x_continuous(breaks = NULL, position = 'bottom') +
    scale_fill_manual(values = ZT_colors) +
    labs(x = 'ZT') +
    my_theme

# UMI plot
pUMI <- ggplot(meta) +
  aes(y = subclass_name, x = nCount_RNA) +
  geom_violin(fill = 'gray50', width = 2) +
  geom_violin(aes(fill = subclass_name), width = 2) +
  scale_fill_manual(values = subclass_colors) +
  geom_boxplot(fill = 'gray80', width = 0.3, outliers = FALSE) +
  scale_x_log10(position = 'bottom') +
  labs(x='UMI') +
  my_theme

# gene plot
pGenes <- ggplot(meta) +
  aes(y = subclass_name, x = nFeature_RNA) +
  geom_violin(aes(fill = subclass_name), width = 2) +
  scale_fill_manual(values = subclass_colors) +
  geom_boxplot(fill = 'gray80', width = 0.3, outliers = FALSE) +
  scale_x_continuous(position = 'bottom') +
  labs(x = 'genes') +
  my_theme

# pct.mt plot
pMt <- ggplot(meta) +
  aes(y = subclass_name, x = percent.mt) +
  geom_violin(aes(fill = subclass_name), width = 2) +
  scale_fill_manual(values = subclass_colors) +
  geom_boxplot(fill = 'gray80', width = 0.3, outliers = FALSE) +
  scale_x_continuous(position = 'bottom') +
  labs(y = '', x = '% mt reads') +
  my_theme

p <- pNum + pActivity + pZT + pSex + pSublib + plot_layout(widths = c(30, 3, 3, 3, 3))
print(p)

# subclass QC for figure
save_plots_if_requested(list(
  list(
    plot = p,
    filename = 'quality_by_subclass_short.png',
    path = save_path,
    width = 15,
    height = 10,
    dpi = 300
  ),
  list(
    plot = p * LoadBarebonesTheme(ticks = 'x'),
    filename = 'quality_by_subclass_short.svg',
    path = save_path,
    width = 15,
    height = 10,
    dpi = 300
  )
))

p <- pNum + pUMI + pGenes + pMt + pActivity + pSex + pSublib + pZT + 
  plot_layout(widths = c(10, 10, 10, 10, 3, 3, 3, 3, 3))
print(p)
save_plots_if_requested(list(
  list(
    plot = p,
    filename = 'quality_by_subclass.png',
    path = save_path,
    width = 27,
    height = 9,
    dpi = 600
  )
))


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
# pReads <- ggplot(meta) +
#   aes(y = sample, x = mread_count, fill = activity_condition) +
#   geom_violin(width = 2) +
#   geom_boxplot(width = 0.2, outliers = FALSE) +
#   scale_fill_manual(values = activity_colors) +
#   scale_x_log10(position = 'top') +
#   scale_y_discrete(limits = rev) +
#   labs(y = '', x = 'read count') +
#   my_theme +
#   theme(
#     axis.text.x.top = element_text(angle = 270, hjust = 1, vjust = 0.5),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#   )

# UMI plot
pUMI <- ggplot(meta) +
  aes(y = sample, x = nCount_RNA, fill = activity_condition) +
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
  aes(y = sample, x = nFeature_RNA, fill = activity_condition) +
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
p <- pNum + pSex + pUMI + pGenes + pMt + plot_layout(widths = c(5, 3, 25, 25, 10), nrow = 1)
print(p)

save_plots_if_requested(list(
  list(
    plot = p,
    filename = 'quality_by_sample.png',
    path = save_path,
    width = 27,
    height = 14
  )
))


# plot sample N ----------------------------------------------------------------------
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


# # distro for reads ----------------------------------------------------------------------
# # plot reads by class
# pReads <- ggplot(meta) +
#   aes(x = sample, y = mread_count, fill = class_name) +
#   geom_violin() +
#   geom_boxplot(width = 0.2, outliers = FALSE) +
#   scale_y_log10(labels = scales::label_comma()) +
#   scale_x_discrete() +
#   facet_wrap(~class_name, 
#              nrow = length(unique(meta$class_name)),
#              scales = 'free_y') +
#   my_theme + 
#   theme(strip.text = element_text(size = 25),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank())
# 
# pReads <- pReads / pN + plot_layout(heights = c(10, 1))
# print(pReads)
# ggsave('class_reads.png', 
#        path = save_path,
#        width = 9, height = 18, dpi = 600)


# distro for UMIs ----------------------------------------------------------------------
# plot reads by class
pUMI <- ggplot(meta) +
  aes(x = sample, y = nCount_RNA, fill = class_name) +
  geom_violin() +
  geom_boxplot(width = 0.2, outliers = FALSE) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_x_discrete() +
  facet_wrap(~class_name, 
             nrow = length(unique(meta$class_name)),
             scales = 'free_y') +
  my_theme + 
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pUMI <- pUMI / pN + plot_layout(heights = c(10, 1))
print(pUMI)


# distro for genes ----------------------------------------------------------------------
# plot reads by class
pGenes <- ggplot(meta) +
  aes(x = sample, y = nFeature_RNA, fill = class_name) +
  geom_violin() +
  geom_boxplot(width = 0.2, outliers = FALSE) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_x_discrete() +
  facet_wrap(~class_name, 
             nrow = length(unique(meta$class_name)),
             scales = 'free_y') +
  my_theme + 
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pGenes <- pGenes / pN + plot_layout(heights = c(10, 1))
print(pGenes)
save_plots_if_requested(list(
  list(
    plot = pUMI,
    filename = 'class_UMI.png',
    path = save_path,
    width = 9,
    height = 18,
    dpi = 600
  ),
  list(
    plot = pGenes,
    filename = 'class_genes.png',
    path = save_path,
    width = 9,
    height = 18,
    dpi = 600
  )
))
## ----