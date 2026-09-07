## ---- SuppFig2
source("03-scripts/R/seq_functions.R")
library(paletteer)

if (!exists("nuclei")) {nuclei <- LoadDataset("Dec2024")}
class_colors <- LoadAllenColors('class')
subclass_colors <- LoadAllenColors('subclass')
activity_colors <- LoadActivityColors()
ZT_colors <- LoadZTColors()
sex_colors <- LoadSexColors()

subclasses_to_use <- LoadSubclassesToUse(nuclei)

meta <- nuclei@meta.data |> 
  filter(subclass_name %in% subclasses_to_use) |>
  mutate(count = n(), .by = subclass_name) |>
  mutate(subclass_name = fct_reorder(subclass_name, count, .desc = FALSE)) |> 
  select(-count)
      
      
pZT <- meta |> 
  group_by(sample) |> 
  summarize(
    ZT = dplyr::first(ZT),
    x = 1
  ) |> 
ggplot() +
  aes(y = sample, x = x, color = ZT) +
  geom_point(shape=16, size = 3) +
  scale_color_manual(values = ZT_colors) +
  theme_void() +
  theme(
    legend.position = 'none',
    plot.margin = margin(t=5, r=0, b=5, l=0),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8, hjust = 1),
  )


pN <- meta |> 
  group_by(sample) |> 
  summarize(n = n(),
            ZT = dplyr::first(ZT),
            activity_condition = dplyr::first(activity_condition)) |> 
ggplot() +
  aes(y = sample, x = n, fill = activity_condition) +
  # aes(y = sample, x = n) +
  geom_col() +
  scale_fill_manual(values = activity_colors) +
  scale_x_log10(breaks = c(100, 1000, 3000),
                labels = scales::label_comma()) +
  labs(x = 'Nuclei') +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 8, angle = 60, hjust = 1),
    legend.position = 'none'
  )

pUMI <- ggplot(meta) +
  aes(y = sample, x = nCount_RNA, fill = class_name) +
  geom_violin() +
  geom_boxplot(width = 0.2, outliers = FALSE) +
  scale_fill_manual(values = class_colors) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_discrete() +
  facet_wrap(~class_name, 
             ncol = length(unique(meta$class_name)),
             scales = 'free_y') +
  labs(title='UMI', x='UMI count') +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    strip.text = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 8, angle = 60, hjust = 1),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = 'none'
  )

pUMI_alone <- pZT + plot_spacer() + pN + pUMI + plot_layout(widths = c(0.3, -0.28, 1.5, 10), axes='collect')
  

# pGene --------------------------------------------------
pGene <- ggplot(meta) +
  aes(y = sample, x = nFeature_RNA, fill = class_name) +
  geom_violin() +
  geom_boxplot(width = 0.2, outliers = FALSE) +
  scale_fill_manual(values = class_colors) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_discrete() +
  facet_wrap(
    ~class_name, 
    ncol = length(unique(meta$class_name)),
    scales = 'free_x') +
  labs(title='Gene', x='Gene count') +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    strip.text = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 8, angle = 60, hjust = 1),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = 'none'
  )


pGene_alone <- pZT + plot_spacer() + pN + pGene + plot_layout(widths = c(0.3, -0.28, 1.5, 10), axes='collect')

p_layout <- c("
ABCD
EFGH
")
p <- pZT + plot_spacer() + pN + pUMI +
     pZT + plot_spacer() + pN + pGene +
     plot_layout(design = p_layout, widths = c(0.3, -0.28, 1.5, 10))
# p

pUMI_alone
pGene_alone

if (exists("SAVE_PLOTS") & SAVE_PLOTS==T) {
  # UMIs are Figure 1-supp 3, genes are Figure 1-supp 4; each supplement owns
  # its own results directory.
  save_path_UMI  <- "05-results/Figure1_SuppFig3/raw_R_plots"
  save_path_gene <- "05-results/Figure1_SuppFig4/raw_R_plots"

  # png
  ggsave(file.path(save_path_UMI,  "Figure1_SuppFig3.png"), pUMI_alone,  width = 6.5, height = 7, dpi = 300)
  ggsave(file.path(save_path_gene, "Figure1_SuppFig4.png"), pGene_alone, width = 6.5, height = 7, dpi = 300)

  # svg
  ggsave(file.path(save_path_UMI,  "Figure1_SuppFig3.svg"), pUMI_alone,  width = 6.5, height = 7, dpi = 300)
  ggsave(file.path(save_path_gene, "Figure1_SuppFig4.svg"), pGene_alone, width = 6.5, height = 7, dpi = 300)
}
## ----