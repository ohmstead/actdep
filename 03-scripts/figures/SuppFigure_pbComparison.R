library(tidyverse)
library(ggpp)

source('03-scripts/R/seq_functions.R')
activity_colors <- LoadActivityColors()
subclass_colors <- LoadAllenColors()

# rename colors for plotting
contrast_list <- c('EE30m_vs_SE', 'EE6h_vs_SE', 'KA30m_vs_SE', 'KA6h_vs_SE')
names(activity_colors[c(2,3,4,6)]) <- contrast_list

subclass_list <- LoadSubclassesToUse()[c(1,4:9,11,17:20)]

# read jaccard data
df <- read.csv("04-analysis/reviewer_response/0_method_concordance.csv") |> 
  filter(contrast %in% contrast_list) |> 
  filter(subclass %in% subclass_list)


# jaccard plot ----
p1 <- ggplot(df) +
  aes(x = subclass, y = jaccard_DESeq2_voom, fill = contrast) +
  geom_hline(yintercept = seq(0,1,by=0.25), linetype = "dashed", color = "grey80") +
  geom_bar(stat = 'identity', position = position_dodge(0.9), width = 0.9) +
  geom_point(shape = 21, position = position_dodge(0.9), size = 5) +
  scale_fill_manual(values = activity_colors) +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
  )


# table plots ----
dfm <- LoadMeta() |> 
  filter(subclass_name %in% subclass_list, activity_condition != 'SE')

p2 <- dfm |> 
  group_by(subclass_name, activity_condition) |> 
  summarize(n = n()) |> 
  mutate(n = as.character(n)) |> 
ggplot() +
  aes(x = subclass_name, y = 1, label = n, color = activity_condition, group = activity_condition) +
  geom_text(angle = 90, position = position_dodge(width = 0.8), size = 4) +
  geom_tile(fill = 'grey80', color = 'white') +
  scale_color_manual(values = LoadActivityColors()) +
  scale_y_continuous(limits = c(0.98, 1.02), expand = c(0, 0)) +
  labs(y = 'N cells') +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = 'none'
  )

p1 <- p1 + theme(plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5))
p2 <- p2 + theme(plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5))

(p1 / p2) + plot_layout(heights = c(4, 0.5, 1))
