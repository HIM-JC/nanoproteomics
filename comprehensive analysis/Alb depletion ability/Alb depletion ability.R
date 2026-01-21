# 加载必要的包
library(tidyverse)
library(ggplot2)
library(patchwork)

# 读取数据
data <- read.csv("E:/nanoparticl proteomics/Figures/Protein depletion/Gen name-PC-RAP-intensity-Materials order.csv", 
                 header = TRUE, 
                 row.names = 1,
                 check.names = FALSE)

# 材料分类函数
get_material <- function(sample_name) {
  case_when(
    str_detect(sample_name, "^AgNPs") ~ "Silver",
    str_detect(sample_name, "^AuNPs|^GNR") ~ "Gold",
    str_detect(sample_name, "^GO|^MWCNT|^SWCNT") ~ "Carbon",
    str_detect(sample_name, "^LNP[1-5]") ~ "LNP",
    str_detect(sample_name, "^CeO2|^Fe3O4|^MnO2|^ZnO") ~ "MOx",
    str_detect(sample_name, "^MIL|^PB|^UiO|^ZIF") ~ "MOF",
    str_detect(sample_name, "^MSN|^SiO2|^Si\\.NH2|^Si-NH2") ~ "Silica",
    str_detect(sample_name, "^PDA|^PLGA|^PS") ~ "Polymer",
    str_detect(sample_name, "^LDH|^TiCNS") ~ "Others",
    str_detect(sample_name, "^Serum") ~ "Serum",
    TRUE ~ "Unknown"
  )
}

# 定义统一的组顺序（基于原始数据列名）
original_order <- colnames(data)
group_order <- str_replace(original_order, "\\.\\d+$", "") %>% unique()
material_order <- c("Silver", "Gold", "Carbon", "LNP", "MOx", "MOF", 
                    "Silica", "Polymer", "Others", "Serum")

# 颜色映射（确保所有材料都有定义）
material_colors <- c(
  "Silver" = "#1f77b4",
  "Gold" = "#ff7f0e",
  "Carbon" = "#2ca02c",
  "LNP" = "#d62728",
  "MOx" = "#9467bd",
  "MOF" = "#8c564b",
  "Silica" = "#e377c2",
  "Polymer" = "#7f7f7f",
  "Others" = "#bcbd22",
  "Serum" = "#17becf",
  "Unknown" = "#333333"
)

# 1. 计算每列非零数值的个数
nonzero_counts <- data %>%
  summarise(across(everything(), ~sum(.x > 0))) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "Nonzero_Count") %>%
  mutate(
    Group = str_replace(Sample, "\\.\\d+$", ""),
    Group = factor(Group, levels = group_order),  # 保持组顺序一致
    Material = factor(get_material(Group), levels = material_order),  # 保持材料顺序一致
    Replicate = str_extract(Sample, "\\d+$")
  )

# 非零计数统计量
nonzero_stats <- nonzero_counts %>%
  group_by(Group, Material) %>%
  summarise(
    Mean_Count = mean(Nonzero_Count),
    SD = sd(Nonzero_Count),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    SE = SD / sqrt(n),
    CI = qt(0.975, df = n-1) * SE
  )

# 2. 计算Alb在每个样本中的排名
rank_data <- data %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
  group_by(Sample) %>%
  mutate(Rank = rank(-Expression)) %>%
  filter(Gene == "Alb") %>%
  mutate(
    Group = str_replace(Sample, "\\.\\d+$", ""),
    Group = factor(Group, levels = group_order),  # 保持组顺序一致
    Material = factor(get_material(Group), levels = material_order),  # 保持材料顺序一致
    Replicate = str_extract(Sample, "\\d+$")
  )

# 排名统计量
stats_data <- rank_data %>%
  group_by(Group, Material) %>%
  summarise(
    Mean_Rank = mean(Rank),
    SD = sd(Rank),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    SE = SD / sqrt(n),
    CI = qt(0.975, df = n-1) * SE
  )

# 3. 直接提取Alb的相对百分比含量
alb_data <- data %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  filter(Gene == "Alb") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Percentage") %>%
  mutate(
    Group = str_replace(Sample, "\\.\\d+$", ""),
    Group = factor(Group, levels = group_order),
    Material = factor(get_material(Group), levels = material_order),
    Replicate = str_extract(Sample, "\\d+$")
  )

# Alb百分比统计量
alb_stats <- alb_data %>%
  group_by(Group, Material) %>%
  summarise(
    Mean_Percentage = mean(Percentage),
    SD = sd(Percentage),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    SE = SD / sqrt(n),
    CI = qt(0.975, df = n-1) * SE
  )

# 创建材料标签数据 - 为所有图使用，排除Serum
material_labels <- alb_stats %>%
  distinct(Group, Material) %>%
  filter(Material != "Serum") %>%
  mutate(Group_index = as.numeric(Group)) %>%
  group_by(Material) %>%
  summarise(
    xmin = min(Group_index) - 0.4,
    xmax = max(Group_index) + 0.4,
    xmid = mean(c(xmin, xmax)),
    .groups = "drop"
  )

# 创建三个图形
# 1. 非零计数图（最上方）- 添加浅色背景条
nonzero_plot <- ggplot(nonzero_stats, aes(x = Group, y = Mean_Count)) +
  # 添加材料背景条（排除Serum）
  geom_rect(data = material_labels, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Material),
            alpha = 0.2, inherit.aes = FALSE) +
  geom_point(aes(color = Material), size = 2) +
  geom_errorbar(aes(ymin = Mean_Count - SE, ymax = Mean_Count + SE, color = Material), 
                width = 0.2, linewidth = 0.8) +
  geom_jitter(data = nonzero_counts, aes(y = Nonzero_Count, color = Material),
              position = position_jitter(width = 0.2), size = 2, alpha = 0.6) +
  scale_color_manual(values = material_colors) +
  scale_fill_manual(values = material_colors) +
  labs(y = "Number of\ndetected proteins", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 20, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90"),
    text = element_text(face = "bold")
  )

# 2. Alb排名图（中间）- 添加浅色背景条
rank_plot <- ggplot(stats_data, aes(x = Group, y = Mean_Rank)) +
  # 添加材料背景条（排除Serum）
  geom_rect(data = material_labels, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Material),
            alpha = 0.2, inherit.aes = FALSE) +
  geom_segment(aes(xend = Group, yend = 0, color = Material), 
               size = 1.2, show.legend = FALSE) +
  geom_point(aes(color = Material), size = 3) +
  scale_color_manual(values = material_colors) +
  scale_fill_manual(values = material_colors) +
  scale_y_reverse() +
  labs(y = "Albumin\nrank", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 20, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90"),
    text = element_text(face = "bold")
  )

# 3. Alb表达量图（最下方）- 添加浅色背景条和材料名称
expression_plot <- ggplot(alb_stats, aes(x = Group, y = Mean_Percentage)) +
  # 添加材料背景条（排除Serum）
  geom_rect(data = material_labels, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Material),
            alpha = 0.2, inherit.aes = FALSE) +
  geom_col(aes(fill = Material), width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_Percentage - SE, ymax = Mean_Percentage + SE), 
                width = 0.2, color = "black", linewidth = 0.8) +
  geom_jitter(data = alb_data, aes(y = Percentage, color = Material),
              position = position_jitter(width = 0.2), size = 2, alpha = 0.7) +
  scale_fill_manual(values = material_colors) +
  scale_color_manual(values = material_colors) +
  labs(y = "Albumin\nnormalized abundance", x = "") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 20),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey90"),
    text = element_text(face = "bold")
  ) +
  # 添加材料名称标签（排除Serum）- 调整y位置到更低的负值
  geom_text(data = material_labels, 
            aes(x = xmid, y = -max(alb_stats$Mean_Percentage)*0.08,  # 进一步降低y位置
                label = Material),
            vjust = 1, size = 6, fontface = "bold") +
  # 扩展y轴范围以确保标签可见
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))  # 底部扩展更多空间

# 组合图形
combined_plot <- nonzero_plot / rank_plot / expression_plot +
  plot_layout(heights = c(1, 1, 1.8))  # 增加第三个图的高度比例

# 显示图形
print(combined_plot)

# 保存结果 - 增加高度以确保标签完全显示
ggsave("E:/nanoparticl proteomics/Figures/Protein depletion/Alb_plot.pdf", combined_plot, 
       width = 25, height = 10, dpi = 300)  # 进一步增加高度到12