# 加载必要的包
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(readr)
library(purrr)

# ----------------------------
# 1. 读取CSV文件（替换为实际路径）
# ----------------------------
data <- read_csv("E:/nanoparticl proteomics/Figures/pathway enhancement/Nanoparticles Protein.csv")  # 确保文件路径正确
np_genes <- na.omit(unique(data$NP))  # 提取NP列并去重

# ----------------------------
# 2. 基因名转换为Entrez ID
# ----------------------------
np_entrez <- bitr(np_genes, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Mm.eg.db,
                  drop = FALSE) %>%  
  na.omit()  # 移除未匹配的基因

# ----------------------------
# 3. 富集分析（不设置阈值）
# ----------------------------
run_enrichment_no_threshold <- function(gene_list, db_type) {
  if (db_type == "KEGG") {
    enrichKEGG(
      gene = gene_list,
      organism = "mmu",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,  # 不设置p值阈值
      qvalueCutoff = 1   # 不设置q值阈值
    )
  } else {
    enrichGO(
      gene = gene_list,
      OrgDb = org.Mm.eg.db,
      ont = db_type,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,  # 不设置p值阈值
      qvalueCutoff = 1,  # 不设置q值阈值
      readable = TRUE
    )
  }
}

# 并行运行富集分析
db_types <- c("KEGG", "BP", "CC", "MF")
enrich_results <- map(db_types, ~ run_enrichment_no_threshold(np_entrez$ENTREZID, .x))
names(enrich_results) <- db_types

# 合并结果
enrich_df <- map_dfr(enrich_results, ~ {
  if (!is.null(.x)) {
    as.data.frame(.x) %>% 
      mutate(Database = ifelse("ONTOLOGY" %in% colnames(.), ONTOLOGY, "KEGG"))
  }
}, .id = "DB_Source") %>% 
  mutate(
    log2_odds = log2((as.numeric(gsub("/.*", "", GeneRatio)) / 
                        as.numeric(gsub(".*/", "", GeneRatio))) /
                       (as.numeric(gsub("/.*", "", BgRatio)) / 
                          as.numeric(gsub(".*/", "", BgRatio)))),
    log10_p = -log10(pvalue)
  )

# ----------------------------
# 4. 输出合并后的表格
# ----------------------------
# 打印前几行查看
head(enrich_df)

# 保存为CSV文件
write_csv(enrich_df, "E:/nanoparticl proteomics/Figures/pathway enhancement/nanoparticles_all_pathways_results.csv")


# 加载必要的包
library(ggplot2)
library(dplyr)
library(ggrepel)
library(stringr)

# 读取数据并处理Description列
enrichment_data <- read.csv("E:/nanoparticl proteomics/Figures/pathway enhancement/nanoparticles_all_pathways_results-1.csv") %>%
  # 确保列名唯一且小写
  `colnames<-`(tolower(make.unique(names(.)))) %>%
  # 精确删除 "- Mus musculus (house mouse)"（包括前面的空格和后面的可能空格）
  mutate(description = str_remove(description, " - Mus musculus \\(house mouse\\)\\s*")) %>%
  # 计算 -log10(p-value) 和方向
  mutate(log10_p = -log10(p.adjust)) %>%
  mutate(direction = ifelse(log2_odds > 0, "Enriched", "Depleted"))

# 获取y轴的最大值，用于确定标签位置
y_max <- max(enrichment_data$log10_p, na.rm = TRUE)

# 计算标签的y位置（在扩展区域）
label_y_position <- y_max * 1.15

# 绘制火山图 - 增加所有字体大小
volcano_plot <- ggplot(enrichment_data, 
                       aes(x = log2_odds, 
                           y = log10_p, 
                           color = direction,
                           size = count)) +
  geom_point(alpha = 0.6) +
  
  # 添加竖直虚线（x=0）
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  
  # 右侧富集通路标注（Enriched Pathways）- 按 log10_p 选择最显著的通路
  geom_label_repel(
    data = filter(enrichment_data, 
                  log2_odds > 0) %>%
      mutate(category = case_when(
        grepl("lipoprotein", description, ignore.case = TRUE) ~ "lipoprotein",
        grepl("cell adhesion", description, ignore.case = TRUE) ~ "cell adhesion",
        grepl("vesicles", description, ignore.case = TRUE) ~ "vesicles",
        grepl("complement", description, ignore.case = TRUE) ~ "complement",
        grepl("extracellular matrix", description, ignore.case = TRUE) ~ "extracellular matrix",
        grepl("immune response", description, ignore.case = TRUE) ~ "immune response",
        grepl("blood coagulation", description, ignore.case = TRUE) ~ "blood coagulation",
        TRUE ~ "other"
      )) %>%
      filter(category != "other") %>%
      group_by(category) %>%
      arrange(desc(log10_p)) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      mutate(description = str_to_sentence(description)),
    aes(label = description),
    color = "darkgreen", 
    size = 10,
    box.padding = 0.5,
    point.padding = 0.3,
    nudge_x = 0.5,
    nudge_y = 0.5,
    max.overlaps = 20,
    segment.color = "darkgreen",
    direction = "both",
    min.segment.length = 0.2,
    force = 1
  ) +
  
  # 左侧耗尽通路标注（Depleted Pathways）
  geom_label_repel(
    data = filter(enrichment_data, 
                  grepl("glutamate receptor activity|ligand-gated calcium channel activity", 
                        description, ignore.case = TRUE) &
                    log2_odds < 0 ) %>%
      mutate(description = str_to_sentence(description)),
    aes(label = description),
    color = "black", 
    size = 10,
    box.padding = 0.5,
    point.padding = 0.3,
    nudge_x = -0.5,
    nudge_y = 0.5,
    max.overlaps = 20,
    segment.color = "black",
    direction = "both",
    min.segment.length = 0.2,
    force = 1
  ) +
  
  # 添加区域标注 - 左侧标注"Depleted Pathways"
  annotate("text", x = -Inf, 
           y = label_y_position, 
           label = "Depleted Pathways", 
           hjust = -0.1, 
           vjust = 0,
           size = 8, color = "#4DBBD5", fontface = "bold") +
  
  # 添加区域标注 - 右侧标注"Enriched Pathways"
  annotate("text", x = Inf, 
           y = label_y_position, 
           label = "Enriched Pathways", 
           hjust = 1.1, 
           vjust = 0,
           size = 8, color = "#F39B7F", fontface = "bold") +
  
  # 调整颜色和大小比例
  scale_color_manual(values = c("Enriched" = "#F39B7F", "Depleted" = "#4DBBD5")) +
  scale_size_continuous(range = c(3, 10)) +
  
  # 使用更大的基础字体大小并调整主题
  theme_minimal(base_size = 18) +
  labs(
    x = "Log2 Odds Ratio", 
    y = "-log10(adjusted p-value)"
  ) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    # 确保绘图区域外的内容可见
    plot.margin = unit(c(1, 1, 3, 1), "cm"),  # 增加上边距为3cm
    # 添加边框
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  guides(
    color = "none",
    size = guide_legend(title = "Gene Count")
  ) +
  # 扩展y轴上限，为标签留出足够空间
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # 增加到0.2
  # 调整坐标轴限制以包含标签
  coord_cartesian(ylim = c(0, y_max * 1.1), clip = "off")  # 增加上限到1.25倍

# 打印图形
print(volcano_plot)

# 保存为PDF（A4尺寸）并添加边框
ggsave("E:/nanoparticl proteomics/Figures/pathway enhancement/pathway_enrichment_volcano.pdf", 
       plot = volcano_plot, 
       device = "pdf", 
       width = 561, 
       height = 317, 
       units = "mm", 
       dpi = 300)

# 或者使用另一种方法添加边框 - 使用plot.background
volcano_plot_with_border <- volcano_plot + 
  theme(plot.background = element_rect(colour = "black", fill = NA, linewidth = 2))

# 打印带边框的图形
print(volcano_plot_with_border)

# 保存带边框的PDF
ggsave("E:/nanoparticl proteomics/Figures/pathway enhancement/pathway_enrichment_volcano_with_border.pdf", 
       plot = volcano_plot_with_border, 
       device = "pdf", 
       width = 561, 
       height = 317, 
       units = "mm", 
       dpi = 300)

# 统计通路总数及Enriched/Depleted数量
pathway_stats <- enrichment_data %>%
  summarise(
    Total_Pathways = n(),
    Enriched = sum(direction == "Enriched"),
    Depleted = sum(direction == "Depleted")
  )

# 打印统计结果
cat(sprintf(
  "Total pathways: %d\nEnriched pathways: %d\nDepleted pathways: %d",
  pathway_stats$Total_Pathways,
  pathway_stats$Enriched,
  pathway_stats$Depleted
))