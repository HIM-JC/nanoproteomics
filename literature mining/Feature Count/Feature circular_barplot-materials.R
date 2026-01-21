# 读取CSV文件
data <- read.csv("E:/nanoparticl proteomics/Figures/data mining/feature/morphology-count.csv")

# 计算百分比
data$percentage <- round(data$count / sum(data$count) * 100, 1)

# 使用您指定的颜色
colors <- c("#64B5F6", "#81C784", "#FFD54F", "#BA68C8", 
            "#BBDEFB", "#C8E6C9", "#FFECB3", "#E1BEE7",
            "#FFB74D", "#80DEEA")

# 加载ggplot2包
library(ggplot2)
library(dplyr)

# 将组名转换为首字母大写格式
data$group <- tools::toTitleCase(tolower(data$group))

# 创建PDF文件
pdf("E:/nanoparticl proteomics/Figures/data mining/feature/circular_barplot_morphology.pdf", width = 14, height = 14)

# 按计数排序数据
data <- data %>% 
  arrange(desc(count)) %>%
  mutate(
    id = seq(1, n()),
    angle = 90 - 360 * (id - 0.5) / n(),
    hjust = ifelse(angle < -90, 1, 0),
    angle = ifelse(angle < -90, angle + 180, angle),
    # 计算文字位置 - 放在条形高度的中间位置
    label_y = count * 0.5
  )

# 绘制增强版圆形条形图
p <- ggplot(data, aes(x = factor(id), y = count, fill = group)) +
  geom_bar(stat = "identity", alpha = 0.9, width = 0.9, color = "white", size = 0.5) +
  # 注释掉或删除下面这行，去掉条形中间的线段
  # geom_segment(aes(y = 0, yend = count, x = id, xend = id), 
  #              color = "gray70", alpha = 0.4, size = 0.3) +
  coord_polar(start = 0) +
  ylim(-max(data$count)*0.3, max(data$count)*1.2) +
  scale_fill_manual(values = colors[1:nrow(data)]) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.text = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 13, face = "bold"),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.key.size = unit(1, "cm")
  ) +
  # 修改这里：将百分比文字放在条形内部
  geom_text(
    aes(y = label_y, 
        label = paste0(percentage, "%"), 
        angle = angle, hjust = hjust), 
    color = "white",  # 使用白色文字，在深色条形上更清晰
    size = 4.5,
    fontface = "bold"
  ) +
  labs(
    fill = "Material"
  ) +
  theme(
    plot.title = element_text(
      color = "black", 
      hjust = 0.5, 
      size = 22,
      face = "bold",
      margin = margin(b = 10)
    ),
    plot.subtitle = element_text(
      color = "black", 
      hjust = 0.5, 
      size = 16,
      face = "italic",
      margin = margin(b = 20)
    )
  ) +
  # 添加中心装饰文字
  annotate("text", x = 0, y = -max(data$count)*0.2, 
           label = "Materials", 
           color = "black", size = 5, fontface = "bold", 
           alpha = 0.8, lineheight = 0.9)

# 打印图形
print(p)

# 关闭PDF设备
dev.off()

cat("\n增强版圆形条形图已保存到: E:/nanoparticl proteomics/Figures/data mining/feature/circular_barplot_enhanced.pdf\n")
cat("百分比文字现在显示在条形内部\n")
cat("已移除条形中间的线段\n")