# 0. 安装和加载必要的包
# install.packages(c("uwot", "ks", "dplyr", "tidyr", "ggplot2", "raster", "RColorBrewer", "parallel"))
library(uwot)
library(ks)
library(dplyr)
library(tidyr)
library(ggplot2)
library(raster)
library(RColorBrewer)
library(parallel)

# 1. 加载原始数据
proteomics_data_original <- read.csv("E:/nanoparticl proteomics/Figures/Shape/Shape-Protein genes-on tube.csv", check.names = FALSE)

# 2. 使用调试找到的最佳参数（请替换为你的实际最佳参数）
best_params <- list(
  n_top_proteins_percentile = 0.90,           # 替换为你的最佳值
  single_label_fraction_threshold = 0.5,      # 替换为你的最佳值  
  min_value_multilabel_threshold = 0.25,      # 替换为你的最佳值
  no_class_threshold_sum_signal = 0.0001      # 替换为你的最佳值
)

# 3. 数据预处理函数
preprocess_data <- function(proteomics_data_original) {
  labels_original <- proteomics_data_original$Label
  protein_features_original <- proteomics_data_original[, -1]
  
  # 标准化
  protein_features_scaled <- as.data.frame(scale(protein_features_original))
  
  # 处理NA值
  if(any(is.na(protein_features_scaled))){
    for(i in 1:ncol(protein_features_scaled)){
      if(any(is.na(protein_features_scaled[,i]))) {
        mean_val <- mean(protein_features_scaled[,i], na.rm = TRUE)
        if (is.nan(mean_val)) mean_val <- 0
        protein_features_scaled[is.na(protein_features_scaled[,i]), i] <- mean_val
      }
    }
  }
  
  # 移除零方差列
  zero_var_cols <- which(apply(protein_features_scaled, 2, var, na.rm = TRUE) < 1e-6)
  if(length(zero_var_cols) > 0) {
    protein_features_scaled <- protein_features_scaled[, -zero_var_cols]
    cat(paste("移除了", length(zero_var_cols), "个零方差列\n"))
  }
  
  return(list(
    labels = labels_original,
    features = protein_features_scaled
  ))
}

# 4. 执行UMAP降维
run_umap <- function(protein_features_scaled, labels_original) {
  # 准备蛋白质表达谱用于UMAP (蛋白质 x 样本)
  protein_profiles_for_umap <- t(protein_features_scaled)
  if (!is.null(rownames(proteomics_data_original))) {
    colnames(protein_profiles_for_umap) <- rownames(proteomics_data_original)
  } else {
    colnames(protein_profiles_for_umap) <- paste0("Sample_", 1:ncol(protein_profiles_for_umap))
  }
  
  # 对蛋白质进行UMAP降维
  set.seed(123)  # 固定随机种子
  umap_proteins_result <- umap(protein_profiles_for_umap,
                               n_neighbors = min(15, nrow(protein_profiles_for_umap)-1, na.rm=TRUE),
                               min_dist = 0.05,
                               metric = "correlation",
                               n_components = 2,
                               verbose = FALSE,
                               n_threads = detectCores() - 1,
                               ret_model = TRUE)
  
  umap_proteins_df <- as.data.frame(umap_proteins_result$embedding)
  colnames(umap_proteins_df) <- c("Prot_UMAP1", "Prot_UMAP2")
  umap_proteins_df$ProteinID <- rownames(protein_profiles_for_umap)
  
  return(umap_proteins_df)
}

# 5. 使用最佳参数进行分类
run_classification_with_best_params <- function(protein_features_scaled, labels_original, umap_proteins_df, best_params) {
  
  categories <- unique(labels_original)
  protein_kde_list <- list()
  
  kde_grid_points <- 100
  umap_x_range <- range(umap_proteins_df$Prot_UMAP1, na.rm = TRUE)
  umap_y_range <- range(umap_proteins_df$Prot_UMAP2, na.rm = TRUE)
  
  min_proteins_for_kde <- 10
  
  # 从best_params中提取参数
  n_top_proteins_percentile <- best_params$n_top_proteins_percentile
  single_label_fraction_threshold <- best_params$single_label_fraction_threshold
  min_value_multilabel_threshold <- best_params$min_value_multilabel_threshold
  no_class_threshold_sum_signal <- best_params$no_class_threshold_sum_signal
  
  # KDE计算
  for (category in categories) {
    samples_indices_this_category <- which(labels_original == category)
    if (length(samples_indices_this_category) > 0) {
      avg_expression_this_category <- colMeans(protein_features_scaled[samples_indices_this_category, , drop = FALSE], na.rm = TRUE)
      threshold_expr <- quantile(avg_expression_this_category, probs = n_top_proteins_percentile, na.rm = TRUE)
      selected_proteins_names <- names(avg_expression_this_category[avg_expression_this_category >= threshold_expr])
      if (length(selected_proteins_names) >= min_proteins_for_kde) {
        protein_coords_for_kde <- umap_proteins_df[umap_proteins_df$ProteinID %in% selected_proteins_names, c("Prot_UMAP1", "Prot_UMAP2")]
        protein_coords_for_kde <- na.omit(protein_coords_for_kde)
        if (nrow(protein_coords_for_kde) >= min_proteins_for_kde) {
          fhat_prot <- kde(x = protein_coords_for_kde,
                           eval.points = expand.grid(seq(umap_x_range[1], umap_x_range[2], length.out = kde_grid_points),
                                                     seq(umap_y_range[1], umap_y_range[2], length.out = kde_grid_points)),
                           xmin = umap_x_range, xmax = umap_y_range)
          protein_kde_list[[as.character(category)]] <- fhat_prot
        }
      }
    }
  }
  
  # 像素级标签分配
  valid_kdes <- names(which(!sapply(protein_kde_list, is.null)))
  if(length(valid_kdes) == 0) {
    stop("没有有效的KDE模型")
  }
  
  first_valid_kde_name <- valid_kdes[1]
  ref_kde <- protein_kde_list[[first_valid_kde_name]]
  x_coords_grid <- unique(ref_kde$eval.points[,1])
  y_coords_grid <- unique(ref_kde$eval.points[,2])
  
  pixel_signal_contributions <- array(0, dim = c(length(y_coords_grid), length(x_coords_grid), length(categories)))
  for (i_cat in 1:length(categories)) {
    category_name <- categories[i_cat]
    if (!is.null(protein_kde_list[[as.character(category_name)]])) {
      current_kde_estimate <- protein_kde_list[[as.character(category_name)]]$estimate
      if(all(dim(current_kde_estimate) == c(length(x_coords_grid), length(y_coords_grid)))){
        pixel_signal_contributions[,,i_cat] <- t(current_kde_estimate)
      } else if (all(dim(current_kde_estimate) == c(length(y_coords_grid), length(x_coords_grid)))){
        pixel_signal_contributions[,,i_cat] <- current_kde_estimate
      }
    }
  }
  
  pixel_final_label <- matrix("Unassigned", nrow = length(y_coords_grid), ncol = length(x_coords_grid))
  for (r_idx in 1:length(y_coords_grid)) {
    for (c_idx in 1:length(x_coords_grid)) {
      signals_at_pixel <- pixel_signal_contributions[r_idx, c_idx, ]
      total_signal_at_pixel <- sum(signals_at_pixel, na.rm = TRUE)
      if (is.na(total_signal_at_pixel) || total_signal_at_pixel < no_class_threshold_sum_signal || total_signal_at_pixel == 0) {
        pixel_final_label[r_idx, c_idx] <- "BelowThreshold"
        next
      }
      relative_contributions <- signals_at_pixel / total_signal_at_pixel
      relative_contributions[is.na(relative_contributions) | is.infinite(relative_contributions)] <- 0
      max_contribution <- max(relative_contributions, na.rm = TRUE)
      if (max_contribution >= single_label_fraction_threshold) {
        assigned_category_index <- which.max(relative_contributions)
        pixel_final_label[r_idx, c_idx] <- categories[assigned_category_index]
      } else {
        multilabel_candidates_indices <- which(relative_contributions >= min_value_multilabel_threshold)
        if (length(multilabel_candidates_indices) > 0 && length(multilabel_candidates_indices) <= 1) {
          pixel_final_label[r_idx, c_idx] <- paste(sort(categories[multilabel_candidates_indices]), collapse = "_&_")
        } else if (length(multilabel_candidates_indices) > 1) {
          pixel_final_label[r_idx, c_idx] <- "Common"
        } else {
          pixel_final_label[r_idx, c_idx] <- "Ambiguous"
        }
      }
    }
  }
  
  # 蛋白质分类分配
  protein_assigned_categories_df <- data.frame(
    ProteinID = umap_proteins_df$ProteinID,
    Prot_UMAP1 = umap_proteins_df$Prot_UMAP1,
    Prot_UMAP2 = umap_proteins_df$Prot_UMAP2,
    Assigned_Category = NA_character_
  )
  
  for(i in 1:nrow(protein_assigned_categories_df)) {
    umap1 <- protein_assigned_categories_df$Prot_UMAP1[i]
    umap2 <- protein_assigned_categories_df$Prot_UMAP2[i]
    if(is.na(umap1) || is.na(umap2)) {
      protein_assigned_categories_df$Assigned_Category[i] <- NA
      next
    }
    closest_x_index <- which.min(abs(x_coords_grid - umap1))
    closest_y_index <- which.min(abs(y_coords_grid - umap2))
    protein_assigned_categories_df$Assigned_Category[i] <- pixel_final_label[closest_y_index, closest_x_index]
  }
  
  return(protein_assigned_categories_df)
}

# 6. 主运行流程
cat("开始运行分类流程（使用最佳参数）...\n")

# 数据预处理
preprocessed_data <- preprocess_data(proteomics_data_original)
labels_original <- preprocessed_data$labels
protein_features_scaled <- preprocessed_data$features

# UMAP降维
cat("执行UMAP降维...\n")
umap_proteins_df <- run_umap(protein_features_scaled, labels_original)

# 使用最佳参数进行分类
cat("使用最佳参数进行分类...\n")
final_result <- run_classification_with_best_params(
  protein_features_scaled, 
  labels_original, 
  umap_proteins_df, 
  best_params
)

# 7. 保存结果
final_categories <- split(final_result$ProteinID, final_result$Assigned_Category)
final_output <- do.call(rbind, lapply(names(final_categories), function(cat) {
  data.frame(Category = cat, ProteinID = final_categories[[cat]])
}))

write.csv(
  final_output,
  "E:/nanoparticl proteomics/Figures/Shape/Shape classification-Protein_Final.csv",
  row.names = FALSE
)

cat("分类结果已保存到: Shape classification-Protein_Final.csv\n")

# 8. 显示分类统计
cat("\n=== 分类结果统计 ===\n")
category_counts <- table(final_result$Assigned_Category)
print(category_counts)

# 9. 生成可视化图表
cat("生成可视化图表...\n")

# 自定义颜色
custom_colors <- c(
  'Spherical' = '#D55640',
  'Sheet' = '#479D88',
  'Rod' = '#6CB8D2',
  'Tube' = '#415284',
  'Ambiguous' = '#89B780',
  'BelowThreshold' = '#F5D8D0',
  'Common' = '#E69F84',
  'Unassigned' = '#8ECFF8'
)

# 过滤掉BelowThreshold点
protein_assigned_categories_filtered <- final_result %>%
  filter(Assigned_Category != "BelowThreshold")

# 获取实际类别
actual_categories_in_plot <- levels(factor(protein_assigned_categories_filtered$Assigned_Category))

# 创建颜色映射
plot_colors_for_umap <- sapply(actual_categories_in_plot, function(cat_level) {
  if (cat_level %in% names(custom_colors)) {
    return(custom_colors[cat_level])
  } else if (grepl("_&_", cat_level)) {
    return("#B2DF8A")
  } else {
    return("black")
  }
}, USE.NAMES = FALSE)
names(plot_colors_for_umap) <- actual_categories_in_plot

# 绘制UMAP图
protein_umap_plot <- ggplot(protein_assigned_categories_filtered, 
                            aes(x = Prot_UMAP1, y = Prot_UMAP2, color = Assigned_Category)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = plot_colors_for_umap, name = NULL) +
  theme_bw(base_size = 14) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.grid.major = element_line(colour = "grey95"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

print(protein_umap_plot)

# 保存图形
ggsave("E:/nanoparticl proteomics/Figures/Shape/protein_umap_Shape_Final.pdf", 
       plot = protein_umap_plot, 
       width = 5, height = 4, dpi = 300)

cat("UMAP图已保存到: protein_umap_Shape_Final.pdf\n")
cat("运行完成！\n")

# 10. 保存使用的参数信息
cat("\n=== 使用的参数 ===\n")
cat("n_top_proteins_percentile:", best_params$n_top_proteins_percentile, "\n")
cat("single_label_fraction_threshold:", best_params$single_label_fraction_threshold, "\n")
cat("min_value_multilabel_threshold:", best_params$min_value_multilabel_threshold, "\n")
cat("no_class_threshold_sum_signal:", best_params$no_class_threshold_sum_signal, "\n")


# --- 步骤11: 统计每个分类的蛋白数目并绘制柱状图 ---

# 统计每个类别的蛋白数量（排除BelowThreshold）
category_counts <- protein_assigned_categories_filtered %>%
  count(Assigned_Category, name = "Count") %>%
  arrange(desc(Count))

# 打印统计结果
cat("\nProtein counts by category (excluding BelowThreshold):\n")
print(category_counts)

# 定义特定的颜色映射（只使用指定的三个颜色）
specific_colors <- c(
  'Spherical' = '#D55640',  # 深蓝色
  'Sheet' = '#479D88',      # 橙色
  'Rod' = '#6CB8D2'        # 绿色
)

# 确保Assigned_Category是因子，并按计数排序
category_counts$Assigned_Category <- factor(category_counts$Assigned_Category, 
                                            levels = category_counts$Assigned_Category)

# 为每个因子级别分配颜色
fill_colors <- ifelse(levels(category_counts$Assigned_Category) %in% names(specific_colors),
                      specific_colors[levels(category_counts$Assigned_Category)],
                      "#CCCCCC")

names(fill_colors) <- levels(category_counts$Assigned_Category)

# 绘制柱状图（带背景网格线和指定颜色）
category_barplot <- ggplot(category_counts, aes(x = Assigned_Category, y = Count, fill = Assigned_Category)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  scale_fill_manual(values = fill_colors) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4, fontface = "bold") + # 在柱子上方添加数量标签
  theme_minimal(base_size = 14) +
  labs(
    title = "Protein Count by Shape Category",
    x = "Shape Category",
    y = "Number of Proteins"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5), # 添加主要网格线
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25), # 添加次要网格线
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5) # 添加面板边框
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) # 调整y轴扩展，为标签留出空间

# 显示柱状图
print(category_barplot)

# 保存柱状图
ggsave("E:/nanoparticl proteomics/Figures/Shape/protein_counts_by_category.pdf", 
       plot = category_barplot, 
       width = 6, height = 8, dpi = 300)

# 检查颜色分配
cat("\nColor assignment check:\n")
for(i in 1:nrow(category_counts)) {
  cat(as.character(category_counts$Assigned_Category[i]), ":", 
      ifelse(as.character(category_counts$Assigned_Category[i]) %in% names(specific_colors),
             "Specific color", "Grey"), "\n")
}
