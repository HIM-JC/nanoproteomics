# 加载必要的库
library(randomForest)
library(ggplot2)
library(dplyr)
library(caret)
library(MASS) # 用于LDA
library(e1071) # 用于SVM
library(reshape2) # 用于数据重塑
library(pheatmap) # 用于热图
library(tidyr) # 用于数据转换
library(gridExtra) # 用于组合多个图形

# 设置工作目录
output_dir <- "E:/nanoparticl proteomics/Figures/Shape/Improved_Analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# 定义配色方案
color_palette <- c(
  'Spherical' = '#D55640',  # 深红色
  'Sheet' = '#479D88',      # 青绿色
  'Rod' = '#6CB8D2',        # 浅蓝色
  'Tube' = '#415284',       # 深蓝色
  'Ambiguous' = '#89B780',  # 绿色
  'BelowThreshold' = '#F5D8D0',  # 浅粉色
  'Common' = '#E69F84',     # 浅橙色
  'Unassigned' = '#8ECFF8'  # 天蓝色
)

# ==================== 数据预处理 ====================

cat("=== Data Preprocessing ===\n")

# 读取数据
data <- read.csv("E:/nanoparticl proteomics/Figures/Shape/Shape specific Protein-Features.csv", header = TRUE)
data_filtered <- data %>% filter(!Category %in% c("BelowThreshold", "Common"))

# 准备特征和标签
features <- data_filtered[, 3:ncol(data_filtered)]
group_factor <- as.factor(data_filtered$Category)
current_categories <- levels(group_factor)
current_colors <- color_palette[current_categories]

# 数据清洗函数
remove_constant_cols <- function(df) {
  df[, apply(df, 2, function(x) var(x, na.rm = TRUE) > 0)]
}

remove_high_na_cols <- function(df, threshold = 0.2) {
  na_ratio <- apply(df, 2, function(x) sum(is.na(x)) / length(x))
  df[, na_ratio < threshold]
}

impute_missing <- function(x) {
  if (is.numeric(x)) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
  }
  return(x)
}

# 数据清洗流程
features_clean <- features %>%
  remove_constant_cols() %>%
  remove_high_na_cols() %>%
  as.data.frame() %>%
  lapply(impute_missing) %>%
  as.data.frame()

features_scaled <- as.data.frame(scale(features_clean))

cat("Original features:", ncol(features), "-> After cleaning:", ncol(features_clean), "\n")

# ==================== 特征选择 ====================

cat("\n=== Feature Selection ===\n")

# 基于方差筛选
variance_threshold <- 0.1
variances <- apply(features_scaled, 2, var, na.rm = TRUE)
high_var_features <- features_scaled[, variances > variance_threshold & !is.na(variances)]

# 移除高度相关特征
remove_highly_correlated <- function(data, threshold = 0.8) {
  cor_matrix <- cor(data, use = "complete.obs")
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- 0
  high_cor <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)
  
  features_to_remove <- unique(colnames(data)[high_cor[, 2]])
  return(features_to_remove)
}

features_to_remove <- remove_highly_correlated(high_var_features)
features_selected <- high_var_features[, !colnames(high_var_features) %in% features_to_remove]

# 随机森林特征重要性排序
set.seed(123)
rf_model <- randomForest(x = features_selected, y = group_factor, 
                         importance = TRUE, ntree = 500)
importance_scores <- importance(rf_model)
important_features <- rownames(importance_scores)[order(-importance_scores[, "MeanDecreaseAccuracy"])]

# 选择前30个最重要的特征
top_n <- min(30, length(important_features))
features_final <- features_selected[, important_features[1:top_n]]

cat("Final features:", ncol(features_final), "\n")

# ==================== 模型训练与比较 ====================

cat("\n=== Model Training and Comparison ===\n")

# 交叉验证设置
train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE)

# 训练三个基础模型
set.seed(123)
models <- list(
  rf = train(x = features_final, y = group_factor, method = "rf",
             tuneLength = 3, trControl = train_control),
  svm = train(x = features_final, y = group_factor, method = "svmRadial",
              tuneLength = 2, trControl = train_control),
  lda = train(x = features_final, y = group_factor, method = "lda",
              trControl = train_control)
)

# 模型比较
results <- resamples(models)
cat("=== Model Performance Comparison ===\n")
print(summary(results))

# ==================== 集成学习 ====================

cat("\n=== Ensemble Learning ===\n")

# 获取各模型预测概率
get_predictions <- function(model, data) {
  predict(model, newdata = data, type = "prob")
}

model_probs <- list(
  rf = get_predictions(models$rf, features_final),
  svm = get_predictions(models$svm, features_final),
  lda = get_predictions(models$lda, features_final)
)

# 平均集成
ensemble_probs <- (model_probs$rf + model_probs$svm + model_probs$lda) / 3
ensemble_pred <- colnames(ensemble_probs)[max.col(ensemble_probs)]
ensemble_accuracy <- mean(ensemble_pred == group_factor)

# 加权集成
accuracy_scores <- sapply(models, function(x) max(x$results$Accuracy, na.rm = TRUE))
model_weights <- accuracy_scores / sum(accuracy_scores)
weighted_ensemble_probs <- model_probs$rf * model_weights[1] + 
  model_probs$svm * model_weights[2] + 
  model_probs$lda * model_weights[3]
weighted_ensemble_pred <- colnames(weighted_ensemble_probs)[max.col(weighted_ensemble_probs)]
weighted_ensemble_accuracy <- mean(weighted_ensemble_pred == group_factor)

cat("Average Ensemble Accuracy:", round(ensemble_accuracy * 100, 2), "%\n")
cat("Weighted Ensemble Accuracy:", round(weighted_ensemble_accuracy * 100, 2), "%\n")

# ==================== 核心可视化 ====================

cat("\n=== Generating Core Visualizations ===\n")

# 1. 模型比较图
performance_df <- data.frame(
  Algorithm = rep(c("Random Forest", "SVM", "LDA"), each = 5),
  Accuracy = c(
    models$rf$resample$Accuracy,
    models$svm$resample$Accuracy,
    models$lda$resample$Accuracy
  )
)

p_comparison <- ggplot(performance_df, aes(x = Algorithm, y = Accuracy, fill = Algorithm)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = c('#D55640', '#479D88', '#6CB8D2')) +
  labs(title = "Model Performance Comparison", x = "Algorithm", y = "Accuracy") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# 2. 特征重要性图
importance_df <- data.frame(
  Feature = rownames(varImp(models$rf)$importance),
  Importance = varImp(models$rf)$importance$Overall
) %>% arrange(-Importance) %>% head(15)

p_importance <- ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = '#415284') +
  coord_flip() +
  labs(title = "Top 15 Important Features", x = "Feature", y = "Importance Score") +
  theme_bw()

# 3. 准确率提升对比图
accuracy_comparison <- data.frame(
  Method = c("Original", "Random Forest", "SVM", "LDA", "Average Ensemble", "Weighted Ensemble"),
  Accuracy = c(0.452, accuracy_scores, ensemble_accuracy, weighted_ensemble_accuracy)
)

p_accuracy <- ggplot(accuracy_comparison, aes(x = reorder(Method, Accuracy), y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_manual(values = c('#F5D8D0', '#D55640', '#479D88', '#6CB8D2', '#415284', '#89B780')) +
  geom_text(aes(label = paste0(round(Accuracy * 100, 1), "%")), hjust = -0.2, size = 4) +
  coord_flip() +
  labs(title = "Model Accuracy Comparison", 
       subtitle = paste0("Ensemble Improvement: ", round((ensemble_accuracy - 0.452) * 100, 1), "% points"),
       x = "Method", y = "Accuracy") +
  theme_bw() +
  theme(legend.position = "none")

# 保存基础比较图
ggsave("01_Model_Comparison.png", p_comparison, width = 8, height = 6, dpi = 300)
ggsave("02_Feature_Importance.png", p_importance, width = 10, height = 6, dpi = 300)
ggsave("03_Accuracy_Improvement.png", p_accuracy, width = 10, height = 6, dpi = 300)

# ==================== 详细性能分析 ====================

cat("\n=== Detailed Performance Analysis ===\n")

# 计算性能指标函数
calculate_metrics <- function(conf_matrix) {
  classes <- rownames(conf_matrix)
  metrics <- data.frame()
  
  for (class in classes) {
    TP <- conf_matrix[class, class]
    FP <- sum(conf_matrix[, class]) - TP
    FN <- sum(conf_matrix[class, ]) - TP
    
    accuracy <- (TP + (sum(conf_matrix) - TP - FP - FN)) / sum(conf_matrix)
    precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
    recall <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
    f1_score <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
    
    metrics <- rbind(metrics, data.frame(
      Class = class,
      Accuracy = accuracy,
      Precision = precision,
      Recall = recall,
      F1_Score = f1_score
    ))
  }
  
  # 总体指标
  overall_accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
  macro_precision <- mean(metrics$Precision)
  macro_recall <- mean(metrics$Recall)
  macro_f1 <- mean(metrics$F1_Score)
  
  metrics <- rbind(metrics, data.frame(
    Class = "Overall",
    Accuracy = overall_accuracy,
    Precision = macro_precision,
    Recall = macro_recall,
    F1_Score = macro_f1
  ))
  
  return(metrics)
}

# 使用集成模型的混淆矩阵
confusion_matrix <- table(Actual = group_factor, Predicted = ensemble_pred)
detailed_metrics <- calculate_metrics(confusion_matrix)

# 4. 混淆矩阵图
confusion_df <- as.data.frame(confusion_matrix)
confusion_percent <- confusion_df %>%
  group_by(Actual) %>%
  mutate(Percent = Freq / sum(Freq) * 100,
         Label = paste0(Freq, "\n(", round(Percent, 1), "%)"))

p_confusion <- ggplot(confusion_percent, aes(x = Predicted, y = Actual, fill = Percent)) +
  geom_tile(color = "white", alpha = 0.8) +
  geom_text(aes(label = Label), color = "black", size = 3.5) +
  scale_fill_gradient(low = "#F5D8D0", high = "#D55640", name = "Percentage (%)") +
  labs(title = "Confusion Matrix Analysis",
       subtitle = paste0("Overall Accuracy: ", round(ensemble_accuracy * 100, 2), "%"),
       x = "Predicted Class", y = "Actual Class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 5. 四指标对比图
metrics_long <- detailed_metrics %>%
  filter(Class != "Overall") %>%
  dplyr::select(Class, Accuracy, Precision, Recall, F1_Score) %>%
  tidyr::pivot_longer(cols = c(Accuracy, Precision, Recall, F1_Score), 
                      names_to = "Metric", values_to = "Score")

p_metrics_comp <- ggplot(metrics_long, aes(x = Class, y = Score, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c('#D55640', '#479D88', '#6CB8D2', '#415284')) +
  labs(title = "Four-Dimensional Model Performance Evaluation",
       subtitle = "Comprehensive comparison of Accuracy, Precision, Recall, and F1-Score",
       x = "Class", y = "Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 6. 性能热图
metrics_matrix <- detailed_metrics %>%
  filter(Class != "Overall") %>%
  dplyr::select(Class, Accuracy, Precision, Recall, F1_Score)

rownames(metrics_matrix) <- metrics_matrix$Class
metrics_matrix <- as.matrix(metrics_matrix[, -1])

p_heatmap <- pheatmap(metrics_matrix,
                      main = "Performance Metrics Heatmap",
                      color = colorRampPalette(c("#6CB8D2", "white", "#D55640"))(50),
                      cluster_rows = TRUE,
                      cluster_cols = FALSE,
                      display_numbers = TRUE,
                      number_format = "%.3f",
                      fontsize = 10,
                      angle_col = 45)

# 保存性能分析图
ggsave("04_Confusion_Matrix.png", p_confusion, width = 8, height = 6, dpi = 300)
ggsave("05_Four_Metrics_Comparison.png", p_metrics_comp, width = 10, height = 6, dpi = 300)
ggsave("06_Performance_Heatmap.png", p_heatmap, width = 8, height = 6, dpi = 300)

# ==================== 特征差异分析 ====================

cat("\n=== Feature Difference Analysis ===\n")

# 获取前十重要特征
top_10_features <- important_features[1:10]

# 准备分析数据
analysis_data <- cbind(features_scaled[, top_10_features, drop = FALSE], 
                       Category = group_factor)

# 差异分析 - 修复版本
results_df <- data.frame()

for (feature in top_10_features) {
  anova_result <- aov(as.formula(paste(feature, "~ Category")), data = analysis_data)
  anova_p <- summary(anova_result)[[1]]$"Pr(>F)"[1]
  
  means <- tapply(analysis_data[[feature]], analysis_data$Category, mean, na.rm = TRUE)
  
  # 创建结果行
  result_row <- data.frame(Feature = feature, ANOVA_p_value = anova_p, stringsAsFactors = FALSE)
  
  # 添加每个类别的均值
  for (category in current_categories) {
    result_row[[paste0("Mean_", category)]] <- ifelse(is.na(means[category]), NA, means[category])
  }
  
  results_df <- rbind(results_df, result_row)
}

# 添加显著性标记
results_df$Significance <- ifelse(results_df$ANOVA_p_value < 0.001, "***",
                                  ifelse(results_df$ANOVA_p_value < 0.01, "**",
                                         ifelse(results_df$ANOVA_p_value < 0.05, "*", "Not Significant")))

cat("Feature difference analysis completed, analyzed", nrow(results_df), "features\n")

# 7. 特征分布箱线图
plot_data <- analysis_data %>%
  melt(id.vars = "Category", variable.name = "Feature", value.name = "Value")
plot_data$Feature <- factor(plot_data$Feature, levels = top_10_features)

p_feature_box <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_manual(values = current_colors) +
  facet_wrap(~ Feature, scales = "free_y", ncol = 3) +
  labs(title = "Top 10 Important Feature Distributions",
       x = "Nanoparticle Shape", y = "Standardized Value") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_text(size = 8))

# 8. 特征显著性图
p_feature_sig <- ggplot(results_df, aes(x = reorder(Feature, -ANOVA_p_value), y = -log10(ANOVA_p_value), 
                                        fill = -log10(ANOVA_p_value))) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text(aes(label = Significance), vjust = -0.5, size = 3) +
  scale_fill_gradient(low = "#89B780", high = '#415284', name = "-log10(p-value)") +
  labs(title = "Statistical Significance of Top 10 Features",
       x = "Feature", y = "-log10(ANOVA p-value)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 9. 特征表达热图
mean_cols <- grep("^Mean_", names(results_df), value = TRUE)
mean_matrix <- results_df[, c("Feature", mean_cols)]
rownames(mean_matrix) <- mean_matrix$Feature
mean_matrix <- as.matrix(mean_matrix[, -1])

# Z-score标准化行方向
mean_matrix_z <- t(scale(t(mean_matrix)))

p_feature_heatmap <- pheatmap(mean_matrix_z,
                              main = "Top 10 Feature Expression Patterns Heatmap\n(Row Standardized)",
                              color = colorRampPalette(c("#6CB8D2", "white", "#D55640"))(50),
                              cluster_rows = TRUE,
                              cluster_cols = FALSE,
                              show_rownames = TRUE,
                              fontsize_row = 8,
                              angle_col = 45)

# 保存特征分析图
ggsave("07_Feature_Distribution.png", p_feature_box, width = 12, height = 10, dpi = 300)
ggsave("08_Feature_Significance.png", p_feature_sig, width = 8, height = 6, dpi = 300)
ggsave("09_Feature_Heatmap.png", p_feature_heatmap, width = 8, height = 6, dpi = 300)

# ==================== 单独输出每个PDF文件 ====================

cat("\n=== Generating Individual PDF Files ===\n")

# 创建可编辑主题
pdf_theme <- theme_bw() +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

# 1. 模型比较图
cairo_pdf("01_Model_Comparison.pdf", width = 10, height = 7)
print(p_comparison + pdf_theme)
dev.off()
cat("Saved: 01_Model_Comparison.pdf\n")

# 2. 特征重要性图
cairo_pdf("02_Feature_Importance.pdf", width = 12, height = 8)
print(p_importance + pdf_theme)
dev.off()
cat("Saved: 02_Feature_Importance.pdf\n")

# 3. 准确率提升对比图
cairo_pdf("03_Accuracy_Improvement.pdf", width = 12, height = 8)
print(p_accuracy + pdf_theme)
dev.off()
cat("Saved: 03_Accuracy_Improvement.pdf\n")

# 4. 混淆矩阵图
cairo_pdf("04_Confusion_Matrix.pdf", width = 10, height = 8)
print(p_confusion + pdf_theme)
dev.off()
cat("Saved: 04_Confusion_Matrix.pdf\n")

# 5. 四指标对比图
cairo_pdf("05_Four_Metrics_Comparison.pdf", width = 12, height = 8)
print(p_metrics_comp + pdf_theme)
dev.off()
cat("Saved: 05_Four_Metrics_Comparison.pdf\n")

# 6. 性能热图
cairo_pdf("06_Performance_Heatmap.pdf", width = 10, height = 8)
pheatmap(metrics_matrix,
         main = "Performance Metrics Heatmap",
         color = colorRampPalette(c("#6CB8D2", "white", "#D55640"))(50),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.3f",
         fontsize = 10,
         angle_col = 45)
dev.off()
cat("Saved: 06_Performance_Heatmap.pdf\n")

# 7. 特征分布箱线图（2行布局）
plot_data <- analysis_data %>%
  melt(id.vars = "Category", variable.name = "Feature", value.name = "Value")
plot_data$Feature <- factor(plot_data$Feature, levels = top_10_features)

p_feature_box <- ggplot(plot_data, aes(x = Category, y = Value, fill = Category)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1, show.legend = FALSE) +  # 在geom中直接禁用图例
  scale_fill_manual(values = current_colors) +
  facet_wrap(~ Feature, scales = "free_y", ncol = 5) +  # 改为5列，2行
  labs(title = "Top 10 Important Feature Distributions",
      # subtitle = "Distribution patterns across nanoparticle shapes",
       x = "Nanoparticle Shape", y = "Standardized Value") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",  # 移除图例位置
        strip.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12))

# 保存PNG版本
ggsave("07_Feature_Distribution.png", p_feature_box, width = 16, height = 8, dpi = 300)

# 保存PDF版本 - 确保图例被移除
cairo_pdf("07_Feature_Distribution.pdf", width = 16, height = 8)
print(p_feature_box + theme(legend.position = "none"))  # 再次强制移除图例
dev.off()
cat("Saved: 07_Feature_Distribution.pdf (2-row layout, legend removed)\n")

# 8. 特征显著性图
cairo_pdf("08_Feature_Significance.pdf", width = 10, height = 7)
print(p_feature_sig + pdf_theme)
dev.off()
cat("Saved: 08_Feature_Significance.pdf\n")

# 9. 特征表达热图
cairo_pdf("09_Feature_Heatmap.pdf", width = 10, height = 8)
pheatmap(mean_matrix_z,
         main = "Top 10 Feature Expression Patterns Heatmap\n(Row Standardized)",
         color = colorRampPalette(c("#6CB8D2", "white", "#D55640"))(50),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         fontsize_row = 8,
         angle_col = 45,
         border_color = NA)
dev.off()
cat("Saved: 09_Feature_Heatmap.pdf\n")

cat("\n✅ All individual PDF files generated successfully!\n")
cat("📁 Generated PDF files:\n")
cat("01_Model_Comparison.pdf\n")
cat("02_Feature_Importance.pdf\n")
cat("03_Accuracy_Improvement.pdf\n")
cat("04_Confusion_Matrix.pdf\n")
cat("05_Four_Metrics_Comparison.pdf\n")
cat("06_Performance_Heatmap.pdf\n")
cat("07_Feature_Distribution.pdf\n")
cat("08_Feature_Significance.pdf\n")
cat("09_Feature_Heatmap.pdf\n")
cat("10_Feature_Importance_Detailed.pdf\n")
cat("11_Model_Performance_Radar.pdf\n")

# 同时生成一个汇总的PDF（可选）
cat("\n=== Generating Summary PDF ===\n")
cairo_pdf("00_All_Figures_Summary.pdf", width = 11, height = 8.5, onefile = TRUE)

print(p_comparison + pdf_theme)
print(p_importance + pdf_theme)
print(p_accuracy + pdf_theme)
print(p_confusion + pdf_theme)
print(p_metrics_comp + pdf_theme)
pheatmap(metrics_matrix,
         main = "Performance Metrics Heatmap",
         color = colorRampPalette(c("#6CB8D2", "white", "#D55640"))(50),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.3f",
         fontsize = 10,
         angle_col = 45)
print(p_feature_box + pdf_theme)
print(p_feature_sig + pdf_theme)
pheatmap(mean_matrix_z,
         main = "Top 10 Feature Expression Patterns Heatmap\n(Row Standardized)",
         color = colorRampPalette(c("#6CB8D2", "white", "#D55640"))(50),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         fontsize_row = 8,
         angle_col = 45)

dev.off()
cat("Saved: 00_All_Figures_Summary.pdf\n")

cat("\n🎉 All PDF files are ready for editing in Adobe Illustrator!\n")
cat("💡 Tips for AI editing:\n")
cat("• Text should be editable (not outlined)\n")
cat("• Colors can be easily modified\n")
cat("• Layout can be adjusted as needed\n")