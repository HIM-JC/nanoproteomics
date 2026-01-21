library(tabplot)
data = read.csv("C:\\Users\\huang\\OneDrive\\Desktop\\HIM\\Project\\NP-PC data mining\\formal experiment\\Data analysis\\Tabplot\\top_features_counts_per_pdf-biological identity-tabplot1.csv")

# 设置PDF保存路径
pdf("C:\\Users\\huang\\OneDrive\\Desktop\\HIM\\Project\\NP-PC data mining\\formal experiment\\Data analysis\\Tabplot\\tabplot_output.pdf", width = 15, height = 8)  # 设置宽度和高度


# 绘制 tableplot 并保存到 PDF
tableplot(data,
          pals = list(
            # 为每一列自定义颜色范围
            `NPmaterials` = colorRampPalette(c("#1E88E5", "#E3F2FD"))(15),
            NPshape = colorRampPalette(c("#1E88E5", "#E3F2FD"))(8),
            Coated = colorRampPalette(c("#64B5F6", "#E3F2FD"))(3),
            Charge = colorRampPalette(c("#64B5F6", "#E3F2FD"))(4),
            Pegylation = c("#64B5F6", "#E3F2FD"),
            Toxicity = colorRampPalette(c("#FFC9C9", "#FFF5F5"))(4),
            Uptakepathway = colorRampPalette(c("#FFC9C9", "#FFF5F5"))(5),
            Immuneresponse = colorRampPalette(c("#F44336", "#FFF5F5"))(8),
            PCprotein = colorRampPalette(c("#FFC9C9", "#FFF5F5"))(5),
            Incubationcondition = colorRampPalette(c("#F44336", "#FFF5F5"))(5),
            Otherbiologicaleffect = colorRampPalette(c("#F44336", "#FFF5F5"))(10),
            PCmethods = colorRampPalette(c("#F44336", "#FFF5F5"))(10)
          ))

# 关闭 PDF 输出
dev.off()

