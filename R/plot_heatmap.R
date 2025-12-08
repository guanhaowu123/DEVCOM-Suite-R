#heatmap#

plot_comm_heatmap <- function(
    input_file,
    output_file,
    title = "IVF",
    legend_title = "MAX LRI",
    text_threshold = 100,   # 大于该值时数字用白色
    width = 8,
    height = 6
) {
  # 加载需要的包
  suppressPackageStartupMessages({
    library(ggplot2)
    library(reshape2)
    library(grid)  # 为 unit() 函数
  })
  
  # 读取数据（第一列作为行名）
  data <- read.csv(input_file, row.names = 1, check.names = FALSE)
  
  # 四舍五入到整数
  data <- round(data)
  
  # 按行、列均值进行排序
  row_order <- order(rowMeans(data, na.rm = TRUE))
  col_order <- order(colMeans(data, na.rm = TRUE))
  data <- data[row_order, col_order, drop = FALSE]
  
  # 行名作为新列
  data$Row <- rownames(data)
  
  # 转为长表
  data_melted <- melt(
    data,
    id.vars = "Row",
    variable.name = "Column",
    value.name = "Value"
  )
  
  # 颜色渐变
  gradient_colors <- scale_fill_gradientn(
    colors = c("white", "lightyellow", "#F5834A", "#E25E55", "#C43F5E"),
    name   = legend_title
  )
  
  # 数值颜色：根据阈值动态调整
  dynamic_number_color <- function(value) {
    ifelse(value > text_threshold, "white", "black")
  }
  data_melted$number_color <- dynamic_number_color(data_melted$Value)
  
  # 因子顺序保持排序后顺序
  data_melted$Row    <- factor(data_melted$Row,    levels = rownames(data))
  data_melted$Column <- factor(data_melted$Column, levels = colnames(data))
  
  # 画图
  p <- ggplot(data_melted, aes(x = Column, y = Row, fill = Value)) +
    geom_tile(color = "white") +
    gradient_colors +
    geom_text(
      aes(label = Value, color = number_color),
      size = 4,
      fontface = "bold"
    ) +
    scale_color_identity() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 90, hjust = 1, size = 14, face = "bold"
      ),
      axis.text.y = element_text(size = 14, face = "bold"),
      plot.title  = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "right"
    ) +
    labs(title = title, x = NULL, y = NULL) +
    guides(fill = guide_colorbar(
      barheight = unit(0.7, "npc"),
      barwidth  = unit(0.03, "npc")
    ))
  
  # 保存 PDF
  ggsave(output_file, plot = p, width = width, height = height, device = cairo_pdf)
  
  invisible(p)
}
