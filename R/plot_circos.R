
plot_global_comm_circos <- function(
    cell_types,
    base_path,
    modes = c("activation", "inhibition"),
    output_dir = file.path(base_path, "circos_total"),
    custom_palette = c("#f57e74", "#ee86d4", "#c7c1de",
                       "#91ccc0", "#fbd54a", "#7fabd1", "#8ce1fb"),
    max_per_receiver = 5,         # 每个 receiver 选 top N 配体-受体对
    score_col   = "score",
    ligand_col  = "ligand",
    receptor_col = "receptor",
    file_suffix = "_filtered.csv",
    pdf_width  = 36,
    pdf_height = 34,
    seed = 123
) {
  suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(circlize)
  })
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # -------- 颜色设置 --------
  set.seed(seed)
  if (length(custom_palette) < length(cell_types)) {
    custom_palette <- rep(custom_palette, length.out = length(cell_types))
  }
  cell_colors <- setNames(sample(custom_palette, length(cell_types)),
                          cell_types)
  
  # -------- 内部函数：提取某个 mode 的全局 link --------
  get_global_links <- function(mode) {
    final_links <- dplyr::tibble()
    
    for (receiver in cell_types) {
      for (sender in cell_types) {
        if (sender == receiver) next
        
        f <- file.path(
          base_path,
          paste0(sender, "_vs_", receiver),
          paste0(mode, file_suffix)
        )
        if (!file.exists(f)) next
        
        df <- readr::read_csv(f, show_col_types = FALSE)
        if (!nrow(df)) next
        
        df_sub <- df %>%
          dplyr::filter(.data[[score_col]] > 0) %>%
          dplyr::mutate(
            pair = paste0(.data[[ligand_col]], "-", .data[[receptor_col]]),
            to   = receiver,
            from = paste0(pair, "@", sender)
          ) %>%
          dplyr::arrange(dplyr::desc(.data[[score_col]])) %>%
          dplyr::slice_head(n = max_per_receiver) %>%
          dplyr::select(from, to, score = dplyr::all_of(score_col))
        
        final_links <- dplyr::bind_rows(final_links, df_sub)
      }
    }
    final_links
  }
  
  # -------- 内部函数：画一个 mode 的 circos --------
  plot_one_mode <- function(links_df, mode) {
    if (!nrow(links_df)) {
      message("No links for mode = ", mode, ", skipped.")
      return(invisible(NULL))
    }
    
    # 扇区排序：按照 sender 聚合 LR，再接收细胞
    from_df <- links_df %>%
      dplyr::mutate(
        sender = sub(".*@", "", .data$from)
      ) %>%
      dplyr::arrange(factor(.data$sender, levels = cell_types))
    
    from_order <- unique(from_df$from)
    to_order   <- cell_types[cell_types %in% links_df$to]
    sector_order <- c(from_order, to_order)
    
    # 扇区颜色：发送 LR 按 sender 颜色、接收按 cell_types 颜色
    sender_color_map <- setNames(
      cell_colors[sub(".*@", "", from_order)],
      from_order
    )
    receiver_color_map <- setNames(
      cell_colors[to_order],
      to_order
    )
    sector_colors <- c(sender_color_map, receiver_color_map)
    
    # 扇区间隙
    gap_vec <- rep(1, length(sector_order))
    names(gap_vec) <- sector_order
    gap_vec[to_order]  <- 3
    gap_vec[from_order] <- 6
    
    # 输出 PDF
    pdf(file.path(output_dir, paste0("total_", mode, "_circos.pdf")),
        width = pdf_width, height = pdf_height)
    layout(matrix(c(1, 2), ncol = 2), widths = c(9, 1))
    par(mar = c(1, 1, 1, 1), cex = 1.5)
    
    # 主图
    circos.par(
      canvas.xlim   = c(-1.3, 1.3),
      canvas.ylim   = c(-1.3, 1.3),
      gap.after     = gap_vec,
      start.degree  = 180,
      track.margin  = c(0.02, 0.02),
      track.height  = 0.25,
      cell.padding  = c(0, 0, 0, 0)
    )
    
    chordDiagram(
      x              = links_df,
      order          = sector_order,
      grid.col       = sector_colors,
      directional    = 1,
      direction.type = "arrows",
      link.arr.type  = "big.arrow",
      link.sort      = TRUE,
      transparency   = 0.2,
      annotationTrack = "grid",
      preAllocateTracks = 1,
      link.lwd = pmax(links_df$score * 10 / max(links_df$score, na.rm = TRUE),
                      0.5)
    )
    
    circos.trackPlotRegion(
      track.index = 1,
      ylim = c(0, 1),
      panel.fun = function(x, y) {
        sector <- CELL_META$sector.index
        label  <- ifelse(grepl("@", sector),
                         strsplit(sector, "@")[[1]][1],
                         sector)
        circos.text(
          CELL_META$xcenter,
          CELL_META$ylim[1],
          label,
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0, 0.5),
          cex = 4
        )
      },
      bg.border = NA
    )
    
    circos.clear()
    
    # 图例
    par(mar = c(5, 0, 5, 0))
    plot.new()
    legend("center",
           legend = names(cell_colors),
           fill   = unname(cell_colors),
           border = "black",
           box.lwd = 0.5,
           cex     = 1.2,
           title   = "Cell Types",
           text.font = 1)
    
    dev.off()
    message("Saved: ", file.path(output_dir, paste0("total_", mode, "_circos.pdf")))
  }
  
  # -------- 主循环：对每个 mode 出图 --------
  for (md in modes) {
    links <- get_global_links(md)
    plot_one_mode(links, md)
  }
  
  invisible(TRUE)
}

