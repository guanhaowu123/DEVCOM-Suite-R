plot_lr_network <- function(
    in_dir,
    cell_types,
    ligand_of_interest,
    receptor_of_interest,
    winsor_q        = 0.95,
    drop_zero_edges = TRUE,
    lwd_min_raw     = 0.6,
    lwd_max_raw     = 4.0,
    lwd_min_plot    = 0.35,
    lwd_max_plot    = 1.6,
    arrow_size_mm   = 1.2,
    loop_strength   = 0.08,
    node_palette = c(
      "#f57e74", "#c2e4d0", "#c3aaf7", "#ffb5e8", "#c7ceea",
      "#faa9a0", "#bae1ff", "#aaead7", "#badac1", "#d8b4f8",
      "#ffae99", "#ffcd45", "#8ce1fb", "#faaab5", "#eee0cb"
    ),
    mutual_edge_color = "#c59d94",
    pdf_width  = 8,
    pdf_height = 7
) {
  # --- 参数检查：强制用户显式指定 L / R ---
  if (missing(ligand_of_interest) || missing(receptor_of_interest)) {
    stop("请显式传入 ligand_of_interest 和 receptor_of_interest，例如:\n",
         'plot_lr_network(..., ligand_of_interest = "MIF", receptor_of_interest = "ITGA4")')
  }
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(stringr)
    library(igraph)
    library(ggraph)
    library(scales)
  })
  if (!requireNamespace("tidygraph", quietly = TRUE)) {
    install.packages("tidygraph", repos = "https://cloud.r-project.org")
  }
  library(tidygraph)
  
  # 颜色映射
  if (length(node_palette) < length(cell_types)) {
    node_palette <- rep(node_palette, length.out = length(cell_types))
  }
  color_map <- setNames(node_palette[seq_along(cell_types)], cell_types)
  arrow_size <- grid::unit(arrow_size_mm, "mm")
  
  # 读取 sender→receiver 上该 LR 的总分
  read_edge_score_for_pair <- function(sender, receiver, dir_path) {
    fp <- if (sender == receiver) {
      file.path(dir_path, sprintf("%s_ligand_receptor_with_score.csv", sender))
    } else {
      file.path(dir_path, sprintf("%s_%s_ligand_receptor_unique_with_score.csv", sender, receiver))
    }
    if (!file.exists(fp)) return(0)
    
    df <- suppressMessages(readr::read_csv(fp, guess_max = 1e6))
    
    nms <- tolower(names(df))
    score_col <- names(df)[nms == "activity_score"]
    if (length(score_col) == 0) {
      cand <- names(df)[str_detect(nms, "activity") & str_detect(nms, "score")]
      if (length(cand) == 0) cand <- names(df)[str_detect(nms, "score")]
      if (length(cand) == 0) return(0)
      score_col <- cand[1]
    }
    
    L_cols <- grep("^L", names(df), value = TRUE, ignore.case = TRUE)
    R_cols <- grep("^R", names(df), value = TRUE, ignore.case = TRUE)
    if (length(L_cols) == 0 || length(R_cols) == 0) return(0)
    
    to_up <- function(x) toupper(trimws(as.character(x)))
    L_hit <- Reduce(`|`, lapply(L_cols, function(c) to_up(df[[c]]) == toupper(ligand_of_interest)))
    R_hit <- Reduce(`|`, lapply(R_cols, function(c) to_up(df[[c]]) == toupper(receptor_of_interest)))
    L_hit[is.na(L_hit)] <- FALSE
    R_hit[is.na(R_hit)] <- FALSE
    
    hit <- L_hit & R_hit
    if (!any(hit)) return(0)
    
    sc <- df[[score_col]][hit]
    sc <- sc[is.finite(sc)]
    if (!length(sc)) return(0)
    
    sum(sc, na.rm = TRUE)
  }
  
  # 构建边
  edges <- expand.grid(sender = cell_types,
                       receiver = cell_types,
                       stringsAsFactors = FALSE) |>
    as_tibble() |>
    mutate(weight_raw = purrr::pmap_dbl(
      list(sender, receiver),
      ~ read_edge_score_for_pair(..1, ..2, in_dir)
    ))
  
  if (drop_zero_edges) {
    edges <- edges |> dplyr::filter(weight_raw > 0)
  }
  if (!nrow(edges)) {
    warning("No non-zero edges found for ",
            ligand_of_interest, "-", receptor_of_interest,
            " in directory: ", in_dir)
    return(invisible(NULL))
  }
  
  # 动态范围
  edges <- edges |> mutate(w_log = log1p(weight_raw))
  cap <- quantile(edges$w_log, winsor_q, na.rm = TRUE)
  edges <- edges |> mutate(w_bal = pmin(w_log, cap))
  edges$lwd <- if (diff(range(edges$w_bal)) < 1e-12) {
    rep(mean(c(lwd_min_raw, lwd_max_raw)), nrow(edges))
  } else {
    rescale(edges$w_bal, to = c(lwd_min_raw, lwd_max_raw))
  }
  
  # 识别互向边
  edges_g <- edges |>
    transmute(
      from = sender, to = receiver,
      sender, receiver, weight_raw, w_bal, lwd,
      pair_key = if_else(sender < receiver,
                         paste(sender, receiver, sep = "|"),
                         paste(receiver, sender, sep = "|")),
      is_loop = sender == receiver
    ) |>
    group_by(pair_key) |>
    mutate(n_in_pair = dplyr::n()) |>
    ungroup()
  
  dom_w <- edges_g |>
    dplyr::filter(!is_loop) |>
    group_by(pair_key) |>
    summarise(lwd_pair = max(lwd), .groups = "drop")
  
  edges_g <- edges_g |>
    left_join(dom_w, by = "pair_key") |>
    mutate(
      is_mutual      = (!is_loop & n_in_pair == 2),
      plot_lwd_raw   = if_else(is_mutual, lwd_pair, lwd),
      plot_color_key = if_else(is_mutual, "Mutual", as.character(sender))
    )
  
  edges_g$plot_lwd <- rescale(edges_g$plot_lwd_raw,
                              to = c(lwd_min_plot, lwd_max_plot))
  
  # 建图
  g <- tidygraph::tbl_graph(
    nodes = tibble(name = cell_types),
    edges = edges_g |> select(from, to, plot_lwd, plot_color_key),
    directed = TRUE
  )
  
  # 画图
  p <- ggraph(g, layout = "circle") +
    geom_edge_link(
      aes(width = plot_lwd, color = plot_color_key),
      lineend  = "round",
      arrow    = arrow(length = arrow_size, type = "closed"),
      end_cap  = circle(1.6, "mm"),
      start_cap = circle(1.6, "mm"),
      filter   = ~ !edge_is_loop()
    ) +
    geom_edge_loop(
      aes(width = plot_lwd, color = plot_color_key),
      arrow     = arrow(length = arrow_size, type = "closed"),
      strength  = loop_strength,
      start_cap = circle(1.2, "mm"),
      end_cap   = circle(1.2, "mm"),
      filter    = ~ edge_is_loop()
    ) +
    geom_node_point(aes(fill = name),
                    size = 10, shape = 21,
                    color = "grey20", stroke = 0.6) +
    geom_node_text(aes(label = name),
                   size = 3.5, fontface = "bold") +
    scale_fill_manual(values = color_map, guide = "none") +
    scale_edge_colour_manual(
      values = c(color_map, Mutual = mutual_edge_color),
      name   = "Edge color"
    ) +
    scale_edge_width(range = c(lwd_min_plot, lwd_max_plot),
                     guide = "none") +
    coord_fixed() +
    theme_void(base_size = 12) +
    labs(
      title    = sprintf("%s − %s communication network",
                         ligand_of_interest, receptor_of_interest),
      subtitle = "Mutual pairs use a distinct color; width = max of two directions\nSelf loops are smaller"
    ) +
    theme(plot.title = element_text(face = "bold"))
  
  base_name <- sprintf("network_%s_%s",
                       ligand_of_interest, receptor_of_interest)
  out_pdf <- file.path(in_dir, paste0(base_name, ".pdf"))
  out_csv <- file.path(in_dir, paste0(base_name, "_edges.csv"))
  
  ggsave(filename = out_pdf, plot = p,
         width = pdf_width, height = pdf_height,
         limitsize = FALSE)
  message("✅ 已输出网络图：", normalizePath(out_pdf, winslash = "/"))
  
  write.csv(
    edges_g |> select(from, to, sender, receiver,
                      weight_raw, w_bal,
                      plot_color_key, plot_lwd),
    out_csv,
    row.names = FALSE
  )
  message("✅ 已导出边表：", normalizePath(out_csv, winslash = "/"))
  
  invisible(list(plot = p, edges = edges_g))
}
