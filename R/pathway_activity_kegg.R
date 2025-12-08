# R language

## KEGG ##


# kegg_pathway_activity.R
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
})

# ------------------------------------------------------------
# Helper: load gene importance scores -> two columns: Gene, Score
# score_col can be a column name or index; by default it first looks
# for a column named "Score", otherwise it uses the 3rd column.
# ------------------------------------------------------------
load_gene_scores <- function(file, gene_col = "Gene", score_col = NULL) {
  stopifnot(file.exists(file))
  df <- readr::read_csv(file, show_col_types = FALSE)
  
  # Gene column
  if (is.numeric(gene_col)) {
    gname <- names(df)[[gene_col]]
  } else if (is.character(gene_col) && gene_col %in% names(df)) {
    gname <- gene_col
  } else {
    cand <- names(df)[tolower(names(df)) %in% c("gene","symbol","genes")]
    if (!length(cand)) stop("无法在 ", file, " 中识别基因列（gene_col）。")
    gname <- cand[[1]]
  }
  
  # Score column
  if (is.null(score_col)) {
    sname <- if ("Score" %in% names(df)) "Score" else {
      if (ncol(df) < 3) stop("无法定位得分列（score_col），且文件列数 < 3：", file)
      names(df)[[3]]
    }
  } else if (is.numeric(score_col)) {
    sname <- names(df)[[score_col]]
  } else {
    sname <- score_col
    if (!sname %in% names(df)) stop("指定的 score_col 在文件中不存在：", score_col)
  }
  
  out <- df %>%
    dplyr::select(Gene = all_of(gname), Score = all_of(sname)) %>%
    dplyr::mutate(Gene = as.character(Gene)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(Score = mean(Score, na.rm = TRUE), .groups = "drop")
  
  return(out)
}

# ------------------------------------------------------------
# Helper: get score for a single gene (return 0 if missing)
# ------------------------------------------------------------

.get_score_simple <- function(gene, score_map) {
  if (is.na(gene) || gene == "") return(0)
  
  # Force to character (avoid issues from factors or other types)
  gene_chr <- as.character(gene)
  
  # Use [ instead of [[; when the name does not exist it returns NA instead of error
  vals <- suppressWarnings(as.numeric(score_map[gene_chr]))
  
  # If nothing valid is retrieved, return 0
  if (length(vals) == 0 || all(is.na(vals))) {
    return(0)
  }
  
  v <- vals[1]
  ifelse(is.na(v), 0, v)
}


# ------------------------------------------------------------
# Helper: score for multi-gene strings (split by sep), combine = sum/mean/max
# Only used when split_complex = TRUE.
# ------------------------------------------------------------
.sum_scores_for_string <- function(gene_str, score_map, sep = ",", combine = c("sum","mean","max")) {
  combine <- match.arg(combine)
  if (is.na(gene_str) || gene_str == "") return(0)
  genes <- str_split(gene_str, sep, simplify = TRUE)
  genes <- str_trim(genes[genes != ""])
  if (!length(genes)) return(0)
  vals <- as.numeric(score_map[genes]); vals[is.na(vals)] <- 0
  switch(combine,
         sum  = sum(vals),
         mean = mean(vals),
         max  = ifelse(length(vals)==0, 0, max(vals)))
}

# ------------------------------------------------------------
# Helper: pick the first existing column from a list of candidate names
# ------------------------------------------------------------
.pick_col <- function(df, candidates) {
  for (c in candidates) if (c %in% names(df)) return(c)
  return(NULL)
}

# Provide a non-dot alias so later code calling pick_col() works
pick_col <- .pick_col

# ------------------------------------------------------------
# Core function: compute pathway activity for a KEGG-annotated LR table
#
# Arguments:
#   lr_df            : data.frame containing at least
#                      {Complex_Ligands, Complex_Receptors, TF, Target, KEGG_Pathways_Names}
#                      If your column names differ (e.g. L/R/Name), specify via col_map.
#   ligand_scores_df : gene scores for the sender cell (Gene, Score)
#   target_scores_df : gene scores for the receiver cell (Gene, Score)
#   kegg_name_col    : column name for KEGG pathway names; default candidates:
#                      c("KEGG_Pathways_Names","KEGG_Pathways","Name")
#   split_complex    : whether to split complexes (e.g. "ITGAV,ITGB3") into genes
#                      and aggregate scores (default FALSE = consistent with original script)
#   multi_gene_sep   : separator within complexes (default ",")
#   combine_multi    : aggregation mode for multi-gene scores: sum/mean/max
#                      (only used when split_complex = TRUE)
#   weighting        : "count" => Final_Score = Pathway_Score * count
#                      "none"  => Final_Score = Pathway_Score
#   weights          : component weights list(ligand=1, receptor=1, tf=1, target=1)
#   kegg_sep_regex   : regex for splitting multiple pathways (default ";\\s*")
#   return_lr_scored : whether to return per-row scored LR table
#
# Returns: list(pathway_scores = ..., lr_scored = optional)
# ------------------------------------------------------------
compute_kegg_pathway_activity <- function(
    lr_df,
    ligand_scores_df,
    target_scores_df,
    kegg_name_col = NULL,
    split_complex = FALSE,
    multi_gene_sep = ",",
    combine_multi = c("sum","mean","max"),
    weighting = c("count","none"),
    weights = list(ligand=1, receptor=1, tf=1, target=1),
    kegg_sep_regex = ";\\s*",
    return_lr_scored = FALSE,
    col_map = list(L = c("Complex_Ligands","L","Ligand","Ligands"),
                   R = c("Complex_Receptors","R","Receptor","Receptors"),
                   TF = c("TF","tf"),
                   Target = c("Target","target"))
) {
  combine_multi <- match.arg(combine_multi)
  weighting     <- match.arg(weighting)
  
  # Column mapping
  L_col  <- pick_col(lr_df, col_map$L);       if (is.null(L_col))  stop("未找到配体列（L/Complex_Ligands）")
  R_col  <- pick_col(lr_df, col_map$R);       if (is.null(R_col))  stop("未找到受体列（R/Complex_Receptors）")
  TF_col <- pick_col(lr_df, col_map$TF);      if (is.null(TF_col)) stop("未找到 TF 列（TF/tf）")
  T_col  <- pick_col(lr_df, col_map$Target);  if (is.null(T_col))  stop("未找到 Target 列（Target）")
  
  if (is.null(kegg_name_col)) {
    kegg_name_col <- pick_col(lr_df, c("KEGG_Pathways_Names","KEGG_Pathways","Name","Pathway"))
    if (is.null(kegg_name_col)) stop("未找到 KEGG 通路列（如 KEGG_Pathways_Names/KEGG_Pathways/Name）")
  }
  
  # Named vectors for fast score lookups
  lig_map <- setNames(ligand_scores_df$Score, ligand_scores_df$Gene)
  tar_map <- setNames(target_scores_df$Score, target_scores_df$Gene)
  
  # Scoring strategy for ligand/receptor
  score_L <- if (split_complex) {
    function(s) .sum_scores_for_string(s, lig_map, sep = multi_gene_sep, combine = combine_multi)
  } else {
    function(s) .get_score_simple(s, lig_map)
  }
  score_R <- if (split_complex) {
    function(s) .sum_scores_for_string(s, tar_map, sep = multi_gene_sep, combine = combine_multi)
  } else {
    function(s) .get_score_simple(s, tar_map)
  }
  
  # Per-row scoring (consistent with original formula + optional weights)
  lr_scored <- lr_df %>%
    dplyr::mutate(
      !!L_col  := as.character(.data[[L_col]]),
      !!R_col  := as.character(.data[[R_col]]),
      !!TF_col := as.character(.data[[TF_col]]),
      !!T_col  := as.character(.data[[T_col]])
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      ligand_score   = score_L(.data[[L_col]]),
      receptor_score = score_R(.data[[R_col]]),
      tf_score       = .get_score_simple(.data[[TF_col]], tar_map),
      target_score   = .get_score_simple(.data[[T_col]], tar_map),
      total_score    = weights$ligand  * ligand_score +
        weights$receptor* receptor_score +
        weights$tf      * tf_score +
        weights$target  * target_score,
      CCI_Pair = paste(
        .data[[L_col]],
        paste0(na.omit(c(.data[[R_col]])), collapse = "+"),
        .data[[TF_col]],
        .data[[T_col]],
        sep = "->"
      )
    ) %>%
    dplyr::ungroup()
  
  # Split multi-pathway entries and aggregate pathway-level scores
  pathway_scores <- lr_scored %>%
    tidyr::separate_rows(dplyr::all_of(kegg_name_col), sep = kegg_sep_regex) %>%
    dplyr::rename(KEGG_Name = dplyr::all_of(kegg_name_col)) %>%
    dplyr::filter(KEGG_Name != "") %>%
    dplyr::group_by(KEGG_Name) %>%
    dplyr::summarise(
      Pathway_Score = sum(.data$total_score),
      Weight        = dplyr::n(),
      CCI_Pairs     = paste(.data$CCI_Pair, collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Final_Score = dplyr::case_when(
        weighting == "count" ~ .data$Pathway_Score * .data$Weight,
        TRUE                 ~ .data$Pathway_Score
      )
    ) %>%
    dplyr::arrange(dplyr::desc(.data$Final_Score))
  
  out <- list(pathway_scores = pathway_scores)
  if (isTRUE(return_lr_scored)) out$lr_scored <- lr_scored
  return(out)
}

# ------------------------------------------------------------
# Single-file wrapper:
#   Automatically infer A_to_B or A autocrine from file name,
#   and write outputs with names consistent with your original script.
#   Input files are usually the integrated "_merged_LR_RTF_TFT.csv"
#   that already contain KEGG_Pathways_Names.
# ------------------------------------------------------------
kegg_pathway_activity_for_file <- function(
    comm_file,
    ligand_score_file,
    target_score_file,
    output_dir,
    kegg_name_col = NULL,
    split_complex = FALSE,
    multi_gene_sep = ",",
    combine_multi = c("sum","mean","max"),
    weighting = c("count","none"),
    weights = list(ligand=1, receptor=1, tf=1, target=1),
    kegg_sep_regex = ";\\s*",
    score_gene_col = "Gene",
    score_val_col  = NULL,   # NULL = auto; equivalent to “3rd column” logic you used
    save = TRUE,
    return_lr_scored = FALSE
) {
  stopifnot(file.exists(comm_file))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  fname <- basename(comm_file)
  if (str_detect(fname, "_to_")) {
    cell1 <- str_extract(fname, "^[^_]+")
    cell2 <- str_extract(fname, "(?<=_to_)[^_]+")
  } else {
    cell1 <- str_remove(fname, "_merged_LR_RTF_TFT\\.csv$")
    cell2 <- cell1
  }
  
  lr_df         <- readr::read_csv(comm_file, show_col_types = FALSE)
  ligand_scores <- load_gene_scores(ligand_score_file, gene_col = score_gene_col, score_col = score_val_col)
  target_scores <- load_gene_scores(target_score_file, gene_col = score_gene_col, score_col = score_val_col)
  
  res <- compute_kegg_pathway_activity(
    lr_df, ligand_scores, target_scores,
    kegg_name_col = kegg_name_col,
    split_complex = split_complex,
    multi_gene_sep = multi_gene_sep,
    combine_multi = combine_multi,
    weighting = weighting,
    weights = weights,
    kegg_sep_regex = kegg_sep_regex,
    return_lr_scored = return_lr_scored
  )
  
  if (isTRUE(save)) {
    out_file <- if (cell1 == cell2) {
      file.path(output_dir, paste0(cell1, "_autocrine_pathway_activity.csv"))
    } else {
      file.path(output_dir, paste0(cell1, "_to_", cell2, "_pathway_activity.csv"))
    }
    readr::write_csv(res$pathway_scores, out_file)
    message("✅ 已保存：", out_file)
  }
  return(res)
}

# ------------------------------------------------------------
# Directory-level batch wrapper:
#   Automatically match communication files and score files and
#   compute pathway activities in batch.
#   - pattern: by default matches "_merged_LR_RTF_TFT.csv"
#   - score files are assumed to be in score_dir (default = input_dir)
#   - suffix for score files: "_gene_importance_scores_transformed.csv"
# ------------------------------------------------------------
kegg_pathway_activity_dir <- function(
    input_dir,
    output_dir,
    pattern = "_merged_LR_RTF_TFT\\.csv$",
    score_dir = NULL,
    ligand_suffix = "_gene_importance_scores_transformed.csv",
    target_suffix = ligand_suffix,
    kegg_name_col = NULL,
    split_complex = FALSE,
    multi_gene_sep = ",",
    combine_multi = c("sum","mean","max"),
    weighting = c("count","none"),
    weights = list(ligand=1, receptor=1, tf=1, target=1),
    kegg_sep_regex = ";\\s*",
    score_gene_col = "Gene",
    score_val_col  = NULL,
    skip_missing = TRUE,
    return_lr_scored = FALSE
) {
  input_dir  <- normalizePath(input_dir, mustWork = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  score_dir  <- if (is.null(score_dir)) input_dir else normalizePath(score_dir, mustWork = TRUE)
  
  comm_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (!length(comm_files)) {
    warning("在 ", input_dir, " 下未找到匹配通讯文件：", pattern)
    return(invisible(list()))
  }
  
  results <- purrr::map(comm_files, function(f) {
    fname <- basename(f)
    if (str_detect(fname, "_to_")) {
      cell1 <- str_extract(fname, "^[^_]+")
      cell2 <- str_extract(fname, "(?<=_to_)[^_]+")
    } else {
      cell1 <- str_remove(fname, "_merged_LR_RTF_TFT\\.csv$")
      cell2 <- cell1
    }
    lig_file <- file.path(score_dir, paste0(cell1, ligand_suffix))
    tar_file <- file.path(score_dir, paste0(cell2, target_suffix))
    
    if (!file.exists(lig_file) || !file.exists(tar_file)) {
      msg <- paste0("得分文件缺失：",
                    ifelse(!file.exists(lig_file), basename(lig_file), ""),
                    " ",
                    ifelse(!file.exists(tar_file), basename(tar_file), ""))
      if (skip_missing) {
        message("⏭️ 跳过 ", basename(f), " —— ", msg)
        return(NULL)
      } else stop("❌ ", msg)
    }
    
    kegg_pathway_activity_for_file(
      comm_file         = f,
      ligand_score_file = lig_file,
      target_score_file = tar_file,
      output_dir        = output_dir,
      kegg_name_col     = kegg_name_col,
      split_complex     = split_complex,
      multi_gene_sep    = multi_gene_sep,
      combine_multi     = combine_multi,
      weighting         = weighting,
      weights           = weights,
      kegg_sep_regex    = kegg_sep_regex,
      score_gene_col    = score_gene_col,
      score_val_col     = score_val_col,
      save              = TRUE,
      return_lr_scored  = return_lr_scored
    )
  })
  
  names(results) <- basename(comm_files)
  results <- results[!vapply(results, is.null, logical(1))]
  invisible(results)
}

# ------------------------------------------------------------
# Filter pathway activity results by keyword(s)
#
# - You can pass 'keywords' as a character vector, or
#   'keywords_file' where each line is one keyword.
# - Matching is case-insensitive by default.
# - If name_col = NULL, the first column is treated as the pathway name.
# ------------------------------------------------------------
filter_kegg_activity_by_keywords <- function(
    input_dir,
    output_dir,
    pattern = "_pathway_activity\\.csv$",
    keywords = NULL,
    keywords_file = NULL,
    case_insensitive = TRUE,
    name_col = NULL
) {
  input_dir  <- normalizePath(input_dir, mustWork = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Read keywords
  if (is.null(keywords)) {
    if (is.null(keywords_file)) stop("请提供 keywords 或 keywords_file。")
    kw <- readr::read_lines(keywords_file)
  } else kw <- keywords
  kw <- kw %>% stringr::str_trim() %>% discard(~ .x == "") %>% unique()
  if (!length(kw)) stop("关键词为空。")
  
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) {
    warning("在 ", input_dir, " 下未找到匹配文件：", pattern)
    return(invisible(NULL))
  }
  
  for (f in files) {
    df <- readr::read_csv(f, show_col_types = FALSE)
    nm_col <- if (is.null(name_col)) names(df)[1] else name_col
    if (!nm_col %in% names(df)) {
      message("⚠️ 跳过（找不到通路名列）: ", basename(f)); next
    }
    pat <- stringr::str_c(if (case_insensitive) tolower(kw) else kw, collapse = "|")
    
    matched <- df %>%
      dplyr::filter(
        stringr::str_detect(
          if (case_insensitive) tolower(.data[[nm_col]]) else .data[[nm_col]],
          pat
        )
      )
    
    if (nrow(matched) > 0) {
      out_f <- file.path(output_dir, basename(f))
      readr::write_csv(matched, out_f)
      message("✅ 已筛出并保存：", out_f)
    } else {
      message("— 无匹配：", basename(f))
    }
  }
  invisible(TRUE)
}
