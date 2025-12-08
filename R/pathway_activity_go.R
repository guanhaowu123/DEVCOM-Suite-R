# R language
## Pathway score ##
## GO ##

# go_pathway_activity.R
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

# ------------------------------------------------------------
# Utility: read gene importance score file -> unify two columns: Gene, Score
# gene_col can be a column name or index; score_col can be a column name or
# index. By default, it first looks for a column named "Score"; otherwise it
# uses the 3rd column.
# ------------------------------------------------------------
load_gene_scores <- function(file, gene_col = "Gene", score_col = NULL) {
  stopifnot(file.exists(file))
  df <- readr::read_csv(file, show_col_types = FALSE)
  
  # handle gene column
  if (is.numeric(gene_col)) {
    gname <- names(df)[[gene_col]]
  } else if (is.character(gene_col) && gene_col %in% names(df)) {
    gname <- gene_col
  } else {
    # tolerant fallback: try lower-case variants
    cand <- names(df)[tolower(names(df)) %in% c("gene", "symbol", "genes")]
    if (length(cand) == 0) stop("无法在 ", file, " 中识别基因列（gene_col）。")
    gname <- cand[[1]]
  }
  
  # handle score column
  if (is.null(score_col)) {
    if ("Score" %in% names(df)) {
      sname <- "Score"
    } else {
      if (ncol(df) < 3) stop("无法定位得分列（score_col），且文件列数 < 3：", file)
      sname <- names(df)[[3]]
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
    dplyr::group_by(Gene) %>%                 # if duplicated genes, use mean (can be changed to sum/max)
    dplyr::summarise(Score = mean(Score, na.rm = TRUE), .groups = "drop")
  
  return(out)
}

# ------------------------------------------------------------
# Utility: for a string like "A,B,C", look up scores in named vector score_map
# and aggregate them.
# combine = "sum" / "mean" / "max"
# ------------------------------------------------------------
.sum_scores_for_string <- function(gene_str, score_map, sep = ",", combine = c("sum","mean","max")) {
  if (is.na(gene_str) || gene_str == "") return(0)
  combine <- match.arg(combine)
  genes <- str_split(gene_str, sep, simplify = TRUE)
  genes <- str_trim(genes[genes != ""])
  if (length(genes) == 0) return(0)
  vals <- as.numeric(score_map[genes])
  vals[is.na(vals)] <- 0
  switch(combine,
         sum  = sum(vals),
         mean = mean(vals),
         max  = ifelse(length(vals)==0, 0, max(vals)))
}

# ------------------------------------------------------------
# Core function: compute GO pathway activity for a single communication table
#
# Arguments:
#   lr_df            : data.frame with columns L, R, TF, Target, Name, Ontology
#   ligand_scores_df : sender gene scores (two columns: Gene, Score)
#   target_scores_df : receiver gene scores (two columns: Gene, Score)
#   ontology         : "BP" (default) / "CC" / "MF" or a vector like c("BP","MF");
#                      "ALL" means no ontology filtering
#   multi_gene_sep   : separator for multiple genes in L/R (default ",")
#   combine_multi    : aggregation method for multi-gene complexes ("sum" default)
#   weighting        : pathway weighting: "count" = multiply by number of entries;
#                      "none" = no extra weighting
#   return_lr_scored : whether to return the scored row-level table as well
#
# Returns: list(pathway_scores = ..., lr_scored = optional)
# ------------------------------------------------------------
compute_go_pathway_activity <- function(
    lr_df,
    ligand_scores_df,
    target_scores_df,
    ontology = "BP",
    multi_gene_sep = ",",
    combine_multi = c("sum","mean","max"),
    weighting = c("count","none"),
    return_lr_scored = FALSE
) {
  combine_multi <- match.arg(combine_multi)
  weighting     <- match.arg(weighting)
  
  needed_cols <- c("L","R","TF","Target","Name","Ontology")
  missing_cols <- setdiff(needed_cols, names(lr_df))
  if (length(missing_cols) > 0)
    stop("通讯表缺少必要列：", paste(missing_cols, collapse = ", "))
  
  # filter by ontology
  if (!identical(ontology, "ALL")) {
    lr_df <- lr_df %>% dplyr::filter(.data$Ontology %in% ontology)
  }
  if (nrow(lr_df) == 0) {
    return(list(pathway_scores = tibble::tibble(
      Name = character(), Pathway_Score = numeric(), Weight = numeric(),
      Final_Score = numeric(), CCI_Pairs = character()
    )))
  }
  
  # named vectors for fast lookup
  lig_map <- setNames(ligand_scores_df$Score, ligand_scores_df$Gene)
  tar_map <- setNames(target_scores_df$Score, target_scores_df$Gene)
  
  # score each L–R–TF–Target row
  lr_scored <- lr_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      ligand_score   = .sum_scores_for_string(.data$L,  lig_map, sep = multi_gene_sep, combine = combine_multi),
      receptor_score = .sum_scores_for_string(.data$R,  tar_map, sep = multi_gene_sep, combine = combine_multi),
      tf_score       = { v <- as.numeric(tar_map[as.character(.data$TF)]); ifelse(is.na(v), 0, v) },
      target_score   = { v <- as.numeric(tar_map[as.character(.data$Target)]); ifelse(is.na(v), 0, v) },
      total_score    = ligand_score + receptor_score + tf_score + target_score,
      CCI_Pair       = paste(
        paste0(na.omit(c(.data$L)), collapse = "+"),
        paste0(na.omit(c(.data$R)), collapse = "+"),
        .data$TF, .data$Target, sep = "->"
      )
    ) %>%
    dplyr::ungroup()
  
  # aggregate at GO pathway level
  pathway_scores <- lr_scored %>%
    dplyr::group_by(.data$Name) %>%
    dplyr::summarise(
      Pathway_Score = sum(.data$total_score),
      Weight        = dplyr::n(),
      CCI_Pairs     = paste(.data$CCI_Pair, collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Final_Score = dplyr::case_when(
        weighting == "count" ~ .data$Pathway_Score * .data$Weight,
        TRUE ~ .data$Pathway_Score
      )
    ) %>%
    dplyr::arrange(dplyr::desc(.data$Final_Score))
  
  out <- list(pathway_scores = pathway_scores)
  if (isTRUE(return_lr_scored)) out$lr_scored <- lr_scored
  return(out)
}

# ------------------------------------------------------------
# File-level wrapper: given three file paths, compute pathway scores
# and (optionally) write results.
#
# Output file naming follows the original logic:
#   autocrine: <cell>_GO_BP_autocrine_pathway_activity.csv
#   pairwise : <cell1>_to_<cell2>_GO_BP_pathway_activity.csv
# ------------------------------------------------------------
go_pathway_activity_for_file <- function(
    comm_file,
    ligand_score_file,
    target_score_file,
    output_dir,
    ontology = "BP",
    multi_gene_sep = ",",
    combine_multi = c("sum","mean","max"),
    weighting = c("count","none"),
    score_gene_col = "Gene",
    score_val_col  = NULL,  # NULL = auto detect; same as using 3rd column
    save = TRUE,
    return_lr_scored = FALSE
) {
  stopifnot(file.exists(comm_file))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # infer cell names from file name (same rule as your previous script)
  fname <- basename(comm_file)
  if (str_detect(fname, "_to_")) {
    cell1 <- str_extract(fname, "^[^_]+")
    cell2 <- str_extract(fname, "(?<=_to_)[^_]+")
  } else {
    cell1 <- str_remove(fname, "_GO_L_R_TF_Target\\.csv$")
    cell2 <- cell1
  }
  
  lr_df          <- readr::read_csv(comm_file, show_col_types = FALSE)
  ligand_scores  <- load_gene_scores(ligand_score_file, gene_col = score_gene_col, score_col = score_val_col)
  target_scores  <- load_gene_scores(target_score_file, gene_col = score_gene_col, score_col = score_val_col)
  
  res <- compute_go_pathway_activity(
    lr_df, ligand_scores, target_scores,
    ontology = ontology,
    multi_gene_sep = multi_gene_sep,
    combine_multi = combine_multi,
    weighting = weighting,
    return_lr_scored = return_lr_scored
  )
  
  if (isTRUE(save)) {
    suffix <- if (identical(ontology, "ALL")) "GO_pathway_activity" else paste0("GO_", paste(ontology, collapse="+"), "_pathway_activity")
    out_file <- if (cell1 == cell2) {
      file.path(output_dir, paste0(cell1, "_", suffix, "_autocrine.csv"))
    } else {
      file.path(output_dir, paste0(cell1, "_to_", cell2, "_", suffix, ".csv"))
    }
    readr::write_csv(res$pathway_scores, out_file)
    message("✅ 已保存：", out_file)
  }
  return(res)
}

# ------------------------------------------------------------
# Directory-level batch wrapper: automatically match file names and
# compute pathway scores in bulk.
#
# Arguments:
#   input_dir      : directory containing communication files
#   output_dir     : directory for output files
#   pattern        : regex pattern for communication files
#                    (default matches your naming scheme)
#   score_dir      : directory containing gene score files
#                    (default: same as input_dir)
#   ligand_suffix  : suffix for sender score files
#                    (default "_gene_importance_scores_transformed.csv")
#   target_suffix  : suffix for receiver score files (default same as above)
#   skip_missing   : if TRUE, skip when any score file is missing;
#                    if FALSE, stop with error
#
# Returns:
#   A named list of results (one per file).
# ------------------------------------------------------------
go_pathway_activity_dir <- function(
    input_dir,
    output_dir,
    ontology = "BP",
    pattern = "_GO_L_R_TF_Target\\.csv$",
    score_dir = NULL,
    ligand_suffix = "_gene_importance_scores_transformed.csv",
    target_suffix = ligand_suffix,
    multi_gene_sep = ",",
    combine_multi = c("sum","mean","max"),
    weighting = c("count","none"),
    score_gene_col = "Gene",
    score_val_col  = NULL,
    skip_missing = TRUE,
    return_lr_scored = FALSE
) {
  input_dir  <- normalizePath(input_dir, mustWork = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  score_dir  <- if (is.null(score_dir)) input_dir else normalizePath(score_dir, mustWork = TRUE)
  
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    warning("在 ", input_dir, " 下未找到匹配通讯文件：", pattern)
    return(invisible(list()))
  }
  
  results <- purrr::map(files, function(f) {
    fname <- basename(f)
    if (str_detect(fname, "_to_")) {
      cell1 <- str_extract(fname, "^[^_]+")
      cell2 <- str_extract(fname, "(?<=_to_)[^_]+")
    } else {
      cell1 <- str_remove(fname, "_GO_L_R_TF_Target\\.csv$")
      cell2 <- cell1
    }
    ligand_score_file <- file.path(score_dir, paste0(cell1, ligand_suffix))
    target_score_file <- file.path(score_dir, paste0(cell2, target_suffix))
    
    if (!file.exists(ligand_score_file) || !file.exists(target_score_file)) {
      msg <- paste0("得分文件缺失：", ifelse(!file.exists(ligand_score_file), basename(ligand_score_file), ""),
                    " ", ifelse(!file.exists(target_score_file), basename(target_score_file), ""))
      if (skip_missing) {
        message("⏭️ 跳过 ", basename(f), " —— ", msg)
        return(NULL)
      } else {
        stop("❌ ", msg)
      }
    }
    
    go_pathway_activity_for_file(
      comm_file = f,
      ligand_score_file = ligand_score_file,
      target_score_file = target_score_file,
      output_dir = output_dir,
      ontology = ontology,
      multi_gene_sep = multi_gene_sep,
      combine_multi = combine_multi,
      weighting = weighting,
      score_gene_col = score_gene_col,
      score_val_col  = score_val_col,
      save = TRUE,
      return_lr_scored = return_lr_scored
    )
  })
  
  # drop NULLs (files that were skipped due to missing scores)
  names(results) <- basename(files)
  results <- results[!vapply(results, is.null, logical(1))]
  invisible(results)
}
