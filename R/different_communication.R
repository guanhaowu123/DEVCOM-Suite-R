#' Differential communication analysis between two conditions
#'
#' Performs differential cell–cell communication analysis between two
#' conditions (e.g. IVF as control and SCNT as case), comparing
#' ligand–receptor (LR) activities across all sender–receiver cell-type pairs.
#'
#' The output tables only contain the following columns:
#'   - mode
#'   - sender
#'   - receiver
#'   - ligand_subunits
#'   - receptor_subunits
#'   - delta_log10_score
#'   - delta_log2_score
#'   - p_perm_fisher_BH
#'
#' where p_perm_fisher_BH is the BH-adjusted Fisher combined p-value
#' based on permutation tests on ligand and receptor expression.
#'
#' @param cell_types Character vector of cell-type names. Each name must match
#'   the prefix of the expression/score CSV files.
#' @param ctrl_base Directory of the control group expression/score files,
#'   e.g. \code{"/path/IVF/filtered_gene_expression"}.
#' @param case_base Directory of the case group expression/score files,
#'   e.g. \code{"/path/SCNT/filtered_gene_expression"}.
#' @param lr_pairs_path Path to the LR definition CSV. Column names may contain
#'   L1/L2/... and R1/R2/... for ligand and receptor subunits.
#' @param out_dir Output directory. If \code{NULL}, it will be constructed as
#'   \code{"<case_label>_vs_<ctrl_label>_diff_comm"}.
#' @param n_perm Number of permutations for the permutation test.
#' @param use_perm Logical; whether to perform permutation-based testing.
#' @param case_label Label of the case condition, used in filenames.
#' @param ctrl_label Label of the control condition, used in filenames.
#' @param seed Random seed. Set to \code{NULL} to avoid resetting the seed.
#'
#' @return A list with elements \code{all}, \code{autocrine},
#'   \code{interaction}, \code{significant} and \code{out_dir}. Each table
#'   only contains the 8 columns described above.
#' @export
#'
#' @importFrom stats pchisq p.adjust t.test wilcox.test
#' @importFrom utils read.csv write.csv
#' @importFrom tibble tibble
#' @importFrom dplyr distinct rowwise ungroup mutate select bind_rows arrange filter
#' @importFrom dplyr c_across
#' @importFrom purrr map_dfr pmap_dfr
run_diff_comm <- function(cell_types,
                          ctrl_base,
                          case_base,
                          lr_pairs_path,
                          out_dir   = NULL,
                          n_perm    = 1000L,
                          use_perm  = TRUE,
                          case_label = "case",
                          ctrl_label = "control",
                          seed      = 2025L) {
  
  ## ----- basic options -----
  options(stringsAsFactors = FALSE)
  if (!is.null(seed)) set.seed(seed)
  
  if (is.null(out_dir)) {
    out_dir <- paste0(case_label, "_vs_", ctrl_label, "_diff_comm")
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## ================= helper functions =================
  canonize_names <- function(x){
    x <- sub("^\ufeff", "", x)
    x <- trimws(x)
    x <- tolower(x)
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }
  
  read_expr_mat <- function(file){
    df  <- read.csv(file, check.names = FALSE)
    cn  <- colnames(df)
    cnc <- canonize_names(cn)
    gene_idx <- which(cnc %in% c("gene","symbol","genes","gene_symbol","geneid","gene_id"))
    if (length(gene_idx) == 0) gene_idx <- 1
    colnames(df)[gene_idx[1]] <- "gene"
    
    df$gene <- trimws(df$gene)
    df <- dplyr::distinct(df, gene, .keep_all = TRUE)
    rownames(df) <- df$gene
    df$gene <- NULL
    
    df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
    as.matrix(df)
  }
  
  read_score_vec <- function(file){
    df  <- read.csv(file, check.names = FALSE)
    cn  <- colnames(df)
    cnc <- canonize_names(cn)
    
    gene_idx <- which(cnc %in% c("gene","symbol","genes","gene_symbol","geneid","gene_id"))
    if (length(gene_idx) == 0) gene_idx <- 1
    colnames(df)[gene_idx[1]] <- "gene"
    
    ## priority 1: log10_Score (or equivalent)
    score_idx <- grep("^log10[_-]?score$", cnc, perl = TRUE)
    
    ## priority 2: Importance_Score (or equivalent)
    if (length(score_idx) == 0) {
      score_idx <- grep("^importance[_-]?score$", cnc, perl = TRUE)
    }
    
    if (length(score_idx) == 0) {
      stop(sprintf(
        "Could not find 'log10_Score' or 'Importance_Score' column in score file: %s",
        file
      ))
    }
    
    colnames(df)[score_idx[1]] <- "log10_Score"
    
    tibble::tibble(
      gene        = trimws(df$gene),
      log10_Score = as.numeric(df$log10_Score)
    ) |>
      dplyr::distinct(gene, .keep_all = TRUE)
  }
  
  ## Fisher's method for combining p-values
  fisher_p <- function(pvals){
    pvals <- pvals[is.finite(pvals) & !is.na(pvals)]
    if (!length(pvals)) return(NA_real_)
    stat <- -2 * sum(log(pmax(pvals, .Machine$double.xmin)))
    stats::pchisq(stat, df = 2 * length(pvals), lower.tail = FALSE)
  }
  
  ## Two-sample permutation test on mean difference (log1p transformed)
  perm_p_two_sample <- function(x_case, x_ctrl, B = n_perm){
    x_case <- x_case[is.finite(x_case)]
    x_ctrl <- x_ctrl[is.finite(x_ctrl)]
    if (length(x_case) < 2 || length(x_ctrl) < 2) return(NA_real_)
    xi <- log1p(x_case); xs <- log1p(x_ctrl)
    n1 <- length(xi); z <- c(xi, xs); m <- length(z)
    obs <- mean(xi) - mean(xs)
    cnt <- 0L
    for (b in seq_len(B)){
      idx   <- sample.int(m, n1)
      diffb <- mean(z[idx]) - mean(z[-idx])
      if (abs(diffb) >= abs(obs)) cnt <- cnt + 1L
    }
    (cnt + 1) / (B + 1)
  }
  
  ## For a multi-subunit complex, take max expression across genes per sample
  complex_expr_max <- function(genes, mat){
    g <- intersect(genes, rownames(mat))
    out <- rep(NA_real_, ncol(mat)); names(out) <- colnames(mat)
    if (!length(g)) return(out)
    vals <- mat[g, , drop = FALSE]
    apply(vals, 2, function(v){
      if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    })
  }
  
  ## For a multi-subunit complex, take max log10 score across genes
  complex_score_max_log <- function(genes, score_tbl){
    s <- score_tbl |>
      dplyr::filter(gene %in% genes) |>
      dplyr::pull(log10_Score)
    if (!length(s)) return(NA_real_)
    val <- suppressWarnings(max(s, na.rm = TRUE))
    if (!is.finite(val)) NA_real_ else val
  }
  
  load_condition_data <- function(base_dir, cell_types){
    expr_list  <- vector("list", length(cell_types))
    score_list <- vector("list", length(cell_types))
    names(expr_list)  <- cell_types
    names(score_list) <- cell_types
    
    for (ct in cell_types){
      f_expr  <- file.path(base_dir, sprintf("%s_filtered_gene_expression.csv", ct))
      f_score <- file.path(base_dir, sprintf("%s_gene_importance_scores_transformed.csv", ct))
      if (!file.exists(f_expr))  stop("Expression matrix not found: ", f_expr)
      if (!file.exists(f_score)) stop("Score file not found: ", f_score)
      expr_list[[ct]]  <- read_expr_mat(f_expr)
      score_list[[ct]] <- read_score_vec(f_score)
    }
    list(expr = expr_list, score = score_list)
  }
  
  sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
  
  ## ================= data loading =================
  message("Loading control condition (", ctrl_label, ") ...")
  ctrl <- load_condition_data(ctrl_base, cell_types)
  
  message("Loading case condition (", case_label, ") ...")
  case <- load_condition_data(case_base, cell_types)
  
  ## ---- LR definition table ----
  lr_pairs_raw <- read.csv(lr_pairs_path, check.names = FALSE)
  cn_raw <- colnames(lr_pairs_raw)
  cn_can <- canonize_names(cn_raw)
  
  lig_cols <- cn_raw[grepl("^l[0-9]+$", cn_can)]
  rec_cols <- cn_raw[grepl("^r[0-9]+$", cn_can)]
  
  if (length(lig_cols) == 0 || length(rec_cols) == 0) {
    stop("No L*/R* columns (e.g. L1/L2, R1/R2/R3) were found in the LR table.")
  }
  
  lr_pairs <- lr_pairs_raw |>
    dplyr::rowwise() |>
    dplyr::mutate(
      ligand_subunits = {
        v <- dplyr::c_across(dplyr::all_of(lig_cols))
        v <- v[!is.na(v) & nzchar(trimws(v))]
        if (length(v) == 0) NA_character_ else paste(unique(trimws(v)), collapse = ";")
      },
      receptor_subunits = {
        v <- dplyr::c_across(dplyr::all_of(rec_cols))
        v <- v[!is.na(v) & nzchar(trimws(v))]
        if (length(v) == 0) NA_character_ else paste(unique(trimws(v)), collapse = ";")
      },
      ligand = ifelse(is.na(ligand_subunits), NA_character_,
                      strsplit(ligand_subunits, ";")[[1]][1]),
      receptor = ifelse(is.na(receptor_subunits), NA_character_,
                        strsplit(receptor_subunits, ";")[[1]][1])
    ) |>
    dplyr::ungroup() |>
    dplyr::select(ligand, receptor, ligand_subunits, receptor_subunits) |>
    dplyr::distinct(ligand, receptor, ligand_subunits, receptor_subunits, .keep_all = TRUE)
  
  ## ================= core analysis for one sender/receiver pair =================
  analyze_sender_receiver <- function(sender, receiver,
                                      mode_label = c("autocrine","interaction")){
    mode_label <- match.arg(mode_label)
    
    S_ctrl <- ctrl$expr[[sender]];   S_case <- case$expr[[sender]]
    R_ctrl <- ctrl$expr[[receiver]]; R_case <- case$expr[[receiver]]
    S_ctrl_score <- ctrl$score[[sender]];   S_case_score <- case$score[[sender]]
    R_ctrl_score <- ctrl$score[[receiver]]; R_case_score <- case$score[[receiver]]
    
    results <- vector("list", nrow(lr_pairs))
    
    for (i in seq_len(nrow(lr_pairs))){
      L  <- lr_pairs$ligand[i]
      R  <- lr_pairs$receptor[i]
      Lc <- lr_pairs$ligand_subunits[i]
      Rc <- lr_pairs$receptor_subunits[i]
      
      L_genes <- if (!is.na(Lc) && nzchar(Lc)) strsplit(Lc, ";")[[1]] |> trimws() else L
      R_genes <- if (!is.na(Rc) && nzchar(Rc)) strsplit(Rc, ";")[[1]] |> trimws() else R
      
      ## log10 scores: max across subunits (case/control)
      L_case_log10 <- complex_score_max_log(L_genes, S_case_score)
      L_ctrl_log10 <- complex_score_max_log(L_genes, S_ctrl_score)
      R_case_log10 <- complex_score_max_log(R_genes, R_case_score)
      R_ctrl_log10 <- complex_score_max_log(R_genes, R_ctrl_score)
      
      LR_case_log10 <- mean(c(L_case_log10, R_case_log10), na.rm = TRUE)
      LR_ctrl_log10 <- mean(c(L_ctrl_log10, R_ctrl_log10), na.rm = TRUE)
      if (!is.finite(LR_case_log10)) LR_case_log10 <- NA_real_
      if (!is.finite(LR_ctrl_log10)) LR_ctrl_log10 <- NA_real_
      
      delta_log10 <- if (is.na(LR_case_log10) | is.na(LR_ctrl_log10))
        NA_real_ else LR_case_log10 - LR_ctrl_log10
      delta_log2  <- if (is.na(delta_log10)) NA_real_ else delta_log10 * log2(10)
      
      ## expression: max across subunits for each sample
      eL_case <- complex_expr_max(L_genes, S_case)
      eL_ctrl <- complex_expr_max(L_genes, S_ctrl)
      eR_case <- complex_expr_max(R_genes, R_case)
      eR_ctrl <- complex_expr_max(R_genes, R_ctrl)
      
      ## permutation tests on expression for L and R, then combine via Fisher
      p_perm_L <- if (use_perm) perm_p_two_sample(eL_case, eL_ctrl, B = n_perm) else NA_real_
      p_perm_R <- if (use_perm) perm_p_two_sample(eR_case, eR_ctrl, B = n_perm) else NA_real_
      p_perm_F <- if (all(is.na(c(p_perm_L, p_perm_R)))) NA_real_ else
        fisher_p(c(p_perm_L, p_perm_R))
      
      ## Output subunit labels: use complex string if available, else single gene
      ligand_subunits_out   <- if (!is.na(Lc) && nzchar(Lc)) Lc else L
      receptor_subunits_out <- if (!is.na(Rc) && nzchar(Rc)) Rc else R
      
      results[[i]] <- tibble::tibble(
        mode               = mode_label,
        sender             = sender,
        receiver           = receiver,
        ligand_subunits    = ligand_subunits_out,
        receptor_subunits  = receptor_subunits_out,
        delta_log10_score  = delta_log10,
        delta_log2_score   = delta_log2,
        p_perm_fisher      = p_perm_F
      )
    }
    dplyr::bind_rows(results)
  }
  
  ## ================= full traversal: autocrine + interactions =================
  message("Analyzing autocrine (self) communications ...")
  raw_auto <- purrr::map_dfr(
    cell_types,
    ~ analyze_sender_receiver(.x, .x, mode_label = "autocrine")
  )
  
  message("Analyzing interactions (sender != receiver) ...")
  inter_pairs <- expand.grid(sender = cell_types,
                             receiver = cell_types,
                             stringsAsFactors = FALSE) |>
    dplyr::filter(sender != receiver)
  raw_inter <- purrr::pmap_dfr(
    list(inter_pairs$sender, inter_pairs$receiver),
    ~ analyze_sender_receiver(..1, ..2, mode_label = "interaction")
  )
  
  ## Combine all
  res <- dplyr::bind_rows(raw_auto, raw_inter)
  
  ## ================= multiple testing (BH) =================
  ## Adjust permutation-based Fisher p-values
  pvec <- res$p_perm_fisher
  res$p_perm_fisher_BH <- ifelse(
    is.finite(pvec),
    stats::p.adjust(pvec, method = "BH"),
    NA_real_
  )
  ## Drop raw permutation p if you don't want to keep it
  res$p_perm_fisher <- NULL
  
  res <- res |>
    dplyr::arrange(mode, sender, receiver, ligand_subunits, receptor_subunits)
  
  ## Split back into autocrine / interaction (with BH).
  auto_res <- res |>
    dplyr::filter(mode == "autocrine")
  inter_res <- res |>
    dplyr::filter(mode == "interaction")
  
  ## ================= write output files =================
  case_vs_ctrl <- paste0(case_label, "_vs_", ctrl_label)
  
  outfile_all   <- file.path(
    out_dir,
    paste0(case_vs_ctrl, "_diff_comm_autocrine_interaction.csv")
  )
  outfile_auto  <- file.path(
    out_dir,
    paste0(case_vs_ctrl, "_diff_comm_autocrine_only.csv")
  )
  outfile_inter <- file.path(
    out_dir,
    paste0(case_vs_ctrl, "_diff_comm_interaction_only.csv")
  )
  
  write.csv(res,        outfile_all,   row.names = FALSE)
  write.csv(auto_res,   outfile_auto,  row.names = FALSE)
  write.csv(inter_res,  outfile_inter, row.names = FALSE)
  
  message("Finished. Output files:")
  message("  - ", outfile_all)
  message("  - ", outfile_auto)
  message("  - ", outfile_inter)
  message("Number of permutations (n_perm) = ", n_perm)
  
  ## ---------------- 1) autocrine: one file per cell type ----------------
  dir_auto_split <- file.path(out_dir, "autocrine_by_cell")
  dir.create(dir_auto_split, showWarnings = FALSE, recursive = TRUE)
  for (ct in unique(auto_res$sender)) {
    sub <- auto_res |>
      dplyr::filter(sender == ct) |>
      dplyr::arrange(receptor_subunits, ligand_subunits)
    write.csv(
      sub,
      file.path(dir_auto_split, paste0(sanitize(ct), "_autocrine.csv")),
      row.names = FALSE
    )
  }
  
  ## ---------------- 2) interactions: one file per sender->receiver pair ------
  dir_inter_split <- file.path(out_dir, "interaction_by_pair")
  dir.create(dir_inter_split, showWarnings = FALSE, recursive = TRUE)
  pairs <- inter_res |>
    dplyr::distinct(sender, receiver) |>
    dplyr::arrange(sender, receiver)
  for (k in seq_len(nrow(pairs))) {
    s <- pairs$sender[k]; r <- pairs$receiver[k]
    sub <- inter_res |>
      dplyr::filter(sender == s, receiver == r) |>
      dplyr::arrange(receptor_subunits, ligand_subunits)
    write.csv(
      sub,
      file.path(dir_inter_split, paste0(sanitize(s), "_to_", sanitize(r), ".csv")),
      row.names = FALSE
    )
  }
  
  ## ---------------- 3) one combined table per cell type ----------------------
  ## (sender or receiver = ct); still only 8 columns.
  dir_bycell <- file.path(out_dir, "by_cell_all")
  dir.create(dir_bycell, showWarnings = FALSE, recursive = TRUE)
  all_cells <- sort(unique(c(res$sender, res$receiver)))
  for (ct in all_cells) {
    sub <- res |>
      dplyr::filter(sender == ct | receiver == ct) |>
      dplyr::arrange(mode, sender, receiver, receptor_subunits, ligand_subunits)
    write.csv(
      sub,
      file.path(dir_bycell, paste0(sanitize(ct), "_ALL.csv")),
      row.names = FALSE
    )
  }
  
  ## ---------------- 4) significant results (permutation-based) --------------
  alpha <- 0.05
  res_sig <- res |>
    dplyr::filter(is.finite(p_perm_fisher_BH), p_perm_fisher_BH < alpha)
  
  write.csv(
    res_sig,
    file.path(out_dir,
              paste0(case_vs_ctrl, "_diff_comm_significant.csv")),
    row.names = FALSE
  )
  
  invisible(list(
    all          = res,
    autocrine    = auto_res,
    interaction  = inter_res,
    significant  = res_sig,
    out_dir      = out_dir
  ))
}
