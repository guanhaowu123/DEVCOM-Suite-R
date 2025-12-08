## R language
##
## Active and Inhibitive
##
## Three modes:
##   - no p-value and no log2FC
##   - with p-value
##   - with p-value and log2FC

#' Compute activating / inhibiting cell–cell communication
#'
#' This function computes "activation" and "inhibition" type cell–cell
#' communication between pairs of cell types based on expression matrices
#' and a ligand–receptor prior knowledge base. Optionally, it can perform
#' statistical tests and report log2 fold changes.
#'
#' @param expr_dir   Character. Directory containing expression matrices.
#'   Each cell type should have a file named \code{"<cell>_gene_expression.csv"},
#'   where the first column is \code{gene} and the remaining columns are
#'   sample expression values.
#' @param kb_path    Character. Path to the prior knowledge base CSV.
#'   It should either contain L*/R* columns (e.g. L1/L2, R1/R2) or already
#'   flattened \code{ligand} and \code{receptor} columns (optionally \code{tf}).
#' @param output_dir Character or \code{NULL}. If non-\code{NULL}, results will
#'   be written to this directory (directory structure is similar to the
#'   original scripts: one folder per sender–receiver pair, each with
#'   \code{activation.csv} and \code{inhibition.csv}).
#' @param test       Statistical test to use:
#'   \itemize{
#'     \item \code{"none"}: only determine direction (activation/inhibition),
#'           no statistical testing.
#'     \item \code{"wilcoxon"}: Wilcoxon rank-sum test, returns \code{p_value}.
#'     \item \code{"ttest"}: two-sample t-test, returns \code{p_value} and
#'           \code{log2fc}.
#'   }
#' @param min_expr   Numeric. Expression filter: if the mean expression of
#'   ligand in sender or receptor in receiver is below this value, the pair
#'   is skipped.
#' @param include_pairwise Logical. Whether to compute pairwise interactions
#'   between different cell types (sender != receiver). Typically \code{TRUE}.
#' @param ordered_pairs Logical. If \code{TRUE}, compute all ordered pairs
#'   (A→B and B→A). If \code{FALSE}, iterate over unique unordered pairs but
#'   still output both directions (A→B and B→A) as in the original scripts.
#' @param p_adjust   Character. Method for multiple-testing correction passed
#'   to \code{stats::p.adjust}. Use \code{"none"} for no correction, or e.g.
#'   \code{"BH"} for Benjamini–Hochberg.
#' @param p_min      Numeric. Lower bound for very small p-values to avoid
#'   writing zeros.
#' @param write_coverage Logical. If \code{TRUE}, output ligand/receptor
#'   expression coverage matrices and a summary of missing ligands/receptors
#'   in each cell type.
#' @param verbose    Logical. If \code{TRUE}, print progress messages.
#'
#' @return A list with a single element:
#' \describe{
#'   \item{$pairwise}{Named list. Each name is "<sender>_vs_<receiver>" and
#'   each element is a list with two components: \code{activation} and
#'   \code{inhibition}, both are data frames/tibbles with columns:
#'   \code{ligand}, \code{receptor}, \code{tf}, \code{p_value}, \code{log2fc},
#'   and optionally \code{p_adj} if multiple-testing correction is applied.}
#' }
#' If \code{output_dir} is not \code{NULL}, results are also saved as CSV
#' files.
#'
#' @export
compute_active_inhibit <- function(
    expr_dir,
    kb_path,
    output_dir = NULL,
    test = c("none", "wilcoxon", "ttest"),
    min_expr = 0,
    include_pairwise  = TRUE,
    ordered_pairs     = FALSE,
    p_adjust = c("none", "BH"),
    p_min = 1e-300,
    write_coverage = FALSE,
    verbose = TRUE
){
  suppressPackageStartupMessages(library(tidyverse))
  
  test     <- match.arg(test)
  p_adjust <- match.arg(p_adjust)
  
  .log <- function(...) if (isTRUE(verbose)) cat(sprintf(...), "\n")
  
  dir_create_safe <- function(path){
    if (!is.null(path)) dir.create(path, showWarnings = FALSE, recursive = TRUE)
  }
  if (!is.null(output_dir)) dir_create_safe(output_dir)
  
  # ===== 1) Read KB (support L*/R* or already flattened ligand/receptor) =====
  kb_raw <- readr::read_csv(kb_path, show_col_types = FALSE)
  if (!all(c("ligand","receptor") %in% names(kb_raw))){
    kb <- kb_raw %>%
      tidyr::pivot_longer(cols = tidyselect::starts_with("L"),
                          names_to = "L_col", values_to = "ligand") %>%
      tidyr::pivot_longer(cols = tidyselect::starts_with("R"),
                          names_to = "R_col", values_to = "receptor") %>%
      dplyr::filter(!is.na(ligand), !is.na(receptor)) %>%
      dplyr::select(ligand, receptor, dplyr::any_of("tf")) %>%
      dplyr::distinct()
  } else {
    kb <- kb_raw %>%
      dplyr::select(ligand, receptor, dplyr::any_of("tf")) %>%
      dplyr::mutate(tf = if ("tf" %in% names(.)) tf else NA_character_) %>%
      dplyr::distinct()
  }
  
  # ===== 2) Read expression matrices as numeric matrices =====
  expr_files <- list.files(expr_dir, pattern = "_gene_expression\\.csv$", full.names = TRUE)
  stopifnot(length(expr_files) > 0)
  cell_names <- gsub("_gene_expression\\.csv$", "", basename(expr_files))
  
  read_expr_as_matrix <- function(fp){
    df <- readr::read_csv(fp, show_col_types = FALSE)
    names(df)[1] <- "gene"
    # force numeric
    num <- suppressWarnings(as.data.frame(lapply(df[,-1], as.numeric)))
    mat <- as.matrix(num)
    rownames(mat) <- df$gene
    # remove rows that are all NA
    keep <- rowSums(!is.na(mat)) > 0
    mat[keep, , drop = FALSE]
  }
  
  expr_list <- purrr::map(expr_files, read_expr_as_matrix) %>%
    stats::setNames(cell_names)
  
  # ===== 3) Optional: write ligand/receptor coverage matrices and summary =====
  if (isTRUE(write_coverage) && !is.null(output_dir)) {
    ligands_all   <- unique(kb$ligand)
    receptors_all <- unique(kb$receptor)
    
    ligand_expr_stat <- matrix(0L, nrow = length(ligands_all), ncol = length(expr_list),
                               dimnames = list(ligands_all, names(expr_list)))
    receptor_expr_stat <- matrix(0L, nrow = length(receptors_all), ncol = length(expr_list),
                                 dimnames = list(receptors_all, names(expr_list)))
    for (cell in names(expr_list)) {
      genes <- rownames(expr_list[[cell]])
      ligand_expr_stat[, cell]   <- as.integer(ligands_all %in% genes)
      receptor_expr_stat[, cell] <- as.integer(receptors_all %in% genes)
    }
    readr::write_csv(
      as.data.frame(ligand_expr_stat) %>%
        tibble::rownames_to_column("ligand"),
      file.path(output_dir, "ligand_expression_matrix.csv")
    )
    readr::write_csv(
      as.data.frame(receptor_expr_stat) %>%
        tibble::rownames_to_column("receptor"),
      file.path(output_dir, "receptor_expression_matrix.csv")
    )
    summary_df <- tibble::tibble(
      cell_type = names(expr_list),
      missing_ligands   = colSums(ligand_expr_stat == 0),
      missing_receptors = colSums(receptor_expr_stat == 0)
    )
    readr::write_csv(summary_df,
                     file.path(output_dir, "expression_missing_summary.csv"))
  }
  
  # ===== 4) Statistical test helpers =====
  safe_wilcox_p <- function(x, y){
    if (length(x) >= 2 && length(y) >= 2 && all(is.finite(c(x,y)))) {
      p <- tryCatch(stats::wilcox.test(x, y)$p.value, error = function(e) NA_real_)
      p
    } else NA_real_
  }
  safe_ttest_p <- function(x, y){
    if (length(x) >= 2 && length(y) >= 2 && all(is.finite(c(x,y)))) {
      p <- tryCatch(stats::t.test(x, y)$p.value, error = function(e) NA_real_)
      p
    } else NA_real_
  }
  
  # ===== 5) Compare one direction sender->receiver using "classic rule" =====
  # Classic rule:
  #   activation: mean(lig_s) > mean(lig_r)  &  mean(rec_r) > mean(rec_s)
  #   inhibition: mean(lig_s) > mean(lig_r)  &  mean(rec_r) < mean(rec_s)
  compare_direction <- function(sender_mat, receiver_mat, kb, min_expr, test, p_min){
    genes_s <- rownames(sender_mat)
    genes_r <- rownames(receiver_mat)
    
    res_act <- vector("list", length = 256)
    res_inh <- vector("list", length = 256)
    ia <- 0L; ii <- 0L
    
    for (k in seq_len(nrow(kb))) {
      lig <- kb$ligand[k]; rec <- kb$receptor[k]
      tf  <- if ("tf" %in% names(kb)) kb$tf[k] else NA_character_
      
      # require all four vectors (lig_s, lig_r, rec_s, rec_r) to exist
      if (!(lig %in% genes_s && lig %in% genes_r && rec %in% genes_s && rec %in% genes_r)) next
      
      lig_s <- as.numeric(sender_mat[lig, , drop = TRUE])
      lig_r <- as.numeric(receiver_mat[lig, , drop = TRUE])
      rec_s <- as.numeric(sender_mat[rec, , drop = TRUE])
      rec_r <- as.numeric(receiver_mat[rec, , drop = TRUE])
      
      # expression filter
      if (mean(lig_s, na.rm = TRUE) < min_expr ||
          mean(rec_r, na.rm = TRUE) < min_expr) next
      
      lig_s_avg <- mean(lig_s, na.rm = TRUE)
      lig_r_avg <- mean(lig_r, na.rm = TRUE)
      rec_s_avg <- mean(rec_s, na.rm = TRUE)
      rec_r_avg <- mean(rec_r, na.rm = TRUE)
      
      direction <- NA_character_
      if (is.finite(lig_s_avg) && is.finite(lig_r_avg) &&
          is.finite(rec_s_avg) && is.finite(rec_r_avg)) {
        if (lig_s_avg > lig_r_avg && rec_r_avg > rec_s_avg) {
          direction <- "activation"
        } else if (lig_s_avg > lig_r_avg && rec_r_avg < rec_s_avg) {
          direction <- "inhibition"
        } else {
          next  # skip if no clear classification
        }
      } else next
      
      # statistical tests and log2FC (depending on user choice)
      p_value <- NA_real_
      log2fc  <- NA_real_
      
      if (test == "wilcoxon") {
        p_lig <- safe_wilcox_p(lig_s, lig_r)
        p_rec <- safe_wilcox_p(rec_r, rec_s) # rec_r vs rec_s
        p_value <- suppressWarnings(max(p_lig, p_rec, na.rm = TRUE))
      } else if (test == "ttest") {
        p_lig <- safe_ttest_p(lig_s, lig_r)
        p_rec <- safe_ttest_p(rec_r, rec_s)
        p_value <- suppressWarnings(max(p_lig, p_rec, na.rm = TRUE))
        # log2FC consistent with your third script
        log2fc <- log2(
          (mean(c(lig_s, rec_s), na.rm = TRUE) + 1e-6) /
            (mean(c(lig_r, rec_r), na.rm = TRUE) + 1e-6)
        )
      }
      
      if (is.finite(p_value)) p_value <- max(p_value, p_min, na.rm = TRUE)
      
      row <- tibble::tibble(
        ligand   = lig,
        receptor = rec,
        tf       = tf,
        p_value  = if (test == "none") NA_real_ else p_value,
        log2fc   = if (test == "ttest") log2fc else NA_real_
      )
      
      if (direction == "activation") {
        ia <- ia + 1L; res_act[[ia]] <- row
      } else {
        ii <- ii + 1L; res_inh[[ii]] <- row
      }
    }
    
    act <- dplyr::bind_rows(res_act[seq_len(ia)])
    inh <- dplyr::bind_rows(res_inh[seq_len(ii)])
    
    # multiple-testing correction (applied separately to activation and inhibition)
    if (p_adjust != "none" && test != "none") {
      if (nrow(act)) act$p_adj <- stats::p.adjust(act$p_value, method = p_adjust)
      if (nrow(inh)) inh$p_adj <- stats::p.adjust(inh$p_value, method = p_adjust)
    }
    
    list(activation = act, inhibition = inh)
  }
  
  # ===== 6) Compute pairwise interactions and optionally write to disk =====
  results <- list(pairwise = list())
  
  if (isTRUE(include_pairwise)) {
    if (isTRUE(ordered_pairs)) {
      # all ordered pairs (sender != receiver)
      for (c1 in names(expr_list)) {
        for (c2 in names(expr_list)) {
          if (c1 == c2) next
          .log("[pair] %s -> %s ...", c1, c2)
          res <- compare_direction(expr_list[[c1]], expr_list[[c2]],
                                   kb, min_expr, test, p_min)
          key <- paste0(c1, "_vs_", c2)
          results$pairwise[[key]] <- res
          if (!is.null(output_dir)) {
            dir_pair <- file.path(output_dir, key)
            dir_create_safe(dir_pair)
            readr::write_csv(res$activation, file.path(dir_pair, "activation.csv"))
            readr::write_csv(res$inhibition, file.path(dir_pair, "inhibition.csv"))
          }
        }
      }
    } else {
      # unique unordered pairs, but still output both directions (A->B and B->A)
      cn <- names(expr_list)
      if (length(cn) >= 2) {
        for (i in 1:(length(cn) - 1)) {
          for (j in (i + 1):length(cn)) {
            c1 <- cn[i]; c2 <- cn[j]
            
            .log("[pair] %s -> %s ...", c1, c2)
            res1 <- compare_direction(expr_list[[c1]], expr_list[[c2]],
                                      kb, min_expr, test, p_min)
            key1 <- paste0(c1, "_vs_", c2)
            results$pairwise[[key1]] <- res1
            if (!is.null(output_dir)) {
              dir1 <- file.path(output_dir, key1); dir_create_safe(dir1)
              readr::write_csv(res1$activation, file.path(dir1, "activation.csv"))
              readr::write_csv(res1$inhibition, file.path(dir1, "inhibition.csv"))
            }
            
            .log("[pair] %s -> %s ...", c2, c1)
            res2 <- compare_direction(expr_list[[c2]], expr_list[[c1]],
                                      kb, min_expr, test, p_min)
            key2 <- paste0(c2, "_vs_", c1)
            results$pairwise[[key2]] <- res2
            if (!is.null(output_dir)) {
              dir2 <- file.path(output_dir, key2); dir_create_safe(dir2)
              readr::write_csv(res2$activation, file.path(dir2, "activation.csv"))
              readr::write_csv(res2$inhibition, file.path(dir2, "inhibition.csv"))
            }
          }
        }
      }
    }
  }
  
  invisible(results)
}
