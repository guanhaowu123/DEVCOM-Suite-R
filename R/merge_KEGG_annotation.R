# ==========================================
# Integrate KEGG results from LR / RTF / TFT
# ==========================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

#' Integrate KEGG-annotated LR / RTF / TFT blocks
#'
#' This function takes three directories containing KEGG-annotated
#' GO block files for LR, RTF and TFT, and merges them into unified
#' L–R–TF–Target–KEGG tables per sample.
#'
#' Expected file name patterns:
#'   <key>_LR_block_complex_alias_replaced_GO_annotation_KEGG_with_names.csv
#'   <key>_RTF_block_complex_alias_replaced_GO_annotation_KEGG_with_names.csv
#'   <key>_TFT_block_TF_Target_GO_annotation_KEGG_with_names.csv
#'
#' For each shared <key>, the function:
#'   1) Reads LR, RTF and TFT tables
#'   2) Standardises RTF columns so that:
#'        - receptor column is "R1" or "Complex_Ligands"
#'        - TF column is "TF" or "Complex_Receptors"
#'   3) Filters out rows with empty genes / KEGG in each table
#'   4) Joins LR and RTF on receptor (LR$Complex_Receptors == RTF$Receptor)
#'   5) Joins TFT on TF
#'   6) Outputs <key>_merged_LR_RTF_TFT.csv with columns:
#'        Complex_Ligands, Complex_Receptors, TF, Target, KEGG_Pathways_Names
#'
#' KEGG_Pathways_Names is taken from the TFT block (TF–Target level).
#'
#' @param lr_dir  Directory containing LR KEGG-with-names files.
#' @param rtf_dir Directory containing RTF KEGG-with-names files.
#' @param tft_dir Directory containing TFT KEGG-with-names files.
#' @param out_dir Directory to write integrated results. Will be created if needed.
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A named character vector of output file paths (names are sample keys).
#' @export
integrate_kegg_lr_rtf_tft <- function(
    lr_dir,
    rtf_dir,
    tft_dir,
    out_dir,
    verbose = TRUE
) {
  # basic checks
  stopifnot(dir.exists(lr_dir), dir.exists(rtf_dir), dir.exists(tft_dir))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  .log <- function(...) if (isTRUE(verbose)) cat(sprintf(...), "\n")
  
  # file name suffixes (literal)
  lr_suffix  <- "_LR_block_complex_alias_replaced_GO_annotation_KEGG_with_names.csv"
  rtf_suffix <- "_RTF_block_complex_alias_replaced_GO_annotation_KEGG_with_names.csv"
  tft_suffix <- "_TFT_block_TF_Target_GO_annotation_KEGG_with_names.csv"
  
  # list file names (not full paths)
  lr_names  <- list.files(lr_dir,  pattern = paste0(gsub("\\.", "\\\\.", lr_suffix), "$"))
  rtf_names <- list.files(rtf_dir, pattern = paste0(gsub("\\.", "\\\\.", rtf_suffix), "$"))
  tft_names <- list.files(tft_dir, pattern = paste0(gsub("\\.", "\\\\.", tft_suffix), "$"))
  
  # derive sample keys by stripping suffix
  lr_keys  <- sub(lr_suffix,  "", lr_names,  fixed = TRUE)
  rtf_keys <- sub(rtf_suffix, "", rtf_names, fixed = TRUE)
  tft_keys <- sub(tft_suffix, "", tft_names, fixed = TRUE)
  
  # intersection of keys
  common_keys <- Reduce(intersect, list(lr_keys, rtf_keys, tft_keys))
  .log("Found %d common keys across LR/RTF/TFT.", length(common_keys))
  
  if (!length(common_keys)) {
    warning("No common keys across LR/RTF/TFT. Nothing to integrate.")
    return(invisible(character(0)))
  }
  
  out_files <- character(0)
  
  for (key in common_keys) {
    .log("---- Integrating KEGG for key: %s ----", key)
    
    lr_file  <- file.path(lr_dir,  paste0(key, lr_suffix))
    rtf_file <- file.path(rtf_dir, paste0(key, rtf_suffix))
    tft_file <- file.path(tft_dir, paste0(key, tft_suffix))
    
    # ---------- LR ----------
    df_lr <- readr::read_csv(lr_file, show_col_types = FALSE) %>%
      dplyr::select(Complex_Ligands, Complex_Receptors, KEGG_Pathways_Names) %>%
      dplyr::filter(
        !is.na(Complex_Ligands)    & Complex_Ligands    != "",
        !is.na(Complex_Receptors)  & Complex_Receptors  != "",
        !is.na(KEGG_Pathways_Names) & KEGG_Pathways_Names != ""
      ) %>%
      dplyr::mutate(
        Complex_Ligands   = trimws(as.character(Complex_Ligands)),
        Complex_Receptors = trimws(as.character(Complex_Receptors))
      ) %>%
      dplyr::distinct() %>%
      # we only keep L/R for merging; KEGG from TFT will be used later
      dplyr::select(Complex_Ligands, Complex_Receptors)
    
    # ---------- RTF ----------
    df_rtf_raw <- readr::read_csv(rtf_file, show_col_types = FALSE)
    
    # detect receptor + TF columns for RTF
    receptor_col <- NULL
    tf_col       <- NULL
    
    if ("R1" %in% names(df_rtf_raw)) {
      receptor_col <- "R1"
    } else if ("Complex_Ligands" %in% names(df_rtf_raw)) {
      receptor_col <- "Complex_Ligands"
    }
    
    if ("TF" %in% names(df_rtf_raw)) {
      tf_col <- "TF"
    } else if ("Complex_Receptors" %in% names(df_rtf_raw)) {
      tf_col <- "Complex_Receptors"
    }
    
    if (is.null(receptor_col) || is.null(tf_col)) {
      warning("Key ", key, ": cannot find receptor/TF columns in RTF file; skip.")
      next
    }
    
    df_rtf <- df_rtf_raw %>%
      dplyr::select(
        Receptor = dplyr::all_of(receptor_col),
        TF       = dplyr::all_of(tf_col),
        KEGG_Pathways_Names
      ) %>%
      dplyr::filter(
        !is.na(Receptor) & Receptor != "",
        !is.na(TF)       & TF != "",
        !is.na(KEGG_Pathways_Names) & KEGG_Pathways_Names != ""
      ) %>%
      dplyr::mutate(
        Receptor = trimws(as.character(Receptor)),
        TF       = trimws(as.character(TF))
      ) %>%
      dplyr::distinct() %>%
      # we only keep receptor + TF for merging
      dplyr::select(Receptor, TF)
    
    # ---------- TFT ----------
    df_tft <- readr::read_csv(tft_file, show_col_types = FALSE) %>%
      dplyr::select(TF, Target, KEGG_Pathways_Names) %>%
      dplyr::filter(
        !is.na(TF)       & TF != "",
        !is.na(Target)   & Target != "",
        !is.na(KEGG_Pathways_Names) & KEGG_Pathways_Names != ""
      ) %>%
      dplyr::mutate(
        TF     = trimws(as.character(TF)),
        Target = trimws(as.character(Target))
      ) %>%
      dplyr::distinct()
    
    # ---------- Merge LR + RTF ----------
    df_merge1 <- df_rtf %>%
      dplyr::inner_join(
        df_lr,
        by = c("Receptor" = "Complex_Receptors")
      ) %>%
      dplyr::transmute(
        Complex_Ligands   = Complex_Ligands,
        Complex_Receptors = Receptor,
        TF
      ) %>%
      dplyr::distinct()
    
    .log("  Matched LR–RTF rows: %d", nrow(df_merge1))
    
    if (!nrow(df_merge1)) {
      .log("  No LR–RTF matches for key %s; skip.", key)
      next
    }
    
    # ---------- Merge with TFT (TF + Target + KEGG) ----------
    df_merge2 <- df_tft %>%
      dplyr::inner_join(df_merge1, by = "TF") %>%
      dplyr::transmute(
        Complex_Ligands,
        Complex_Receptors,
        TF,
        Target,
        KEGG_Pathways_Names
      ) %>%
      dplyr::distinct()
    
    .log("  After merging TFT: %d rows.", nrow(df_merge2))
    
    # ---------- Write output ----------
    out_file <- file.path(out_dir, paste0(key, "_merged_LR_RTF_TFT.csv"))
    readr::write_csv(df_merge2, out_file)
    .log("  -> Saved: %s\n", out_file)
    
    out_files[key] <- out_file
  }
  
  invisible(out_files)
}
