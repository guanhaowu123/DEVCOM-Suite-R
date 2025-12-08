#' Integrate GO-annotated LR / RTF / TFT blocks into L–R–TF–Target tables
#'
#' @description
#' This function integrates GO-annotated results from three types of blocks:
#' \itemize{
#'   \item LR:  files named like \code{*_LR_block_complex_alias_replaced_GO_annotation.csv}
#'   \item RTF: files named like \code{*_RTF_block_complex_alias_replaced_GO_annotation.csv}
#'   \item TFT: files named like \code{*_TFT_block_TF_Target_GO_annotation.csv}
#' }
#' For each common prefix, the function:
#' \enumerate{
#'   \item Reads LR / RTF / TFT GO annotation tables.
#'   \item Renames RTF complex columns so that receptors are stored in column \code{R}
#'         and TFs in column \code{TF}.
#'   \item Joins LR and RTF on \code{GO} and receptor complexes.
#'   \item Joins the result with TFT on \code{GO} and \code{TF}.
#'   \item Outputs an integrated table with columns:
#'         \code{GO, Name, Definition, Ontology, L, R, TF, Target}.
#' }
#'
#' @param lr_dir   Directory containing LR GO-annotated files
#'   (pattern: \code{"_LR_block_complex_alias_replaced_GO_annotation.csv$"}).
#' @param rtf_dir  Directory containing RTF GO-annotated files
#'   (pattern: \code{"_RTF_block_complex_alias_replaced_GO_annotation.csv$"}).
#' @param tft_dir  Directory containing TFT GO-annotated files
#'   (pattern: \code{"_TFT_block_TF_Target_GO_annotation.csv$"}).
#' @param out_dir  Output directory for integrated
#'   \code{*_GO_L_R_TF_Target.csv} files.
#' @param lr_pattern  Optional regex pattern for LR files
#'   (default: \code{"_LR_block_complex_alias_replaced_GO_annotation.csv$"}).
#' @param rtf_pattern Optional regex pattern for RTF files
#'   (default: \code{"_RTF_block_complex_alias_replaced_GO_annotation.csv$"}).
#' @param tft_pattern Optional regex pattern for TFT files
#'   (default: \code{"_TFT_block_TF_Target_GO_annotation.csv$"}).
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#'
#' @return
#' A character vector of output file paths (one per integrated prefix),
#' returned invisibly.
#'
#' @importFrom tibble tibble
#' @importFrom dplyr inner_join rename transmute distinct
#' @importFrom readr read_csv write_csv
#' @importFrom stringr str_replace str_detect
#' @export
integrate_go_lr_rtf_tft <- function(
    lr_dir,
    rtf_dir,
    tft_dir,
    out_dir,
    lr_pattern  = "_LR_block_complex_alias_replaced_GO_annotation.csv$",
    rtf_pattern = "_RTF_block_complex_alias_replaced_GO_annotation.csv$",
    tft_pattern = "_TFT_block_TF_Target_GO_annotation.csv$",
    verbose = TRUE
){
  # ---- basic checks ----
  stopifnot(dir.exists(lr_dir), dir.exists(rtf_dir), dir.exists(tft_dir))
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  .log <- function(...) if (isTRUE(verbose)) cat(sprintf(...), "\n")
  
  # ---- list files ----
  lr_files  <- list.files(lr_dir,  pattern = lr_pattern,  full.names = TRUE)
  rtf_files <- list.files(rtf_dir, pattern = rtf_pattern, full.names = TRUE)
  tft_files <- list.files(tft_dir, pattern = tft_pattern, full.names = TRUE)
  
  if (!length(lr_files) || !length(rtf_files) || !length(tft_files)) {
    warning("No LR / RTF / TFT files found with the given patterns.")
  }
  
  # ---- extract prefix from filename ----
  get_prefix <- function(path) {
    stringr::str_replace(basename(path), "_(LR_block|RTF_block|TFT_block).*", "")
  }
  
  lr_map  <- tibble::tibble(prefix = vapply(lr_files,  get_prefix, character(1)),
                            lr_file  = lr_files)
  rtf_map <- tibble::tibble(prefix = vapply(rtf_files, get_prefix, character(1)),
                            rtf_file = rtf_files)
  tft_map <- tibble::tibble(prefix = vapply(tft_files, get_prefix, character(1)),
                            tft_file = tft_files)
  
  file_pairs <- lr_map %>%
    dplyr::inner_join(rtf_map, by = "prefix") %>%
    dplyr::inner_join(tft_map, by = "prefix")
  
  if (!nrow(file_pairs)) {
    warning("No common prefixes found across LR / RTF / TFT GO files.")
    return(invisible(character(0)))
  }
  
  # helper: ensure the first column matching `pattern` is renamed to `new`
  rename_first <- function(df, pattern, new) {
    if (new %in% names(df)) return(df)
    cand <- grep(pattern, names(df), value = TRUE)
    if (length(cand) == 0) return(df)
    names(df)[match(cand[1], names(df))] <- new
    df
  }
  
  out_files <- character(0)
  
  # ---- main loop over prefixes ----
  for (i in seq_len(nrow(file_pairs))) {
    prefix <- file_pairs$prefix[i]
    .log("Integrating GO blocks for prefix: %s", prefix)
    
    # 1) read three tables
    lr  <- readr::read_csv(file_pairs$lr_file[i],  show_col_types = FALSE)
    rtf <- readr::read_csv(file_pairs$rtf_file[i], show_col_types = FALSE)
    tft <- readr::read_csv(file_pairs$tft_file[i], show_col_types = FALSE)
    
    if (!nrow(lr) || !nrow(rtf) || !nrow(tft)) {
      .log("  Skipped (one of LR/RTF/TFT is empty): %s", prefix)
      next
    }
    
    # 2) standardise key columns
    lr  <- rename_first(lr,  "^GO",             "GO")
    lr  <- rename_first(lr,  "Complex_Ligands", "Complex_Ligands")
    lr  <- rename_first(lr,  "Complex_Receptors","Complex_Receptors")
    
    rtf <- rename_first(rtf, "^GO",             "GO")
    rtf <- rename_first(rtf, "Complex_Ligands", "Complex_Ligands")
    rtf <- rename_first(rtf, "Complex_Receptors","Complex_Receptors")
    
    tft <- rename_first(tft, "^GO",             "GO")
    tft <- rename_first(tft, "^TF$",            "TF")
    tft <- rename_first(tft, "^Target$",        "Target")
    
    # 3) rename RTF complex columns: R (receptor complex), TF (TF complex)
    rtf <- rtf %>%
      dplyr::rename(
        R  = Complex_Ligands,
        TF = Complex_Receptors
      )
    
    # 4) join LR and RTF by GO and receptor complex
    merged1 <- lr %>%
      dplyr::inner_join(
        rtf,
        by = c("GO", "Complex_Receptors" = "R"),
        suffix = c("", "")
      )
    
    if (!nrow(merged1)) {
      .log("  No LR–RTF overlap for prefix: %s", prefix)
      next
    }
    
    # 5) join with TFT by GO and TF
    merged2 <- merged1 %>%
      dplyr::inner_join(
        tft,
        by = c("GO", "TF" = "TF"),
        suffix = c("", "")
      )
    
    if (!nrow(merged2)) {
      .log("  No LR–RTF–TFT overlap for prefix: %s", prefix)
      next
    }
    
    # 6) pick Name / Definition / Ontology columns safely
    all_cols <- colnames(merged2)
    
    name_col <- all_cols[stringr::str_detect(all_cols, "^Name")][1]
    def_col  <- all_cols[stringr::str_detect(all_cols, "^Definition")][1]
    ont_col  <- all_cols[stringr::str_detect(all_cols, "Ontology")][1]
    
    if (is.na(name_col) || is.na(def_col) || is.na(ont_col)) {
      .log("  Missing Name / Definition / Ontology columns, skip prefix: %s", prefix)
      next
    }
    
    # 7) ensure Complex_Ligands / Complex_Receptors survived joins
    merged2 <- rename_first(merged2, "Complex_Ligands",   "Complex_Ligands")
    merged2 <- rename_first(merged2, "Complex_Receptors", "Complex_Receptors")
    merged2 <- rename_first(merged2, "^TF$",              "TF")
    merged2 <- rename_first(merged2, "^Target$",          "Target")
    merged2 <- rename_first(merged2, "^GO",               "GO")
    
    required <- c("GO","Complex_Ligands","Complex_Receptors","TF","Target")
    if (!all(required %in% names(merged2))) {
      .log("  Required columns missing after join, skip prefix: %s", prefix)
      next
    }
    
    # 8) construct final integrated table
    final_df <- merged2 %>%
      dplyr::transmute(
        GO,
        Name       = .data[[name_col]],
        Definition = .data[[def_col]],
        Ontology   = .data[[ont_col]],
        L          = Complex_Ligands,
        R          = Complex_Receptors,
        TF         = TF,
        Target     = Target
      ) %>%
      dplyr::distinct()
    
    if (!nrow(final_df)) {
      .log("  Empty integrated table, skip prefix: %s", prefix)
      next
    }
    
    out_file <- file.path(out_dir, paste0(prefix, "_GO_L_R_TF_Target.csv"))
    readr::write_csv(final_df, out_file)
    out_files <- c(out_files, out_file)
    .log("  Saved: %s", out_file)
  }
  
  invisible(out_files)
}
