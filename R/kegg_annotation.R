# ==========================================
# KEGG annotation for LR / RTF / TFT blocks
# ==========================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(KEGGREST)
  library(AnnotationDbi)
  library(readr)
})

#' KEGG annotation for LR / RTF / TFT GO blocks
#'
#' Given a root directory containing subfolders "LR", "RTF", "TFT" with
#' GO-annotated block files, this function:
#'   1) annotates each file with KEGG pathway IDs (intersection across genes)
#'   2) adds human-readable pathway names for each KEGG ID
#'
#' Output files:
#'   *_KEGG.csv            – adds column KEGG_Pathways
#'   *_KEGG_with_names.csv – adds column KEGG_Pathways_Names
#'   If `output_root` is NULL, they are written into the same LR/RTF/TFT
#'   folders under `root_dir`. Otherwise they are written into
#'   file.path(output_root, "LR"/"RTF"/"TFT").
#'
#' @param root_dir   Root directory that contains subfolders "LR", "RTF", "TFT".
#' @param species    Species to use. One of:
#'                   "human", "mouse", "cow", "rhesus", "rat".
#' @param org_db     Optional OrgDb object (e.g. org.Hs.eg.db). If NULL,
#'                   it will be inferred from `species`.
#' @param kegg_org   Optional KEGG organism code ("hsa", "mmu", "bta", "mcc", "rno").
#'                   If NULL, it will be inferred from `species`.
#' @param output_root Optional directory for result files. If NULL, write to
#'                   the same LR/RTF/TFT directories under `root_dir`.
#' @param blocks     Which blocks to process: subset of c("LR","RTF","TFT").
#' @param lr_pattern  Regex pattern for LR GO files.
#' @param rtf_pattern Regex pattern for RTF GO files.
#' @param tft_pattern Regex pattern for TFT GO files.
#' @param retries    Number of retries for each KEGG REST call.
#' @param sleep_sec  Delay (seconds) between retries.
#' @param verbose    Whether to print progress.
#'
#' @return A list with components LR/RTF/TFT, each containing:
#'         $kegg_files and $with_name_files.
#' @export
compute_kegg_annotations_blocks <- function(
    root_dir,
    species    = c("human","mouse","cow","rhesus","rat"),
    org_db     = NULL,
    kegg_org   = NULL,
    output_root = NULL,
    blocks     = c("LR","RTF","TFT"),
    lr_pattern  = "_LR_block_complex_alias_replaced_GO_annotation\\.csv$",
    rtf_pattern = "_RTF_block_complex_alias_replaced_GO_annotation\\.csv$",
    tft_pattern = "_TFT_block_TF_Target_GO_annotation\\.csv$",
    retries   = 3,
    sleep_sec = 0.3,
    verbose   = TRUE
){
  root_dir <- normalizePath(root_dir, mustWork = TRUE)
  
  # ---------- species -> OrgDb & KEGG code ----------
  species <- match.arg(species)
  
  if (is.null(org_db)) {
    org_pkg <- switch(
      species,
      human  = "org.Hs.eg.db",
      mouse  = "org.Mm.eg.db",
      cow    = "org.Bt.eg.db",
      rhesus = "org.Mmu.eg.db",
      rat    = "org.Rn.eg.db"
    )
    if (!requireNamespace(org_pkg, quietly = TRUE)) {
      stop("Package ", org_pkg, " must be installed for species '",
           species, "'. Please install it or pass `org_db` explicitly.")
    }
    org_db <- get(org_pkg, envir = asNamespace(org_pkg))
  }
  
  if (!inherits(org_db, "OrgDb"))
    stop("`org_db` must be an OrgDb object.")
  
  if (is.null(kegg_org)) {
    kegg_org <- switch(
      species,
      human  = "hsa",
      mouse  = "mmu",
      cow    = "bta",
      rhesus = "mcc",
      rat    = "rno"
    )
  }
  
  blocks <- intersect(blocks, c("LR","RTF","TFT"))
  if (!length(blocks))
    stop("`blocks` must contain at least one of 'LR','RTF','TFT'.")
  
  # 输出根目录：如果没给，就用 root_dir
  if (is.null(output_root)) {
    output_root <- root_dir
  }
  output_root <- normalizePath(output_root, mustWork = FALSE)
  dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
  
  .log <- function(...) if (isTRUE(verbose)) cat(sprintf(...), "\n")
  
  # ---------- SYMBOL -> ENTREZ ----------
  symbol_to_entrez <- function(symbols) {
    symbols <- unique(na.omit(symbols))
    if (!length(symbols)) return(character(0))
    suppressMessages({
      ids <- AnnotationDbi::mapIds(
        org_db,
        keys      = symbols,
        column    = "ENTREZID",
        keytype   = "SYMBOL",
        multiVals = "first"
      )
    })
    ids <- ids[!is.na(ids)]
    as.character(ids)
  }
  
  # ---------- KEGG intersection with cache ----------
  kegg_cache <- new.env(parent = emptyenv())
  
  get_kegg_paths_intersect <- function(gene_symbols) {
    entrez_ids <- symbol_to_entrez(gene_symbols)
    if (!length(entrez_ids)) return(character(0))
    
    kegg_lists <- lapply(entrez_ids, function(eid) {
      if (exists(eid, envir = kegg_cache, inherits = FALSE)) {
        return(get(eid, envir = kegg_cache, inherits = FALSE))
      } else {
        res <- NULL
        for (i in seq_len(retries)) {
          Sys.sleep(sleep_sec)
          res <- tryCatch({
            r  <- KEGGREST::keggGet(paste0(kegg_org, ":", eid))[[1]]
            pw <- if (!is.null(r$PATHWAY)) names(r$PATHWAY) else character(0)
            assign(eid, pw, envir = kegg_cache)
            pw
          }, error = function(e) NULL)
          if (!is.null(res)) break
        }
        if (is.null(res)) {
          assign(eid, character(0), envir = kegg_cache)
          return(character(0))
        }
        res
      }
    })
    
    Reduce(intersect, kegg_lists)
  }
  
  # ---------- KEGG ID -> name ----------
  .log("Fetching KEGG pathway list for organism '%s' ...", kegg_org)
  kegg_pathways <- KEGGREST::keggList("pathway", kegg_org)
  kegg_id_to_name <- as.character(kegg_pathways)
  names(kegg_id_to_name) <- sub("^path:", "", names(kegg_pathways))
  
  convert_kegg_ids_to_names <- function(kegg_str) {
    if (is.na(kegg_str) || kegg_str == "") return("")
    ids <- strsplit(kegg_str, ";")[[1]]
    ids <- stringr::str_trim(ids)
    if (!length(ids)) return("")
    nm <- vapply(ids, function(id) {
      if (id %in% names(kegg_id_to_name)) {
        paste0(id, ": ", kegg_id_to_name[[id]])
      } else {
        paste0(id, ": Unknown")
      }
    }, character(1))
    paste(nm, collapse = "; ")
  }
  
  # ---------- LR / RTF annotator ----------
  annotate_lr_rtf_dir <- function(in_dir, out_dir, pattern, gene_cols) {
    if (!dir.exists(in_dir)) {
      .log("Directory not found, skip: %s", in_dir)
      return(character(0))
    }
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    files <- list.files(in_dir, pattern = pattern, full.names = TRUE)
    if (!length(files)) {
      .log("No GO files found in %s (pattern: %s)", in_dir, pattern)
      return(character(0))
    }
    
    out_files <- character(0)
    for (f in files) {
      .log("  [KEGG] %s", basename(f))
      df <- tryCatch(
        readr::read_csv(f, show_col_types = FALSE),
        error = function(e) { .log("    ! read failed: %s", e$message); NULL }
      )
      if (is.null(df)) next
      if (!all(gene_cols %in% colnames(df))) {
        .log("    ! missing columns %s, skip",
             paste(gene_cols, collapse = ","))
        next
      }
      
      unique_combo <- df %>%
        dplyr::select(dplyr::all_of(gene_cols)) %>%
        dplyr::distinct()
      
      annotated <- unique_combo %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          Genes_used = list({
            g1 <- unlist(strsplit(as.character(.data[[gene_cols[1]]]),
                                  ",\\s*"))
            g2 <- unlist(strsplit(as.character(.data[[gene_cols[2]]]),
                                  ",\\s*"))
            unique(c(g1, g2))
          }),
          KEGG_Pathways = list(get_kegg_paths_intersect(Genes_used))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          KEGG_Pathways = vapply(
            KEGG_Pathways,
            function(x) paste(x, collapse = ";"),
            character(1)
          )
        ) %>%
        dplyr::select(-Genes_used)
      
      df_final <- dplyr::left_join(
        df,
        annotated,
        by = setNames(gene_cols, gene_cols)
      )
      
      out_file <- file.path(
        out_dir,
        paste0(tools::file_path_sans_ext(basename(f)), "_KEGG.csv")
      )
      readr::write_csv(df_final, out_file)
      out_files <- c(out_files, out_file)
      .log("    -> wrote %s", basename(out_file))
    }
    out_files
  }
  
  # ---------- TFT annotator ----------
  annotate_tft_dir <- function(in_dir, out_dir, pattern) {
    if (!dir.exists(in_dir)) {
      .log("Directory not found, skip: %s", in_dir)
      return(character(0))
    }
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    files <- list.files(in_dir, pattern = pattern, full.names = TRUE)
    if (!length(files)) {
      .log("No GO files found in %s (pattern: %s)", in_dir, pattern)
      return(character(0))
    }
    
    out_files <- character(0)
    for (f in files) {
      .log("  [KEGG] %s", basename(f))
      df <- tryCatch(
        readr::read_csv(f, show_col_types = FALSE),
        error = function(e) { .log("    ! read failed: %s", e$message); NULL }
      )
      if (is.null(df)) next
      
      if (!all(c("TF","Target") %in% colnames(df))) {
        .log("    ! missing TF/Target, skip")
        next
      }
      
      unique_combo <- df %>%
        dplyr::select(TF, Target) %>%
        dplyr::distinct()
      
      annotated <- unique_combo %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          Genes_used    = list(na.omit(c(.data$TF, .data$Target))),
          KEGG_Pathways = list(get_kegg_paths_intersect(Genes_used))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          KEGG_Pathways = vapply(
            KEGG_Pathways,
            function(x) paste(x, collapse = ";"),
            character(1)
          )
        ) %>%
        dplyr::select(-Genes_used)
      
      df_final <- df %>%
        dplyr::left_join(annotated, by = c("TF","Target"))
      
      out_file <- file.path(
        out_dir,
        paste0(tools::file_path_sans_ext(basename(f)), "_KEGG.csv")
      )
      readr::write_csv(df_final, out_file)
      out_files <- c(out_files, out_file)
      .log("    -> wrote %s", basename(out_file))
    }
    out_files
  }
  
  # ---------- add KEGG names ----------
  add_names_dir <- function(dir_ke, pattern = "_KEGG\\.csv$") {
    if (!dir.exists(dir_ke)) return(character(0))
    files <- list.files(dir_ke, pattern = pattern, full.names = TRUE)
    if (!length(files)) return(character(0))
    
    out_files <- character(0)
    for (f in files) {
      df <- tryCatch(
        readr::read_csv(f, show_col_types = FALSE),
        error = function(e) { .log("    ! read failed: %s", e$message); NULL }
      )
      if (is.null(df) || !"KEGG_Pathways" %in% colnames(df)) next
      
      df2 <- df %>%
        dplyr::mutate(
          KEGG_Pathways_Names = vapply(
            KEGG_Pathways,
            convert_kegg_ids_to_names,
            character(1)
          )
        )
      
      out_file <- file.path(
        dir_ke,
        paste0(tools::file_path_sans_ext(basename(f)), "_with_names.csv")
      )
      readr::write_csv(df2, out_file)
      out_files <- c(out_files, out_file)
      .log("    -> wrote %s", basename(out_file))
    }
    out_files
  }
  
  # ---------- run for each block ----------
  res <- list(
    LR  = list(kegg_files = character(0), with_name_files = character(0)),
    RTF = list(kegg_files = character(0), with_name_files = character(0)),
    TFT = list(kegg_files = character(0), with_name_files = character(0))
  )
  
  if ("LR" %in% blocks) {
    lr_in_dir  <- file.path(root_dir,  "LR")
    lr_out_dir <- file.path(output_root, "LR")
    .log("== LR block: in %s -> out %s ==", lr_in_dir, lr_out_dir)
    res$LR$kegg_files      <- annotate_lr_rtf_dir(
      lr_in_dir, lr_out_dir, lr_pattern,
      gene_cols = c("Complex_Ligands","Complex_Receptors")
    )
    res$LR$with_name_files <- add_names_dir(lr_out_dir)
  }
  
  if ("RTF" %in% blocks) {
    rtf_in_dir  <- file.path(root_dir,  "RTF")
    rtf_out_dir <- file.path(output_root, "RTF")
    .log("== RTF block: in %s -> out %s ==", rtf_in_dir, rtf_out_dir)
    res$RTF$kegg_files      <- annotate_lr_rtf_dir(
      rtf_in_dir, rtf_out_dir, rtf_pattern,
      gene_cols = c("Complex_Ligands","Complex_Receptors")
    )
    res$RTF$with_name_files <- add_names_dir(rtf_out_dir)
  }
  
  if ("TFT" %in% blocks) {
    tft_in_dir  <- file.path(root_dir,  "TFT")
    tft_out_dir <- file.path(output_root, "TFT")
    .log("== TFT block: in %s -> out %s ==", tft_in_dir, tft_out_dir)
    res$TFT$kegg_files      <- annotate_tft_dir(tft_in_dir, tft_out_dir, tft_pattern)
    res$TFT$with_name_files <- add_names_dir(tft_out_dir)
  }
  
  invisible(res)
}
