## R language
##
## GO annotation pipeline (block-based, LR / RTF / TFT)

#' Block-based GO annotation and integration for LR / RTF / TFT
#'
#' @description
#' This function reads L–R–TF–Target network CSV files
#' (autocrine: `*_L-R-TF-T-network-auto.csv`,
#'  pairwise: `*_L-R-TF-T_network.csv`) and performs:
#' \enumerate{
#'   \item Block extraction for LR / RTF / TFT.
#'   \item Building complex pairs (collapse multi-column L/R into comma-separated strings).
#'   \item Gene alias normalisation (ALIAS → SYMBOL).
#'   \item GO annotation (optionally using GO OBO to intersect ancestor sets).
#'   \item Integration of LR / RTF / TFT GO results into L–R–TF–Target–GO tables.
#' }
#'
#' @param input_dir Directory containing network CSVs
#'   (`*_L-R-TF-T-network-auto.csv` and/or `*_L-R-TF-T_network.csv`).
#' @param output_dir Output root directory. Subdirectories `LR/`, `RTF/`, `TFT/`,
#'   and `integrated/` will be created.
#' @param org_db Annotation package object, e.g. \code{org.Hs.eg.db}.
#' @param stages Character vector of stages to run. Subset of
#'   \code{c("extract_blocks","build_complex",
#'          "alias_normalize","go_annot","integrate")}.
#' @param which_files Which network types to process. Subset of
#'   \code{c("autocrine","pairwise")}.
#' @param use_obo Logical. If \code{TRUE}, download/use GO OBO and intersect
#'   ancestor GO sets (recommended for LR/RTF).
#' @param obo_url URL for GO OBO file.
#' @param verbose Logical. If \code{TRUE}, print progress messages.
#'
#' @return A list with vectors of output paths for each stage:
#'   \code{extracted_blocks}, \code{complex_files},
#'   \code{alias_replaced}, \code{go_outputs}, \code{integrated}.
#' @export
compute_go_blocks_and_integrate <- function(
    input_dir,
    output_dir,
    org_db = org.Hs.eg.db,
    stages = c("extract_blocks","build_complex",
               "alias_normalize","go_annot","integrate"),
    which_files = c("autocrine","pairwise"),
    use_obo = TRUE,
    obo_url = "http://purl.obolibrary.org/obo/go.obo",
    verbose = TRUE
){
  stopifnot(dir.exists(input_dir))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(AnnotationDbi)
    library(GO.db)
  })
  
  .log <- function(...) if (isTRUE(verbose)) cat(sprintf(...), "\n")
  
  ## --------------------------- dirs & patterns ---------------------------
  dir_lr        <- file.path(output_dir, "LR")
  dir_rtf       <- file.path(output_dir, "RTF")
  dir_tft       <- file.path(output_dir, "TFT")
  dir_integrated<- file.path(output_dir, "integrated")
  
  dir.create(dir_lr,        showWarnings = FALSE, recursive = TRUE)
  dir.create(dir_rtf,       showWarnings = FALSE, recursive = TRUE)
  dir.create(dir_tft,       showWarnings = FALSE, recursive = TRUE)
  dir.create(dir_integrated,showWarnings = FALSE, recursive = TRUE)
  
  pats <- c()
  if ("autocrine" %in% which_files) {
    pats <- c(pats, "_L-R-TF-T-network-auto\\.csv$")
  }
  if ("pairwise"  %in% which_files) {
    pats <- c(pats, "_L-R-TF-T_network\\.csv$")
  }
  if (length(pats) == 0) {
    stop("`which_files` must contain at least 'autocrine' or 'pairwise'.")
  }
  
  list_input <- function() {
    unlist(lapply(
      pats,
      function(pt) list.files(input_dir, pattern = pt, full.names = TRUE)
    ))
  }
  
  ## ------------------------ helper: prefix ------------------------
  get_prefix <- function(path) {
    tools::file_path_sans_ext(basename(path)) %>%
      sub("_L-R-TF-T-network-auto$", "", .) %>%
      sub("_L-R-TF-T_network$", "", .)
  }
  
  ## ------------------------ Stage 1: extract blocks ----------------------
  out_extract <- character(0)
  if ("extract_blocks" %in% stages) {
    .log("Stage 1: extracting LR / RTF / TFT blocks ...")
    in_files <- list_input()
    
    for (fp in in_files) {
      prefix <- get_prefix(fp)
      df <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
      
      l_cols <- grep("^L", names(df), value = TRUE)
      r_cols <- grep("^R", names(df), value = TRUE)
      tf_col <- "tf"
      target_col <- "target"
      
      ## LR block
      if (length(l_cols) && length(r_cols)) {
        lr_block <- df[, c(l_cols, r_cols), drop = FALSE] %>%
          dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !is.na(.) & . != "")) %>%
          dplyr::distinct() %>%
          dplyr::mutate(dplyr::across(dplyr::everything(),
                                      ~ ifelse(is.na(.), "", .)))
        fp_out <- file.path(dir_lr, paste0(prefix, "_LR_block.csv"))
        write.csv(lr_block, fp_out, row.names = FALSE)
        out_extract <- c(out_extract, fp_out)
      }
      
      ## RTF block
      if (length(r_cols) && tf_col %in% names(df)) {
        rtf_block <- df[, c(r_cols, tf_col), drop = FALSE] %>%
          dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !is.na(.) & . != "")) %>%
          dplyr::distinct() %>%
          dplyr::mutate(dplyr::across(dplyr::everything(),
                                      ~ ifelse(is.na(.), "", .)))
        fp_out <- file.path(dir_rtf, paste0(prefix, "_RTF_block.csv"))
        write.csv(rtf_block, fp_out, row.names = FALSE)
        out_extract <- c(out_extract, fp_out)
      }
      
      ## TFT block
      if (all(c(tf_col, target_col) %in% names(df))) {
        tft_block <- df[, c(tf_col, target_col), drop = FALSE] %>%
          dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !is.na(.) & . != "")) %>%
          dplyr::distinct() %>%
          dplyr::mutate(dplyr::across(dplyr::everything(),
                                      ~ ifelse(is.na(.), "", .)))
        fp_out <- file.path(dir_tft, paste0(prefix, "_TFT_block.csv"))
        write.csv(tft_block, fp_out, row.names = FALSE)
        out_extract <- c(out_extract, fp_out)
      }
    }
    .log("✓ Stage 1 done, %d block files written.", length(out_extract))
  }
  
  ## ---------------------- Stage 2: build complex pairs -------------------
  out_complex <- character(0)
  
  build_complex_pairs <- function(df, left_cols, right_cols,
                                  left_lab = "Complex_Ligands",
                                  right_lab = "Complex_Receptors"){
    tibble::tibble(
      !!left_lab := apply(df[, left_cols, drop = FALSE], 1, function(x){
        x <- na.omit(x)
        paste(unique(x[x != ""]), collapse = ",")
      }),
      !!right_lab := apply(df[, right_cols, drop = FALSE], 1, function(x){
        x <- na.omit(x)
        paste(unique(x[x != ""]), collapse = ",")
      })
    ) %>%
      dplyr::filter(.data[[left_lab]] != "", .data[[right_lab]] != "") %>%
      dplyr::distinct()
  }
  
  if ("build_complex" %in% stages) {
    .log("Stage 2: building complex (multi-subunit) pairs ...")
    
    ## LR complex
    lr_files <- list.files(dir_lr, pattern = "_LR_block\\.csv$", full.names = TRUE)
    for (fp in lr_files) {
      df <- read.csv(fp, check.names = FALSE)
      l_cols <- grep("^L", names(df), value = TRUE)
      r_cols <- grep("^R", names(df), value = TRUE)
      if (length(l_cols) && length(r_cols)) {
        out <- build_complex_pairs(df, l_cols, r_cols)
        fp_out <- sub("_LR_block\\.csv$", "_LR_block_complex.csv", fp)
        write.csv(out, fp_out, row.names = FALSE)
        out_complex <- c(out_complex, fp_out)
      }
    }
    
    ## RTF complex: treat R as "left", TF as "right"
    rtf_files <- list.files(dir_rtf, pattern = "_RTF_block\\.csv$", full.names = TRUE)
    for (fp in rtf_files) {
      df <- read.csv(fp, check.names = FALSE)
      r_cols <- grep("^R", names(df), value = TRUE)
      if (length(r_cols) && "tf" %in% names(df)) {
        out <- build_complex_pairs(df, r_cols, "tf")
        fp_out <- sub("_RTF_block\\.csv$", "_RTF_block_complex.csv", fp)
        write.csv(out, fp_out, row.names = FALSE)
        out_complex <- c(out_complex, fp_out)
      }
    }
    
    .log("✓ Stage 2 done, %d complex files written.", length(out_complex))
  }
  
  ## --------------------- Stage 3: alias normalisation --------------------
  out_alias <- character(0)
  if ("alias_normalize" %in% stages) {
    .log("Stage 3: alias normalisation (ALIAS → SYMBOL) ...")
    
    ## Build alias→symbol mapping
    symbol_alias_map <- AnnotationDbi::select(
      org_db,
      keys    = keys(org_db, keytype = "SYMBOL"),
      columns = c("SYMBOL","ALIAS"),
      keytype = "SYMBOL"
    ) %>%
      dplyr::filter(!is.na(ALIAS))
    
    alias_to_symbol <- setNames(symbol_alias_map$SYMBOL,
                                symbol_alias_map$ALIAS)
    
    replace_aliases <- function(glist){
      genes <- unlist(strsplit(glist, ","))
      rep <- vapply(
        genes,
        function(g) {
          if (!is.na(alias_to_symbol[g])) alias_to_symbol[[g]] else g
        },
        character(1)
      )
      paste(unique(rep), collapse = ",")
    }
    
    complex_files <- c(
      list.files(dir_lr,  pattern = "_LR_block_complex\\.csv$",  full.names = TRUE),
      list.files(dir_rtf, pattern = "_RTF_block_complex\\.csv$", full.names = TRUE)
    )
    
    for (fp in complex_files) {
      df <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
      if (!all(c("Complex_Ligands","Complex_Receptors") %in% names(df))) next
      df$Complex_Ligands   <- vapply(df$Complex_Ligands,
                                     replace_aliases, character(1))
      df$Complex_Receptors <- vapply(df$Complex_Receptors,
                                     replace_aliases, character(1))
      fp_out <- sub("\\.csv$", "_alias_replaced.csv", fp)
      write.csv(df, fp_out, row.names = FALSE)
      out_alias <- c(out_alias, fp_out)
    }
    
    .log("✓ Stage 3 done, %d alias-normalised files written.", length(out_alias))
  }
  
  ## --------------------- Stage 4: GO annotation --------------------------
  out_go <- character(0)
  if ("go_annot" %in% stages) {
    .log("Stage 4: GO annotation ...")
    
    ## optionally load GO OBO
    go_obo <- NULL
    if (isTRUE(use_obo)) {
      suppressPackageStartupMessages(library(ontologyIndex))
      go_obo <- tryCatch(
        ontologyIndex::get_OBO(obo_url),
        error = function(e) NULL
      )
      if (is.null(go_obo)) {
        .log("⚠ OBO download failed, falling back to direct GO intersection.")
      }
    }
    
    ## gene → GO mapping
    valid_symbols <- keys(org_db, keytype = "SYMBOL")
    select_go <- function(symbols){
      symbols <- unique(symbols[symbols %in% valid_symbols])
      if (!length(symbols)) return(tibble())
      AnnotationDbi::select(
        org_db,
        keys    = symbols,
        columns = c("GO","SYMBOL","ONTOLOGY"),
        keytype = "SYMBOL"
      ) %>%
        dplyr::filter(!is.na(GO)) %>%
        dplyr::distinct()
    }
    
    get_ancestors_or_self <- function(ids){
      if (is.null(go_obo)) return(ids)
      anc <- unique(unlist(lapply(
        ids,
        function(go) ontologyIndex::get_ancestors(go_obo, go)
      )))
      unique(anc[!is.na(anc)])
    }
    
    go_term_df <- function(go_ids){
      go_ids <- unique(go_ids[go_ids %in% keys(GO.db)])
      if (!length(go_ids)) return(tibble())
      tibble::tibble(
        GO   = go_ids,
        Name = vapply(go_ids,
                      function(x) Term(GO.db::GOTERM[[x]]),
                      character(1)),
        Definition = vapply(go_ids,
                            function(x) Definition(GO.db::GOTERM[[x]]),
                            character(1)),
        Ontology   = vapply(go_ids,
                            function(x) Ontology(GO.db::GOTERM[[x]]),
                            character(1))
      ) %>%
        dplyr::filter(!is.na(Name))
    }
    
    ## ---- 4.1 LR: intersect ancestor GO sets of L and R ----
    lr_src <- list.files(
      dir_lr,
      pattern = "_LR_block_complex(_alias_replaced)?\\.csv$",
      full.names = TRUE
    )
    
    for (fp in lr_src) {
      df <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
      if (!nrow(df)) next
      res_list <- vector("list", nrow(df))
      
      for (i in seq_len(nrow(df))) {
        Ls <- unlist(strsplit(df$Complex_Ligands[i],  ","))
        Rs <- unlist(strsplit(df$Complex_Receptors[i],","))
        
        go_L <- select_go(Ls)$GO %>% unique()
        go_R <- select_go(Rs)$GO %>% unique()
        if (!length(go_L) || !length(go_R)) next
        
        go_L_full <- get_ancestors_or_self(go_L)
        go_R_full <- get_ancestors_or_self(go_R)
        
        cross_go <- intersect(go_L_full, go_R_full)
        if (length(cross_go)) {
          info <- go_term_df(cross_go)
          if (nrow(info)) {
            info$Complex_Ligands   <- df$Complex_Ligands[i]
            info$Complex_Receptors <- df$Complex_Receptors[i]
            res_list[[i]] <- info
          }
        }
      }
      
      res <- dplyr::bind_rows(res_list) %>% dplyr::distinct()
      if (nrow(res)) {
        fp_out <- sub("\\.csv$", "_GO_annotation.csv", fp)
        write.csv(res, fp_out, row.names = FALSE)
        out_go <- c(out_go, fp_out)
      }
    }
    
    ## ---- 4.2 RTF: intersect ancestor GO sets of R and TF ----
    rtf_src <- list.files(
      dir_rtf,
      pattern = "_RTF_block_complex(_alias_replaced)?\\.csv$",
      full.names = TRUE
    )
    
    for (fp in rtf_src) {
      df <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
      if (!nrow(df)) next
      res_list <- vector("list", nrow(df))
      
      for (i in seq_len(nrow(df))) {
        Rs  <- unlist(strsplit(df$Complex_Ligands[i],  ","))  # left = R
        TFs <- unlist(strsplit(df$Complex_Receptors[i],","))  # right = TF
        
        go_R  <- select_go(Rs)$GO  %>% unique()
        go_TF <- select_go(TFs)$GO %>% unique()
        if (!length(go_R) || !length(go_TF)) next
        
        go_R_full  <- get_ancestors_or_self(go_R)
        go_TF_full <- get_ancestors_or_self(go_TF)
        
        cross_go <- intersect(go_R_full, go_TF_full)
        if (length(cross_go)) {
          info <- go_term_df(cross_go)
          if (nrow(info)) {
            info$Complex_Ligands   <- df$Complex_Ligands[i]
            info$Complex_Receptors <- df$Complex_Receptors[i]
            res_list[[i]] <- info
          }
        }
      }
      
      res <- dplyr::bind_rows(res_list) %>% dplyr::distinct()
      if (nrow(res)) {
        fp_out <- sub("\\.csv$", "_GO_annotation.csv", fp)
        write.csv(res, fp_out, row.names = FALSE)
        out_go <- c(out_go, fp_out)
      }
    }
    
    ## ---- 4.3 TFT: direct GO intersection of TF and Target ----
    tft_src <- list.files(
      dir_tft,
      pattern = "_TFT_block\\.csv$",
      full.names = TRUE
    )
    
    for (fp in tft_src) {
      df <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE) %>%
        tidyr::drop_na(tf, target) %>%
        dplyr::distinct()
      
      if (!nrow(df)) next
      
      go_map <- AnnotationDbi::select(
        org_db,
        keys    = unique(c(df$tf, df$target)),
        columns = c("GO","SYMBOL"),
        keytype = "SYMBOL"
      ) %>%
        dplyr::filter(!is.na(GO)) %>%
        dplyr::distinct()
      
      gene2go <- split(go_map$GO, go_map$SYMBOL)
      
      res_list <- vector("list", nrow(df))
      for (i in seq_len(nrow(df))) {
        tf_g <- df$tf[i]
        tg   <- df$target[i]
        go_tf <- gene2go[[tf_g]]
        go_tg <- gene2go[[tg]]
        
        if (!is.null(go_tf) && !is.null(go_tg)) {
          shared <- intersect(unique(go_tf), unique(go_tg))
          if (length(shared)) {
            info <- go_term_df(shared)
            if (nrow(info)) {
              info$TF     <- tf_g
              info$Target <- tg
              res_list[[i]] <- info
            }
          }
        }
      }
      
      res <- dplyr::bind_rows(res_list) %>% dplyr::distinct()
      if (nrow(res)) {
        fp_out <- sub("\\.csv$", "_TF_Target_GO_annotation.csv", fp)
        write.csv(res, fp_out, row.names = FALSE)
        out_go <- c(out_go, fp_out)
      }
    }
    
    .log("✓ Stage 4 done, %d GO-annotated files written.", length(out_go))
  }
  
  ## --------------------- Stage 5: integrate LR / RTF / TFT ---------------
  out_integrated <- character(0)
  if ("integrate" %in% stages) {
    .log("Stage 5: integrating GO annotations from LR / RTF / TFT ...")
    
    lr_files <- list.files(
      dir_lr,
      pattern = "_LR_block_complex(_alias_replaced)?_GO_annotation\\.csv$",
      full.names = TRUE
    )
    rtf_files <- list.files(
      dir_rtf,
      pattern = "_RTF_block_complex(_alias_replaced)?_GO_annotation\\.csv$",
      full.names = TRUE
    )
    tft_files <- list.files(
      dir_tft,
      pattern = "_TFT_block_TF_Target_GO_annotation\\.csv$",
      full.names = TRUE
    )
    
    get_prefix2 <- function(path) {
      nm <- basename(path)
      nm <- sub("_LR_block.*$",  "", nm)
      nm <- sub("_RTF_block.*$", "", nm)
      nm <- sub("_TFT_block.*$", "", nm)
      nm
    }
    
    lr_map  <- tibble::tibble(
      prefix = vapply(lr_files,  get_prefix2, character(1)),
      lr_file  = lr_files
    )
    rtf_map <- tibble::tibble(
      prefix = vapply(rtf_files, get_prefix2, character(1)),
      rtf_file = rtf_files
    )
    tft_map <- tibble::tibble(
      prefix = vapply(tft_files, get_prefix2, character(1)),
      tft_file = tft_files
    )
    
    pairs <- lr_map %>%
      dplyr::inner_join(rtf_map, by = "prefix") %>%
      dplyr::inner_join(tft_map, by = "prefix")
    
    if (!nrow(pairs)) {
      .log("No matching prefixes with LR + RTF + TFT GO annotations; skipping integration.")
    } else {
      ## helper: rename first column matching pattern to new name
      rename_first <- function(df, pattern, new) {
        if (new %in% names(df)) return(df)
        cand <- grep(pattern, names(df), value = TRUE)
        if (length(cand) == 0) return(df)
        names(df)[match(cand[1], names(df))] <- new
        df
      }
      
      for (i in seq_len(nrow(pairs))) {
        prefix <- pairs$prefix[i]
        
        lr  <- readr::read_csv(pairs$lr_file[i],  show_col_types = FALSE)
        rtf <- readr::read_csv(pairs$rtf_file[i], show_col_types = FALSE)
        tft <- readr::read_csv(pairs$tft_file[i], show_col_types = FALSE)
        
        if (!nrow(lr) || !nrow(rtf) || !nrow(tft)) next
        
        ## Standardise key columns in each table
        lr  <- rename_first(lr,  "^GO",              "GO")
        lr  <- rename_first(lr,  "Complex_Ligands",  "Complex_Ligands")
        lr  <- rename_first(lr,  "Complex_Receptors","Complex_Receptors")
        
        rtf <- rename_first(rtf, "^GO",              "GO")
        rtf <- rename_first(rtf, "Complex_Ligands",  "Complex_Ligands")
        rtf <- rename_first(rtf, "Complex_Receptors","Complex_Receptors")
        
        tft <- rename_first(tft, "^GO",              "GO")
        tft <- rename_first(tft, "^TF$",             "TF")
        tft <- rename_first(tft, "^Target$",         "Target")
        
        ## Rename RTF complexes to R / TF for clarity
        rtf_std <- rtf %>%
          dplyr::rename(
            R  = Complex_Ligands,
            TF = Complex_Receptors
          )
        
        ## Join: LR (L,R) + RTF (R,TF) on GO and receptor-complex,
        ## then join TFT (TF,Target) on GO only (keeps original behaviour).
        merged1 <- lr %>%
          dplyr::inner_join(
            rtf_std,
            by = c("GO", "Complex_Receptors" = "R")
          )
        if (!nrow(merged1)) next
        
        merged2 <- merged1 %>%
          dplyr::inner_join(tft, by = c("GO"))
        if (!nrow(merged2)) next
        
        ## After joins, column names like Name/Definition may be duplicated
        ## (Name.x, Name.y, ...). Pick the first matching ones.
        name_candidates <- grep("^Name", names(merged2), value = TRUE)
        def_candidates  <- grep("^Definition", names(merged2), value = TRUE)
        ont_candidates  <- grep("Ontology", names(merged2), value = TRUE)
        
        name_col <- if (length(name_candidates)) name_candidates[1] else NA_character_
        def_col  <- if (length(def_candidates))  def_candidates[1]  else NA_character_
        ont_col  <- if (length(ont_candidates))  ont_candidates[1]  else NA_character_
        
        ## Ensure key columns exist (may be .x/.y after joins)
        merged2 <- rename_first(merged2, "Complex_Ligands",   "Complex_Ligands")
        merged2 <- rename_first(merged2, "Complex_Receptors", "Complex_Receptors")
        merged2 <- rename_first(merged2, "^TF$",              "TF")
        merged2 <- rename_first(merged2, "^Target$",          "Target")
        merged2 <- rename_first(merged2, "^GO",               "GO")
        
        if (!all(c("GO","Complex_Ligands","Complex_Receptors","TF","Target") %in% names(merged2))) {
          .log("Skipping prefix %s: required columns missing after merging.", prefix)
          next
        }
        
        out_df <- merged2 %>%
          dplyr::transmute(
            GO,
            Name = if (!is.na(name_col)) .data[[name_col]] else NA_character_,
            Definition = if (!is.na(def_col)) .data[[def_col]] else NA_character_,
            Ontology   = if (!is.na(ont_col)) .data[[ont_col]] else NA_character_,
            L  = Complex_Ligands,
            R  = Complex_Receptors,
            TF = TF,
            Target = Target
          ) %>%
          dplyr::distinct()
        
        if (!nrow(out_df)) next
        
        fp_out <- file.path(dir_integrated,
                            paste0(prefix, "_GO_L_R_TF_Target.csv"))
        readr::write_csv(out_df, fp_out)
        out_integrated <- c(out_integrated, fp_out)
      }
    }
    
    .log("✓ Stage 5 done, %d integrated L–R–TF–Target GO files written.",
         length(out_integrated))
  }
  
  list(
    extracted_blocks = out_extract,
    complex_files    = out_complex,
    alias_replaced   = out_alias,
    go_outputs       = out_go,
    integrated       = out_integrated
  )
}
