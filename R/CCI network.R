## R language
##
## CCI network: L–R–TF–Target (DevcomDB-based)

#' Construct CCI networks (L–R–TF–Target) based on DevcomDB
#'
#' @description
#' This function performs the following steps in one shot:
#' \enumerate{
#'   \item Load DevcomDB interaction data for the specified species
#'         (ligand–receptor, receptor–TF, TF–target).
#'   \item For each cell type, read its hub-gene list
#'         (\code{<cell>_filtered_genes.csv} containing a \code{Gene} column)
#'         and filter L–R–TF interactions that are compatible with
#'         sender/receiver gene sets (autocrine and/or pairwise).
#'   \item Merge TF–target relations to obtain full L–R–TF–Target networks.
#'   \item Optionally export per-network tables, L/R unique pairs and a
#'         count matrix summarising how many LR pairs exist between each
#'         sender/receiver pair.
#' }
#'
#' Ligand and receptor columns in the DevcomDB ligand–receptor table are
#' detected automatically as all columns whose names start with \code{"L"}
#' (for ligands) or \code{"R"} (for receptors), case-insensitive.
#'
#' @param species Species key. One of
#'   \code{"human"}, \code{"mouse"}, \code{"bos"},
#'   \code{"sheep"}, \code{"mmulatta"}, \code{"mfascicularis"}.
#' @param cell_types Character vector of cell-type names. Each must match the
#'   prefix of the hub-gene files (see \code{filtered_dir}).
#' @param filtered_dir Directory containing hub-gene files. For each cell type
#'   \code{c}, there must be a file
#'   \code{file.path(filtered_dir, paste0(c, "_filtered_genes.csv"))}
#'   with at least a \code{Gene} column.
#' @param output_dir Output directory. If \code{NULL}, nothing is written
#'   to disk and only in-memory objects are returned.
#' @param require_all_ligands Logical. If \code{TRUE} (default), all
#'   non-empty ligand subunits in a LR complex must be present in the sender
#'   gene set. If \code{FALSE}, at least one ligand subunit in the complex
#'   must be present.
#' @param require_all_receptors Logical. Same as \code{require_all_ligands}
#'   but for receptor subunits and the receiver gene set.
#' @param build_autocrine Logical. Whether to build autocrine networks
#'   (sender == receiver).
#' @param build_pairwise Logical. Whether to build pairwise sender→receiver
#'   networks for all ordered pairs of distinct cell types.
#' @param write_lr_unique Logical. If \code{TRUE}, export ligand–receptor
#'   unique pairs for every network.
#' @param write_counts Logical. If \code{TRUE}, export a LR count matrix
#'   (file name: \code{"cell_communication_summary-LR.csv"}).
#' @param verbose Logical. If \code{TRUE}, print progress messages.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{kb$lr_tf}}{Assembled L–R–TF data frame.}
#'   \item{\code{kb$tf_target}}{TF–target data frame (\code{tf}, \code{target}).}
#'   \item{\code{autocrine}}{Named list; each element is a cell type and is a
#'     list with elements \code{network} (L–R–TF–Target) and
#'     \code{lr_unique} (unique L/R combos).}
#'   \item{\code{pairwise}}{Named list; each element is "<A>_to_<B>" and is a
#'     list with elements \code{network} and \code{lr_unique}.}
#'   \item{\code{counts}}{Data frame of LR pair counts; rows/columns are
#'     \code{cell_types}.}
#' }
#'
#' If \code{output_dir} is non-\code{NULL}, multiple CSV files are also
#' written to disk (KB tables, per-network tables, LR-unique tables and the
#' summary LR count matrix).
#'
#' @export
compute_cc_networks_devcom <- function(
    species = c("human","mouse","bos","sheep","mmulatta","mfascicularis"),
    cell_types,
    filtered_dir,
    output_dir = NULL,
    require_all_ligands = TRUE,
    require_all_receptors = TRUE,
    build_autocrine = TRUE,
    build_pairwise = TRUE,
    write_lr_unique = TRUE,
    write_counts = TRUE,
    verbose = TRUE
){
  species <- match.arg(species)
  stopifnot(length(cell_types) >= 1)
  stopifnot(dir.exists(filtered_dir))
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(DevcomDB)
  })
  
  .log <- function(...) if (isTRUE(verbose)) cat(sprintf(...), "\n")
  .dir_create <- function(p) if (!is.null(p)) dir.create(p, showWarnings = FALSE, recursive = TRUE)
  
  if (!is.null(output_dir)) .dir_create(output_dir)
  
  ## ===== 1) Load DevcomDB priors and assemble L–R–TF KB =====
  data_list <- DevcomDB::load_interaction_data()
  
  key_lr <- paste0(species, "_ligand_receptor")
  key_rt <- paste0(species, "_receptor_tf")
  key_tt <- paste0(species, "_tf_target_gene")
  
  stopifnot(all(c(key_lr, key_rt, key_tt) %in% names(data_list)))
  
  ligand_receptor_df <- data_list[[key_lr]]
  receptor_tf_df     <- data_list[[key_rt]]
  tf_target_df       <- data_list[[key_tt]]
  
  # Standardise column names: receptor_tf should have columns receptor, tf;
  # tf_target should have columns tf, target (first two columns).
  if (ncol(receptor_tf_df) < 2)
    stop("receptor_tf table must contain at least two columns: receptor, tf")
  names(receptor_tf_df)[1:2] <- c("receptor", "tf")
  
  if (ncol(tf_target_df) < 2)
    stop("tf_target table must contain at least two columns: tf, target")
  names(tf_target_df)[1:2] <- c("tf", "target")
  
  # Automatically detect ligand and receptor columns:
  # all columns whose names start with 'L' or 'R' (case-insensitive).
  all_names <- names(ligand_receptor_df)
  L_cols <- all_names[grepl("^L", all_names, ignore.case = TRUE)]
  R_cols <- all_names[grepl("^R", all_names, ignore.case = TRUE)]
  
  if (length(L_cols) == 0 || length(R_cols) == 0) {
    stop(
      "No ligand/receptor columns detected. ",
      "Columns starting with 'L' (ligands) and 'R' (receptors) are required."
    )
  }
  
  # Attach an internal row id so we can pivot receptors and then re-attach
  # the full original ligand-receptor complex.
  ligand_receptor_df <- ligand_receptor_df %>%
    dplyr::mutate(.row_id = dplyr::row_number())
  
  # Map each LR row (.row_id) to one or more TFs via receptor_tf:
  tf_map <- ligand_receptor_df %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(R_cols),
      names_to  = ".Rcol",
      values_to = "receptor"
    ) %>%
    dplyr::filter(!is.na(receptor), receptor != "") %>%
    dplyr::left_join(
      receptor_tf_df %>% dplyr::distinct(receptor, tf),
      by = "receptor",
      relationship = "many-to-many"  # many-to-many is expected here
    ) %>%
    dplyr::filter(!is.na(tf), tf != "") %>%
    dplyr::distinct(.row_id, tf)
  
  # Build L–R–TF table by joining back to the original LR table,
  # keeping all L*/R* columns plus TF.
  lr_tf <- tf_map %>%
    dplyr::left_join(ligand_receptor_df, by = ".row_id") %>%
    dplyr::select(dplyr::all_of(L_cols), dplyr::all_of(R_cols), tf) %>%
    dplyr::distinct()
  
  # Optionally write KB tables
  if (!is.null(output_dir)) {
    readr::write_csv(
      lr_tf,
      file.path(output_dir, paste0(species, "_L-R-TF_develop.csv"))
    )
    readr::write_csv(
      tf_target_df %>% dplyr::select(tf, target),
      file.path(output_dir, paste0(species, "_tf_target_develop.csv"))
    )
  }
  
  ## ===== 2) Read hub genes for each cell type =====
  read_hub <- function(cell){
    fp <- file.path(filtered_dir, paste0(cell, "_filtered_genes.csv"))
    if (!file.exists(fp)) stop("Hub-gene file not found: ", fp)
    df <- readr::read_csv(fp, show_col_types = FALSE)
    cn <- names(df)
    gene_col <- if ("Gene" %in% cn) "Gene" else cn[1]
    unique(stats::na.omit(df[[gene_col]]))
  }
  hub_list <- setNames(lapply(cell_types, read_hub), cell_types)
  
  ## ===== 3) Helper: filter LR–TF rows by sender/receiver gene sets =====
  filter_lr_tf_by_sets <- function(df, sender_genes, receiver_genes){
    # L columns correspond to sender; R columns correspond to receiver.
    Ldf <- df %>% dplyr::select(dplyr::all_of(L_cols))
    Rdf <- df %>% dplyr::select(dplyr::all_of(R_cols))
    
    # Non-empty entries ("" or NA are considered empty).
    nonempty_L <- as.data.frame(lapply(Ldf, function(x) !is.na(x) & x != ""))
    nonempty_R <- as.data.frame(lapply(Rdf, function(x) !is.na(x) & x != ""))
    
    in_sender <- as.data.frame(lapply(Ldf, function(x) x %in% sender_genes))
    in_recv   <- as.data.frame(lapply(Rdf, function(x) x %in% receiver_genes))
    
    L_hit  <- rowSums(as.matrix(nonempty_L) & as.matrix(in_sender))
    L_need <- rowSums(as.matrix(nonempty_L))
    
    R_hit  <- rowSums(as.matrix(nonempty_R) & as.matrix(in_recv))
    R_need <- rowSums(as.matrix(nonempty_R))
    
    ok_L <- if (require_all_ligands) L_hit == L_need else L_hit >= pmin(1L, L_need)
    ok_R <- if (require_all_receptors) R_hit == R_need else R_hit >= pmin(1L, R_need)
    
    df[ok_L & ok_R, , drop = FALSE]
  }
  
  ## ===== 4) Build L–R–TF–Target networks and LR-unique helpers =====
  make_network <- function(lr_tf_sub, tf_target_df){
    if (nrow(lr_tf_sub) == 0) return(tibble())
    lr_tf_sub %>%
      dplyr::inner_join(tf_target_df %>% dplyr::select(tf, target), by = "tf") %>%
      dplyr::select(dplyr::all_of(L_cols), dplyr::all_of(R_cols), tf, target) %>%
      dplyr::distinct()
  }
  
  lr_unique <- function(df){
    if (nrow(df) == 0) return(tibble())
    df %>%
      dplyr::select(dplyr::all_of(L_cols), dplyr::all_of(R_cols)) %>%
      dplyr::distinct()
  }
  
  ## ===== 5) Build autocrine and pairwise networks =====
  results <- list(
    kb = list(
      lr_tf     = lr_tf,
      tf_target = tf_target_df %>% dplyr::select(tf, target)
    ),
    autocrine = list(),
    pairwise  = list(),
    counts    = NULL
  )
  
  ## --- autocrine (sender == receiver) ---
  if (isTRUE(build_autocrine)) {
    for (cell in cell_types) {
      .log("[autocrine] %s ...", cell)
      sender <- hub_list[[cell]]
      recv   <- sender
      lr_tf_sub <- filter_lr_tf_by_sets(lr_tf, sender, recv)
      net  <- make_network(lr_tf_sub, tf_target_df)
      uniq <- lr_unique(net)
      results$autocrine[[cell]] <- list(network = net, lr_unique = uniq)
      
      if (!is.null(output_dir)) {
        readr::write_csv(
          net,
          file.path(output_dir, paste0(cell, "_L-R-TF-T_network-auto.csv"))
        )
        if (isTRUE(write_lr_unique)) {
          readr::write_csv(
            uniq,
            file.path(output_dir, paste0(cell, "_ligand_receptor_unique_auto.csv"))
          )
        }
      }
    }
  }
  
  ## --- pairwise (ordered pairs: i -> j, i != j) ---
  if (isTRUE(build_pairwise)) {
    for (i in seq_along(cell_types)) {
      for (j in seq_along(cell_types)) {
        if (i == j) next
        c1 <- cell_types[i]
        c2 <- cell_types[j]
        
        .log("[pairwise] %s -> %s ...", c1, c2)
        lr_tf_sub <- filter_lr_tf_by_sets(lr_tf, hub_list[[c1]], hub_list[[c2]])
        net  <- make_network(lr_tf_sub, tf_target_df)
        uniq <- lr_unique(net)
        key  <- paste0(c1, "_to_", c2)
        results$pairwise[[key]] <- list(network = net, lr_unique = uniq)
        
        if (!is.null(output_dir)) {
          readr::write_csv(
            net,
            file.path(output_dir, paste0(key, "_L-R-TF-T_network.csv"))
          )
          if (isTRUE(write_lr_unique)) {
            readr::write_csv(
              uniq,
              file.path(output_dir, paste0(key, "_ligand_receptor_unique.csv"))
            )
          }
        }
      }
    }
  }
  
  ## ===== 6) LR count matrix across cell types =====
  if (isTRUE(write_counts) || !is.null(output_dir)) {
    n <- length(cell_types)
    mat <- matrix(0L, nrow = n, ncol = n,
                  dimnames = list(cell_types, cell_types))
    
    # autocrine counts
    if (isTRUE(build_autocrine)) {
      for (cell in cell_types) {
        uniq <- results$autocrine[[cell]]$lr_unique
        mat[cell, cell] <- if (is.null(uniq)) 0L else nrow(uniq)
      }
    }
    
    # pairwise counts
    if (isTRUE(build_pairwise)) {
      for (i in seq_along(cell_types)) {
        for (j in seq_along(cell_types)) {
          if (i == j) next
          key  <- paste0(cell_types[i], "_to_", cell_types[j])
          uniq <- results$pairwise[[key]]$lr_unique
          mat[i, j] <- if (is.null(uniq)) 0L else nrow(uniq)
        }
      }
    }
    
    results$counts <- as.data.frame(mat)
    if (!is.null(output_dir) && isTRUE(write_counts)) {
      write.csv(
        results$counts,
        file.path(output_dir, "cell_communication_summary-LR.csv"),
        row.names = TRUE
      )
    }
  }
  
  results
}
