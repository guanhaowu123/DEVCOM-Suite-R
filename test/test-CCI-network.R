
##CCI network##


res <- compute_cc_networks_devcom(
  species      = "human",
  cell_types   = c("EVT","CTB"),
  filtered_dir = "E:/DEVCOM Suite/tests/testthat/CCI network",
  output_dir   = "E:/DEVCOM Suite/tests/testthat/CCI network",
  require_all_ligands   = TRUE,
  require_all_receptors = TRUE,
  build_autocrine = TRUE,
  build_pairwise  = TRUE,
  write_lr_unique = TRUE,
  write_counts    = TRUE,
  verbose = TRUE
)

