###Active and Inhibitive##


res <- compute_active_inhibit(
  expr_dir   = "E:/DEVCOM Suite/tests/testthat/IVF",
  kb_path    = "E:/DEVCOM Suite/tests/testthat/mouse_L-R-TF_develop.csv",
  output_dir = "E:/DEVCOM Suite/tests/testthat/active_inhibit/",
  test = "none",          # "none" | "wilcoxon" | "ttest"
  min_expr = 0,            # low-expression filtering threshold
  include_pairwise  = TRUE,
  ordered_pairs     = FALSE,
  p_adjust = "none",        # or "none"
  p_min = 1e-300,
  write_coverage = TRUE,  # Whether to additionally output the gene coverage matrix and missing value statistics?
  verbose = TRUE
)

