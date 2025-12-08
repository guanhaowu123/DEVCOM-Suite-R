##go-annotation#

library("org.Hs.eg.db")
res <- compute_go_blocks_and_integrate(
  input_dir  = "E:/DEVCOM Suite/tests/testthat/go-annotation",
  output_dir = "E:/DEVCOM Suite/tests/testthat/go-annotation",
  org_db     = org.Hs.eg.db, #human
  stages     = c("extract_blocks","build_complex","alias_normalize","go_annot","integrate"),
  which_files = c("autocrine","pairwise"),   # Only autocrine: c("autocrine"); Only intercellular interaction: c("pairwise")
  use_obo    = TRUE,                         # Perform GO term intersection analysis using OBO ancestor sets (recommended for LR/RTF analysis)
  verbose    = TRUE
)




