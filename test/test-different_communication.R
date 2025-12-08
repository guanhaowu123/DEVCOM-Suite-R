#Differential communication analysis between two conditions
#install DevcomDB#

library(DevcomDB)

res <- run_diff_comm(
  cell_types = c("Erythrocyte","Endothelial"),
  ctrl_base  = "E:/DEVCOM Suite/tests/testthat/IVF/",
  case_base  = "E:/DEVCOM Suite/tests/testthat/SCNT/",
  lr_pairs_path = "E:/DEVCOM Suite/tests/testthat/mouse_L-R-TF_develop.csv",
  out_dir    = "E:/DEVCOM Suite/tests/testthat/diff_comm_IVF_SCNT",
  case_label = "SCNT",
  ctrl_label = "IVF"
)


