#merge GO files

res_files <- integrate_go_lr_rtf_tft(
  lr_dir  = "E:/DEVCOM Suite/tests/testthat/go-annotation/LR",
  rtf_dir = "E:/DEVCOM Suite/tests/testthat/go-annotation/RTF",
  tft_dir = "E:/DEVCOM Suite/tests/testthat/go-annotation/TFT",
  out_dir = "E:/DEVCOM Suite/tests/testthat/go-annotation/integrated",
  verbose = TRUE
)

