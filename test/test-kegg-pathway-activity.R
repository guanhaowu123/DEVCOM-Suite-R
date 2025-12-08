# 载入
#source("kegg_pathway_activity.R")

kegg_pathway_activity_dir(
  input_dir   = "E:\\DEVCOM Suite\\tests\\testthat\\go-annotation_kegg_result\\integrated\\",
  output_dir  = "E:\\DEVCOM Suite\\tests\\testthat\\go-annotation_kegg_result\\integrated\\keggpathway score",
  pattern     = "_merged_LR_RTF_TFT\\.csv$",
  ligand_suffix = "_gene_importance_scores_transformed.csv",
  target_suffix = "_gene_importance_scores_transformed.csv",
  split_complex = FALSE,
  weighting     = "count"
)


#  seclect developmental pathway
filter_kegg_activity_by_keywords(
  input_dir   = "E:\\DEVCOM Suite\\tests\\testthat\\go-annotation_kegg_result\\integrated\\keggpathway score",
  output_dir  = "E:\\DEVCOM Suite\\tests\\testthat\\go-annotation_kegg_result\\integrated\\keggpathway score\\development pathway",
  keywords_file = "E:/DEVCOM Suite/tests/testthat/development.txt"  # key words
)
