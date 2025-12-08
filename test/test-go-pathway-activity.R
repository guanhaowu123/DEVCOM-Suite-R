


#source("go_pathway_activity.R")

# 1) only file
res1 <- go_pathway_activity_for_file(
  comm_file         = "E:\\DEVCOM Suite\\tests\\testthat\\go-annotation\\integrated\\CTB_to_EVT_GO_L_R_TF_Target.csv",
  ligand_score_file = "E:\\DEVCOM Suite\\tests\\testthat\\go-annotation\\CTB_gene_importance_scores_transformed.csv",
  target_score_file = "E:\\DEVCOM Suite\\tests\\testthat\\go-annotation\\EVT_gene_importance_scores_transformed.csv",
  output_dir        = "E:\\DEVCOM Suite\\tests\\testthat\\go-annotation\\Go_out\\",
  ontology          = "BP"          # or c("BP","MF") or "ALL"
)

# 2) Batch processing of multiple directories (auto-match X/Y by filename)）
go_pathway_activity_dir(
  input_dir  = "E:\\DEVCOM Suite\\tests\\testthat\\go-annotation\\integrated\\",
  output_dir = "E:\\DEVCOM Suite\\tests\\testthat\\go-annotation\\paythway score",
  ontology   = "BP",
  ligand_suffix = "_gene_importance_scores_transformed.csv"
)
