
cell_types <- c("Bergmann_glial_cell", "B_cell", "Club_cell", "Erythroid_cell",
                "Endothelial_cell","GABAergic_neuron",
                "Interneuron","Neural_stem_cell","Stem_cell","Tanycyte")

plot_lr_network(
  in_dir   = "E:\\DEVCOM Suite\\tests\\testthat\\plot_network\\",
  cell_types = cell_types,
  ligand_of_interest   = "MIF",
  receptor_of_interest = "ITGA4"
)
