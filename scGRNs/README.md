# scGRN code and supporting files
1. Initial preprocessing is done in R (scGRN1.Preprocess.Batch.R)
2. GENIE3/SCENIC is implemented through arboreto and scPlant (scGRN2.run_scenic_all_clusters.sh)
3. Shuffled weights are identified and used as cutoffs for significance (scGRN3.collect_adj.tsv.sh)

Network visualization was done manually in Cytoscape, and the .cys files are provided here
