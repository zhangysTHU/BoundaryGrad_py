args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
  stop("Usage: run_infercnv_external.R raw_counts.mtx genes.tsv cells.tsv annotation.txt reference_cluster.txt gene_order.txt out_dir cnv_calls.tsv threads k")
}

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_file_arg) > 0) sub("^--file=", "", script_file_arg[[1]]) else NA_character_
if (!is.na(script_file)) {
  repo_root <- normalizePath(file.path(dirname(script_file), "..", "..", ".."), winslash = "/", mustWork = FALSE)
  renv_activate <- file.path(repo_root, "scripts_format_R", "renv", "activate.R")
  if (file.exists(renv_activate)) {
    source(renv_activate)
  }
}

raw_counts_file <- args[[1]]
genes_file <- args[[2]]
cells_file <- args[[3]]
annotation_file <- args[[4]]
reference_cluster_file <- args[[5]]
gene_order_file <- args[[6]]
out_dir <- args[[7]]
calls_file <- args[[8]]
threads <- as.integer(args[[9]])
k <- as.integer(args[[10]])

suppressPackageStartupMessages({
  library(Matrix)
  library(infercnv)
  library(ape)
  library(dendextend)
})

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
mat <- Matrix::readMM(raw_counts_file)
genes <- read.delim(genes_file, stringsAsFactors = FALSE)[[1]]
cells <- read.delim(cells_file, stringsAsFactors = FALSE)[[1]]
rownames(mat) <- genes
colnames(mat) <- cells

normal_cluster <- readLines(reference_cluster_file, warn = FALSE)[1]

infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = mat,
  annotations_file = annotation_file,
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = normal_cluster
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = out_dir,
  cluster_by_groups = FALSE,
  analysis_mode = "subclusters",
  denoise = TRUE,
  HMM = TRUE,
  tumor_subcluster_partition_method = "random_trees",
  HMM_type = "i6",
  BayesMaxPNormal = 0,
  num_threads = threads
)

tree_file <- file.path(out_dir, "infercnv.17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations_dendrogram.txt")
obs_file <- file.path(out_dir, "infercnv.17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations.txt")

cell_groupings <- ape::read.tree(file = tree_file)
labels <- as.data.frame(dendextend::cutree(cell_groupings, k = k))
colnames(labels) <- "CNVLabel"
missing_cells <- cells[!cells %in% rownames(labels)]
labels <- rbind(labels, data.frame(row.names = missing_cells, CNVLabel = rep("Normal", length(missing_cells))))

cnv_table <- read.table(obs_file, header = TRUE)
score_table <- abs(as.matrix(cnv_table) - 3)
scores <- as.data.frame(colSums(score_table))
colnames(scores) <- "cnv_score"
rownames(scores) <- gsub("\\.", "-", rownames(scores))

calls <- data.frame(cell_ID = cells)
calls$CNVLabel <- labels$CNVLabel[match(calls$cell_ID, rownames(labels))]
calls$cnv_score <- scores$cnv_score[match(calls$cell_ID, rownames(scores))]
calls$cnv_score[is.na(calls$cnv_score)] <- 0
calls$cnv_score[calls$CNVLabel == "Normal"] <- 0
write.table(calls, calls_file, sep = "\t", row.names = FALSE, quote = FALSE)
