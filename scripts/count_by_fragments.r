# command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(interactive() | args[1] == "test") {
  print("running in test mode")
  fragments_file = "10X_PBMC/01_raw_data/ATAC/head_fragments.tsv.gz"
  project_name = "PBMC_10X"
  unique_enhancers_path = "GSE126074_SNARE_seq/05_non_gene_body_peaks/uniqe_enhancers.bed"
  whitelist = "10X_PBMC/01_raw_data/gex_matrix/filtered_feature_bc_matrix/barcodes.tsv.gz"
  output_dir = "10X_PBMC/ATAC/01_count_fragments/"
  threads = 1
} else {
  # running from command line
  # Rscript scripts/count_by_fragments.r /path/to/fragments.tsv.gz MyProject
  fragments_file = args[1]
  project_name = args[2]
  unique_enhancers_path = args[3]
  whitelist = args[4]
  output_dir = args[5]
  out_file = args[6]
  threads = as.numeric(args[7])

}
whitelist = read.table(gzfile(whitelist), header = FALSE, sep = "\t")[[1]] 
# log files
createArrowFiles_log = file.path(output_dir, "ArchR_createArrowFiles.log")
addFeatureMatrix_log = file.path(output_dir, "ArchR_addFeatureMatrix.log")
getMatrixFromProject_log = file.path(output_dir, "ArchR_getMatrixFromProject.log")
library(ArchR)
library(magrittr)
addArchRThreads(threads)
addArchRGenome("hg38")

inputFiles = c(fragments_file) %>% setNames(project_name)

# Create Arrow Files and ArchR Project
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  validBarcodes = whitelist,
  sampleNames = names(inputFiles),
  minTSS = 0, #keep all cells 
  minFrags = 0, #keep all cells 
  maxFrags = Inf, #keep all cells 
  minFragSize  = 10, #default
  maxFragSize  = 2000,  #default
  addTileMat = F,
  addGeneScoreMat = F,
  logFile = createArrowFiles_log,
  QCDir = file.path(output_dir, "QCDir"),
  threads = threads
  )

# Create ArchR Project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = output_dir,
  copyArrows = T #This is recommened so that you maintain an unaltered copy for later usage.
)

#remove origin arrow file
file.remove(ArrowFiles)

# read unique enhancers and convert to GRanges
library(GenomicRanges)
unique_enhancers <- read.table(unique_enhancers_path, header = F) 

unique_enhancers_gr <- GRanges(
  seqnames = paste0("chr", unique_enhancers$V1),
  ranges = IRanges(start = unique_enhancers$V2,
                   end = unique_enhancers$V3),
  name = unique_enhancers$V4
)


# Add Feature Matrix for unique enhancers
addFeatureMatrix( 
  input = proj,
  features = unique_enhancers_gr,
  matrixName = "EnhancersMatrix",
  ceiling = 10^9,
  binarize = FALSE,
  verbose = TRUE,
  threads = threads,
  force = TRUE,
  logFile = addFeatureMatrix_log
)
# Get the counts matrix
enhancer_counts <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "EnhancersMatrix",
  logFile = getMatrixFromProject_log
)
# save counts matrix as mtx
writeMM(assays(enhancer_counts)$EnhancersMatrix, file = file.path(output_dir, paste0(project_name, "_enhancer_counts.mtx")))


# save features and barcodes
features <- rowData(enhancer_counts)$name
barcodes <- colnames(enhancer_counts)
# remove prefix from barcodes i.e  "pbmc_granulocyte_sorted_10k#TGAGCAAAGGCCTGGT-1"
barcodes <- sub( "^.*#", "", barcodes)

write.table(as.data.frame(features), file = file.path(output_dir, paste0(project_name, "_features.tsv")), sep = "\t", quote = F, row.names = F, col.names = F)
write.table(as.data.frame(barcodes), file = file.path(output_dir, paste0(project_name, "_barcodes.tsv")), sep = "\t", quote = F, row.names = F, col.names = F)


