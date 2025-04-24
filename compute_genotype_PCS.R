library(tidyverse)
library(data.table)
library(optparse)

###### PARSER ############
# Define options
option_list <- list(
  make_option(c("-v", "--vcf_path"), type = "character", default = "",
              help = "gsutil path to VCF file", metavar = "VCF"),
  make_option(c("-p", "--prefix"), type = "character", default = "output",
              help = "prefix for outputfile", metavar = "PREFIX")
)

# Create the parser object
opt_parser <- OptionParser(option_list = option_list)

# Parse the arguments
opt <- parse_args(opt_parser)


###### FUNCTIONS ##########

# function to calculate genetic PCs from VCF file
compute_genetic_PCs <- function(VCF_path,label='converted'){
require(SNPRelate)
set.seed(1000)


# convert VCF to GDS file and prune variants
gds_file <- paste0(label,'.gds')
snpgdsVCF2GDS(VCF_path, gds_file, method="biallelic.only")
genofile <- snpgdsOpen(gds_file)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, maf=0.01)

# run PCA analysis 
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, maf=0.01, eigen.cnt=0, eigen.method = "DSPEV")
    
# run Gavish Donoho to determine number of significant PCs to include
genetic_n_pcs <- chooseGavishDonoho(.dim = c(length(pca$eigenval), length(pca$snp.id)), var.explained = pca$eigenval, noise = 1)
print(paste0(genetic_n_pcs, ' genetic PCs selected'))

# take output PCA data and subset to number of PCs chosen 
# by gavish donoho and create final dataframe
subset_PCA_data <- pca$eigenvect %>% data.frame() %>% 
    select(1:genetic_n_pcs) %>% 
    bind_cols(sample_id = as.numeric(pca$sample.id)) %>% 
    dplyr::rename_with(.cols =starts_with('X'),~str_replace(.,'X','GENETICPC') )
subset_PCA_data   
}


########## START ANALYSIS #########
# load args from opt
input_VCF_path <- opt$VCF
prefix <- opt$PREFIX

# get VCF with path removed. This is used to
# to load the VCF after it is copied to the
# bucket
base_VCF_name <- basename(input_VCF_path)

print('Copying VCF file to bucket')
system(paste0('gsutil cp ',input_VCF_path,' .'))

# create output file name
ouput_file_name <- paste0(prefix,'_genetic_PCs.tsv')


# compute genetic PCs on VCF
genetic_PCs <- base_VCF_name %>% compute_genetic_PCs

# write data to output
genetic_PCs %>% write_tsv(ouput_file_name)
