# Bayes Factor Lesion-Deficit Mapping Disconnection
# Main script
####################################################
#
# Within this script, the main analysis parameters and the
# paths to the data have to be defined
####################################################

# empty working space
rm(list = ls())

# load packages
print("Loading dependencies")
library(BayesFactor)

####################################################
# start of input required by user
####################################################

# provide the path to the folder that contains the disconnection matrices in .txt format
# (Attention - R uses "/" (NOT: "\") in path names!)
path = "D:/Folder/Disconnection_matrices"

# provide the path to the csv file that contains the scores
# (you can change the separator by changing the sep argument)
dat_scores <- read.csv("D:/Folder/Deficit_Scores.csv",
header = TRUE,sep = ";")

# define the minimum number of disconnection values >0 for an element to be included
min_disc_threshold <- 5

# provide the path to folder of BLDI main script
script_path <- "D:/Folder/BLDI"

####################################################
# end of input required by user
####################################################

# set working directory to path of the scripts
setwd(script_path)

source(paste(script_path,"/bldi_load_disc_matrices.R",sep = ""))
source(paste(script_path,"/bldi_create_glm_bf_disc.R",sep = ""))


# the inclusion of covariates in GLMs is hard-coded to include maximally 2 covariates,
if (dim(dat_scores)[2]>4){
print("A maximum of 2 covariates can be included.")
print(sprintf("Currently, %d covariates are included, and all after the second one are ignored", dim(dat_scores)[2]-2))
}

####################################################

# load the matrices
disconns_full <- bldi_load_disc_matrices(dat_scores,path)
print("Done - Loading disconnection data")

# get the index of all elements to be included in the analysis
binary_disconns <- as.matrix(disconns_full>0)+0 # the awkward +0 transforms the matrix from logical to numeric
disc_overlap <- colSums(binary_disconns)
elements_included <- as.numeric(disc_overlap>=min_disc_threshold)
idx_elements_included <- which(elements_included == 1)

# report number of tested elements
writeLines(sprintf("%d elements are damaged in at least %d patients \n and will be included in the analysis",
length(idx_elements_included),min_disc_threshold))

print(sprintf("Starting Bayes Factor mapping for score %s", colnames(dat_scores)[2]))
print("The procedure may take up to several hours! Please wait!")

BF_result <- bldi_create_glm_bf_disc(dat_scores,idx_elements_included,disconns_full)

####################################################

# re-load the first matrix as a reference for dimensions
disconn_path <- paste(path, "/", dat_scores[1,1], ".csv", sep="")
disconn <- read.table(disconn_path)
dimension_1st_matrix <- dim(disconn)
disc_vect_length <- dim(disconn)[1]*dim(disconn)[2]


# rebuild a matrix of the BF
BF_full_vect <- matrix(0, disc_vect_length)
BF_full_vect_log <- matrix(0, disc_vect_length)
BF_full_vect[idx_elements_included] <- BF_result
BF_full_vect_log[idx_elements_included] <- log10(BF_result)

# generate a text file of the BF matrix
BF_full_matrix <- array(BF_full_vect, dimension_1st_matrix)
BF_full_matrix_log <- array(BF_full_vect_log, dimension_1st_matrix)

filename <- paste("Output_BLDI_BF_disconnetion_glm_",as.character(colnames(dat_scores)[2]),"_",format(Sys.time(), "%Y-%m-%d_%H_%M"),".txt",sep = "")
write.table(BF_full_matrix, filename, append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)

# additionally, create a log(BF) map
filename <- paste("Output_BLDI_logBF_disconnection_glm_",as.character(colnames(dat_scores)[2]),"_",format(Sys.time(), "%Y-%m-%d_%H_%M"),".txt",sep = "")
write.table(BF_full_matrix_log, filename, append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)


## create a results txt with some additional information
filename <- paste("Output_BLDI_info_disconnection_glm_",as.character(colnames(dat_scores)[2]),"_",format(Sys.time(), "%Y-%m-%d_%H_%M"), ".txt",sep = "")
writeLines("Results of disconnection BLDI with GLMs", filename)
line <- paste("For deficit",as.character(colnames(dat_scores)[2]))
write(line,file = filename,append = TRUE)
write("",file = filename,append = TRUE)
line <- paste("Elements with at least ",as.character(min_disc_threshold),"disconnections were analysed")
write(line,file = filename,append = TRUE)
line <- paste(as.character(length(idx_elements_included)),"elements were analysed (symmetric elements not excluded!)")
write(line,file = filename,append = TRUE)


line <- paste("The maximum BF is ", as.character(max(BF_result)), ", the minimum BF is ", as.character(min(BF_result)))
write(line,file = filename,append = TRUE)

if(dim(dat_scores)[2]>2){
covar_num <- min(c(2,(dim(dat_scores)[2])-2))
line <- paste(as.character(covar_num), "covariate(s) were included as specified in the .csv scores file")
write(line,file = filename,append = TRUE)
}
## end of txt file creation

print("Done - Analyses finished!")
