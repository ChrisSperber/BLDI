# Bayes Factor Lesion-Deficit Mapping 
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
library(RNifti)
library(BayesFactor)

####################################################
# start of input required by user
####################################################

# provide the path to the folder that contains the lesions in NIFTI format
# (Attention - R uses "/" (NOT: "\") in path names!)
path = "D:/Arbeit_Bern/Projekt_Bayesian_LSM/Data"

# provide the path to the csv file that contains the scores
# (you can change the separator by changing the sep argument)
dat_scores <- read.csv("D:/Folder/Deficit_Scores.csv",
header = TRUE,sep = ";")

# define the minimum number of lesions in a voxel for the voxel to be included
min_vox_threshold <- 5

# flag for binary lesion data - 1 yes, 0 no
# (only set this to "1" if you have binary 0 - 1 data!)
binary_data <- 1

# provide the path to folder of BLDI main script
script_path <- "D:/Folder/BLDI"

# flag for adaptive control of lesion size (only available for GLM) - 1 yes, 0 no
# (Attention: lesion size is computed based in a binarised map > 0, and
# lesion size control as computed in the script is intended for binary 
# voxelwise data only)
lesion_size_control <- 0

####################################################
# end of input required by user
####################################################

# set working directory to path of the scripts
setwd(script_path)

source(paste(script_path,"/bldi_create_analysis_mask.R",sep = ""))
source(paste(script_path,"/bldi_load_lesions.R",sep = ""))
source(paste(script_path,"/bldi_create_glm_bf_map.R",sep = ""))
source(paste(script_path,"/bldi_create_glm_bf_map_ls_cont.R",sep = ""))

# create a lesion overlap and get the index of all voxels to be included in the analysis
output_overlap <- bldi_create_analysis_mask(dat_scores,path,min_vox_threshold,lesion_size_control)
idx_voxels_included<-unlist(output_overlap[1])

# report number of voxels
writeLines(sprintf("%d voxles are damaged in at least %d patients \n and will be included in the analysis",
length(idx_voxels_included),min_vox_threshold))

# if lesion size control is desired, lesion size is now added as column
# to the scores data frame
if (lesion_size_control==1){
print("Adding lesion size as a covariate into the GLM; adaptive (!) lesion size control is applied.")
dat_scores[,dim(dat_scores)[2]+1] <- unlist(output_overlap[2])
}

# the inclusion of covariates in GLMs is hard-coded to include maximally 2 covariates,
# of which 1 may be lesion size.
if (dim(dat_scores)[2]>4){
print("A maximum of 2 covariates can be included (of which 1 may be lesion size).")
print(sprintf("Currently, %d covariates are included, and all after the second one are ignored", dim(dat_scores)[2]-2))
if (lesion_size_control==1){
print("Adaptive lesion size control is now removed from the analysis.")
lesion_size_control <- 0
}
}

####################################################

# load the lesions
lesions_full <- bldi_load_lesions(dat_scores,path,idx_voxels_included,binary_data)
print("Done - Loading lesion data")


print(sprintf("Starting Bayes Factor mapping for score %s", colnames(dat_scores)[2]))
print("The procedure may take up to several hours! Please wait!")

if (lesion_size_control==0){
BF_result <- bldi_create_glm_bf_map(dat_scores,idx_voxels_included,lesions_full)
} else if(lesion_size_control==1){
BF_result <- bldi_create_glm_bf_map_ls_cont(dat_scores,idx_voxels_included,lesions_full)
}

####################################################

# re-load the first lesion as a reference for image image
lesion_path <- paste(path, "/", dat_scores[1,1], ".nii", sep="")
lesion <- readNifti(lesion_path)
dimension_1st_lesion <- dim(lesion)
lesion_vect_length <- length(c(lesion))


# rebuild a 3D image of the BF
BF_full_vect <- matrix(0, lesion_vect_length)
BF_full_vect[idx_voxels_included] <- BF_result

# generate a NIFTI of the BF map
BF_full_image <- array(BF_full_vect, dimension_1st_lesion)
# take the first lesion image as reference for header info
BF_image <- asNifti(BF_full_image, reference = lesion_path)

filename <- paste("Output_BLDI_BF_map_glm_",as.character(colnames(dat_scores)[2]),"_",format(Sys.time(), "%Y-%m-%d_%H_%M"),sep = "")
writeNifti(BF_image, filename)

# additionally, create a log(BF) map
filename <- paste("Output_BLDI_logBF_map_glm_",as.character(colnames(dat_scores)[2]),"_",format(Sys.time(), "%Y-%m-%d_%H_%M"),sep = "")
writeNifti(log10(BF_image), filename)

## create a results txt with some additional information
filename <- paste("Output_BLDI_info_glm_",as.character(colnames(dat_scores)[2]),"_",format(Sys.time(), "%Y-%m-%d_%H_%M"), ".txt",sep = "")
writeLines("Results of voxel-wise BLDI with GLMs", filename)
line <- paste("For deficit",as.character(colnames(dat_scores)[2]))
write(line,file = filename,append = TRUE)
write("",file = filename,append = TRUE)
line <- paste("Voxels with at least ",as.character(min_vox_threshold),"lesions were analysed")
write(line,file = filename,append = TRUE)
line <- paste(as.character(length(idx_voxels_included)),"voxels were analysed")
write(line,file = filename,append = TRUE)

if (lesion_size_control==1){
 write("Adaptive lesion size control was included as a covariate",file = filename,append = TRUE)
}

line <- paste("The maximum BF is ", as.character(max(BF_result)), ", the minimum BF is ", as.character(min(BF_result)))
write(line,file = filename,append = TRUE)

if(dim(dat_scores)[2]>2){
covar_num <- min(c(2,(dim(dat_scores)[2])-2))
line <- paste(as.character(covar_num), "covariate(s) were included as specified in the .csv scores file (or lesion size, if included)")
write(line,file = filename,append = TRUE)
}
## end of txt file creation

print("Done - Analyses finished!")
