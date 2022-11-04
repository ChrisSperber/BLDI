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
# (Attention: the symbol - R uses "/" (NOT: "\") in path names!)
path = "D:/Folder/Data"

# provide the path to the csv file that contains the scores
# (you can change the separator by changing the sep argument)
dat_scores <- read.csv("D:/Folder/Deficit_Scores.csv",
header = TRUE,sep = ";")

# define the minimum number of lesions in a voxel for the voxel to be included
min_vox_threshold <- 5

# provide the path to folder of BLDI main script
script_path <- "D:/Folder/BLDI"

####################################################
# end of input required by user
####################################################

# flag for binary lesion data - 1 yes, 0 no
binary_data <- 1

# set working directory to path of the scripts
setwd(script_path)

source(paste(script_path,"/bldi_create_analysis_mask.R",sep = ""))
source(paste(script_path,"/bldi_load_lesions.R",sep = ""))
source(paste(script_path,"/bldi_create_ttest_bf_map.R",sep = ""))

# create a lesion overlap and get the index of all voxels to be included in the analysis
output_overlap <- bldi_create_analysis_mask(dat_scores,path,min_vox_threshold)
idx_voxels_included<-unlist(output_overlap[1])

# report number of voxels
writeLines(sprintf("%d voxles are damaged in at least %d patients \n and will be included in the analysis",
length(idx_voxels_included),min_vox_threshold))

####################################################

# load the lesions
lesions_full <- bldi_load_lesions(dat_scores,path,idx_included_voxels,binary_data)
print("Done - Loading lesion data")


print(sprintf("Starting Bayes Factor mapping for score %s", colnames(dat_scores)[2]))
print("The procedure may take up to several hours! Please wait!")

BF_result <- bldi_create_ttest_bf_map(dat_scores,idx_voxels_included,lesions_full)

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

filename <- paste("Output_BLDI_BF_map_ttest_",as.character(colnames(dat_scores)[2]),"_",format(Sys.time(), "%Y-%m-%d_%H_%M"),sep = "")
writeNifti(BF_image, filename)

# additionally, create a log(BF) map
filename <- paste("Output_BLDI_logBF_map_ttest_",as.character(colnames(dat_scores)[2]),"_",format(Sys.time(), "%Y-%m-%d_%H_%M"),sep = "")
writeNifti(log10(BF_image), filename)

## create a results txt with some additional information
filename <- paste("Output_BLDI_info_ttest_",as.character(colnames(dat_scores)[2]),"_",format(Sys.time(), "%Y-%m-%d_%H_%M"), ".txt",sep = "")
writeLines("Results of voxel-wise BLDI with t-tests", filename)
line <- paste("For deficit",as.character(colnames(dat_scores)[2]))
write(line,file = filename,append = TRUE)
write("",file = filename,append = TRUE)
line <- paste("Voxels with at least ",as.character(min_vox_threshold),"lesions were analysed")
write(line,file = filename,append = TRUE)
line <- paste(as.character(length(idx_voxels_included)),"voxels were analysed")
write(line,file = filename,append = TRUE)

line <- paste("The maximum BF is ", as.character(max(BF_result)), ", the minimum BF is ", as.character(min(BF_result)))
write(line,file = filename,append = TRUE)
## end of txt file creation

print("Done - Analyses finished!")
