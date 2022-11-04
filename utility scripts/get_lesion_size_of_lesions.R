# read the lesion size of a set of lesion

# load package
library(RNifti)

# provide the path to the folder that contains the lesions in NIFTI format
# [make sure that all images have the same resolution to be comparable!]
path = "D:/Arbeit_Bern/Projekt_Bayesian_LSM/Data_Prospective"


files <- list.files(path, pattern = "\\.nii$")

lesion_size_output=matrix(0,length(files),2)
colnames(lesion_size_output) <- c('Image name','Lesion size')

for (x in 1:length(files)){
 lesion_path <- paste(path, "/", files[x],sep="")
 lesion <- readNifti(lesion_path)

 # make sure that lesions are binary
 lesion <- as.numeric(lesion>0)

 lesion_size_output[x,1]=files[x]
 lesion_size_output[x,2]=sum(lesion)
}

### OUTPUT
# check variable 'lesion_size_output' for lesion sizes across all images
# size is in voxels