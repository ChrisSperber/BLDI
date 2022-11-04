bldi_load_lesions <- function(dat_scores,path,idx_included_voxels,binary_data) {

# create blank for full vectorised lesions
lesions_full <- matrix(0, dim(dat_scores)[1], length(idx_voxels_included))

print("Loading lesion data")

# store the relevant voxels of all lesions
for (x in 1:dim(dat_scores)[1]){
lesion_path <- paste(path, "/", dat_scores[x,1], ".nii", sep="")
lesion <- readNifti(lesion_path)

# make sure that lesions are binary if they are supposed to be
if(binary_data==1){
lesion <- as.numeric(lesion>0)
}
lesions_full[x,] <- c(lesion[idx_voxels_included])
}

return(lesions_full)
}