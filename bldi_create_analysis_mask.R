bldi_create_analysis_mask <- function(dat_scores,path,min_vox_threshold,lesion_size_control = 0) {

print("Creating lesion overlap for pre-analysis")

# get first example lesion
lesion_path <- paste(path, "/", dat_scores[1,1], ".nii", sep="")
lesion <- readNifti(lesion_path)
dimension_1st_lesion <- dim(lesion)
lesion_vect <- c(lesion)
lesion_vect_length <- length(lesion_vect)

# blank vector for lesion size
lesion_size <-  rep(0, dim(dat_scores)[1])

# reading a full voxel-wise lesion sample with common resolutions would
# exceed the working memory (for 300 vectorised lesions at 182x218x182 >16GB)
# therefore, an overlap map is generated in a first step and all voxels 
# lesioned in at least x patients are pre-selected by a mask to be 
# stored in the next step
overlap <- array(rep(0, dimension_1st_lesion[1]*dimension_1st_lesion[2]*dimension_1st_lesion[3]),
dim=c(dimension_1st_lesion[1],dimension_1st_lesion[2],dimension_1st_lesion[3])) 

# loop through lesions, add each binarised lesion is added to the overlap
for (x in 1:dim(dat_scores)[1]){
lesion_path <- paste(path, "/", dat_scores[x,1], ".nii", sep="")
lesion <- readNifti(lesion_path)
dimension_lesion <- dim(lesion)

# check if the image size is equal to the first lesion
if (dimension_lesion[1] != dimension_1st_lesion[1] || dimension_lesion[2] != dimension_1st_lesion[2] ||
 dimension_lesion[3] != dimension_1st_lesion[3]){
print(sprintf("WARNING - ERROR! The selected lesions are not all of equal size. The first lesion had the dimensions: %d x %d x %d",
 dimension_1st_lesion[1],dimension_1st_lesion[2],dimension_1st_lesion[3]))
print(sprintf("The lesion in line %s (and potentially further lesions) has dimensions: %d x %d x %d", lesion_path,
 dimension_lesion[1],dimension_lesion[2],dimension_lesion[3]))
}

if (lesion_size_control == 1){
lesion_size[x]= sum(as.numeric(lesion > 0))
}


overlap <- overlap + as.numeric(lesion > 0)
}

# indicate all voxels that are above the threshold for minimum lesion number
included_voxels=as.numeric(c(overlap>=min_vox_threshold))
idx_voxels_included=which(included_voxels==1)

output_overlap <- list(idx_voxels_included,lesion_size)

return(output_overlap)
}