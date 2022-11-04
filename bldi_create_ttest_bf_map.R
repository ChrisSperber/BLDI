bldi_create_ttest_bf_map <- function(dat_scores,idx_voxels_included,lesions_full) {


# create blank for the BF for each voxel's lesion-deficit association
BF_result <- matrix(0, length(idx_voxels_included))

scores_target <- unlist(dat_scores[2])


### loop through all voxels
for (voxel in 1:length(idx_voxels_included)){
voxel_status <- lesions_full[,voxel]


bf = ttestBF(scores_target[voxel_status==1],scores_target[voxel_status==0])
BF_result[voxel] <- exp(unlist(bf@bayesFactor[1]))
}
###

return(BF_result)
}
