bldi_create_glm_bf_map <- function(dat_scores,idx_voxels_included,lesions_full) {


# create blank for the BF for each voxel's lesion-deficit association
BF_result <- matrix(0, length(idx_voxels_included))

scores_target <- unlist(dat_scores[2])


# for the pontential inclusion of up to 2 covariates,
# if-conditions pick the right formula 


if (dim(dat_scores)[2]==2){
### loop through all voxels
for (voxel in 1:length(idx_voxels_included)){
lesions_temp <- lesions_full[,voxel]
data4bf <- data.frame(scores_target,lesions_temp)

bf = lmBF(scores_target ~ lesions_temp, data = data4bf, progress = FALSE)
BF_result[voxel] <- exp(unlist(bf@bayesFactor[1]))
}
}
###
if (dim(dat_scores)[2]==3){

scores_cov1 <- as.numeric(unlist(dat_scores[3]))

### loop through all voxels
for (voxel in 1:length(idx_voxels_included)){
lesions_temp <- lesions_full[,voxel]
data4bf <- data.frame(scores_target,scores_cov1,lesions_temp)

bf_full = lmBF(scores_target ~ lesions_temp + scores_cov1, data = data4bf, progress = FALSE)
bf_ref =  lmBF(scores_target ~ scores_cov1, data = data4bf, progress = FALSE)


BF_result[voxel] <- (exp(unlist(bf_full@bayesFactor[1]))/exp(unlist(bf_ref@bayesFactor[1])))
}
}
###
if (dim(dat_scores)[2]>=4){

scores_cov1 <- as.numeric(unlist(dat_scores[3]))
scores_cov2 <- as.numeric(unlist(dat_scores[4]))

### loop through all voxels
for (voxel in 1:length(idx_voxels_included)){
lesions_temp <- lesions_full[,voxel]
data4bf <- data.frame(scores_target,scores_cov1,scores_cov2,lesions_temp)

bf_full = lmBF(scores_target ~ lesions_temp + scores_cov1 + scores_cov2, data = data4bf, progress = FALSE)
bf_ref =  lmBF(scores_target ~ scores_cov1 + scores_cov2, data = data4bf, progress = FALSE)


BF_result[voxel] <- (exp(unlist(bf_full@bayesFactor[1]))/exp(unlist(bf_ref@bayesFactor[1])))
}
}
return(BF_result)
}
