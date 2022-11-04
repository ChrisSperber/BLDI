bldi_create_glm_bf_map_ls_cont <- function(dat_scores,idx_voxels_included,lesions_full) {
# this BLDI glm script includes adaptive lesion size control

# create blank for the BF for each voxel's lesion-deficit association
BF_result <- matrix(0, length(idx_voxels_included))

scores_target <- unlist(dat_scores[2])

# for the pontential inclusion of up to 2 covariates,
# if-conditions pick the right formula 


if (dim(dat_scores)[2]==3){

scores_cov_lesionsize <- unlist(dat_scores[3])

### loop through all voxels
for (voxel in 1:length(idx_voxels_included)){
lesions_temp <- lesions_full[,voxel]
data4bf <- data.frame(scores_target,scores_cov_lesionsize,lesions_temp)

bf = lmBF(scores_target ~ lesions_temp, data = data4bf, progress = FALSE)

if (exp(unlist(bf@bayesFactor[1]))<3){
BF_result[voxel] <- exp(unlist(bf@bayesFactor[1]))
} else {
bf_full = lmBF(scores_target ~ lesions_temp + scores_cov_lesionsize, data = data4bf, progress = FALSE)
bf_ref =  lmBF(scores_target ~ scores_cov_lesionsize, data = data4bf, progress = FALSE)

BF_result[voxel] <- (exp(unlist(bf_full@bayesFactor[1]))/exp(unlist(bf_ref@bayesFactor[1])))
}
}
}
###
if (dim(dat_scores)[2]==4){

scores_cov1 <- unlist(dat_scores[3])
scores_cov_lesionsize  <- unlist(dat_scores[4])

### loop through all voxels
for (voxel in 1:length(idx_voxels_included)){
lesions_temp <- lesions_full[,voxel]
data4bf <- data.frame(scores_target,scores_cov1,scores_cov_lesionsize,lesions_temp)

bf_full = lmBF(scores_target ~ lesions_temp + scores_cov1, data = data4bf, progress = FALSE)
bf_ref =  lmBF(scores_target ~ scores_cov1, data = data4bf, progress = FALSE)
bf <- (exp(unlist(bf_full@bayesFactor[1]))/exp(unlist(bf_ref@bayesFactor[1])))

if (bf<3){
BF_result[voxel] <- bf
} else {
bf_full = lmBF(scores_target ~ lesions_temp + scores_cov1 + scores_cov_lesionsize, data = data4bf, progress = FALSE)
bf_ref =  lmBF(scores_target ~ scores_cov1 + scores_cov_lesionsize, data = data4bf, progress = FALSE)

BF_result[voxel] <- (exp(unlist(bf_full@bayesFactor[1]))/exp(unlist(bf_ref@bayesFactor[1])))
}
}
}
return(BF_result)
}
