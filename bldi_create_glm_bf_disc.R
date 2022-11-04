bldi_create_glm_bf_disc <- function(dat_scores,idx_elements_included,disconns_full) {


# create blank for the BF for each voxel's lesion-deficit association
BF_result <- matrix(0, length(idx_elements_included))

scores_target <- unlist(dat_scores[2])


# for the pontential inclusion of up to 2 covariates,
# if-conditions pick the right formula 


if (dim(dat_scores)[2]==2){
### loop through all elements
for (element in 1:length(idx_elements_included)){
disconns_temp <- disconns_full[,idx_elements_included[element]]
data4bf <- data.frame(scores_target,disconns_temp)

bf = lmBF(scores_target ~ disconns_temp, data = data4bf, progress = FALSE)
BF_result[element] <- exp(unlist(bf@bayesFactor[1]))
}
}
###
if (dim(dat_scores)[2]==3){

scores_cov1 <- as.numeric(unlist(dat_scores[3]))

### loop through all elements
for (element in 1:length(idx_elements_included)){
disconns_temp <- disconns_full[,idx_elements_included[element]]
data4bf <- data.frame(scores_target,scores_cov1,disconns_temp)

bf_full = lmBF(scores_target ~ disconns_temp + scores_cov1, data = data4bf, progress = FALSE)
bf_ref =  lmBF(scores_target ~ scores_cov1, data = data4bf, progress = FALSE)


BF_result[element] <- (exp(unlist(bf_full@bayesFactor[1]))/exp(unlist(bf_ref@bayesFactor[1])))
}
}
###
if (dim(dat_scores)[2]>=4){

scores_cov1 <- as.numeric(unlist(dat_scores[3]))
scores_cov2 <- as.numeric(unlist(dat_scores[4]))

### loop through all elements
for (element in 1:length(idx_elements_included)){
disconns_temp <- disconns_full[,idx_elements_included[element]]
data4bf <- data.frame(scores_target,scores_cov1,scores_cov2,disconns_temp)

bf_full = lmBF(scores_target ~ disconns_temp + scores_cov1 + scores_cov2, data = data4bf, progress = FALSE)
bf_ref =  lmBF(scores_target ~ scores_cov1 + scores_cov2, data = data4bf, progress = FALSE)


BF_result[element] <- (exp(unlist(bf_full@bayesFactor[1]))/exp(unlist(bf_ref@bayesFactor[1])))
}
}
return(BF_result)
}
