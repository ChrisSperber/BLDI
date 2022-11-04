bldi_load_disc_matrices <- function(dat_scores,path) {

# get first matrix to get the dimensions
disconn_path <- paste(path, "/", dat_scores[1,1], ".csv", sep="")
disconn <- read.table(disconn_path)

element_num <- dim(disconn)[1]*dim(disconn)[2]

# create blank for full vectorised matrices
disconns_full <- matrix(0, dim(dat_scores)[1], element_num)

print("Loading disconnection data")

# store the data of all files
for (x in 1:dim(dat_scores)[1]){
disconn_path <- paste(path, "/", dat_scores[x,1], ".csv", sep="")
disconn <- as.matrix(read.table(disconn_path))


disconns_full[x,] <- c(disconn)
}

return(disconns_full)
}