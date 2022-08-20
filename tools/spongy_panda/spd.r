### spd.R ###

library(Matrix)
library(PTXQC)
library(dplyr)
library(hash)


get_name <- function(path_array) {
    name <- rev(unlist(strsplit(LCSn(path_array))))[1]
    return(name)
    }

form_matrix <- function(path_array){
    barcode.name <- read.delim(
        path_array[1], header = FALSE, stringAsFactors = FALSE
        )
    feat.name <- read.delim(
        path_array[2], header = FALSE, stringAsFactors = FALSE
        )
    mat <- readMM(file = path_array[3])
    colnames(mat)
    ronames(mat)
    return(mat)
}

read_as_sparse <- function(datadir, idx){
    barcodes <- Sys.glob(file.path(datadir, "*barcodes.tsv.gz"))
    features <- Sys.glob(file.path(datadir, "*features.tsv.gz"))
    matrices <- Sys.glob(file.path(datadir, "*matrix.mtx.gz"))
    h <- hash()
    for (i in Map(c, barcodes, features, matrices)){
        h[[get_name(i)]] <- form_matrix(i)
    }
    return(h)
    }