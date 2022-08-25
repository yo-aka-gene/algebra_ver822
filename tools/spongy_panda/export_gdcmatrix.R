### export_gdcmatrix.R ###

library(Matrix)
library(jsonlite)

read_as_sparse <- function(matrix, filenames){

    col <- colnames(matrix)
    idx <- rownames(matrix)
    max.len <- max(c(length(idx), length(col)))

    index <- idx[1:max.len]
    columns <- col[1:max.len]
    lst <- as.list(data.frame(index, columns))

    if (length(col) > length(idx)) {
        lst$index <- lst$index[1:length(idx)]
    } else if (length(idx) > length(col)) {
        lst$columns <- lst$columns[1:length(col)]
    }

    writeMM(matrix, file = paste0(filenames, ".mtx"))
    write_json(toJSON(lst, pretty=TRUE), paste0(filenames, ".json"))
    }
