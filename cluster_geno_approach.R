# read files, do calcs, and name columns

setwd("/Users/mfc/git.repos/clustering_genotypes/samples/genotyped/")

# function to read multiple files at once adapted from:
# http://stackoverflow.com/questions/2104483/how-to-read-table-multiple-files-into-a-single-table-in-r/2104532#2104532
ReadMultiTables <- function(file.names, ...) {
    require(plyr)
    ldply(file.names, function(fn) data.frame(Filename=fn, read.table(fn, ...)))
}

CalcScore <- function(sample.df) {
    colnames(sample.df) <- c("filename", "chr", "pos", "par1", "par2", "tot")
    sample.df <- cbind(sample.df, (sample.df$par2 - sample.df$par1) / sample.df$tot)
    colnames(sample.df)[7] <- "ratio"
    sample.df$score <- sample.df$ratio
    threshold <- 0.75
    sample.df$score[sample.df$score      >  threshold] <- 1
    sample.df$score[sample.df$score      < -threshold] <- -1
    sample.df$score[abs(sample.df$score) <= threshold] <- 0
    return(sample.df)
}

GenoCor <- function(files) {
    files <- files[file.info(files)$size > 0]
    ids <- unique(gsub(pattern = "(.*)\\.[^.]+\\.genotyped.*", replacement = "\\1", files))

    merged <- data.frame(chr = NA, pos = NA)
    for (id in ids) {
        print(id)
        data <- ReadMultiTables(c(list.files(pattern = id)))
        data <- CalcScore(data)
        merged <- merge(merged, data[, c(2, 3, 8)], by = c("pos", "chr"), all = TRUE)
        colnames(merged)[ncol(merged)] <- id
    }

    genocor.mat <- cor(merged[, 3:ncol(merged)], use = "pairwise.complete.obs")
    return(genocor.mat)
}

files <- list.files(pattern = "\\.genotyped\\.nr$")
genocor.mat <- GenoCor(files)
