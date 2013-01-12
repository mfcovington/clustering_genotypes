# read files, do calcs, and name columns

setwd("/Users/mfc/git.repos/clustering_genotypes/samples/genotyped/")

# function to read multiple files at once adapted from:
# http://stackoverflow.com/questions/2104483/how-to-read-table-multiple-files-into-a-single-table-in-r/2104532#2104532
read.multi.tables <- function(file.names, ...) {
    require(plyr)
    ldply(file.names, function(fn) data.frame(Filename=fn, read.table(fn, ...)))
}

calc.score <- function(sample.df) {
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


files <- list.files(pattern = "\\.genotyped\\.nr$")
files <- files[file.info(files)$size > 0]
ids <- unique(gsub(pattern = "(.*)\\.[^.]+\\.genotyped.*", replacement = "\\1", files))

merged <- data.frame(chr = NA, pos = NA)

for (id in ids) {
    print(id)
    data <- read.tables(c(list.files(pattern = id)))
    data <- calc.score(data)
    merged <- merge(merged, data[, c(2, 3, 8)], by = c("pos", "chr"), all = TRUE)
    colnames(merged)[ncol(merged)] <- id
}

cor(merged[, 3:ncol(merged)], use = "pairwise.complete.obs")


