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

blk.03 <- read.tables(c(list.files(pattern = "RIL_1.03.")))
blk.03 <- calc.score(blk.03)
blk.04 <- read.tables(c(list.files(pattern = "RIL_1.04.")))
blk.04 <- calc.score(blk.04)
blk.11 <- read.tables(c(list.files(pattern = "RIL_1.11.")))
blk.11 <- calc.score(blk.11)
blk.12 <- read.tables(c(list.files(pattern = "RIL_1.12.")))
blk.12 <- calc.score(blk.12)

merged <- merge(blk.03, blk.04, by = c("pos", "chr"), all = TRUE)
cor(merged$score.x, merged$score.y, use = "pairwise.complete.obs")
# [1] 0.009969097

merged <- merge(blk.03, blk.11, by = c("pos", "chr"), all = TRUE)
cor(merged$score.x, merged$score.y, use = "pairwise.complete.obs")
# [1] 0.007874977

merged <- merge(blk.03, blk.12, by = c("pos", "chr"), all = TRUE)
cor(merged$score.x, merged$score.y, use = "pairwise.complete.obs")
# [1] 0.005195858

merged <- merge(blk.04, blk.11, by = c("pos", "chr"), all = TRUE)
cor(merged$score.x, merged$score.y, use = "pairwise.complete.obs")
# [1] 0.9670474

merged <- merge(blk.04, blk.12, by = c("pos", "chr"), all = TRUE)
cor(merged$score.x, merged$score.y, use = "pairwise.complete.obs")
# [1] 0.9849189

merged <- merge(blk.11, blk.12, by = c("pos", "chr"), all = TRUE)
cor(merged$score.x, merged$score.y, use = "pairwise.complete.obs")
# [1] 0.9681997

merged <- merge(
    blk.11[, c(2, 3, 8)],
    blk.12[, c(2, 3, 8)],
    by = c("pos", "chr"),
    all = TRUE,
    suffix = c("a", "b")
)



library("reshape")


df1 <- data.frame(x=rnorm(1000), y=rnorm(1000), z=factor(letters[2:5]))

ldf <- lapply(seq(1, 100), function(.) df1)

merge.rec <- function(.list, ...){
    if(length(.list)==1) return(.list[[1]])
    Recall(c(list(merge(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
}

system.time(Reduce(function(x, y) merge(x, y, all=T), ldf, accumulate=F))
system.time(merge_all(ldf))
system.time(merge.rec(ldf, all=T))

list.df <- lapply()


all.files <- list.files(pattern = "RIL_.*\\.genotyped\\.nr$")
ril.ids <- unique(gsub(pattern = "^([^.]+)\\..*", replacement = "\\1", all.files))

for (ril in ril.ids) {
    print(ril)
}

ril.files <- list.files(pattern = paste(ril.ids, "\\.genotyped\\.nr$", sep = '')
block.ids <- unique(gsub(pattern = "^[^.]+\\.([^.]+)\\..*", replacement = "\\1", all.files))




# from here

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

# to here looks good



df.list <- vector("list", 3) # create list

for(i in 1:3) {
    df.list[[i]] <- matrix(data = rnorm(3),
                                     nrow = i,
                                     ncol = 3,
                                     byrow = FALSE,
                                     dimnames = NULL)
}

do.call(rbind, df.list) # rbind list elements

