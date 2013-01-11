# read files

setwd("/Users/mfc/git.repos/clustering_genotypes/samples/genotyped/")

# function to read multiple files at once adapted from:
# http://stackoverflow.com/questions/2104483/how-to-read-table-multiple-files-into-a-single-table-in-r/2104532#2104532
read.multi.tables <- function(file.names, ...) {
    require(plyr)
    ldply(file.names, function(fn) data.frame(Filename=fn, read.table(fn, ...)))
}

blk.03 <- read.tables(c(list.files(pattern = "RIL_1.03.")))
blk.04 <- read.tables(c(list.files(pattern = "RIL_1.04.")))
blk.12 <- read.tables(c(list.files(pattern = "RIL_1.12.")))


# do calcs and name columns

colnames(blk.03) <- c("filename", "chr", "pos", "par1", "par2", "tot")
head(blk.03)
#                    filename chr     pos par1 par2 tot ratio score
# 1 RIL_1.03.A01.genotyped.nr A01  669437    0    1   1     1     1
# 2 RIL_1.03.A01.genotyped.nr A01 1512215    0    0   0   NaN   NaN
# 3 RIL_1.03.A01.genotyped.nr A01 1843646    0    0   0   NaN   NaN
# 4 RIL_1.03.A01.genotyped.nr A01 1843647    0    0   0   NaN   NaN
# 5 RIL_1.03.A01.genotyped.nr A01 2414553    1    0   1    -1    -1
# 6 RIL_1.03.A01.genotyped.nr A01 2430867    4    0   4    -1    -1

blk.03 <- cbind(blk.03, (blk.03$par2 - blk.03$par1) / blk.03$tot)
colnames(blk.03)[7] <- "ratio"

summary(blk.03)
#                      filename          chr             pos                par1             par2             tot            ratio            score
# RIL_1.03.A03.genotyped.nr:18710   A03    :18710   Min.   :     997   Min.   : 0.000   Min.   : 0.000   Min.   : 0.00   Min.   :-1.000   Min.   :-1.00
# RIL_1.03.A09.genotyped.nr:15458   A09    :15458   1st Qu.: 7310387   1st Qu.: 0.000   1st Qu.: 0.000   1st Qu.: 1.00   1st Qu.:-1.000   1st Qu.:-1.00
# RIL_1.03.A01.genotyped.nr:11445   A01    :11445   Median :14423124   Median : 1.000   Median : 0.000   Median : 2.00   Median :-1.000   Median :-1.00
# RIL_1.03.A06.genotyped.nr:11298   A06    :11298   Mean   :14848865   Mean   : 2.265   Mean   : 1.707   Mean   : 3.98   Mean   :-0.171   Mean   :-0.17
# RIL_1.03.A10.genotyped.nr: 9197   A10    : 9197   3rd Qu.:21443262   3rd Qu.: 2.000   3rd Qu.: 1.000   3rd Qu.: 4.00   3rd Qu.: 1.000   3rd Qu.: 1.00
# RIL_1.03.A07.genotyped.nr: 8029   A07    : 8029   Max.   :37110566   Max.   :83.000   Max.   :89.000   Max.   :89.00   Max.   : 1.000   Max.   : 1.00
# (Other)                  :24304   (Other):24304                                                                        NA's   :10770    NA's   :10770

table(blk.03$ratio == 1)["TRUE"]
#  TRUE
# 35557
table(blk.03$ratio == -1)["TRUE"]
#  TRUE
# 50374
table(abs(blk.03$ratio) < 1)["TRUE"]
# TRUE
# 1740
table(is.na(blk.03$ratio))["TRUE"]
#  TRUE
# 10770
table(abs(blk.03$ratio) > 0.9)["TRUE"]
#  TRUE
# 86210
table(abs(blk.03$ratio) < 0.75)["TRUE"]
# TRUE
# 1019

blk.03$score <- blk.03$ratio
threshold <- 0.75
blk.03$score[blk.03$score      >  threshold] <- 1
blk.03$score[blk.03$score      < -threshold] <- -1
blk.03$score[abs(blk.03$score) <= threshold] <- 0
summary(blk.03$score)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#   -1.00   -1.00   -1.00   -0.17    1.00    1.00   10770


colnames(blk.04) <- c("filename", "chr", "pos", "par1", "par2", "tot")
blk.04 <- cbind(blk.04, (blk.04$par2 - blk.04$par1) / blk.04$tot)
colnames(blk.04)[7] <- "ratio"
blk.04$score <- blk.04$ratio
threshold <- 0.75
blk.04$score[blk.04$score      >  threshold] <- 1
blk.04$score[blk.04$score      < -threshold] <- -1
blk.04$score[abs(blk.04$score) <= threshold] <- 0


colnames(blk.12) <- c("filename", "chr", "pos", "par1", "par2", "tot")
blk.12 <- cbind(blk.12, (blk.12$par2 - blk.12$par1) / blk.12$tot)
colnames(blk.12)[7] <- "ratio"
blk.12$score <- blk.12$ratio
threshold <- 0.75
blk.12$score[blk.12$score      >  threshold] <- 1
blk.12$score[blk.12$score      < -threshold] <- -1
blk.12$score[abs(blk.12$score) <= threshold] <- 0


# merge and calc correlations

dim(blk.03)
# [1] 98441     8
dim(blk.04)
# [1] 82920     8
dim(blk.12)
# [1] 98716     8


dim(merge(blk.03, blk.04, by = "pos"))
[1] 11464    13

merged <- merge(blk.03, blk.04, by = "pos")
cor(merged$score.x, merged$score.y, use = "pairwise.complete.obs")
# [1] 0.01005245

merged <- merge(blk.03, blk.12, by = "pos")
cor(merged$score.x, merged$score.y, use = "pairwise.complete.obs")
# [1] 0.005213639

merged <- merge(blk.04, blk.12, by = "pos")
cor(merged$score.x, merged$score.y, use = "pairwise.complete.obs")
# [1] 0.9819545


