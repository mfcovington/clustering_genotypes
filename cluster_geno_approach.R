# read files

setwd("/Users/mfc/git.repos/clustering_genotypes/samples/genotyped/")

blk.03 <- read.table("RIL_1.03.A03.genotyped.nr")
blk.04 <- read.table("RIL_1.04.A03.genotyped.nr")
blk.12 <- read.table("RIL_1.12.A03.genotyped.nr")


# do calcs and name columns

colnames(blk.03) <- c("chr", "pos", "par1", "par2", "tot")
head(blk.03)
#   chr   pos par1 par2 tot
# 1 A03 16662    1    0   1
# 2 A03 35417    1    0   1
# 3 A03 70578    3    0   3
# 4 A03 93765    2    0   2
# 5 A03 93884    3    0   3
# 6 A03 93886    3    0   3

blk.03 <- cbind(blk.03, (blk.03$par2 - blk.03$par1) / blk.03$tot)
colnames(blk.03)[6] <- "ratio"
#   chr   pos par1 par2 tot ratio
# 1 A03 16662    1    0   1    -1
# 2 A03 35417    1    0   1    -1
# 3 A03 70578    3    0   3    -1
# 4 A03 93765    2    0   2    -1
# 5 A03 93884    3    0   3    -1
# 6 A03 93886    3    0   3    -1

summary(blk.03)
 #  chr             pos                par1            par2             tot             ratio
 # A03:18710   Min.   :   16662   Min.   : 0.00   Min.   : 0.000   Min.   : 0.000   Min.   :-1.0000
 #             1st Qu.: 6106632   1st Qu.: 0.00   1st Qu.: 0.000   1st Qu.: 1.000   1st Qu.:-1.0000
 #             Median :13972974   Median : 1.00   Median : 0.000   Median : 2.000   Median :-1.0000
 #             Mean   :14341319   Mean   : 2.67   Mean   : 1.458   Mean   : 4.136   Mean   :-0.3185
 #             3rd Qu.:21931606   3rd Qu.: 3.00   3rd Qu.: 1.000   3rd Qu.: 5.000   3rd Qu.: 1.0000
 #             Max.   :31665469   Max.   :75.00   Max.   :88.000   Max.   :88.000   Max.   : 1.0000
 #                                                                                  NA's   :2007

table(blk.03$ratio == 1)["TRUE"]
# TRUE
# 5551
table(blk.03$ratio == -1)["TRUE"]
#  TRUE
# 10785
table(abs(blk.03$ratio) < 1)["TRUE"]
# TRUE
#  367
table(is.na(blk.03$ratio))["TRUE"]
# TRUE
# 2007
table(abs(blk.03$ratio) > 0.9)["TRUE"]
#  TRUE
# 16401
table(abs(blk.03$ratio) < 0.75)["TRUE"]
# TRUE
#  205

blk.03$score <- blk.03$ratio
threshold <- 0.75
blk.03$score[blk.03$score      >  threshold] <- 1
blk.03$score[blk.03$score      < -threshold] <- -1
blk.03$score[abs(blk.03$score) <= threshold] <- 0
summary(blk.03$score)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# -1.0000 -1.0000 -1.0000 -0.3167  1.0000  1.0000    2007


colnames(blk.04) <- c("chr", "pos", "par1", "par2", "tot")
blk.04 <- cbind(blk.04, (blk.04$par2 - blk.04$par1) / blk.04$tot)
colnames(blk.04)[6] <- "ratio"
blk.04$score <- blk.04$ratio
threshold <- 0.75
blk.04$score[blk.04$score      >  threshold] <- 1
blk.04$score[blk.04$score      < -threshold] <- -1
blk.04$score[abs(blk.04$score) <= threshold] <- 0


colnames(blk.12) <- c("chr", "pos", "par1", "par2", "tot")
blk.12 <- cbind(blk.12, (blk.12$par2 - blk.12$par1) / blk.12$tot)
colnames(blk.12)[6] <- "ratio"
blk.12$score <- blk.12$ratio
threshold <- 0.75
blk.12$score[blk.12$score      >  threshold] <- 1
blk.12$score[blk.12$score      < -threshold] <- -1
blk.12$score[abs(blk.12$score) <= threshold] <- 0


# merge and calc correlations

dim(blk.03)
# [1] 18710     7
dim(blk.04)
# [1] 15481     7
dim(blk.12)
# [1] 18683     7


dim(merge(blk.03, blk.04, by = "pos"))
[1] 11464    13

merged <- merge(blk.03, blk.04, by = "pos")
cor(merged$score.x, merged$score.y, use = "pairwise.complete.obs")
# [1] -0.1862907

merged <- merge(blk.03, blk.12, by = "pos")
cor(merged$score.x, merged$score.y, use = "pairwise.complete.obs")
# [1] -0.1979788

merged <- merge(blk.04, blk.12, by = "pos")
cor(merged$score.x, merged$score.y, use = "pairwise.complete.obs")
# [1] 0.9827633


