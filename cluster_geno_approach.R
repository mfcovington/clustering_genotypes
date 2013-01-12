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

    # filter out empty files
    files <- files[file.info(files)$size > 0]

    # extract sample ids from input files
    ids <- unique(gsub(pattern = "(.*)\\.[^.]+\\.genotyped.*", replacement = "\\1", files))

    # initiate data frame w/ merge-by columns
    merged <- data.frame(chr = NA, pos = NA)

    # read in all files for each sample,
    # calculate genotype scores, and
    # merge resulting scores
    # (columns are named after sample id)
    for (id in ids) {
        print(id)
        data <- ReadMultiTables(c(list.files(pattern = id)))
        data <- CalcScore(data)
        merged <- merge(merged, data[, c(2, 3, 8)], by = c("pos", "chr"), all = TRUE)
        colnames(merged)[ncol(merged)] <- id
    }

    # generate correlation matrix
    genocor.mat <- cor(merged[, 3:ncol(merged)], use = "pairwise.complete.obs")

    return(genocor.mat)
}

CorPval <- function(x, alternative="two-sided", ...) {

    # from: http://tolstoy.newcastle.edu.au/R/help/05/04/2659.html
    corMat <- cor(x, ...)
    n <- nrow(x)
    df <- n - 2
    STATISTIC <- sqrt(df) * corMat / sqrt(1 - corMat^2)
    p <- pt(STATISTIC, df)
    p <- if (alternative == "less") {
        p
    }
    else if (alternative == "greater") {
        1 - p
    }
    else 2 * pmin(p, 1 - p)
    p
}

CorPlot <- function(cor.mat) {

    # adapted from: http://theatavism.blogspot.com/2009/05/plotting-correlation-matrix-with.html
    library("ggplot2")
    cor.pval <- CorPval(cor.mat)
    stars <- as.character(
        symnum(
            cor.pval,
            cutpoints = c( 0,     0.001,    0.01,   0.05,  1),
            symbols   = c(   '***',     '**',    '*',    '' ), legend = FALSE))

    cor.df <- cbind(melt(cor.mat), stars)
    names(cor.df) <- c("id.a", "id.b", "cor.val", "p.sig")
    sample.order <- as.character(unique(cor.df$id.b))

    cor.df <- subset(cor.df[upper.tri(cor.mat), ], id.a != id.b)

    cor.plot <- ggplot(cor.df, aes(id.a, id.b, fill = cor.val)) +
      geom_tile() +
      geom_text(aes(label = p.sig), color = "green") +
      theme(axis.text.x = element_text(angle = 270, hjust = 0)) +
      scale_fill_gradientn(
        colours = c("red", "white", "blue"),
        limits  = c(-1, 1)) +
      scale_x_discrete(limits = sample.order[1:length(sample.order) - 1]) +
      scale_y_discrete(limits = sample.order[length(sample.order):2]) +
      labs(x = '', y = '')
}


files <- list.files(pattern = "\\.genotyped\\.nr$")
genocor.mat <- GenoCor(files)
CorPlot(genocor.mat)





