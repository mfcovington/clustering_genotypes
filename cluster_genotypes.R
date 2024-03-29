GenoCor <- function(files, verbose = TRUE) {
  # function to generate a correlation matrix for
  # genotype data from multiple samples
  #
  # requires functions:
  # - ReadMultiTables
  # - CalcScore

  # filter out empty files
  files <- files[file.info(files)$size > 0]

  # extract sample ids from input files
  ids <- unique(gsub(pattern     = "(.*)\\.[^.]+\\.genotyped.*",
                     replacement = "\\1",
                     files))

  # initiate data frame w/ merge-by columns
  merged <- data.frame(chr = NA, pos = NA)

  # read in all files for each sample,
  # calculate genotype scores, and
  # merge resulting scores
  # (columns are named after sample id)
  for (id in ids) {
    if (verbose == TRUE)
      print(id)
    data <- ReadMultiTables(c(list.files(pattern = paste("^", id, "\\.", sep = ''))))
    data <- CalcScore(data)
    merged <- merge(x   = merged,
                    y   = data[, c(2, 3, 8)],
                    by  = c("pos", "chr"),
                    all = TRUE)
    colnames(merged)[ncol(merged)] <- id
  }

  # generate correlation matrix
  if (verbose == TRUE)
    print("Generating correlation matrix")
  genocor.mat <- cor(merged[, 3:ncol(merged)], use = "pairwise.complete.obs")

  return(genocor.mat)
}

ReadMultiTables <- function(file.names, ...) {
  # function to read multiple files at once
  # adapted from: http://stackoverflow.com/questions/2104483/how-to-read-table-multiple-files-into-a-single-table-in-r/2104532#2104532
  #
  # called by function:
  # - GenoCor

  require(plyr)
  ldply(file.names,
        function(fn) data.frame(filename = fn,
                                read.table(fn, ...)))
}

CalcScore <- function(sample.df) {
  # function to calculate genotype scores
  # par1 ==  1
  # het  ==  0
  # par2 == -1
  #
  # called by function:
  # - GenoCor

  colnames(sample.df) <- c("filename", "chr", "pos", "par1", "par2", "tot")
  sample.df <- cbind(sample.df,
                     (sample.df$par2 - sample.df$par1) / sample.df$tot)
  colnames(sample.df)[7] <- "ratio"
  sample.df$score <- sample.df$ratio
  threshold <- 0.75
  sample.df$score[sample.df$score      >  threshold] <- 1
  sample.df$score[sample.df$score      < -threshold] <- -1
  sample.df$score[abs(sample.df$score) <= threshold] <- 0
  return(sample.df)
}

CorPlot <- function(cor.mat,
                    cor.colors = c("red", "white", "blue"),
                    star.color = "green",
                    plot       = TRUE,
                    save       = FALSE,
                    filename   = "corplot.png") {
  # function to plot correlation matrix
  # adapted from: http://theatavism.blogspot.com/2009/05/plotting-correlation-matrix-with.html
  #
  # requires function:
  # - CorPval

  require("ggplot2")
  require("reshape")
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
                geom_text(aes(label = p.sig),
                          color = star.color) +
                theme(axis.text.x = element_text(angle = 270,
                                                 hjust = 0)) +
                scale_fill_gradientn(colours = cor.colors,
                                     limits  = c(-1, 1)) +
                scale_x_discrete(limits = sample.order[1:length(sample.order) - 1]) +
                scale_y_discrete(limits = sample.order[length(sample.order):2]) +
                labs(x = '', y = '')
  if (plot == TRUE)
    return(cor.plot)
  if (save == TRUE)
    ggsave(cor.plot, filename = filename)
}

CorPval <- function(x, alternative="two-sided", ...) {
  # function to calculate p-values from correlation matrix
  # from: http://tolstoy.newcastle.edu.au/R/help/05/04/2659.html
  #
  # called by function:
  # - CorPlot

  corMat <- cor(x, ...)
  n <- nrow(x)
  df <- n - 2
  statistic <- sqrt(df) * corMat / sqrt(1 - corMat^2)
  p <- pt(statistic, df)
  p <- if (alternative == "less")
         p
       else if (alternative == "greater")
         1 - p
       else
         2 * pmin(p, 1 - p)
  return(p)
}

setwd("/Users/mfc/git.repos/clustering_genotypes/samples/genotyped/")
files <- list.files(pattern = "\\.genotyped\\.nr$")
genocor.mat <- GenoCor(files)
CorPlot(genocor.mat, plot = F, save = T, filename = "test.png")





