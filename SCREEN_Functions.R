rmNoQuant <- function(df) {
  # removes all rows from data.set without quantitation
  # df: sample data.frame, eg, x227.df
  ix <- grep('Abundances\\.\\.Grouped', colnames(df))
  xdf <- df[,ix]
  xl <- apply(xdf, 1, function(x) all(is.na(x)))
  df[!xl,]
}

# function from MSnbase library
normalise_MSnSet <- function(object, method, ...) {
  if (method == "vsn") {
    e <- exprs(vsn2(exprs(object), ...))
  } else if (method == "quantiles") {
    e <- preprocessCore::normalize.quantiles(exprs(object), ...)
  } else if (method == "quantiles.robust") {
    e <- preprocessCore::normalize.quantiles.robust(exprs(object), ...)
  } else if (method == "center.mean") {
    e <- exprs(object)
    center <- colMeans(e, na.rm = TRUE)
    e <- sweep(e, 2L, center, check.margin = FALSE, ...)
  } else if (method == "center.median") {
    e <- exprs(object)
    center <- apply(e, 2L, median, na.rm = TRUE)
    e <- sweep(e, 2L, center, check.margin = FALSE, ...)
  } else if (method == "diff.median") {
    e <- exprs(object)
    med <- median(as.numeric(e), na.rm = TRUE)
    cmeds <- apply(e, 2L, median, na.rm = TRUE)
    e <- sweep(e, 2L, cmeds - med)
  } else {
    switch(method,
           max = div <- .rowMaxs(exprs(object), na.rm = TRUE),
           sum = div <- rowSums(exprs(object), na.rm = TRUE))
    e <- exprs(object)/div
  }
  rownames(e) <- rownames(exprs(object))
  colnames(e) <- colnames(exprs(object))
  exprs(object) <- e
  object@processingData@processing <-
    c(object@processingData@processing,
      paste("Normalised (", method ,"): ",
            date(), sep = ""))
  object@processingData@normalised <- TRUE
  if (validObject(object))
    return(object)
}

.plot <- function(x,ttl=NULL) {
  boxplot(exprs(x),
          main=ifelse(is.null(ttl),processingData(x)@processing[2],ttl),
          cex.main=.8,
          cex.lab=0.5,
          cex.axis=.9,
          cex=0.7, las=2)
  grid()
}

plotPCA <- function(m) {
  pca <- prcomp(m, scale=TRUE)
  df <- as.data.frame(pca$rotation[, 1:2])
  df <- namerows(df, col.name='Samples')
  
  ##df$Samples <- rep(c('A', 'B', 'C', 'D', 'E', 'F'), each=3)
  #df$Samples <- gsub('[1-3]', '', df$Samples)
  
  p <- ggplot(df, aes(PC1, PC2, colour=Samples)) + geom_point(size=2)
  p <- p + theme(legend.position='right', legend.title=element_blank())
  return(p)
}


prepAnnot <- function(df) {
  aggMaxVect <- function(vect) {
    return(ifelse(!all(is.na(vect)), max(vect[!is.na(vect)], na.rm=TRUE), NA))
  }
  
  adf <- aggregate(df[-1], df[1], aggMaxVect)
  adf <- adf[order(adf$PepSeq, decreasing=FALSE), ]
  
  fdf <- data.frame(ID=adf$PepSeq, Acc=adf$PepSeq)
  rownames(fdf) <- fdf$ID
  
  rownames(adf) <- adf$PepSeq
  bm <- as.matrix(adf[-1])
  
  bm.cnames <- sapply(colnames(bm), function(x) {
    if (grepl('128|129', x)) {
      x <- paste(x, '0', sep=',')
    } else if (grepl('130|131', x)) {
      x <- paste(x, '1', sep=',')
    }
  })
  write.table(bm.cnames, file='pData.txt', col.names='Treatment', row.names=FALSE, quote=FALSE, sep='\t')
  
  return(list(bm, fdf))
}

rmAllMiss <- function(ms) {
  miss <- apply(exprs(ms), 1, function(x) all(is.na(x)))
  ms <- ms[!miss, ]
  
  return(ms)
}

plotPercMissing <- function(ms, lessthan=100) {
  # lessthan: percent value, less than lessthan missing!
  lessthan <- as.numeric(lessthan)
  percmiss <- apply(exprs(ms), 2, function(x) round((sum(is.na(x))/length(x))*100, 0))
  percmiss.df <- as.data.frame(percmiss)
  percmiss.df <- namerows(percmiss.df, col.name='Sample')
  
  percmiss.df <- percmiss.df[percmiss.df$percmiss <= lessthan, ]
  keepspls <- percmiss.df$Sample
  keepspls <- paste(keepspls, collapse='|')
  
  plot(barchart(percmiss ~ Sample, percmiss.df, horizontal=FALSE, scales=list(rot=90, cex=0.7), ylab='Percent Missing Abundance Values', xlab='Samples'))
  
  return(ms[, grep(keepspls, sampleNames(ms))])
}

hmapMissing <- function(ms) {
  require(gplots)
  show.m <- exprs(ms)
  show.m[is.na(show.m)] <- 0
  show.m[show.m != 0] <- 1
  
  heatmap.2(show.m, col=c("lightgray", "black"), scale="none", dendrogram="none", trace="none", keysize = 0.5, key=FALSE)
}

logTransformWithMissing <- function(ms) {
  e <- exprs(ms)
  e <- apply(e, 2, function(x) ifelse(is.na(x), NA, log(x)))
  exprs(ms) <- e
  
  return(ms)
}

plotPCA_sc_v1 <- function(m, pdat) {
  pca <- prcomp(m, scale=TRUE)
  df <- as.data.frame(pca$rotation[, 1:2])
  df <- namerows(df, col.name='Samples')
  
  spl <- df$Samples
  cl <- pdat[match(spl, names(pdat))]
  spl <- ifelse(cl==1, 'Group1', 'Group2')
  df$Samples <- spl
  
  p <- ggplot(df, aes(PC1, PC2, colour=Samples)) + geom_point(size=2)
  p <- p + theme(legend.position='right', legend.title=element_blank())
  p <- p + labs(title='Group1 vs Group2 Cells')
  return(p)
}

plotPCA_sc<- function(m) {
  pca <- prcomp(m, scale=TRUE)
  df <- as.data.frame(pca$rotation[, 1:2])
  df <- namerows(df, col.name='Samples')
  
  df$Samples <- c(rep('N', 2), rep('A', 5), rep(c(rep('N', 3), rep('A', 5)), 3), rep('N', 3), rep('A', 4), rep(c(rep('N', 3), rep('A', 5)), 3))
  
  p <- ggplot(df, aes(PC1, PC2, colour=Samples)) + geom_point(size=2)
  p <- p + theme(legend.position='right', legend.title=element_blank())
  p <- p + labs(title='Single Cell Group1 vs. Group2')
  return(p)
}

plotSIMLRclusters <- function(obj) {
  # obj: SIMLR object
  col <- ifelse(obj$y$cluster==1, 'red', ifelse(obj$y$cluster==2, 'blue',
                                                ifelse(obj$y$cluster==3, 'green',
                                                       ifelse(obj$y$cluster==4, 'yellow',
                                                              ifelse(obj$y$cluster==5, 'magenta',
                                                                     ifelse(obj$y$cluster==6, 'turquoise',
                                                                            ifelse(obj$y$cluster==7, 'sienna',
                                                                                   ifelse(obj$y$cluster==8, 'olivedrab',
                                                                                          ifelse(obj$y$cluster==9, 'lightsteelblue1', 'black')))))))))
  #'orchid', 'lavender', 'lawngreen', 'sienna'
  #opar <- par(bg='grey80')
  opar <- par(bg='grey90')
  plot(obj$ydata, col=col, xlab = "SIMLR component 1", ylab = "SIMLR component 2", pch=20,  cex=0.7)
  grid(col='blue', nx=12, ny=12)
  par(opar)
}