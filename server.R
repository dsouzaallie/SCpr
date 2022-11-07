#Created by Alexandria D'Souza and Juerg Straubhaar

options(shiny.maxRequestSize=300*1024^2) 

library(MSnbase)
library(reshape)
library(pheatmap)
library(ggplot2)
library(lattice)
library(limma)
library(gplots)
library(MASS)
library(corrplot)
library(SIMLR)

shinyServer(function(input, output, session) {
  observeEvent(input$runPCA, {
    xJU.df <- read.csv(input$filename$datapath, header=TRUE)
    # annotation
    xJU.df <- JU.df
    # annotation
    annot.df <- xJU.df[, c(5, 8)]; colnames(annot.df) <- c('PepSeq', 'UniProt')
    annot.df$UniProt <- as.character(sapply(annot.df$UniProt, function(x) unlist(strsplit(x, split=';'))[1]))
    annot.df <- annot.df[!is.na(annot.df$UniProt), ]
    write.table(annot.df$UniProt, file='uniprot.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
    accannot.df <- read.table('uniprot.txt', quote='', sep='\t')
    
    # data
    xJU.df <- xJU.df[, c(5, 8, 31:38, 41)]
    xJU.df <- xJU.df[!xJU.df$Contaminant, ]
    xJU.df <- xJU.df[, -2]
    colnames(xJU.df) <- gsub('\\.\\.', '_', gsub('Abundances\\.\\.Grouped\\.\\.\\.', '', colnames(xJU.df)))
    colnames(xJU.df)[1] <- 'PepSeq'
    xJU.df <- xJU.df[1:9]
    
    # lookup
    pep2uniprot <- new.env(hash=TRUE)
    apply(annot.df, 1, function(x) {
      pep2uniprot[[x[1]]] <- x[2]
    })
    uniprot2sym <- new.env(hash=TRUE)
    apply(annot.df, 1, function(x) {
      uniprot2sym[[x[1]]] <- x[2]
    })
    
    xJU.lst <- prepAnnot(xJU.df)
    xJU.pd <- read.csv('pData.txt')
    xJU.mss <- MSnSet(xJU.lst[[1]], xJU.lst[[2]], xJU.pd)
    xJU.mss <- xJU.mss[, grep('126|127', sampleNames(xJU.mss), invert=TRUE)]
    
    ju.mss <- rmAllMiss(xJU.mss)
    
    aju.mss <- plotPercMissing(ju.mss, lessthan=70)
    aju.mss <- rmAllMiss(aju.mss)
    # 1. log transform
    xju.mss <- logTransformWithMissing(aju.mss) # omit missing values when taking log
    # 2. impute
    bju.mss <- impute(xju.mss, method='bpca') # takes a day!!!
    #save(bju.mss, file='bju.mss.rda')
    # 3. outliers
    # 4. normalize
    nbju.mss <- normalise_MSnSet(bju.mss, method='quantiles') # 5838x231; 1410(uniprot)
    pd <- phenoData(nbju.mss)$Treatment
    names(pd) <- sampleNames(nbju.mss)
    nbju.pd <- as.data.frame(pd)
    colnames(nbju.pd) <- 'Treatment'
    phenoData <- new('AnnotatedDataFrame', data=nbju.pd)
    nbju.mses <- ExpressionSet(assayData=exprs(nbju.mss), phenoData=phenoData)
    
    dsgn <- model.matrix(~ 0+factor(pData(nbju.mses)$Treatment))
    colnames(dsgn) <- c('Group1', 'Group2')
    
    fit <- lmFit(nbju.mses, dsgn)
    residuals.m <- residuals.MArrayLM(fit, exprs(nbju.mses))
    e <- exprs(nbju.mses)
    te <- t(e)
    E <- t(residuals.m)
    svdWa <- svd(E)
    W <- svdWa$u[, 1:2]
    alpha <- solve(t(W) %*% W) %*% t(W) %*% te
    cY <- te - W %*% alpha
    cY <- t(cY)
    
    p <- plotPCA_sc_v1(cY, pd)
    png('pca.png')
    plot(p)
    dev.off()
  })
  observeEvent( input$normMatrix, {
    cY.mses <- nbju.mses
    exprs(cY.mses) <- cY
    
    dsgn <- model.matrix(~ 0+factor(pData(cY.mses)$Treatment))
    colnames(dsgn) <- c('Group1', 'Group2')
    
    cYfit <- lmFit(cY.mses, dsgn)
    contrast.matrix <- makeContrasts(Group1-Group2, levels=dsgn)
    cYfit2 <- contrasts.fit(cYfit, contrast.matrix)
    cYfit2 <- eBayes(cYfit2)
    res.tt <- topTable(cYfit2, number=Inf, p.value=0.1, lfc=0.59) 
    
    deAcc <- unique(as.character(unlist(mget(rownames(res.tt), pep2uniprot))))
    write.table(deAcc, file='normalizedMatrix.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
    
  })
  observeEvent(input$createHeat, {
    cy.m <- cY
    rownames(cy.m) <- as.character(unlist(mget(rownames(cy.m), pep2uniprot, ifnotfound=NA)))
    cy.df <- aggregate(cy.m, list(rownames(cy.m)), mean)
    rownames(cy.df) <- cy.df[, 1]
    cy.df <- cy.df[-1]
    cy.df <- cy.df[rownames(cy.df) %in% deAcc, ] # 56x231
    
    ord <- pd
    ord <- sort(ord)
    cy_col <- as.data.frame(ifelse(ord == 1, 'Group2', 'Group1'))
    colnames(cy_col) <- 'Tissue'
    cy.df <- cy.df[, match(names(ord), colnames(cy.df))]
    cy.m <- as.matrix(cy.df)
    
    #Testing From here
    
    htMap <- pheatmap(cy.m,
                      color = inferno(10),
                      kmeans_k=NA,
                      scale = 'row',
                      breaks = NA,
                      clustering_distance_rows = 'euclidean',
                      #clustering_method= 'ward.D',
                      cluster_rows = TRUE,
                      cluster_cols = FALSE,
                      show_rownames = FALSE, 
                      show_colnames = FALSE,
                      annotation_col = cy_col,
                      drop_levels = TRUE,
                      main = 'Single Cell Group 1 vs Group 2'
    )
    
    png('heatmap.png')
    plot(htMap)
    dev.off()
  })
  observeEvent(input$createFrequency, {
    jx <- grep('128|129', colnames(cy.df))
    ux <- grep('130|131', colnames(cy.df))
    xmeancy.df <- as.data.frame(t(apply(cy.df, 1, function(x) {
      jmn <- mean(x[jx])
      umn <- mean(x[ux])
      return(c(jmn, umn))
    })))
    
    meancy.df <- exp(xmeancy.df)
    meancy.df <- namerows(meancy.df, col.name='Prot')
    meltcy.df <- melt(meancy.df)
    meltcy.df$variable <- ifelse(meltcy.df$variable=='V1', 'Group 1', 'Group 2')
    
    pAvg <- ggplot(meltcy.df)
    pAvg <- pAvg + geom_bar(mapping=aes(x=Prot, y=value, fill=variable), position='dodge', stat='identity')
    pAvg + theme(axis.text.x=element_text(angle=90, hjust=1))
    
    png('protAvg.png')
    plot(pAvg)
    dev.off()
  })
  observeEvent(input$createCorrelation, {
    gr <- rep(c(rep('N', 4), rep('A', 4)), 1)
    sn <- sampleNames(nbju.mses)
    names(sn) <- gr
    
    xte <- data.frame(te, Treatment=gr)
    myz <- lda(Treatment ~ ., xte) # warning: variables are collinear
    
    # look at  correlations
    te.cor <- cor(xte[-dim(xte)[2]])
    corrplot.mixed(te.cor)
    
  })
  
  observeEvent(input$createCorrelation, {
    gr <- rep(c(rep('N', 4), rep('A', 4)), 1)
    sn <- sampleNames(nbju.mses)
    names(sn) <- gr
    
    xte <- data.frame(te, Treatment=gr)
    myz <- lda(Treatment ~ ., xte) # warning: variables are collinear
    
    # look at  correlations
    te.cor <- cor(xte[-dim(xte)[2]])
    corrplot.mixed(te.cor)
    
  })
  
  observeEvent(input$SIMLR, {
    
    fx <- factanal(exprs(nbju.mss),factors=2,scores="regression")
    
    set.seed(123410)
    x12346.simlr <- SIMLR(exprs(nbju.mss), 5)
    x12348.simlr <- SIMLR(exprs(nbju.mss), 40)
    x123410.simlr <- SIMLR(exprs(nbju.mss), 4)
    
    pdf('simlr.pdf')
    plotSIMLRclusters(x12348.simlr)
    dev.off()
  })
  

})