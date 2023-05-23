#Created by Alexandria D'Souza and Juerg.

options(shiny.maxRequestSize=300*1024^2) 
source('SCREEN_Functions.R')
library(MSnbase)
library(reshape)
library(pheatmap)
library(ggplot2)
library(lattice)
library(limma)
library(gplots)
library(MASS)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(SIMLR)
library(dplyr)
library("tidyverse")

shinyServer(function(input, output, session){
  
  observeEvent(input$runimputation, {
    impute.df <- read.csv(input$psmfilename, header=TRUE)
    startAbundance <- as.integer(grep(paste(input$abundancecolumn,"$",sep=''), colnames(impute.df)))-1
    
    system(paste("/Library/Frameworks/Python.framework/Versions/3.10/bin/python3 ImputationMethods/missForest_model.py ", input$psmfilename, " ", input$replicatenum1, " ", startAbundance, " ", input$replicatenum2,  " ", input$psmfilename, wait=FALSE))
    
    
  })
  observeEvent(input$downloadPCA, {
    


    screenFile.df <- read.csv(input$ImputedFile, header=TRUE)
    xJU.df <- screenFile.df
    peptidePCA <- grep("Annotated.Sequence", colnames(screenFile.df))
    proteinPCA <- grep("Master.Protein.Accessions", colnames(screenFile.df))
    idPCA <- grep("File.ID", colnames(screenFile.df))
    # annot.df <- xJU.df[, c(peptidePCA, proteinPCA)]; colnames(annot.df) <- c('PepSeq', 'UniProt')
    # write.table(annot.df$UniProt, file='uniprot.txt', quote=FALSE, sep='\t')
    annot.df <- data.frame(xJU.df$Annotated.Sequence, xJU.df$Master.Protein.Accession)
    annot2.df <- data.frame(paste(xJU.df$Annotated.Sequence, xJU.df$Master.Protein.Accessions), xJU.df$Master.Protein.Accession)
    colnames(annot.df)[colnames(annot.df) == "Annotated.Sequence"] ="PepSeq"
    colnames(annot.df)[colnames(annot.df) == "Master.Protein.Accessions"] ="UniProt"
    write.csv(file='atRef.csv',annot2.df)
    
    # data
    startAbundancePCA <- as.integer(grep(paste(input$abundancecolumn,"$",sep=''), colnames(screenFile.df)))
    endAbundancePCA <- startAbundancePCA + as.integer(input$tmtchannels) -1
    
    #Select columns for Data Frame
    xJU.df <- xJU.df[, c(peptidePCA, proteinPCA, idPCA, startAbundancePCA:endAbundancePCA)]
    xJU_.df=melt(xJU.df, id=c('Annotated.Sequence','Master.Protein.Accessions', 'File.ID'))
    colnames(xJU_.df)[colnames(xJU_.df) == "value"] ="Abundance"
    colnames(xJU_.df)[colnames(xJU_.df) == "variable"] ="ChannelSpecs"
    xJU_.df$Specs.Label <- paste(xJU_.df$File.ID, xJU_.df$ChannelSpecs)
    xJU_.df$Specs.Protein <- paste(xJU_.df$Annotated.Sequence, xJU_.df$Master.Protein.Accessions)
    xJU_.df <- xJU_.df[,!names(xJU_.df) %in% c("ChannelSpecs","File.ID","Annotated.Sequence","Master.Protein.Accessions")]
    xJU_P.df <- reshape(xJU_.df, idvar="Specs.Protein", timevar="Specs.Label", direction="wide")
    write.csv(file='PSM_.csv',xJU_P.df)
    # lookup
    pep2uniprot <- new.env(hash=TRUE)
    apply(annot.df, 1, function(x) {
      pep2uniprot[[x[1]]] <- x[2]
    })
    pep2uniprot_ <- new.env(hash=TRUE)
    apply(annot2.df, 1, function(x) {
      pep2uniprot_[[x[1]]] <- x[2]
    })
    uniprot2sym <- new.env(hash=TRUE)
    apply(annot.df, 1, function(x) {
      uniprot2sym[[x[1]]] <- x[2]
    })
    prepAnnot <- function(df) {
      aggMaxVect <- function(vect) {
        return(ifelse(!all(is.na(vect)), max(vect[!is.na(vect)], na.rm=TRUE), NA))
      }
      
      adf <- aggregate(df[-1], df[1], aggMaxVect)
      adf <- adf[order(adf$Specs.Protein, decreasing=FALSE), ]
      
      fdf <- data.frame(ID=adf$Specs.Protein, Acc=adf$Specs.Protein)
      rownames(fdf) <- fdf$ID
      
      rownames(adf) <- adf$Specs.Protein
      bm <- as.matrix(adf[-1])
      
      bm.cnames <- sapply(colnames(bm), function(x) {
        if (grepl(input$controls, x)) {
          x <- paste(x, '0', sep=',')
        } else if (grepl(input$treatments, x)) {
          x <- paste(x, '1', sep=',')
        }
      })
      write.table(bm.cnames, file='pData.txt', col.names='Treatment', row.names=FALSE, quote=FALSE, sep='\t')
      
      return(list(bm, fdf))
    }
    #Creating Mass Spec Data Frame
    xJU.lst <- prepAnnot(xJU_P.df)
    xJU.pd <- read.csv('pData.txt')
    xJU.mss <- MSnSet(xJU.lst[[1]], xJU.lst[[2]], xJU.pd)
    SIMLR_.mss <- logTransformWithMissing(xJU.mss)
    write.csv(SIMLR_.mss, 'SIMLR_.csv')
    
    ju.mss <- rmAllMiss(xJU.mss)
    
    
    aju.mss <- plotPercMissing(ju.mss, lessthan=70)
    aju.mss <- rmAllMiss(aju.mss)
    xju.mss <- logTransformWithMissing(aju.mss)
    bju.mss <- impute(xju.mss, method='bpca')
    nbju.mss <- normalise_MSnSet(bju.mss, method='quantiles')
    pd <- phenoData(nbju.mss)$Treatment
    names(pd) <- sampleNames(nbju.mss)
    nbju.pd <- as.data.frame(pd)
    colnames(nbju.pd) <- 'Treatment'
    phenoData <- new('AnnotatedDataFrame', data=nbju.pd)
    nbju.mses <- ExpressionSet(assayData=exprs(nbju.mss), phenoData=phenoData)
    dsgn <- model.matrix(~ 0+factor(pData(nbju.mses)$Treatment))
    colnames(dsgn) <- c(input$grp1label, input$grp2label)
    
    #creating the Matrix, Normalized and Bias Corrected
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
    
    plotPCA_sc_v1 <- function(m, pdat) {
      pca <- prcomp(m, scale=TRUE)
      df <- as.data.frame(pca$rotation[, 1:2])
      df <- namerows(df, col.name='Samples')
      
      spl <- df$Samples
      cl <- pdat[match(spl, names(pdat))]
      spl <- ifelse(cl==1, input$grp1label, input$grp2label)
      df$Samples <- spl
      
      p <- ggplot(df, aes(PC1, PC2, colour=Samples)) + geom_point(size=2)
      p <- p + theme(legend.position='right', legend.title=element_blank())
      p <- p + labs(title=input$pcatitle)
      return(p)
    }
    #plot and output pca files
    p <- plotPCA_sc_v1(cY, pd)
    pdf('InputFiles/pca.pdf')
    plot(p)
    dev.off()
  })
  
  observeEvent(input$downloadHeatmap, {
    
    
    screenFile.df <- read.csv(input$ImputedFileHM, header=TRUE)
    xJU.df <- screenFile.df
    peptideHM <- grep("Annotated.Sequence", colnames(screenFile.df))
    proteinHM <- grep("Master.Protein.Accessions", colnames(screenFile.df))
    idHM <- grep("File.ID", colnames(screenFile.df))
    annot.df <- data.frame(xJU.df$Annotated.Sequence, xJU.df$Master.Protein.Accession)
    annot2.df <- data.frame(paste(xJU.df$Annotated.Sequence, xJU.df$Master.Protein.Accessions), xJU.df$Master.Protein.Accession)
    colnames(annot.df)[colnames(annot.df) == "Annotated.Sequence"] ="PepSeq"
    colnames(annot.df)[colnames(annot.df) == "Master.Protein.Accessions"] ="UniProt"
    write.csv(file='atRef.csv',annot2.df)
    
    # data
    startAbundanceHM <- as.integer(grep(paste(input$abundancecolumnHM,"$",sep=''), colnames(screenFile.df)))
    endAbundanceHM <- startAbundanceHM + as.integer(input$tmtchannelsHM) -1
    
    #Select columns for Data Frame
    xJU.df <- xJU.df[, c(peptideHM, proteinHM, idHM, startAbundanceHM:endAbundanceHM)]
    xJU_.df=melt(xJU.df, id=c('Annotated.Sequence','Master.Protein.Accessions', 'File.ID'))
    colnames(xJU_.df)[colnames(xJU_.df) == "value"] ="Abundance"
    colnames(xJU_.df)[colnames(xJU_.df) == "variable"] ="ChannelSpecs"
    xJU_.df$Specs.Label <- paste(xJU_.df$File.ID, xJU_.df$ChannelSpecs)
    xJU_.df$Specs.Protein <- paste(xJU_.df$Annotated.Sequence, xJU_.df$Master.Protein.Accessions)
    xJU_.df <- xJU_.df[,!names(xJU_.df) %in% c("ChannelSpecs","File.ID","Annotated.Sequence","Master.Protein.Accessions")]
    xJU_P.df <- reshape(xJU_.df, idvar="Specs.Protein", timevar="Specs.Label", direction="wide")
    write.csv(file='PSM_.csv',xJU_P.df)
    # lookup
    pep2uniprot <- new.env(hash=TRUE)
    apply(annot.df, 1, function(x) {
      pep2uniprot[[x[1]]] <- x[2]
    })
    pep2uniprot_ <- new.env(hash=TRUE)
    apply(annot2.df, 1, function(x) {
      pep2uniprot_[[x[1]]] <- x[2]
    })
    uniprot2sym <- new.env(hash=TRUE)
    apply(annot.df, 1, function(x) {
      uniprot2sym[[x[1]]] <- x[2]
    })
    prepAnnot <- function(df) {
      aggMaxVect <- function(vect) {
        return(ifelse(!all(is.na(vect)), max(vect[!is.na(vect)], na.rm=TRUE), NA))
      }
      
      adf <- aggregate(df[-1], df[1], aggMaxVect)
      adf <- adf[order(adf$Specs.Protein, decreasing=FALSE), ]
      
      fdf <- data.frame(ID=adf$Specs.Protein, Acc=adf$Specs.Protein)
      rownames(fdf) <- fdf$ID
      
      rownames(adf) <- adf$Specs.Protein
      bm <- as.matrix(adf[-1])
      
      bm.cnames <- sapply(colnames(bm), function(x) {
        if (grepl(input$controlsHM, x)) {
          x <- paste(x, '0', sep=',')
        } else if (grepl(input$treatmentsHM, x)) {
          x <- paste(x, '1', sep=',')
        }
      })
      write.table(bm.cnames, file='pData.txt', col.names='Treatment', row.names=FALSE, quote=FALSE, sep='\t')
      
      return(list(bm, fdf))
    }
    #Creating Mass Spec Data Frame
    xJU.lst <- prepAnnot(xJU_P.df)
    xJU.pd <- read.csv('pData.txt')
    xJU.mss <- MSnSet(xJU.lst[[1]], xJU.lst[[2]], xJU.pd)
    SIMLR_.mss <- logTransformWithMissing(xJU.mss)
    write.csv(SIMLR_.mss, 'SIMLR_.csv')
    
    ju.mss <- rmAllMiss(xJU.mss)
    
    
    aju.mss <- plotPercMissing(ju.mss, lessthan=70)
    aju.mss <- rmAllMiss(aju.mss)
    xju.mss <- logTransformWithMissing(aju.mss)
    bju.mss <- impute(xju.mss, method='bpca')
    nbju.mss <- normalise_MSnSet(bju.mss, method='quantiles')
    pd <- phenoData(nbju.mss)$Treatment
    names(pd) <- sampleNames(nbju.mss)
    nbju.pd <- as.data.frame(pd)
    colnames(nbju.pd) <- 'Treatment'
    phenoData <- new('AnnotatedDataFrame', data=nbju.pd)
    nbju.mses <- ExpressionSet(assayData=exprs(nbju.mss), phenoData=phenoData)
    dsgn <- model.matrix(~ 0+factor(pData(nbju.mses)$Treatment))
    colnames(dsgn) <- c(input$grp1labelHM, input$grp2labelHM)
    
    #creating the Matrix, Normalized and Bias Corrected
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
    
    cY.mses <- nbju.mses
    
    exprs(cY.mses) <- cY
    
    #Create Matrix for Heatmap
    dsgn <- model.matrix(~ 0+factor(pData(cY.mses)$Treatment))
    colnames(dsgn) <- c(input$grp1labelHM, input$grp2labelHM)
    
    #Create the Contrast Matrix with Control and Treatment Design
    cYfit <- lmFit(cY.mses, dsgn)
    contrast.matrix <- makeContrasts(Group1-Group2, levels=dsgn)
    cYfit2 <- contrasts.fit(cYfit, contrast.matrix)
    cYfit2 <- eBayes(cYfit2)
    res.tt <- topTable(cYfit2, number=Inf, p.value=0.1, lfc=0.59)
    
    deAcc <- unique(as.character(unlist(mget(rownames(cYfit2), pep2uniprot_, ifnotfound=NA))))
    write.table(deAcc, file='normalizedMatrix.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
    
    
    cy.m <- cY
    rownames(cy.m) <- as.character(unlist(mget(rownames(cy.m), pep2uniprot_, ifnotfound=NA)))
    cy.df <- aggregate(cy.m, list(rownames(cy.m)), mean)
    rownames(cy.df) <- cy.df[, 1]
    cy.df <- cy.df[-1]
    cy.df <- cy.df[rownames(cy.df) %in% deAcc, ] # 56x231
    
    ord <- pd
    ord <- sort(ord)
    cy_col <- as.data.frame(ifelse(ord == 1, input$grp2labelHM, input$grp1labelHM))
    colnames(cy_col) <- 'Cells'
    cy.df <- cy.df[, match(names(ord), colnames(cy.df))]
    cy.m <- as.matrix(cy.df)
    
    #Create Heatmap
    htMap <- pheatmap(cy.m,
                      color = inferno(10),
                      kmeans_k=NA,
                      scale = 'row',
                      breaks = NA,
                      clustering_distance_rows = 'euclidean',
                      #clustering_method= 'ward.D',
                      cluster_rows = FALSE,
                      cluster_cols = FALSE,
                      show_rownames = FALSE,
                      show_colnames = TRUE,
                      annotation_col = cy_col,
                      drop_levels = TRUE,
                      main = 'Heatmap.'
    )
    
    save_pheatmap_png <- function(x, filename, width=12000, height=10000, res = 300) {
      png(filename, width = width, height = height, res = res)
      grid::grid.newpage()
      grid::grid.draw(x$gtable)
      dev.off()
    }
    try(plot(htMap), silent=TRUE)
    save_pheatmap_png(htMap, 'InputFiles/heatmap.png')
  })
  
  observeEvent(input$downloadSIMLR, {
    #SIMLR

    #Read in the SIMLR Matrix
    SIMLR_.df <- read.table(input$SIMLRFile, quote='', sep=',', header=TRUE)
    SIMLR_transpose.df <- setNames(data.frame(t(SIMLR_.df[,-1])), SIMLR_.df[,1])
    
    #Create Computation from User Identified Values
    set.seed(as.integer(input$seed))
    x12348.simlr <- SIMLR(SIMLR_.df, as.integer(input$c))
    
    pdf('simlr.pdf')
    plotSIMLRclusters(x12348.simlr)
    dev.off()
  })

})