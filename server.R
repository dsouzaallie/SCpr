#Created by Alexandria D'Souza

options(shiny.maxRequestSize=300*1024^2) 
options(warn=-1)

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

shinyServer(function(input, output, session){
  observeEvent(input$runSCREEN, {
    source('SCREEN_Functions.R')
    

    screenFile.df <- read.csv(input$ImputedFile$datapath, header=TRUE)
    xJU.df <- screenFile.df
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

    xju.mss <- logTransformWithMissing(aju.mss) # omit missing values when taking log

    bju.mss <- impute(xju.mss, method='bpca') # takes a day!!!

    nbju.mss <- normalise_MSnSet(bju.mss, method='quantiles') # 5838x231; 1410(uniprot)

    fx <- factanal(exprs(nbju.mss),factors=2,scores="regression")
    write.table(exprs(nbju.mss), 'SIMLR_.txt')
    SIMLR_.df <- read.table('SIMLR_.txt', quote='', sep='\t', header=TRUE)
    
    set.seed(123410)
    x12348.simlr <- SIMLR(SIMLR_.df, 40)
    pdf('simlr.pdf')
    plotSIMLRclusters(x12348.simlr)
    dev.off()
  })
  
})