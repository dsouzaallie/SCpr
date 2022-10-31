source('SCREEN_Functions.R')

JU.df <- read.csv('072018X_SAM00093_AS_Megakaryocytes_TMT11plexPlate1_1-12_UniprotMouse2_PSM.csv', header=TRUE)
xJU.df <- JU.df
# annotation
annot.df <- xJU.df[, c(5, 8)]; colnames(annot.df) <- c('PepSeq', 'UniProt')
annot.df$UniProt <- as.character(sapply(annot.df$UniProt, function(x) unlist(strsplit(x, split=';'))[1]))
annot.df <- annot.df[!is.na(annot.df$UniProt), ]
write.table(annot.df$UniProt, file='uniprot.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
accannot.df <- read.table('uniprot.txt', quote='', sep='\t')

# data
xJU.df <- xJU.df[, c(5, 8, 31:41, 45)]
xJU.df <- xJU.df[!xJU.df$Contaminant, ]
xJU.df <- xJU.df[, -2]
colnames(xJU.df) <- gsub('\\.\\.', '_', gsub('Abundances\\.\\.Grouped\\.\\.\\.', '', colnames(xJU.df)))
colnames(xJU.df)[1] <- 'PepSeq'

# lookup
pep2uniprot <- new.env(hash=TRUE)
apply(annot.df, 1, function(x) {
  pep2uniprot[[x[1]]] <- x[2]
})
uniprot2sym <- new.env(hash=TRUE)
apply(annot.df, 1, function(x) {
  uniprot2sym[[x[1]]] <- x[2]
})

#Stopped Testing Here
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
colnames(dsgn) <- c('Jurkat', 'U937')

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
pdf('pca.pdf')
plot(p)
dev.off()