library(edgeR)            #3.34.1
library(ggfortify)        #0.4.13
library(ggplot2)          #3.3.5
library(stats)            #4.1.1
library(factoextra)       #1.0.7
library(ggConvexHull)     #0.1.0
library(DOSE)             #3.18.3
library(clusterProfiler)  #4.0.5
library(org.Hs.eg.db)     #3.13.0
library(BiocParallel)     #1.26.2
library(gridExtra)

# sources----
path <- 'C:/Users/tehil/Dropbox/Projects/APRT/paper/RNAseq/'
expressionPath <- 'inputFiles/expression_matrix_RNAseq_star2206.csv' #0621.csv' 2204
sampleInfoPath <- 'inputFiles/information_table_RNAseq.csv'

#function
createDEGobject <- function(path, expPath, infoPath){
  # take two table- 1. RNA seq with genes on the columns and samples in row 2. information table- age, genotype, feed and sex for each sample
  # create DEG object, order the samples(first - young WT male full)
  # return the DE-object
  
  # read the expression (counts) table and order the names of genes
  exMatrixT <- read.csv(paste0(path, expPath), row.names = 1)
  exMatrix <- data.frame(t(exMatrixT))
  colnames(exMatrix) <- row.names(exMatrixT)

  # read the samples information table, make the parameter as factor and re-factor them(wt, male, full, young as default)
  infoTable <- read.csv(paste0(path, infoPath), row.names = 1, stringsAsFactors=T)
  infoTable$genotype <- relevel(infoTable$genotype, ref = "WT")
  infoTable$sex <- relevel(infoTable$sex, ref = "male")
  infoTable$age <- relevel(infoTable$age, ref = "6.5weeks")
  infoTable$feed <- relevel(infoTable$feed, ref = "full")
  infoTable$sample <- row.names(infoTable)
  
  infoTable <- infoTable[, names(infoTable) %in% c("age", "genotype", "feed", "sex", "sample")]
  infoTable <- infoTable[colnames(exMatrix),] # very important line!! order the info according the RNA-seq data; the first column in count table will be corresponded to the first row in sample table
  
  # create DGE object
  DEobj <- DGEList(exMatrix, samples = infoTable, group=paste(infoTable$age, infoTable$genotype, infoTable$feed, infoTable$sex),
                    genes = row.names(exMatrix))
  
  DEobj <- DEobj[,order(DEobj$samples$genotype)]
  DEobj <- DEobj[,order(DEobj$samples$age)]
  DEobj <- DEobj[,order(DEobj$samples$feed)]
  DEobj <- DEobj[,order(DEobj$samples$sex)]
  DEobj$samples$group <- factor(DEobj$samples$group, levels = unique(DEobj$samples$group))
  return(DEobj)
}
filterNormDEGobject <-function(DEobj){
  # filtering low expressed genes
  # normalized each sample by weight average of non-DE genes using TMM
  
  # filtering out lower expressed genes
  keep <- filterByExpr(DEobj, model.matrix(~genotype+age, data=DEobj$samples))
  print(table(keep))
  DEobj <- DEobj[keep,, keep.lib.sizes=FALSE]
  
  #Normalize data using TMM
  # TMM- 
  DEobj <- calcNormFactors(DEobj, method = 'TMM')
  par(mfrow = c(1, 2))
  boxplot(cpm(DEobj, normalized.lib.sizes = F, log = T), outline=F, col='white', main='Before')
  boxplot(cpm(DEobj, normalized.lib.sizes = T, log = T), outline=F, col='white', main='After')
  return(DEobj)
}
dataExplorer <- function(DEobj, titleG='', logTag = T, clusterSamples = F){
  # dataExplorer accept DGE object
  # return PCA (PC1/2 and PC1/3) with and without labels, and figure of variances against the number of dimensions
  
  #using log2(CPM)
  cpm.data <- log2(cpm(DEobj, normalized.lib.sizes = T)+1)
  cpm.data <- cpm.data[apply(cpm.data, 1, sd) != 0,]
  infoTable <- DEobj$samples[,c('age', 'genotype', 'feed', 'sex', 'sample', 'group')]
  
  #*****************PCA************
  #https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
  pca.de <- prcomp(t(cpm.data), scale. = T)
  fz <-fviz_eig(pca.de)
  # PCA with labels PC1/2 and PC1/3
  pcaL12 <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=2,label = TRUE, title=titleG)+
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group), alpha=0.2, show.legend = F)+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  pcaL13 <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=3,label = TRUE, title=titleG)+
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group), alpha=0.2, show.legend = F)+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  #PCA with dots
  pcaD12 <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=2,label =F,
                     size = 2.5, title=titleG)+#, shape = 19
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group, alpha=infoTable$group),
                    show.legend = T)+
    scale_color_manual(name='Group', values = c('#1D71BB', '#E30613', '#1D71BB', '#E30613'))+
    scale_fill_manual(name='Group', values = c('#1D71BB', '#E30613', '#1D71BB', '#E30613'))+
    scale_alpha_manual(name='Group', values = c(0.2,0.2,0,0))+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  pcaD13 <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=3,label =F,
                     size = 2.5, title=titleG)+#, shape = 19
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group, alpha=infoTable$group),
                    show.legend = T)+
    scale_color_manual(name='Group', values = c('#1D71BB', '#E30613', '#1D71BB', '#E30613'))+
    scale_fill_manual(name='Group', values = c('#1D71BB', '#E30613', '#1D71BB', '#E30613'))+
    scale_alpha_manual(name='Group', values = c(0.2,0.2,0,0))+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  return(list(knee= fz, PCAlabels12 = pcaL12, PCAlabels13 = pcaL13, PCAdots12 = pcaD12, PCAdots13 = pcaD13))
}

# create DGE object and separate to three datasets- male fully fed, male fasted and female fasted.
DEobj <- createDEGobject(path, expressionPath, sampleInfoPath)

DEobj = DEobj[,!colnames(DEobj) %in% '35']
DEobjMale <- filterNormDEGobject(DEobj[,DEobj$samples$sex == 'male' & DEobj$samples$feed == 'fasted'])
DEobjFemale <- filterNormDEGobject(DEobj[,DEobj$samples$sex == 'female'])

dev.off()


maleGraphs <- dataExplorer(DEobjMale, 'male fasted')
femaleGraphs <- dataExplorer(DEobjFemale, 'female fasted')

#save the PCA graphs
pdf(paste0(path, '/figures/PCA RNAseq.pdf'), width=10, height=18) 
grid.arrange(maleGraphs$PCAlabels12, femaleGraphs$PCAlabels12, 
             maleGraphs$PCAlabels13, femaleGraphs$PCAlabels13,
             maleGraphs$PCAdots12,   femaleGraphs$PCAdots12,
             maleGraphs$PCAdots13,   femaleGraphs$PCAdots13,
             maleGraphs$knee,        femaleGraphs$knee, nrow = 5)
dev.off()

# To convert NCBI ids to human entrez ids. There are ways to adapt it for Nfur only, but for now I do everything based on human orthologs
humanKfConversion = read.table(paste0(path, "inputFiles/NCBI-Human-orthologs.txt"), head = T, sep = "\t")

# -------DE------
# Differential expression analysis between ages, between genotypes and in the interaction between them in the three datasets

DErankingSaving <- function(DEobj, design, con, test, pathGL, saveRanking=F){
  DEobj <- estimateDisp(DEobj, design)
  fit <- glmQLFit(DEobj, design, robust = T)
  FC <- data.frame(row.names = row.names(DEobj))
  pval <- data.frame(row.names = row.names(DEobj))
  FDR <- data.frame(row.names = row.names(DEobj))
  for (g in 1:length(con)){
    qlf <- glmQLFTest(fit, contrast = con[[g]])
    top <- topTags(qlf, n=Inf)$table
    print(names(con[g]))
    print(paste('+', dim(top[top$FDR < 0.05 & top$logFC > 0,])))
    print(paste('-', dim(top[top$FDR < 0.05 & top$logFC < 0,])))
    print(head(top[top$FDR < 0.05 & top$logFC < 0,]))
    top$mlog10QvalxFC <- -log10(top$FDR) * top$logFC
    FC[[names(con)[g]]] <- top[rownames(FC),]$logFC
    pval[[names(con)[g]]] <- top[rownames(pval),]$PValue
    FDR[[names(con)[g]]] <- top[rownames(FDR),]$FDR
    
    #saving the ranking list
    ranking_genes <- data.frame(Gene = row.names(top), mlog10QvalxFC = top$mlog10QvalxFC)
    if (saveRanking)
      write.csv(ranking_genes, paste0(pathGL, 'ranking_', names(con[g]), ' ', test[g], '.csv'), row.names = F)
  }
  return(list(fc=FC, fdr=FDR, pval=pval))
}

# "6.5weeks WT fasted female" "6.5weeks Het fasted female" "15weeks WT fasted female" "15weeks Het fasted female" 
# we calculated the differential expression between young vs. old instead of old vs. young, therefor in the downstream analysis we multiple by -1 to get the intuitive direction
con <- list(YOUNG = c(1,-1,0,0),
            OLD = c(0,0,1,-1),
            WT = c(1,0,-1,0),
            HET = c(0,1,0,-1),
            c(1,-1,-1,1))

con <- list(YOUNG = c(-1,1,0,0),
           OLD = c(0,0,-1,1),
           WT = c(-1,0,1,0),
           HET = c(0,-1,0,1),
           c(1,-1,-1,1))

tests <- c('genotype', 'genotype', 'age', 'age', 'age+genotype')

experimentGroups <- list(MALE = DEobjMale,FEMALE = DEobjFemale)
FCs <- list()
pvals <- list()
FDRs <- list()

for (i in 1:length(experimentGroups)){
  y <- experimentGroups[[i]]
  design <- model.matrix(~0+group, data=y$samples)
  colnames(design) <- gsub('group', '', colnames(design))
  print(colnames(design))
  
  contrast_new <- con
  names(contrast_new) <- paste(names(experimentGroups)[i], names(contrast_new), sep = '_')
  
  de <- DErankingSaving(y, design, contrast_new, tests, paste0(path, 'GO/GeneSets/'), F)
  FCs[[names(experimentGroups)[i]]] = de$fc
  pvals[[names(experimentGroups)[i]]] = de$pval
  FDRs[[names(experimentGroups)[i]]] = de$fdr
}

#----GSEA ----
GSEAparam <- function(data, savePath, path){ # thanks for Param Priya Singh
  # This input to GSEA is a ranked list of genes. Read the input ranked list.
  
  # Get human ortholog symbols based on the BLAST results file using org.Hs.eg.db package
  # Some Ids will fail to map and will be ignored
  dataH = merge(humanKfConversion, data, by.x = "ncbi", by.y = "Gene") 
  entrezIds = bitr(as.character(dataH[,2]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") # Get entrez ids for annotation
  dataHE = merge(dataH, entrezIds, by.x = "human", by.y = "SYMBOL") # Get human symbols
  head(dataHE)
  
  # There can be duplicate values because of paralogs, I take average of those for quantitative score
  unique = aggregate(dataHE[,3], list(dataHE$human), mean)
  dataHEU = merge(unique, entrezIds, by.x = "Group.1", by.y = "SYMBOL") 
  colnames(dataHEU) = c("human", "mlog10QvalxFC", "entrez")
  head(dataHEU)
  
  geneList = dataHEU[,2]  # gene list for GO 
  names(geneList) = as.character(dataHEU[,1]) # with entrez ids as names
  
  geneListKegg = geneList # gene list for KEGG
  names(geneListKegg) = as.character(dataHEU[,3]) #  with humna symbols as names
  
  # *** Sort the gene list based on quantitative score in decreasing order. This is critical for GSEA  
  geneList = sort(geneList, decreasing = TRUE)
  geneListKegg = sort(geneListKegg, decreasing = TRUE)
  
  head(geneList)
  tail(geneList)
  
  head(geneListKegg)
  tail(geneListKegg)
  
  # *****************  Now do different enrichment analyses *****************************
  
  # Gene Ontology ------------------------------------------------------------------------------------------------------------------------------------
  ego3 <- gseGO(geneList     = geneList,
                OrgDb        = org.Hs.eg.db,
                keyType      = 'SYMBOL',
                ont          = c("ALL"),
                pvalueCutoff = 1)

  write.table(ego3, paste0(path, 'GO/Results/', savePath, "_GOGSEA.csv"), sep = ",", quote = T, row.names = F)
}

# run GSEA on all the ranking lists from the previous step
files <- list.files(path=paste0(path, 'GO/GeneSets/'), pattern="*.csv")
for (i in 1:length(files)){
  f <- gsub('.csv|ranking_', '', files[i])
  data = read.csv(paste0(path, 'GO/GeneSets/', files[i]))
  ego3 <- GSEAparam(data, f, path)
}

# calculate FC matrix for y/o and wt/het----
cpmMale <- log2(cpm(DEobjMale, normalized.lib.sizes = T)+1)
cpmFemale <- log2(cpm(DEobjFemale, normalized.lib.sizes = T)+1)

FCall <- merge(FCs$MALE, FCs$FEMALE, by=0)
rownames(FCall) = FCall$Row.names
#FCall[-c(1,6,11,16)] = -FCall[-c(1,6,11,16)]  # change the direction for age, genotype but no in the interaction
FCall = FCall[2:length(FCall)]

pvalall <- merge(pvals$MALE, pvals$FEMALE, by=0)
rownames(pvalall) = pvalall$Row.names
pvalall = pvalall[2:length(pvalall)]

FDRall <- merge(FDRs$MALE, FDRs$FEMALE, by=0)
rownames(FDRall) = FDRall$Row.names
FDRall = FDRall[2:length(FDRall)]

infoTableFC <- unique(DEobj$samples[DEobj$samples$feed == 'fasted', c('age', 'genotype', 'feed', 'sex', 'group')])

infoTableFC$age <- factor(infoTableFC$age, levels = c('6.5weeks', '15weeks'))
infoTableFC$genotype <- factor(infoTableFC$genotype, levels = c('WT', 'Het'))
infoTableFC$feed <- factor(infoTableFC$feed, levels = c('full', 'fasted'))
infoTableFC$sex <- factor(infoTableFC$sex, levels = c('male', 'female'))

# fold change between WT and Het
FCwtVShet <- FCall[,grep("YOUNG|OLD", colnames(FCall))]
pvalwtVShet <- pvalall[,grep("YOUNG|OLD", colnames(pvalall))]
FDRwtVShet <- FDRall[,grep("YOUNG|OLD", colnames(FDRall))]
infoTableFCwtVShet <- infoTableFC[infoTableFC$genotype == 'WT',]
infoTableFCwtVShet$genotype <- 'het/wt'
colnames(FCwtVShet) <- 1:4
colnames(pvalwtVShet) <- 1:4
colnames(FDRwtVShet) <- 1:4
rownames(infoTableFCwtVShet) <- 1:4
infoTableFCwtVShet <- infoTableFCwtVShet[order(infoTableFCwtVShet$age),]
FCwtVShet <- FCwtVShet[,rownames(infoTableFCwtVShet)]
pvalwtVShet <- pvalwtVShet[,rownames(infoTableFCwtVShet)]
FDRwtVShet <- FDRwtVShet[,rownames(infoTableFCwtVShet)]

# fold change between young and old
FCyVSo <- FCall[,grep("WT|HET", colnames(FCall))]
pvalyVSo <- pvalall[,grep("WT|HET", colnames(pvalall))]
FDRyVSo <- FDRall[,grep("WT|HET", colnames(FDRall))]
infoTableFCyVSo <- infoTableFC[infoTableFC$age == '6.5weeks',]
infoTableFCyVSo$age <- 'o/y'
colnames(FCyVSo) <- 1:4
colnames(FDRyVSo) <- 1:4
colnames(pvalyVSo) <- 1:4
rownames(infoTableFCyVSo) <- 1:4


# saving the object to the visualization script----
save(DEobj, cpmMale, cpmFemale, 
     FCwtVShet, FDRwtVShet, infoTableFCwtVShet, 
     FCyVSo, FDRyVSo, infoTableFCyVSo, file = paste0(path,'outputFiles/data.RData'))


colnames(FCyVSo) <- infoTableFCyVSo$group
colnames(FDRyVSo) <- infoTableFCyVSo$group
colnames(pvalyVSo) <- infoTableFCyVSo$group

write.csv(FCyVSo, paste0(path, 'outputFiles/FC.csv'))
write.csv(FDRyVSo, paste0(path, 'outputFiles/FDR.csv'))
write.csv(pvalyVSo, paste0(path, 'outputFiles/pval.csv'))
