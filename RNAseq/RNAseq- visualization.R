library(ggpubr)           #0.4.0
library(org.Hs.eg.db)     #3.13.0
library(xlsx)             #0.6.5
library(VennDiagram)      #1.7.0
library(ComplexHeatmap)   #2.8.0
library(tidyr)            #1.1.4
library(reshape2)         #1.4.4
library(GSVA)             #1.40.1
library(KEGGREST)         #1.32.0
library(rjson)            #0.2.20
library(curl)             #1.2.0
library(clusterProfiler)  #4.0.5
library(org.Hs.eg.db)     #3.13.0

path <- 'C:/Users/tehil/Dropbox/paper/RNAseq/'
load(paste0(path,'outputFiles/data.RData'))
# To convert NCBI ids to human entrez ids. There are ways to adapt it for nfur only, but for now I do everything based on human orthologs
humanKfConversion = read.table(paste0(path, "inputFiles/NCBI-Human-orthologs.txt"), head = T, sep = "\t")

# APRT expression----
barplotExpression <- function(gene){
  aprt <- data.frame(Exp = c(cpmMaleFull[gene,], cpmMaleFast[gene,], cpmFemale[gene,]), 
                     sample=DEobj$samples$sample, group=DEobj$samples$group, genotype=DEobj$samples$genotype, age=DEobj$samples$age)
  aprt$newExp <- 2^aprt$Exp -1 # CPM without log2
  
  #normalize each Het + wt by correspond average wt.
  averageGroups <- aggregate(aprt$newExp, list(aprt$group), mean)
  averageGroups$x[seq(from=2, to=12, length=6)] <- averageGroups$x[seq(from=2, to=12, length=6)-1] 
  aprt = merge(aprt, averageGroups, by.x="group", by.y="Group.1", sort=F)
  aprt$normExp = aprt$newExp / aprt$x
  
  # we combined young and old together, and did Wilcox test between WT and Het of each group.
  # gruopsWtest- the group within Wilcox test. Wtest- the main group. e.g Wtest- male fully fed, groupWtest- WT male fully-fed, Het male fully-fed
  aprt$groupWtest <- gsub('15weeks |6.5weeks ', '', aprt$group)
  aprt$groupWtest <- factor(aprt$groupWtest, levels = unique(aprt$groupWtest))
  aprt$Wtest <- gsub('WT |Het ', '', aprt$groupWtest)
  aprt$Wtest <- factor(aprt$Wtest, levels = unique(aprt$Wtest))
  
  pvalsdf <- data.frame(matrix(0, ncol=5, nrow = 3))
  colnames(pvalsdf) <- c('group1', 'group2', 'p', 'p.adj', 'y.position')
  for (i in 1:length(unique(aprt$Wtest))){
    pvalsdf[i, c('group1', 'group2')] = unique(aprt$groupWtest)[c(2*i-1,2*i)]
    c1 <- aprt[aprt$Wtest == unique(aprt$Wtest)[i],]
    pvalsdf[i, 'p'] <- wilcox.test(normExp ~ genotype, data=c1)$p.value
  }
  pvalsdf[, 'p.adj'] <- round(p.adjust(pvalsdf$p, method='fdr'), 4) 

  pvalsdf$y.position <- 1.2
  
  #
  aprtPlot <- ggplot(aprt, aes(x=groupWtest, y=normExp)) +
    geom_bar(stat='Summary', fun = 'mean', color="black", fill='white', width=0.75, size=1) +
    geom_errorbar(stat = "Summary", fun.data = mean_se, width=0.3) +
    geom_point(position = position_jitter(width=0.15), aes(color=genotype), size=2.5) +
    scale_color_manual(values = c('black', '#008000')) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    stat_pvalue_manual(pvalsdf, label = 'p.adj', size = 4) +
    ggtitle(toupper(gene)) +
    xlab("") +
    scale_y_continuous(limits = c(0,1.4), expand = c(0, 0)) +
    geom_hline(yintercept=0.5, linetype="dashed")+
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.y = element_text(color="black"), axis.text.x = element_text(color="black"))

  return(aprtPlot)
}


pdf(paste0(path, '/figures/APRT expression.pdf'), width = 5.5, height = 3)
plot(barplotExpression('aprt'))
dev.off()

# functions----
averageSamples <- function(cpm.data, infoTable){
  mean.express <- aggregate(t(cpm.data), list(infoTable$group), mean)
  mean.express <- separate(mean.express, Group.1, c('age', 'genotype', 'feed', 'sex'), sep=" ")
  infoTable_x <- mean.express[,colnames(mean.express) %in% c('age', 'genotype', 'feed', 'sex')]
  infoTable_x$group <- paste(infoTable_x$age, infoTable_x$genotype, infoTable_x$feed, infoTable_x$sex ,sep=' ')
  infoTable_x$age <- factor(infoTable_x$age, levels = c('6.5weeks', '15weeks'))
  infoTable_x$genotype <- factor(infoTable_x$genotype, levels = c('WT', 'Het'))
  infoTable_x$feed <- factor(infoTable_x$feed, levels = c('full','fasted'))
  infoTable_x$sex <- factor(infoTable_x$sex, levels = c('male', 'female'))
  
  infoTable_x <- infoTable_x[order(infoTable_x$genotype),]
  infoTable_x <- infoTable_x[order(infoTable_x$age),]
  infoTable_x <- infoTable_x[order(infoTable_x$feed),]
  infoTable_x <- infoTable_x[order(infoTable_x$sex),]
  
  mean.express <- mean.express[rownames(infoTable_x),]
  cpm.data.core <- t(mean.express[,!colnames(mean.express) %in% c('age', 'genotype', 'feed', 'sex')])
  
  return(list(cpm.mean = cpm.data.core, infoTable = infoTable_x))
}
averageSamplesDEobj <- function(DEobj){
  cpm.data <- log2(cpm(DEobj, normalized.lib.sizes = T, log = F)+1)
  return(averageSamples(cpm.data, DEobj$samples))
}
convertLoc2symbol <- function(oldNames){
  # convert symbol to the final symbol.
  NCBIgeneNames <- read.csv(paste0(path,'inputFiles/genes_names.csv'))
  names(NCBIgeneNames) <- c('NCBI', 'FinalSymbol', 'Human')
  head(NCBIgeneNames)
  newNames <- NCBIgeneNames[NCBIgeneNames$NCBI %in% oldNames,]
  newNames <- newNames[match(oldNames, newNames$NCBI),]$FinalSymbol
  return(newNames)
}

# combined GSEA result files ----
ont <- 'BP'
DBs <- 'GO'

extractSigGenes <- function(path, group, ngroup, test, ont= '', pvalCutoff=0.05) {
  # extract  the list of significant enrichment pathways and table enrichment pathways from certain database, test, ontology option (BP, CC or MF) and condition.
  #  there is option to choose p-value cutoff
  
  file_name <-  paste0(path, 'GO/Results/', group[ngroup], ' ', test, '_GOGSEA.csv')
  print(file_name)
  ExTable <- read.csv(file_name)
  genes <- ExTable[ExTable$qvalues < pvalCutoff & ExTable$ONTOLOGY == ont ,]$ID 
  return(list(table=ExTable, genes=genes))
}

threeGroupsCompare <- function(path, group, test, DB, nameG, ont= '', pvalCutoff=0.05){
  # create table with all the significant enriched pathway in at least one of the given comparisons,
  # drawing Venn diagram with the number of significant pathways in each group (less than 5 groups, if more- without venn diagram)
  # 
  # return table with the go ID, description and for each group:
  # 1. NES (normalized enrichment score) 2.adjust p.value 3.column that contain the core genes in each pathway
  
  tables <- list()
  genes <- list()
  for (i in 1:length(group)){
    GG <- extractSigGenes(path, group, i, test, ont)
    tables[[i]] <- GG$table
    genes[[i]] <- GG$genes
  }
  
  names(genes) <- names(group)
  #https://www.rdocumentation.org/packages/VennDiagram/versions/1.6.20/topics/venn.diagram
  if (length(genes) < 6){
    v <- venn.diagram(genes, category.names = names(group), main= paste(DB, test, ont),
                      filename = NULL)
    ggsave(v, file=paste0(path, 'GO/Results/MergedAnalysis/', test, '_in_', nameG, '_venn.pdf'), device = "pdf", width=5.5, height = 6)
  }
  # one column contain which group have significant enrichment in the pathway
  vennGroups <- get.venn.partitions(genes)
  vennGroups$..set.. <- gsub("\\).*", '', vennGroups$..set..)
  vennGroups$..set.. <- gsub("\\(", '', vennGroups$..set..)
  vennGroups$..set.. <- gsub("n", '#', vennGroups$..set..)
  vennGroups$..set.. <- as.character(lapply(vennGroups$..set.., as.name))
  vennGroups$..set.. <- gsub("n", '+', vennGroups$..set..)
  vennGroups$..set.. <- gsub("#", 'n', vennGroups$..set..)
  
  tableGroups <- data.frame(group=rep('NAN', length(unique(unlist(genes)))), ID=unique(unlist(genes)))
  
  for (gr in 1:nrow(vennGroups)){
    if (length(vennGroups$..values..[[gr]]))
      tableGroups[tableGroups$ID %in% vennGroups$..values..[[gr]],]$group <- vennGroups$..set..[gr]
  }
  print(table(tableGroups$group))
  
  #col names for overlapping tables
  colNamesGroup <- c('ID', 'NES', 'qvalues', 'core_enrichment')
  N <- length(colNamesGroup) - 1
  
  #col names for describe table
  colNamesGeneral <- c('ONTOLOGY', 'ID', 'Description')
  
  
  tablePvalNES <- merge(tables[[1]][colNamesGroup], tables[[2]][colNamesGroup], all=T, by='ID', suffixes=paste0('.', names(group)[1:2]))
  tablePathwaysNames <- merge(tables[[1]], tables[[2]], by=colNamesGeneral,all=T, suffixes=paste0('.', names(group)[1:2]))
  
  if (length(group) > 2){ # more than 2 groups in the comparisons
    for (k in 3:length(group)){
      tablePvalNES <- merge(tablePvalNES, tables[[k]][colNamesGroup], all=T, by='ID')
      
      colnames(tablePvalNES) <- c(colnames(tablePvalNES)[1:(1+(k-1)*N)], paste0(colnames(tablePvalNES)[(2+(k-1)*N):(1+k*N)], '.', names(group[k])))
      tablePathwaysNames <- merge(tablePathwaysNames, tables[[k]], by=colNamesGeneral, all=T)
    }
  }
  
  # merge the column with the groups are significant in each pathways with the pathway id and description
  tableDescriptionValues <- merge(tableGroups, tablePathwaysNames[colNamesGeneral], by='ID', all.x=T)
  
  #merge everything together
  merge(tableDescriptionValues, tablePvalNES, by='ID', all.x=T)
  
}

sendToCompareSaveXlsx <- function(path, group, DBs, nameG, test, ont=''){
  xlsxFileName <- paste0(path, 'GO/Results/MergedAnalysis/', DBs, ont, '_compare_', test, ' in ', nameG, '.xlsx')
  print(xlsxFileName)
  
  write.xlsx(t(data.frame(DB=DBs, ontology=ont, age='positive- decrease with age', genotype='positive-more in WT')),
             file = xlsxFileName, sheetName='summary', append=FALSE)
  
  comparisonsTable <- threeGroupsCompare(path, group, test, DBs, nameG, 'BP')
  if(nrow(comparisonsTable) > 0)
    write.xlsx(comparisonsTable, file = xlsxFileName, sheetName=test, append=TRUE)
  
}

groups <- list(ag = setNames(paste0(c('MALE_FULL', 'MALE_FASTED', 'FEMALE'), '_'), c('MALE_FULL', 'MALE_FASTED', 'FEMALE')),
               age = c('MALE_FULL_WT', 'MALE_FULL_HET', 'MALE_FASTED_WT', 'MALE_FASTED_HET', 'FEMALE_WT', 'FEMALE_HET'),
               genotype = c('MALE_FULL_YOUNG', 'MALE_FULL_OLD', 'MALE_FASTED_YOUNG', 'MALE_FASTED_OLD', 'FEMALE_YOUNG', 'FEMALE_OLD'))
tests <- c('age+genotype', 'age', 'genotype')
sendToCompareSaveXlsx(path, groups$ag, DBs, 'all', 'age+genotype', ont)

for (i in 2:length(groups)){
  # all
  group <- setNames(groups[[i]], groups[[i]])
  sendToCompareSaveXlsx(path, group, DBs, 'all', tests[i], ont)
  
  # males
  group <- setNames(groups[[i]][1:4], gsub('MALE_', '', groups[[i]][1:4]))
  sendToCompareSaveXlsx(path, group, DBs, 'males', tests[i], ont)
  
  # fasted
  group <- setNames(groups[[i]][3:6], gsub('FASTED_', '', groups[[i]][3:6]))
  sendToCompareSaveXlsx(path, group, DBs, 'fasted', tests[i], ont)
}

# after GSEA analysis- draw heatmap on specific core genes in interesting pathway(s)----

# GO to the gene list; the gene list for each GO pathway
goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
go2gene <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys=names(goterms), column='SYMBOL', keytype="GOALL", multiVals='list'))
goAnno <- stack(go2gene) 
colnames(goAnno) <- c('SYMBOL', "GOID")
goAnno <- unique(goAnno[!is.na(goAnno[,1]), ])
goAnno$ONTOLOGYALL <- goterms[goAnno$GOID]


actualHeatmap <- function(cpm.data, infoData, geneList, goID, average=F, convertNames = T, splitCol = 'none', clusterRows=T, clusterSamples = F, geneOrder=NA, scaleT = T,
                          filterGenes=F, fdr=NA, fdrMfc=NA, cutoff=0.05, hirerByAnoutherDS = F, Adataset = c(), infoDS = 'none', ...){
  # cpm.data-   dataframe with expression;    
  # infoData-   sample information: age, genotype, sex, feeding condition
  # geneList-   the genes to plot in the heatmap
  # goID-       name of the GO pathway
  # average-    use the average of each group. default False
  # splitCol-   split the columns by category. could be: 'none'(default), 'age', 'genotype', 'feed', 'sex' and any new col-name you create
  # clusterRows, clusterSamples- cluster the rows or the columns by hierarchical clustering. default- T for the rows, F for the columns
  # geneOrder-  ordering the genes by a list of genes
  # scaleT-     scale each gene by row. Boolean value; default True
  
  # filterGenes-filter the genes by FDR, FC or FDR*log(FC); Boolean value, default False
  # if filterGenes is True:
  # fdr or fdrMfc - should contain dataframe of the FDR/ FDR*log(FC) values. if not will use the FC dataframe
  # cutoff-     the threshold for the filterGenes; default- 0.05
  # 
  # hirerByAnoutherDS- use anouther dataset to order the genes by hierarchical clustering. Boolean value, default False
  # if True- Adataset should contain the dataframe and infoDS the information about the samples in the new dataset.
  # ... others parameters for the Heatmap function
  
  # took the genes are exits in the CMP data
  geneList <- geneList[geneList %in% rownames(cpm.data)]
  
  if (hirerByAnoutherDS){
    geneList <- geneList[geneList %in% rownames(Adataset)  ]
  }
  cpm.data <- cpm.data[geneList,]
  
  if (average){
    av <- averageSamples(cpm.data, infoData)
    cpm.data <- av$cpm.mean
    infoData <- av$infoTable
  }
  
  if (filterGenes){ # at least one of the column are passed the cutoff
    if (!is.na(fdr)){
      fdr <- fdr[geneList,]
      cpm.data <- cpm.data[apply(fdr,1,function(x) sum(x < cutoff) != 0),]
    }
    else if (!is.na(fdrMfc)){
      fdrMfc <- fdrMfc[geneList,]
      cpm.data <- cpm.data[apply(fdrMfc,1,function(x) sum(abs(x) > cutoff) != 0),]
    }
    else{
      cpm.data <- cpm.data[apply(cpm.data,1,function(x) sum(abs(x) > cutoff) != 0),]
    }
  }
  
  
  if (nrow(cpm.data) < 2) # if the matrix too small- draw a random heatmap 
    return(Heatmap(matrix = c(1,2,3,4, name='ERROR!!')))
  
  if (!clusterRows){
    #head(cpm.data)
    geneOrder <- geneOrder[geneOrder %in% row.names(cpm.data)]
    cpm.data <- cpm.data[geneOrder,]
    hr = F
  }
  else
    hr <- hclust(as.dist(1-cor(t(cpm.data), method="pearson"))) # Cluster rows by Pearson correlation.
  #hc <- hclust(as.dist(1-cor(cpm.data, method="pearson"))) # Clusters columns by Pearson correlation.

  if (scaleT)
    data.scaled <- t(scale(t(cpm.data))) 
  else{
    data.scaled <- as.matrix(cpm.data)
  }
  
  # used hierarchical clustering from another dataset
  if(hirerByAnoutherDS){
    geneList <- geneList[geneList %in% rownames(Adataset)]
    geneList <- geneList[geneList %in% row.names(cpm.data)]
    Adataset <- Adataset[geneList,]
    
    if (average){
      av <- averageSamples(Adataset, infoDS)
      Adataset <- av$cpm.mean
    }
    hr <- hclust(as.dist(1-cor(t(Adataset), method="pearson"))) # Cluster rows by Pearson correlation.
  }
  
  if (convertNames)
    rownames(data.scaled) <- convertLoc2symbol(rownames(data.scaled))
  c.annotaion <- HeatmapAnnotation(age=infoData$age, genotype=infoData$genotype, feed=infoData$feed, sex=infoData$sex,
                                   col = list(age=c('6.5weeks'='black', 'o/y'='gray60', '15weeks'='gray'),
                                              genotype=c('WT'= 'blue','het/wt'='darkorchid2', 'Het'='red'),
                                              feed=c('fasted'= 'lightgreen', 'full'= 'forestgreen'),
                                              sex=c('male'= 'darkblue', 'female'='pink')),
                                   simple_anno_size = unit(2, "mm"))
  if (splitCol == 'none')
    infoData['none'] = rep("", nrow(infoData))
  

  hm <- Heatmap(data.scaled, name = goID, top_annotation = c.annotaion, cluster_columns = clusterSamples,
              column_split  = infoData[splitCol], cluster_rows = hr, row_title_rot =0,
              row_title_gp = gpar(fontsize =8), row_names_gp = gpar(fontsize = 5), ...) 
  return(hm)
}
GOdots <- function(GoCompareEnrichmentGroup, chosenPathways, group, title="GO", labelss=c(), cirles=c()){
  # extract the pathways and the group 
  GoCompareEnrichmentGroup <- GoCompareEnrichmentGroup[GoCompareEnrichmentGroup$ID %in% chosenPathways$ID, -grep('core_', colnames(GoCompareEnrichmentGroup))]
  cnames <- c(colnames(GoCompareEnrichmentGroup)[1:5], colnames(GoCompareEnrichmentGroup)[grep(paste(group, collapse = '|'), colnames(GoCompareEnrichmentGroup))])
  GoCompareEnrichmentGroup <- GoCompareEnrichmentGroup[,cnames]
  
  #order the pathways according FDR or the original order
  pathwaysOrder <- chosenPathways$Description
  #
  
  #select the pathways with significant value at least in one group
  if (length(grep('qvalues', colnames(GoCompareEnrichmentGroup))) > 1){
    checkTable <- GoCompareEnrichmentGroup[, grep('qvalues', colnames(GoCompareEnrichmentGroup))]
    row.names(checkTable)
    checkTable[is.na(checkTable)] = 1
    GoCompareEnrichmentGroup <- GoCompareEnrichmentGroup[rowSums(checkTable < 0.05) > 0,]
    #colnames(GoCompareEnrichmentGroup) <- gsub('p.adjust', 'padjust', colnames(GoCompareEnrichmentGroup))
  }
  
  
  
  meltEnrich <- melt(GoCompareEnrichmentGroup, id.vars = c('ID', 'Description'))
  meltEnrich <- meltEnrich[!meltEnrich$variable %in% c('ONTOLOGY', 'NA.'),]
  meltEnrich <- separate(data = meltEnrich, col = variable, into = c("measure", "Group"), sep = "\\.")
  meltEnrich$Group <- factor(meltEnrich$Group, levels = unique(meltEnrich$Group))
  meltEnrich$value <- as.numeric(meltEnrich$value)
  enrichment_score <- meltEnrich[meltEnrich$measure == 'NES', c('ID', 'Group', 'Description', 'value')]
  pval <- meltEnrich[meltEnrich$measure == 'qvalues', c('ID', 'Group', 'Description', 'value')]
  colnames(enrichment_score) <- gsub('value', 'enrichment_score', colnames(enrichment_score))
  colnames(pval) <- gsub('value', 'FDR', colnames(pval))
  goVisDF <- merge(pval, enrichment_score, by=c('ID', 'Group', 'Description'))
  
  goVisDF[is.na(goVisDF)] = 0
  
  
  #goVisDF$p_adjust[goVisDF$p_adjust > 0.05] = 0
  goVisDF$sig <- c(goVisDF$FDR < 0.05)
  goVisDF$sig <- replace(goVisDF$sig, goVisDF$sig==T, 'Sig')
  goVisDF$sig <- replace(goVisDF$sig, goVisDF$sig==F, 'NS')
  goVisDF$sig <- factor(goVisDF$sig, levels = c('Sig', 'NS'))
  
  goVisDF$Description <- paste0(toupper(substr(goVisDF$Description, 1, 1)), substr(goVisDF$Description, 2, nchar(goVisDF$Description)))
  pathwaysOrder <- paste0(toupper(substr(pathwaysOrder, 1, 1)), substr(pathwaysOrder, 2, nchar(pathwaysOrder)))
  goVisDF$Description <- factor(goVisDF$Description, levels = rev(pathwaysOrder))
  
  one <- ggplot(data = goVisDF, aes(x = '1', y = Description, color= enrichment_score)) + 
    geom_point(aes(shape = sig, size = -log10(FDR))) +  
    scale_color_gradient(low = "blue", high = "red") +
    scale_shape_manual(values = c(19, 21), breaks = waiver())+
    theme_bw() +
    ylab("") + 
    xlab("") + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    facet_grid(.~Group) + 
    ggtitle(title)
  
  if (length(labelss) == 0 | length(cirles) == 0)
    return(one)
  
  two <- ggplot(data = goVisDF, aes(x = '1', y = Description, color= enrichment_score, shape = sig, size = -log10(FDR))) + 
    geom_point() +  
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous(name='-log10(FDR)', breaks=labelss, labels=labelss)+
    scale_shape_manual(values = c(19, 21)) +
    guides(size = guide_legend(override.aes = list(shape = cirles)))+
    theme_bw() +
    ylab("") + 
    xlab("") + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    facet_grid(.~Group) + 
    ggtitle(title)
  
  return(two)
}

GSVAheatmap <- function(cpmData, infoTable, geneSets){
  dataGSVA <- gsva(cpmData, geneSets)
  hm <- actualHeatmap(dataGSVA, infoTable, names(geneSets), 'GSVA',
                      average = T, convertNames = F, scaleT = F, clusterRows = F, geneOrder = names(geneSets))
  return(hm)
}
# read the data for dots plots----
GoCompareEnrichmentGroupAge <- read.xlsx(paste0(path, '/GO/Results/MergedAnalysis/GOBP_compare_age in all.xlsx'), sheetName = 'age')
GoCompareEnrichmentGroupAge[,grep('NES', colnames(GoCompareEnrichmentGroupAge))] <- -GoCompareEnrichmentGroupAge[,grep('NES', colnames(GoCompareEnrichmentGroupAge))]

GoCompareEnrichmentGroupGeno <- read.xlsx(paste0(path, '/GO/Results/MergedAnalysis/GOBP_compare_genotype in all.xlsx'), sheetName = 'genotype')
GoCompareEnrichmentGroupGeno[,grep('NES', colnames(GoCompareEnrichmentGroupGeno))] <- -GoCompareEnrichmentGroupGeno[,grep('NES', colnames(GoCompareEnrichmentGroupGeno))]

# age and genotype
GoCompareEnrichmentGroup <- read.xlsx(paste0(path, '/GO/Results/MergedAnalysis/GOBP_compare_age+genotype in all.xlsx'), sheetName = 'age+genotype')

# figure 2b----
chosenPathways <- read.csv(paste0(path, 'inputFiles/pathways to age in males.csv'), row.names=1)
godots <- GOdots(GoCompareEnrichmentGroupAge[,1:8], chosenPathways, c('FULL', 'FASTED'), "GO age in males WT fully-fed")
pdf(paste0(path, 'figures/GO age in males WT fully-fed.pdf'))
plot(godots)
dev.off()
# full vs fasted figure 2D (age + genotype)
chosenPathways <- read.csv(paste0(path, 'inputFiles/pathways to age+genotype in males.csv'), row.names=1)
godots <- GOdots(GoCompareEnrichmentGroup, chosenPathways, c('FULL', 'FASTED'), "GO age+genotype in males", c(1,2,3,4,5,6,7), c(21,19,19,19,19,19,19))
pdf(paste0(path, 'figures/GO age+genotype in males.pdf'))
plot(godots)
dev.off()

# inflammation- Figure S2c ----
goID <-  c('GO:0045088')
names(goID) <-  c('regulation of innate immune response')

for (p in 1:length(goID)){
  pdf(paste0(path, 'figures/Heatmap ', names(goID[p]), '2.pdf'))
  
  geneSetHuman <- unique(goAnno[goAnno$GOID == goID[p],]$SYMBOL)
  geneSetKf <-c(humanKfConversion[humanKfConversion$human %in% geneSetHuman,]$ncbi)
 
  # male full average and sample
  hm <- actualHeatmap(cpmMaleFull, DEobj$samples[DEobj$samples$sex == 'male' & DEobj$samples$feed == 'full',], geneSetKf, names(goID[p])) 
  draw(hm, column_title = "male fully fed samples")
   
  # male fasted average and sample
  hm <- actualHeatmap(cpmMaleFast, DEobj$samples[DEobj$samples$sex == 'male' & DEobj$samples$feed == 'fasted',], geneSetKf, names(goID[p])) 
  draw(hm, column_title = "male fasted samples")

  dev.off()
}

# comparison between ages in fasted (males and females, WT and Het)- Figure 4b and 4c----
chosenPathways <- read.csv(paste0(path, 'inputFiles/pathways age in fasted.csv'), row.names=1)
goID <- chosenPathways$ID
names(goID) <- chosenPathways$Description

godots <- GOdots(GoCompareEnrichmentGroupAge, chosenPathways, c('MALE_FASTED_WT', 'MALE_FASTED_HET', 'FEMALE_WT', 'FEMALE_HET'), "GO age in fasted", c(1,2,3,4,5,6), c(21,19,19,19,19,19))
pdf(paste0(path, 'figures/GO age in fasted.pdf'), width = 6.5, height = 4.5)
plot(godots)
dev.off()

for (p in 1:length(goID)){
  pdf(paste0(path, 'figures/FASTED AGE ', names(goID[p]), '.pdf'))
  
  geneSetHuman <- unique(goAnno[goAnno$GOID == goID[p],]$SYMBOL)
  geneSetKf <-c(humanKfConversion[humanKfConversion$human %in% geneSetHuman,]$ncbi)

  # FC between young and old
  hm <- actualHeatmap(FCyVSo[,3:6], infoTableFCyVSo[3:6,], geneSetKf, "increase with age", scaleT = F, splitCol = 'sex',
                      filterGenes = T, cutoff = log2(1.5), show_column_names =F) 
  draw(hm, column_title = paste(names(goID[p]), "FC > 1.5, old/young"))
  
  dev.off()
}


# male vs female (age and genotype) - Figure 4d and S4c ----
chosenPathways <- read.csv(paste0(path, 'inputFiles/pathways to age+genotype in fasted.csv'), row.names=1)
goID <- chosenPathways$ID
names(goID) <- chosenPathways$Description

godots <- GOdots(GoCompareEnrichmentGroup, chosenPathways, c('MALE_FASTED', 'FEMALE'), "GO age+genotype in fasted", c(1,1.5,2), c(21,19,19))
pdf(paste0(path, 'figures/GO age+genotype in fasted.pdf'), width = 5.5, height = 5.5)
plot(godots)
dev.off()

for (p in 1:length(goID)){
  pdf(paste0(path, 'figures/FASTED ', names(goID[p]), '.pdf'))
  
  geneSetHuman <- unique(goAnno[goAnno$GOID == goID[p],]$SYMBOL)
  geneSetKf <-c(humanKfConversion[humanKfConversion$human %in% geneSetHuman,]$ncbi)
  
  #FC between WT and Het
  cc <- c(2,3,5,6)
  hm <- actualHeatmap(FCwtVShet[,cc], infoTableFCwtVShet[cc,], geneSetKf, "increase in Het", scaleT = F, splitCol = 'sex',
                       filterGenes = T, cutoff = log2(1.5), show_column_names =F)   
  draw(hm, column_title = paste(names(goID[p]), "FC > 1.5, Het/WT"))

  dev.off()
}

# fission and fusion- Figure S4d ----
geneSetKf <- read.csv(paste0(path, 'inputFiles/fission fusion biogenesis.csv'), row.names=1)$symbol
pdf(paste0(path, "figures/Heatmap fusion fission and biogenesis.pdf"))
infoTableFCwtVShet$feedsex <- factor(paste(infoTableFCwtVShet$feed, infoTableFCwtVShet$sex), levels = c('full male', 'fasted male', 'fasted female'))
hm <- actualHeatmap(FCwtVShet[,c(2,3,5,6)], infoTableFCwtVShet[c(2,3,5,6),], geneSetKf, "increase in Het",
                    filterGenes = F, clusterRows = F, geneOrder = geneSetKf, splitCol = 'feedsex', scaleT = F, show_column_names =F)
draw(hm, column_title = "FC het/WT")
dev.off()

# GO age and genotype in all (4 pathways)- Figure S5b ----
chosenPathways <- read.csv(paste0(path, 'inputFiles/pathways to age+genotype in all.csv'), row.names=1)
godots <- GOdots(GoCompareEnrichmentGroup, chosenPathways, c('MALE_FULL', 'MALE_FASTED', 'FEMALE'), "GO age+genotype in all")
pdf(paste0(path, 'figures/GO age+genotype in all.pdf'), width = 7, height = 3)
plot(godots)
dev.off()

# GSVA- Steroid hormones, estrogen and testosterone Figure S6g----
chosenPathways <- read.csv(paste0(path, 'inputFiles/pathway steroid hormone.csv'), row.names=1)
goID <- chosenPathways$ID
names(goID) <- chosenPathways$Description

geneSetList <- list()
for (p in 1:length(goID)){
  geneSetHuman <- unique(goAnno[goAnno$GOID == goID[p],]$SYMBOL)
  geneSetKf <- c(humanKfConversion[humanKfConversion$human %in% geneSetHuman,]$ncbi)
  geneSetList[[names(goID)[p]]] <- geneSetKf
}

pdf(paste0(path, 'figures/SteroidEstrogenTestosterone.pdf'))

GSVAhm <- GSVAheatmap(cpmMaleFast, DEobj$samples[DEobj$samples$sex == 'male' & DEobj$samples$feed == 'fasted',], geneSetList)
draw(GSVAhm, column_title = "male fasted")

GSVAhm <- GSVAheatmap(cpmFemale, DEobj$samples[DEobj$samples$sex == 'female' & DEobj$samples$feed == 'fasted',], geneSetList)
draw(GSVAhm, column_title = "female fasted")
dev.off()

# GSVA- energy pathways- KEGG Figure S7d----
res <- keggList("pathway", "hsa")
names(res) <- gsub('path:', '', names(res))
keggID <- c('hsa04152', 'hsa04150', 'hsa04935', 'hsa04910') # AMPK, mTOR, Growth hormone, insulin 
energy <- list()

for (k in keggID){
  aa = names(unlist(fromJSON(file=curl(paste0("http://togows.dbcls.jp/entry/pathway/",k,"/genes.json")))))
  geneSetHuman = bitr(as.character(aa), fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL
  energy[[res[[k]]]] <- c(humanKfConversion[humanKfConversion$human %in% geneSetHuman,]$ncbi)
}
names(energy) <- gsub(' - Homo sapiens \\(human\\)', '', names(energy))
pdf(paste0(path, 'figures/EnergyKEGG.pdf'))

GSVAhm <- GSVAheatmap(cpmMaleFast, DEobj$samples[DEobj$samples$sex == 'male' & DEobj$samples$feed == 'fasted',], energy)
draw(GSVAhm, column_title = "male fasted")
dev.off()
# AMPK subunits expression ----
geneSetKf <- read.csv(paste0(path, 'inputFiles/ampk subunits.csv'), row.names=1)$symbol
geneSetKf <- geneSetKf[geneSetKf %in% rownames(cpmMaleFull)]

#heatmap
pdf(paste0(path, '/figures/ampk heatmap_FC.pdf'))
hm <- actualHeatmap(FCwtVShet[,c(1,4)], infoTableFCwtVShet[c(1,4),], geneSetKf, "increase in Het",
                    filterGenes = F, clusterRows = F, geneOrder = geneSetKf, scaleT = F, show_column_names =F)
draw(hm, column_title = "FC het/WT")
dev.off()

infoTable <- DEobj$samples[DEobj$samples$feed == 'full' & DEobj$samples$sex == 'male',]
AMPKdata <- 2^cpmMaleFull[c('prkag1','LOC107383943'),]

WTmeans <- rowMeans(AMPKdata[, infoTable[infoTable$age=='6.5weeks' & infoTable$genotype== 'WT',]$sample])

AMPKsu <- merge(melt(AMPKdata[, infoTable$sample] / WTmeans), infoTable, by.x='Var2', by.y='sample')
aprtPlot <- ggplot(AMPKsu[AMPKsu$Var1 == 'LOC107383943',], aes(x=Var1, y=value, xlab=group), fill='') +
  geom_bar(stat='Summary', fun = 'mean', color="black", fill='white',  width=0.8, size=1, position = "dodge") +
  geom_errorbar(stat = "Summary", fun.data = mean_se, width=0.8, size=0.7, position = "dodge") +
  geom_point(position = position_dodge(width = 0.8), aes(x=Var1, color=group), size=2.5) +
  scale_color_manual(values = c('black', '#008000', 'black', '#008000')) +
  xlab("") +
  #scale_y_continuous(limits = c(0,5.5), expand = c(0, 0)) +
  #geom_hline(yintercept=0.5, linetype="dashed")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(color="black"), axis.text.x = element_text(color="black"))

pdf(paste0(path, '/figures/AMPK subunits.pdf'), width = 4, height = 2)
plot(aprtPlot)
dev.off()