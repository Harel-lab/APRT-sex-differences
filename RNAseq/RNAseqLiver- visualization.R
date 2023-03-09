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
library(cowplot)          #1.1.1
library(stringr)

path <- 'RNAseq/'
load(paste0(path,'outputFiles/data.RDataRev'))

# To convert NCBI ids to human entrez ids. There are ways to adapt it for nfur only, but for now I do everything based on human orthologs
humanKfConversion = read.table(paste0(path, "inputFiles/NCBI-Human-orthologs.txt"), head = T, sep = "\t")


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
  NCBIgeneNames <- read.csv("C:/Users/tehil/Dropbox/Projects/genes_names4.csv")
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
  
  file_name <-  paste0(path, 'GO/Results/', group[ngroup], ' ', test, '_ALL.csv') # '_GOGSEA.csv')
  print(file_name)
  ExTable <- read.csv(file_name)
  print(dim(ExTable))
  genes <- ExTable[ExTable$qvalues < pvalCutoff  ,]$ID #& ExTable$ONTOLOGY == ont
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
    GG <- extractSigGenes(path, group, i, test, ont, pvalCutoff)
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
  colNamesGeneral <- c( 'ID', 'Description')#'ONTOLOGY',
  
  
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
  
  comparisonsTable <- threeGroupsCompare(path, group, test, DBs, nameG, 'BP', 1)
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

#rev
group <- setNames(groups[[3]], gsub('FASTED_', '', groups[[3]]))
sendToCompareSaveXlsx(path, group, DBs, 'rev', tests[3], '')

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

  cpm.data <- na.omit(cpm.data)
  
  if (!clusterRows){
    #head(cpm.data)
    geneOrder <- geneOrder[geneOrder %in% row.names(cpm.data)]
    cpm.data <- cpm.data[geneOrder,]
    hr = F
  }
  else
    #hr <- hclust(as.dist(1-cor(t(cpm.data), method="spearman"))) #Cluster rows by Pearson correlation.
    hr <- hclust(as.dist(1-cor(t(cpm.data), method="pearson"))) # Cluster rows by Pearson correlation.

  if (scaleT)
    data.scaled <- t(scale(t(cpm.data))) 
  else{
    data.scaled <- as.matrix(cpm.data)
  }
  print(data.scaled)
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
  GoCompareEnrichmentGroup <- GoCompareEnrichmentGroup[GoCompareEnrichmentGroup$ID %in% chosenPathways$ID,]# -grep('core_', colnames(GoCompareEnrichmentGroup))]
  cnames <- c('ID', 'Description', colnames(GoCompareEnrichmentGroup)[grep(paste(group, collapse = '|'), colnames(GoCompareEnrichmentGroup))])
  GoCompareEnrichmentGroup <- GoCompareEnrichmentGroup[,cnames]
  
  #order the pathways according FDR or the original order
  pathwaysOrder <- chosenPathways$Description
  
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
GoCompareEnrichmentGroupGenotype <- read.xlsx(paste0(path, '/GO/Results/MergedAnalysis/GOBP_compare_genotype in all.xlsx'), sheetName = 'genotype')
GoCompareEnrichmentGroupAG <- read.xlsx(paste0(path, '/GO/Results/MergedAnalysis/GOBP_compare_age+genotype in all.xlsx'), sheetName = 'age+genotype')

# Figure 1F
chosenPathways <- read.csv(paste0(path, 'inputFiles/pathways genotype male full.csv'), row.names=1)
goID <- chosenPathways$ID
names(goID) <- chosenPathways$Description

godots <- GOdots(GoCompareEnrichmentGroupGenotype, chosenPathways, c('MALE_FULL_YOUNG', 'MALE_FULL_OLD'), "GO genotype, Male fed") #, c(1,2,3,4,5,6), c(21,19,19,19,19,19))
pdf(paste0(path, 'figures/GO GSVA genotype in make full.pdf'), width = 4, height = 5)
plot(godots)
dev.off()

# comparison between ages in fasted (males and females, WT and Het)- Figure 2C
chosenPathways <- read.csv(paste0(path, 'inputFiles/pathways age in fasted.csv'), row.names=1)
goID <- chosenPathways$ID
names(goID) <- chosenPathways$Description

godots <- GOdots(GoCompareEnrichmentGroupAge, chosenPathways, c('MALE_FASTED_WT', 'MALE_FASTED_HET', 'FEMALE_WT', 'FEMALE_HET'), "GO age in fasted", c(1,2,3,4,5,6), c(21,19,19,19,19,19))
pdf(paste0(path, 'figures/GO age in fasted.pdf'), width = 6.5, height = 5)
plot(godots)
dev.off()

# barplots ----
barplotExpression <- function(gene){
  aprt <- data.frame(Exp = c(cpmMaleFull[gene,], cpmMaleFast[gene,], cpmFemale[gene,]), 
                     sample=DEobj$samples$sample, group=DEobj$samples$group, genotype=DEobj$samples$genotype, age=DEobj$samples$age)
  aprt$newExp <- 2^aprt$Exp -1 # CPM without log2
  aprt$group <- gsub('6.5weeks |15weeks ', '', aprt$group)
  aprt$group <- factor(aprt$group, levels = unique(aprt$group))
  
  #normalize each Het + wt by correspond average wt.
  averageGroups <- aggregate(aprt$newExp, list(aprt$group), mean)
  averageGroups$x[seq(from=2, to=6, length=3)] <- averageGroups$x[seq(from=2, to=6, length=3)-1] 
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
  pvalsdf[, 'p.adj'] <- round(pvalsdf$p, 4)  #round(p.adjust(pvalsdf$p, method='fdr'), 4) 
  
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
barplotExpression2 <- function(expp, infoData, gene, normBar = T){
  # presenting all 4 groups
  ppar <- data.frame(Exp = expp, 
                     sample=infoData$sample, group=infoData$group, 
                     genotype=infoData$genotype, age=infoData$age)
  ppar$newExp <- 2^ppar$Exp -1 #cpm without log2
  aa <- aggregate(ppar$newExp, list(ppar$group), mean)
  aa$x[1:4] <- aa$x[1] #normalize the males to WT young male.
  ppar = merge(ppar, aa, by.x="group", by.y="Group.1", sort=F)
  ppar$normExp = ppar$newExp / ppar$x
  
  ppar$gg <- gsub('WT |Het ', '', ppar$group)
  ppar$gg <- factor(ppar$gg, levels = unique(ppar$gg))
  
  pvalsdf <- data.frame(matrix(0, ncol=5, nrow = 4))   #####4
  colnames(pvalsdf) <- c('group1', 'group2', 'p', 'p.adj', 'y.position')
  
  pvalsdf$y.position <- 1.2
  
  aov.results <- aov(Exp ~ age * genotype, data = ppar)
  pvals <- round(as.numeric(summary(aov.results)[[1]][["Pr(>F)"]][1:3]), 4)
  tukey <- data.frame(TukeyHSD(aov.results)[[3]])
  tukey[c('group1', 'group2')] <- str_split_fixed(gsub(':', ' ', rownames(tukey)), '-', 2)
  tukey$group1 <- paste(tukey$group1, gsub('6.5weeks ', '', ppar$gg[1])) 
  tukey$group2 <- paste(tukey$group2, gsub('6.5weeks ', '', ppar$gg[1])) 
  tukey$p.adj <- round(tukey$p.adj, 4) 
 
  if (normBar){
    tukey$y.position <- seq(max(ppar$normExp)-0.25,max(ppar$normExp)+1.5, length.out=6)
    pparPlot <- ggplot(ppar, aes(x=group, y=normExp)) +
      geom_bar(stat='Summary', fun = 'mean', color="black", fill='white', width=0.75, size=1) +
      geom_errorbar(stat = "Summary", fun.data = mean_se, width=0.3) +
      geom_point(position = position_jitter(width=0.15), aes(color=genotype), size=2.5) +
      scale_color_manual(values = c('black', '#008000')) +
      scale_x_discrete(guide = guide_axis(n.dodge = 4)) +
      #stat_pvalue_manual(tukey, label = 'p.adj', size = 4, hide.ns = T) +
      ggtitle(toupper(gene), subtitle = paste(pvals, collapse = ' ')) +
      xlab("") +
      scale_y_continuous(limits = c(0,max(ppar$normExp)+0.5), expand = c(0, 0)) +
      theme_classic() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(color="black"), axis.text.x = element_text(color="black"))
  }  
  else {
    tukey$y.position <- seq(max(ppar$Exp)-0.5,max(ppar$Exp)+3.5, length.out=6)
    pparPlot <- ggplot(ppar, aes(x=group, y=Exp)) +
      geom_bar(stat='Summary', fun = 'mean', color="black", fill='white', width=0.75, size=1) +
      geom_errorbar(stat = "Summary", fun.data = mean_se, width=0.3) +
      geom_point(position = position_jitter(width=0.15), aes(color=genotype), size=2.5) +
      scale_color_manual(values = c('black', '#008000')) +
      scale_x_discrete(guide = guide_axis(n.dodge = 4)) +
      stat_pvalue_manual(tukey, label = 'p.adj', size = 4, hide.ns = T) +
      ggtitle(toupper(gene), subtitle = paste(pvals, collapse = ' ')) +
      xlab("") +
      scale_y_continuous(limits = c(0,max(ppar$Exp)+5), expand = c(0, 0)) +
      theme_classic() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(color="black"), axis.text.x = element_text(color="black"))
    
  }
  pparPlot
}
barplotExpression3 <- function(gene){
  #separate by genotype; 4 groups
  infoMalesFast = DEobj$samples[DEobj$samples$feed =='fasted' & DEobj$samples$sex =='male',]
  infoMalesFast$group = factor(infoMalesFast$group, levels = c('6.5weeks WT fasted male', '15weeks WT fasted male', '6.5weeks Het fasted male', '15weeks Het fasted male'))
  ppar <- data.frame(Exp = c(cpmMaleFast[gene,]), 
                     sample=infoMalesFast$sample, group=infoMalesFast$group, 
                     genotype=infoMalesFast$genotype, age=infoMalesFast$age)
  ppar <- ppar[order(ppar$group),]
  ppar$newExp <- 2^ppar$Exp -1 #cpm without log2
  aa <- aggregate(ppar$newExp, list(ppar$group), mean)
  aa$x[1:2] <- aa$x[1] #normalize the WTs to WT young male.
  aa$x[3:4] <- aa$x[3] #normalize the Hets to Het young male.
  
  ppar = merge(ppar, aa, by.x="group", by.y="Group.1", sort=F)
  ppar$normExp = ppar$newExp / ppar$x
  
  ppar$gg <- gsub('6.5weeks |15weeks ', '', ppar$group)
  ppar$gg <- factor(ppar$gg, levels = unique(ppar$gg))
  
  pvalsdf <- data.frame(matrix(0, ncol=5, nrow = 2))   
  colnames(pvalsdf) <- c('group1', 'group2', 'p', 'p.adj', 'y.position')
  
  for (i in 1:length(unique(ppar$gg))){
    pvalsdf[i, c('group1', 'group2')] = unique(ppar$group)[c(2*i-1,2*i)]
    c1 <- ppar[ppar$gg == unique(ppar$gg)[i],]
    pvalsdf[i, 'p'] <- t.test(normExp ~ age, data=c1)$p.value
  }
  pvalsdf[, 'p.adj'] <- round(pvalsdf$p, 4)
  
  pvalsdf$y.position <- 2
  
  pparPlot <- ggplot(ppar, aes(x=group, y=normExp)) +
    geom_bar(stat='Summary', fun = 'mean', color="black", fill='white', width=0.75, size=1) +
    geom_errorbar(stat = "Summary", fun.data = mean_se, width=0.3) +
    geom_point(position = position_jitter(width=0.15), aes(color=genotype), size=2.5) +
    scale_color_manual(values = c('black',  '#008000')) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    stat_pvalue_manual(pvalsdf, label = 'p.adj', size = 4) +
    ggtitle(toupper(gene)) +
    xlab("") +
    scale_y_continuous(limits = c(0,2.5), expand = c(0, 0)) +
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.y = element_text(color="black"), axis.text.x = element_text(color="black"))
  
  pparPlot
}


#figure S1I, S2B
pdf(paste0(path, '/figures/APRT expression.pdf'), width = 5.5, height = 3)
plot(barplotExpression('aprt'))
dev.off()

#Figures 1F, S1J, S1K
genesG <- read.csv(paste0(path, 'inputFiles/ampk subunits.csv'))
genesG <- genesG[genesG$NCBI %in% rownames(cpmMaleFull),]
pdf(paste0(path, '/figures/barplots ampk subunits.pdf'), width = 6, height = 3.5)

genesG <- read.csv(paste0(path, 'inputFiles/starvation full.csv'))
pdf(paste0(path, '/figures/barplots full.pdf'), width = 6, height = 3.5)

genesG <- read.csv(paste0(path, 'inputFiles/glyglu reactom.csv'))
pdf(paste0(path, '/figures/barplots glyglu.pdf'), width = 6, height = 3.5)
for (i in 1:nrow(genesG)){ #
  gene = genesG$NCBI[i]
  names(gene) = genesG$FinalSymbol[i]
  mT = barplotExpression2(cpmMaleFull[gene,], DEobj$samples[DEobj$samples$feed =='full' & DEobj$samples$sex =='male',],paste(names(gene), 'male norm'))
  mF = barplotExpression2(cpmMaleFull[gene,], DEobj$samples[DEobj$samples$feed =='full' & DEobj$samples$sex =='male',],paste(names(gene), 'male'), F)

  print(plot_grid(mT, mF,ncol = 2,align = "v"))
}
dev.off()

#Figures 2C, D, S2C
genesG <- data.frame(NCBI = geneSetKf$`Nucleotide salvage`, FinalSymbol = convertLoc2symbol(geneSetKf$`Nucleotide salvage`))
genesG <- genesG[genesG$NCBI %in% rownames(cpmFemale),]
genesG <- genesG[genesG$NCBI %in% rownames(cpmMaleFast),]
pdf(paste0(path, '/figures/barplots salvage.pdf'), width = 6, height = 5)

genesG <- read.csv(paste0(path, 'inputFiles/genes gly glu.csv'))
pdf(paste0(path, '/figures/barplots.pdf'), width = 6, height = 5)

genesG <- read.csv(paste0(path, 'inputFiles/age.csv'))
pdf(paste0(path, '/figures/barplots age.pdf'), width = 6, height = 5)
for (i in 1:nrow(genesG)){ #
  gene = genesG$NCBI[i]
  names(gene) = genesG$FinalSymbol[i]
  mT = barplotExpression2(cpmMaleFast[gene,], DEobj$samples[DEobj$samples$feed =='fasted' & DEobj$samples$sex =='male',],paste(names(gene), 'male norm'))
  fT = barplotExpression2(cpmFemale[gene,], DEobj$samples[DEobj$samples$feed =='fasted' & DEobj$samples$sex =='female',],paste(names(gene), 'female norm'))
  mF = barplotExpression2(cpmMaleFast[gene,], DEobj$samples[DEobj$samples$feed =='fasted' & DEobj$samples$sex =='male',],paste(names(gene), 'male'), F)
  fF = barplotExpression2(cpmFemale[gene,], DEobj$samples[DEobj$samples$feed =='fasted' & DEobj$samples$sex =='female',],paste(names(gene), 'female'), F)
  
  print(plot_grid(mT, fT, mF, fF, ncol = 2,align = "v"))
}
dev.off()


######Glycolysis and Gluconeogenesis Figure S1K #####

hGenes = read.csv(paste0(path, 'inputFiles/genes glucose lipids.csv'))
goID = unique(hGenes$ï..Des)
names(goID) = unique(hGenes$ï..Des)
geneSetKf = list()
for (it in names(goID)){
  geneSetKf[[it]] <- c(humanKfConversion[humanKfConversion$human %in% hGenes[hGenes$ï..Des == goID[it],]$Hunam.gene,]$ncbi)
}
geneSetKf[2] <- geneSetKf[1]
geneSetKf[4] <- geneSetKf[3]

geneG = merge(hGenes[hGenes$ï..Des %in% c('GlycolysisALL', 'GluconeogenesisALL'), ], humanKfConversion, by.x='Hunam.gene', by.y='human', all.x=T)
geneG$FinalSymbol = convertLoc2symbol(geneG$ncbi)
write.csv(geneG, paste0(path, 'outputFiles/glyglu genes.csv'))

GSVAhmMfedGly <- GSVAheatmap(cpmMaleFull, DEobj$samples[DEobj$samples$sex == 'male' & DEobj$samples$feed == 'full',], geneSetKf[c(1,2)])
GSVAhmMfedglu <- GSVAheatmap(cpmMaleFull, DEobj$samples[DEobj$samples$sex == 'male' & DEobj$samples$feed == 'full',], geneSetKf[c(3,4)])
pdf(paste0(path, '/figures/GSVA glyglu.pdf'), width = 6, height = 4)
draw(GSVAhmMfedGly, column_title = "gly male fed")
draw(GSVAhmMfedglu, column_title = "glu male fed")
dev.off()

#################### compare to mouse data :https://doi.org/10.1016/j.cmet.2019.06.018 #################### 
GoGenotype = read.xlsx(paste0(path, '/GO/Results/MergedAnalysis/GO_compare_genotype in rev.xlsx'), sheetName = 'genotype')
GoGenotype$ID = gsub('_', ' ', GoGenotype$ID)
GoGenotype$ID = gsub('GOBP', 'GO', GoGenotype$ID)
GoGenotype$ID = gsub('GOMF', 'GO', GoGenotype$ID)
GoGenotype = GoGenotype[,colnames(GoGenotype) %in% c("ID", colnames(GoGenotype)[grep('NES|qvalue', colnames(GoGenotype))])]
GoGenotypeFed = GoGenotype[,-grep('MALE_YOUNG|MALE_OLD|FEMALE', colnames(GoGenotype))]
GoGenotypeFedO = GoGenotypeFed[,-grep('YOUNG', colnames(GoGenotypeFed))]

miceData = read.csv(paste0(path, 'mice/mice_data.csv'),row.names = 1)
miceDataFish = merge(miceData, GoGenotypeFedO, by.x=0, by.y='ID', all.x=T)
rownames(miceDataFish) = miceDataFish$Row.names

#dotplot
miceDataFish$ID = miceDataFish$Row.names
miceDataFish$Description = miceDataFish$Row.names
colnames(miceDataFish) = gsub('q.value', 'qvalues', colnames(miceDataFish))
chosenPathways = read.csv(paste0(path, 'mice/selected pathway mice.csv'))
chosenPathways$Description = chosenPathways$ID
godots <- GOdots(miceDataFish, chosenPathways, c('MR.M', 'estradiol.M', 'MALE_FULL_OLD', 'Protandim.M', 'Rapamycin.M', 
                                                 'CR.M', 'Snell.dwarf.M', 'GHRKO.M', 'Acarbose.M'), "GO genotype mice")
pdf(paste0(path, 'mice/dotplot.pdf'), width = 7, height = 4.5)
show(godots)
dev.off()

col_fun = circlize::colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

miceDataFish = miceDataFish[,grep('NES', colnames(miceDataFish))]
miceDataFish = miceDataFish[,grep('.M', colnames(miceDataFish))]
infoMice = DataFrame(group=colnames(miceDataFish), sex=sub('.*(?=.$)', '', colnames(miceDataFish), perl=T))
infoMice$sex = factor(gsub('F', 'female', gsub('M', 'male', infoMice$sex)), levels = c('male', 'female'))  
hrorder = actualHeatmap(miceDataFish, infoMice, row.names(miceDataFish), 'mice', convertNames = F, clusterSamples = T, scaleT = F, show_row_names = FALSE, col=col_fun)
hrorderSelect = actualHeatmap(miceDataFish, infoMice, chosenPathways$ID, 'miceSelected', geneOrder = chosenPathways$ID, convertNames = F, clusterRows = F, col=col_fun,  clusterSamples = T, scaleT = F)
write.csv(hrorder, paste0(path, 'mice/hr order.csv'))

pdf(paste0(path, 'mice/heatmap2.pdf'), width = 6.5, height = 8)
show(hrorder)
dev.off()
pdf(paste0(path, 'mice/heatmap selected.pdf'), width = 6.5, height = 5)
show(hrorderSelect)
dev.off()


