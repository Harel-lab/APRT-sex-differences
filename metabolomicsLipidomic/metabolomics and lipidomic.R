# Metabolomics and lipidomics- 
# read, filtering, normalized.
# PCA and visualization.
# save in different metabolites names.
# calculate statistics for each metabolite or lipid
# box plot for the expression for each metabolite/lipid

library(ggplot2)
library(ggfortify)
library(factoextra)
library(heatmap3)
library(gplots)
library(dplyr)
library(tidyr)
library(devtools)
library(ComplexHeatmap)
library(reshape)
library(data.table)
library(FSA)
library(gridExtra)
library(grid)
library(plyr)
library(ggpubr)
library(ggforce)
library(mutoss)
library(ggConvexHull)

path <- 'C:/Users/tehil/Dropbox/paper/metabolomicsLipidomic/'

# read metabolomics data----
metabolomics <- read.csv(paste0(path, 'inputFiles/Polar compounds_norm_weight.csv'), row.names = 1)
#filter metabolites with less than 70% expression
metabolomics <- metabolomics[apply(metabolomics, 1, function(x){sum(x==0)}) < round(ncol(metabolomics)*0.7),]

#read lipidomics data
lipidomic <- read.csv(paste0(path, 'inputFiles/Lipids_norm_weight.csv'), row.names = 1)
lipidomic <- lipidomic[apply(lipidomic, 1, function(x){sum(x==0)}) < round(ncol(lipidomic)*0.7),]
lipid_groups <- gsub("\\(.*","",rownames(lipidomic)) #lipids speacies


# read information about the samples----
infoTable <- read.csv(paste0(path, 'inputFiles/sample information.csv'), row.names = 1, stringsAsFactors=T)
infoTable$genotype <- relevel(infoTable$genotype, ref = "WT")
infoTable$sex <- relevel(infoTable$sex, ref = "male")
infoTable$age <- relevel(infoTable$age, ref = "6.5weeks")
infoTable$feed <- relevel(infoTable$feed, ref = "full")
infoTable$sample <- row.names(infoTable)
infoTable <- infoTable[order(infoTable$genotype),]
infoTable <- infoTable[order(infoTable$age),]
infoTable <- infoTable[order(infoTable$feed),]
infoTable <- infoTable[infoTable$sex == 'male',]

infoTable$group <- paste(infoTable$age, infoTable$genotype, infoTable$feed, infoTable$sex ,sep=' ')
infoTable$group <- factor(infoTable$group, levels = unique(infoTable$group))

infoTableL <- data.frame(infoTable) # copy infoTable for the lipids

# order the data, the first row in info correspond to the first column in metabolomics
colnames(metabolomics) <- gsub('X', '', colnames(metabolomics))
head(metabolomics)
metabolomics <- metabolomics[,row.names(infoTable)]

#lipids
colnames(lipidomic) <- gsub('X', '', colnames(lipidomic))
head(lipidomic)
lipidomic <- lipidomic[,row.names(infoTableL)]


#log transpose and scaling. adding 1 to avoid log(0)
metabolomics[metabolomics==0] <- 1
metabolomics <- log2(metabolomics)
metabolomics <- metabolomics / apply(metabolomics, 1, mean)

#lipids
lipidomic[lipidomic==0] <- 1
lipidomic <- log2(lipidomic)
lipidomic <- lipidomic / apply(lipidomic, 1, mean)


# data explorer PCA----
PCA_M_L <- function(data.M.L, infoTable, splitM=T, group='all'){
  if (splitM){
    infoTable <- infoTable[infoTable$feed == group,]
    data.M.L <- data.M.L[, row.names(infoTable)]
  }

  pca.de <- prcomp(t(data.M.L), scale. = T)
  knee <- fviz_eig(pca.de)
  ellipseType = "convex"
  pca <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=2,label =F)+
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group),
                    show.legend = F, alpha=0.2)+
    ggtitle(group)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  pcaD <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=2,label =F,
           size = 2.5, title=titleG) + #, shape = 19
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group, alpha=infoTable$group),
                    show.legend = T)+
    scale_color_manual(name='Group', values = c('#1D71BB', '#E30613', '#1D71BB', '#E30613'))+
    scale_fill_manual(name='Group', values = c('#1D71BB', '#E30613', '#1D71BB', '#E30613'))+
    scale_alpha_manual(name='Group', values = c(0.4,0.4,0,0))+
    theme(legend.position="bottom", legend.direction="vertical")+
    ggtitle(group) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  return(list(knee=knee, PCAlabels = pca, PCAdots=pcaD))
}

# metabolomics
all_M <- PCA_M_L(metabolomics, infoTable, F)
full_M <- PCA_M_L(metabolomics, infoTable, group = 'full')
fasted_M <- PCA_M_L(metabolomics, infoTable, group = 'fasted')
pdf(paste0(path, 'figures/PCA Metabolomics.pdf'), width=17, height=4)
grid.arrange(all_M$PCAlabels, full_M$PCAdots, fasted_M$PCAdots, nrow = 1)
dev.off()

#lipids
all_L <- PCA_M_L(lipidomic, infoTableL, F)
full_L <- PCA_M_L(lipidomic, infoTableL, group = 'full')
fasted_L <- PCA_M_L(lipidomic, infoTableL, group = 'fasted')
pdf(paste0(path, 'figures/PCA lipidomic.pdf'), width=17, height=4)
grid.arrange(all_L$PCAlabels, full_L$PCAdots, fasted_L$PCAdots, nrow = 1)
dev.off()

# over all function-----
averageSamples <- function(data.M.L, infoTable){
  # calculate the average of each group
  meanData <- aggregate(t(data.M.L), list(infoTable$group), mean)
  meanData <- separate(meanData, Group.1, c('age', 'genotype', 'feed', 'sex'), sep=" ")
  meanInfoTable <- meanData[,colnames(meanData) %in% c('age', 'genotype', 'feed', 'sex')]
  meanInfoTable$group <- paste(meanInfoTable$age, meanInfoTable$genotype, meanInfoTable$feed, meanInfoTable$sex ,sep=' ')
  meanData <- t(meanData[,!colnames(meanData) %in% c('age', 'genotype', 'feed', 'sex')])
  
  colnames(meanData)<- row.names(meanInfoTable)
  meanData <- meanData[,rownames(meanInfoTable)]
  
  # set the infoTable_f as factors
  meanInfoTable$genotype <- factor(meanInfoTable$genotype, levels = c("WT", "Het"))
  meanInfoTable$age <- factor(meanInfoTable$age, levels = c("6.5weeks", "15weeks"))
  meanInfoTable$sex <- factor(meanInfoTable$sex, levels = c("male", "female"))
  meanInfoTable$feed <- factor(meanInfoTable$feed, levels = c("full", "fasted"))
  
  return(list(meanData = meanData, infoTable = meanInfoTable))
}
heatmapMetbolomics <- function(data.m, infoTable, hmName, splitBy=NA, sampleClustering=F, rowsClustering=T, convertMetaNames=T, ...){
  hc=F
  hr=F
  if (sampleClustering)
    hc <- hclust(as.dist(1-cor(data.m, method="pearson"))) # Clusters columns by Perason correlation.
  if (rowsClustering)
    hr <- hclust(as.dist(1-cor(t(data.m), method="pearson"))) # Cluster rows by Pearson correlation.
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  data.scaled <- t(scale(t(data.m)))
  if(convertMetaNames){
    rownames(data.scaled) <- convertMetabolitesNames(rownames(data.scaled))
  }
  cannotaion <- HeatmapAnnotation(age=infoTable$age, genotype=infoTable$genotype, feed=infoTable$feed,
                                  col = list(age=c('6.5weeks'='gray', '15weeks'='black'),
                                             genotype=c('WT'= 'blue', 'Het'='red'),
                                             feed=c('fasted'= 'lightgreen', 'full'= 'forestgreen')),
                                  simple_anno_size = unit(1, "mm"), annotation_name_gp = gpar(fontsize = 4))
  
  if (!is.na(splitBy)){
    if (rowsClustering)
      hm = Heatmap(data.scaled, split = splitBy, name = hmName, top_annotation = cannotaion,
                   cluster_columns = hc, col=col_fun, row_title_rot =0, row_title_gp = gpar(fontsize =8),
                   row_names_gp = gpar(fontsize = 3), ...)
    else 
      hm = Heatmap(data.scaled, split = splitBy, name = hmName, top_annotation = cannotaion,
             cluster_columns = F, cluster_rows = F, col=col_fun, row_title_rot =0, row_title_gp = gpar(fontsize =7),
             row_names_gp = gpar(fontsize = 3),...)
  }
  else
    hm = Heatmap(data.scaled, name = hmName, top_annotation = cannotaion,
                 cluster_columns = hc, cluster_rows = hr, col=col_fun, row_title_rot =0, row_title_gp = gpar(fontsize =8),
                 row_names_gp = gpar(fontsize = 2),...)
  return(hm)
}
heatmapLipids <- function(data.l, infoTable, hmName, splitBy=NA, sampleClustering=F, rowsClustering=T, ...){
  hc=F
  hr=F
  if (sampleClustering)
    hc <- hclust(as.dist(1-cor(data.l, method="pearson"))) # Clusters columns by Pearson correlation.
  if (rowsClustering)
    hr <- hclust(as.dist(1-cor(t(data.l), method="pearson"))) # Cluster rows by Pearson correlation.
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  data.scaled <- t(scale(t(data.l)))

  cannotaion <- HeatmapAnnotation(age=infoTable$age, genotype=infoTable$genotype, feed=infoTable$feed,
                                  col = list(age=c('6.5weeks'='gray', '15weeks'='black'),
                                             genotype=c('WT'= 'blue', 'Het'='red'),
                                             feed=c('fasted'= 'lightgreen', 'full'= 'forestgreen')),
                                  simple_anno_size = unit(1, "mm"), annotation_name_gp = gpar(fontsize = 4))
  
  if (!is.na(splitBy)){
    rownames(data.scaled) <- gsub(paste(unique(splitBy), collapse = '|'), '', rownames(data.scaled))
    if (rowsClustering)
      hm = Heatmap(data.scaled, split = splitBy, name = hmName, top_annotation = cannotaion,
                   cluster_columns = hc, col=col_fun, row_title_rot =0, row_title_gp = gpar(fontsize =8),
                   row_names_gp = gpar(fontsize = 3),...)
    else 
      hm = Heatmap(data.scaled, split = splitBy, name = hmName, top_annotation = cannotaion,
                   cluster_columns = F, cluster_rows = F, col=col_fun, row_title_rot =0, row_title_gp = gpar(fontsize =7),
                   row_names_gp = gpar(fontsize = 3),...)
    }
  else
    hm = Heatmap(data.scaled, name = hmName, top_annotation = cannotaion,
                 cluster_columns = hc, cluster_rows = hr, col=col_fun, row_title_rot =0, row_title_gp = gpar(fontsize =8),
                 row_names_gp = gpar(fontsize = 2),...)
  return(hm)
}
convertMetabolitesNames <- function(old_names){
  new_names <- converting_M_names[converting_M_names$Biocyc_query %in% old_names,]
  new_names <- new_names[match(old_names, new_names$Biocyc_query),]$final
  return(new_names)
}


# ************* convert the metabolites names to biocyc names, HMDB and KEGG******************----
#create new table to convert metabolites names for BioCyc query
converting_M_names <- read.csv(paste0(path, 'inputFiles/converting metabolites names.csv'))
whole_names <- merge(metabolomics, converting_M_names, by.x = "row.names", by.y = "Metabolite")

# BioCyc
biocycData <- whole_names[whole_names$Biocyc_query != '0', colnames(whole_names) %in% c('Biocyc_query', infoTable$sample)]
row.names(biocycData) <- biocycData$Biocyc_query
biocycData <- biocycData[,!colnames(biocycData) %in% c('Biocyc_query')]
biocycData <- biocycData[,infoTable$sample]
#write.csv(biocycData, paste0(path, 'outputFiles/PC norm weight newnames biocyc.csv'))

#calculate average of each group
biocycData_mean <- aggregate(t(biocycData), list(infoTable$group), mean)
row.names(biocycData_mean) <- biocycData_mean$Group.1
biocycData_mean <- biocycData_mean[,!colnames(biocycData_mean) %in% c("Group.1")]

biocycData_mean <- t(biocycData_mean) # sample in columns
biocycData_mean<- t(scale(t(biocycData_mean))) # biocycData_mean / apply(biocycData_mean, 1, scale)
#write.table(biocycData_mean, paste0(path, 'outputFiles/PC norm weight newnames biocyc_mean.txt'), sep ='\t', quote = FALSE)
# bioCyc: https://biocyc.org/dashboard/dashboard-intro.shtml



# ************* analysis after biocyc- Metabolomics **********************----

# order the metabolomics like infoTable order, the first row in info correspond to the first column in metab and so on
metaboDataBiocycNames <- biocycData[,row.names(infoTable)]
metabColNames <- c("mainGroup", "subPathway", "biocyc_short_names", "biocyc_full_names", "metabolite")

# heatmap for all the metabolites
HM_sample <- heatmapMetbolomics(metaboDataBiocycNames, infoTable, 'metabolomics')
average.metab <- averageSamples(metaboDataBiocycNames, infoTable)
HM_average <- heatmapMetbolomics(average.metab$meanData, average.metab$infoTable, 'metabolomics-average', show_column_names = F)

pdf(paste0(path, 'figures/SUP M heatmap all metabolites samples.pdf'), width = 8, height = 10)
draw(HM_sample)
dev.off()

pdf(paste0(path, 'figures/M heatmap all metabolites average.pdf'), width = 7, height = 8)
draw(HM_average)
dev.off()
# ---- statistic test two-way ANOVA + Tukey----
anovaFULL <- function(df){
  df <- df[df$feed == 'full',]
  df$group <- paste(df$age, df$genotype)
  aov.results <- aov(m ~ group, data=df)
  pvals <- list(ANOVAfull = summary(aov.results)[[1]][["Pr(>F)"]][1])
  return(pvals)
}

anova4groups.tukey <- function(df){
  df$group <- paste(df$feed, df$age, df$genotype)
  df$cluster[df$feed == 'full'] <- 'full'
  df$cluster[df$age == '6.5weeks' & df$feed == 'fasted'] <- 'young'
  df$cluster[df$age == '15weeks' & df$feed == 'fasted' & df$genotype == 'WT'] <- 'oldWT'
  df$cluster[df$age == '15weeks' & df$feed == 'fasted' & df$genotype == 'Het'] <- 'oldHet'
  aov.results <- aov(m ~ cluster, data = df)
  pvals <- list(anova = summary(aov.results)[[1]][["Pr(>F)"]][1])
  tukey <- data.frame(TukeyHSD(aov.results)[[1]])
  pvals[rownames(tukey)] <- tukey$p.adj
  return(pvals)
}


# statistic test metabolomics ----
statsMetabo <- data.frame(row.names = row.names(metaboDataBiocycNames),
                          no = rep(0, nrow(metaboDataBiocycNames)))

metaboData <- metaboDataBiocycNames[,rownames(infoTable)]


for (i in 1:nrow(metaboData)){
  df <- data.frame(m = unlist(metaboData[i,]), feed=infoTable$feed, age=infoTable$age, genotype=infoTable$genotype)
  pvals <- anovaFULL(df)
  statsMetabo[i, names(pvals)] <- pvals
  pvals <- anova4groups.tukey(df)
  statsMetabo[i, names(pvals)] <- pvals
  }

statsMetaboAdjust <- data.frame(statsMetabo)
colnames(statsMetaboAdjust) <- colnames(statsMetabo)
for (a in colnames(statsMetaboAdjust)){
  statsMetaboAdjust[[a]] <- p.adjust(statsMetabo[[a]], method = 'fdr')
}


#------- boxplot of each metabolite-------
boxes <- list()
for (metabolit in rownames(statsMetaboAdjust[statsMetaboAdjust$anova < 0.1 & statsMetaboAdjust$`oldWT-oldHet` < 0.1 & statsMetaboAdjust$`young-full` < 0.1,])){
  df <- data.frame(m= unlist(metaboData[metabolit,]), genotype=infoTable$genotype, age=infoTable$age, feed=infoTable$feed,
                   group=factor(gsub(' male', '', infoTable$group), levels = c("6.5weeks WT full", "6.5weeks Het full", "15weeks WT full", "15weeks Het full", "6.5weeks WT fasted", "6.5weeks Het fasted", "15weeks WT fasted", "15weeks Het fasted")))

  boxp <- ggplot(df, aes(x=group, y=m)) +
    geom_boxplot(color='black', fill = 'white', width=0.5, lwd=0.4, fatten = 0.8, outlier.shape = NA, coef=10) +
    stat_summary(fun = max, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., group = group),
                 width = 0.2, linetype = "solid", size=0.7)+
    stat_summary(fun = min, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., group = group),
                 width = 0.2, linetype = "solid", size=0.7)+
    geom_jitter(width=0.1, size=2.2, aes(color=genotype)) +
    scale_color_manual(values = c('black', '#098036'))+
    xlab("") + ylab("") +
    theme_classic() + 
    #scale_y_continuous(limits = c(0.8,1.1), expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.3),
          axis.ticks = element_line(color = 'black', size=0.3), axis.text.x = element_text(angle = 45, vjust = 0.5, color="black"),
          plot.title = element_text(size=10), plot.subtitle = element_text(size=8), axis.text.y = element_text(color="black", size=13)) +
    ggtitle(convertMetabolitesNames(metabolit),  subtitle = paste('ANOVA: ', round(statsMetaboAdjust[metabolit, 'anova'], 6),
                                                                  ' full-young: ', round(statsMetaboAdjust[metabolit, 'young-full'], 6),
                                                                  '\noldWT-oldHet: ', round(statsMetaboAdjust[metabolit, 'oldWT-oldHet'], 6),
                                                                  ' oldHet-full: ', round(statsMetaboAdjust[metabolit, 'oldHet-full'], 6),
                                                                  '\nyoung-oldHet: ', round(statsMetaboAdjust[metabolit, 'young-oldHet'], 6))) #oldWT  oldWTorFull
  
  boxes[[metabolit]] <- boxp
}


pdf(paste0(path, 'figures/boxplot sig metabolites ANOVA 4clusters.pdf'), width = 14, height = 13)
grid.arrange(boxes[[1]], boxes[[2]], boxes[[3]], boxes[[4]], boxes[[5]], boxes[[6]],
             boxes[[7]], ncol=3)
dev.off()

statsMetaboAdjust <- round(statsMetaboAdjust, 4)
statsMetabo <- round(statsMetabo, 4)
write.csv(statsMetaboAdjust, paste0(path, 'outputFiles/statistics Metabolites FDR.csv'))


# heatmap order by biocyc pathway-----
metaboDataAverage <- averageSamples(metaboDataBiocycNames, infoTable)
selectedMetabolites <- read.csv(paste0(path, 'inputFiles/selected metabolites and pathways pattern.csv'))
BioGly <- read.csv(paste0(path, 'inputFiles/selected metabolites for biosynthesis and glycolosis.csv'))

selectedMetabolites$subPathway <- factor(selectedMetabolites$subPathway, levels = unique(selectedMetabolites$subPathway))
selectedMetabolites$metabolite <- factor(selectedMetabolites$metabolite, levels = unique(selectedMetabolites$metabolite))


orderMandHeatmap <- function(metabolitesInGroups, metaboDataBiocycNames, infoTable,...){
  dataNumbersPathways <- merge(metabolitesInGroups, metaboDataBiocycNames, by.x = "metabolite", by.y = 0)
  dataNumbersPathways <- dataNumbersPathways[,c(metabColNames, row.names(infoTable))]
  
  dataNumbersPathways <- dataNumbersPathways[order(dataNumbersPathways$subPathway),]
  dataNumbersPathways <- dataNumbersPathways[order(dataNumbersPathways$mainGroup),]
  
  df <- as.matrix(dataNumbersPathways[, row.names(infoTable)])
  row.names(df) <- dataNumbersPathways$metabolite
  splitByM <- gsub('Degradation', 'Metabolism', dataNumbersPathways$subPathway)
  splitByM <- gsub('noPathway', 'Others', splitByM)
  hm = heatmapMetbolomics(df, infoTable, hmName = 'from biocyc', splitBy = splitByM, rowsClustering = T, cluster_row_slices = FALSE, ...)
  return(hm)
}

hm_sample <- orderMandHeatmap(selectedMetabolites, metaboDataBiocycNames, infoTable)
hm_average <- orderMandHeatmap(selectedMetabolites, metaboDataAverage$meanData, metaboDataAverage$infoTable, show_column_names = F)

hm_sampleBG <- orderMandHeatmap(BioGly, metaboDataBiocycNames, infoTable)
hm_averageBG <- orderMandHeatmap(BioGly, metaboDataAverage$meanData, metaboDataAverage$infoTable, show_column_names = F)

#save heatmaps
pdf(paste0(path, 'figures/SUP M heatmap patterned metabolites samples.pdf'), width = 12, height = 8)
draw(hm_sample)
dev.off()

pdf(paste0(path, 'figures/M heatmap patterned metabolites average.pdf'), width = 10, height = 8)
draw(hm_average)
dev.off()

pdf(paste0(path, 'figures/SUP M heatmap biosynthsis glycolosis samples.pdf'), width = 8, height = 6)
draw(hm_sampleBG)
dev.off()

pdf(paste0(path, 'figures/M heatmap biosynthsis glycolosis average.pdf'), width = 8, height = 6)
draw(hm_averageBG)
dev.off()
#----




# **************lipidomics- heatmaps and statistics************------
lipid_trend <- read.csv(paste0(path, 'inputFiles/lipids.csv'), row.names = 1)
lipid_groups <- factor(lipid_groups, levels = rownames(lipid_trend))

# heatmap with all lipids clustered in group. two heatmap- samples as is and average of the samples.
HM_sampleLipids <- heatmapLipids(lipidomic, infoTableL, "lipidomic", lipid_groups, cluster_row_slices = FALSE, show_row_names = FALSE)
average.lipid <- averageSamples(lipidomic, infoTableL)
HM_averageLipids <- heatmapLipids(average.lipid$meanData, average.lipid$infoTable, "lipidomic-average",
                                  lipid_groups, cluster_row_slices = FALSE, show_row_names = FALSE, show_column_names = FALSE,
                                  border = F, row_gap = unit(0.5, "mm"))

lipid_trend$lipid.family <- rownames(lipid_trend)
lipids <- c('LPC', 'LPE', 'PI', 'PA', 'PS', 'SM')
hm_sampleL <- heatmapLipids(lipidomic[lipid_groups %in% lipids,], infoTableL, "rescue lipids", lipid_groups[lipid_groups %in% lipids], cluster_row_slices = FALSE,
                            border=T, row_gap = unit(c(0.5, 0.5, 4, 0.5, 4, 4), "mm"))
hm_averageL <- heatmapLipids(average.lipid$meanData[lipid_groups %in% lipids,], average.lipid$infoTable, "lipids-average", lipid_groups[lipid_groups %in% lipids],
                             cluster_row_slices = FALSE, border=T, row_gap = unit(c(0.5, 0.5, 4, 0.5, 4, 4), "mm"), show_column_names = FALSE)


# save heatmaps lipids----
pdf(paste0(path, 'figures/SUP L heatmap all lipids samples.pdf'), width = 13, height = 11)
draw(HM_sampleLipids)
dev.off()

pdf(paste0(path, 'figures/L heatmap all lipids average.pdf'), width = 13, height = 11)
draw(HM_averageLipids)
dev.off()

pdf(paste0(path, 'figures/SUP L heatmap selected families samples.pdf'), width = 10, height = 8)
draw(hm_sampleL)
dev.off()

pdf(paste0(path, 'figures/L heatmap selected families average.pdf'), width = 10, height = 8)
draw(hm_averageL)
dev.off()

# statistics for lipids----
statsLipids <- data.frame(row.names = row.names(lipidomic),
                          no = rep(0, nrow(lipidomic)))

for (i in 1:nrow(lipidomic)){
  df <- data.frame(m = unlist(lipidomic[i,]), feed=infoTableL$feed, age=infoTableL$age, genotype=infoTableL$genotype)
  pvals <- anovaFULL(df)
  statsLipids[i, names(pvals)] <- pvals
  pvals <- anova4groups.tukey(df[df$feed=='fasted',])
  statsLipids[i, names(pvals)] <- pvals
}

statsLipidsAdjust <- data.frame(statsLipids)
colnames(statsLipidsAdjust) <- colnames(statsLipids)
for (a in colnames(statsLipidsAdjust)){
  statsLipidsAdjust[[a]] <- p.adjust(statsLipids[[a]], method = 'fdr')
}


boxesL <- list()
for (lipid in rownames(statsLipidsAdjust[statsLipidsAdjust$anova < 0.1,])){
  df <- data.frame(m= unlist(lipidomic[lipid,]), genotype=infoTableL$genotype, feed=infoTableL$feed,
                   group=factor(gsub(' male', '', infoTableL$group), levels = c("6.5weeks WT full", "6.5weeks Het full", "15weeks WT full", "15weeks Het full", "6.5weeks WT fasted", "6.5weeks Het fasted", "15weeks WT fasted", "15weeks Het fasted")))
  df <- df[df$feed == 'fasted', ]
  
  boxp <- ggplot(df, aes(x=group, y=m)) +
    geom_boxplot(color='black', fill = 'white', width=0.5, lwd=0.4, fatten = 0.8, outlier.shape = NA, coef=10) +
    stat_summary(fun = max, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., group = group),
                 width = 0.2, linetype = "solid", size=0.7)+
    stat_summary(fun = min, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., group = group),
                 width = 0.2, linetype = "solid", size=0.7)+
    geom_jitter(width=0.1, size=2.2, aes(color=genotype)) +
    scale_color_manual(values = c('black', '#098036'))+
    xlab("") + ylab("") +
    theme_classic() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.3),
          axis.ticks = element_line(color = 'black', size=0.3), axis.text.x = element_text(angle = 45, vjust = 0.5, color="black"),
          plot.title = element_text(size=10), plot.subtitle = element_text(size=8), axis.text.y = element_text(color="black", size=13)) +
    ggtitle(lipid,  subtitle = paste('ANOVA: ', round(statsLipidsAdjust[lipid, 'anova'], 6),
                                     ' oldWT-oldHet: ', round(statsLipidsAdjust[lipid, 'oldWT-oldHet'], 6),
                                     '\nyoung-oldWT: ', round(statsLipidsAdjust[lipid, 'young-oldWT'], 6),
                                     '\nyoung-oldHet: ', round(statsLipidsAdjust[lipid, 'young-oldHet'], 6)))
  
  boxesL[[lipid]] <- boxp
}


pdf(paste0(path, 'figures/lipids boxplot fasted.pdf'), width = 9, height = 11)
grid.arrange(boxesL$`SM(d43:1)`, boxesL$`SM(d35:1)`, boxesL$`SM(t38:1)`, boxesL$`PS(18:0_22:6) (1)`,
             boxesL$`PI(17:0_20:4)`, boxesL$`SM(d41:1)`, boxesL$`LPE(16:0)`, boxesL$`LPC(18:2) (1)`, ncol=3)
dev.off()
 
statsLipidsAdjust <- round(statsLipidsAdjust, 4)
statsLipids <- round(statsLipids, 4)
write.csv(statsLipidsAdjust, paste0(path, 'outputFiles/statistics lipids FDR.csv'))

