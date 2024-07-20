rm(list = ls())
library(scRICA)
rds1 <- readRDS('~/Desktop/757_scRNA-seq/RDS/scRICA_demo.rds')
table(rds1$expCond1)
table(rds1$expCond2)
## -------------------------------------------------------------------------------------------------------------------- ##
#### 5.4.1 Cellular compositions summary and visualization
getClusterSummaryReplot(rds = rds1, resDir = 'misc/visualization_results', newAnnotation = F, expCondCheck = 'expCond1', expCondCheckFname = 'donors_summary')
getClusterSummaryReplot(rds = rds1, resDir = 'misc/visualization_results', newAnnotation = T, newAnnotationRscriptName = 'misc/anno_example.R', expCondCheck = 'expCond2', expCondCheckFname = 'donors_summary')
getClusterSummaryReplot(rds = rds1, resDir = 'misc/visualization_results', perPlotVertical = F, perPlotHeight = 3, perPlotWidth = 7, perRemake = T, fpRemake = F, stack =T,
                        expCondNameOrder = c('I', 'A', 'F'),
                        newAnnotation = F, expCondCheck = 'expCond2', expCondCheckFname = 'tissue_sep2')

getClusterSummaryReplot(rds = rds1, resDir = 'misc/visualization_results',
                        perPlotVertical = T, perPlotHeight = 3, perPlotWidth = 7, perRemake = T, fpRemake = F, stack =F,
                        expCondNameOrder = c('I', 'A', 'F'),
                        newAnnotation = F, expCondCheck = 'expCond2', expCondCheckFname = 'tissue_sep1')
## -------------------------------------------------------------------------------------------------------------------- ##
##### 5.4.2 feature genes dot-plot visualization
getGoiDotplot(resDir = 'misc/visualization_results', rds = rds1,
              newAnnotation = F,
              goiFname = 'misc/marker_genes_CC.xlsx',
              expCondCheck = 'expCond2', expCondCheckFname = 'dotplot_tissue',
              # expCondReorderLevels = paste(rep(names(table(Idents(rds1))), each=3), rep(c('I', 'A', 'F'), time = length(names(table(Idents(rds1))))), sep = '_'),
              dotPlotFnamePrefix = 'cc')

getGoiDotplot(resDir = 'misc/visualization_results', rds = rds1,
              newAnnotation = F,
              goiFname = 'misc/marker_genes_CC.xlsx',
              expCondCheck = 'expCond1', expCondCheckFname = 'dotplot_donor',
              expCond = c('D4', 'D5', 'D6'),
              geneTypeOrder = c('S phase', 'G2/M phase'),
              dotPlotFnamePrefix = 'cc_Donor456only',
              dotPlotMinExpCutoff = 0, dotPlotMaxExpCutoff = 4,
              dotPlotWidth = 26, dotPlotHeight = 6,
              legendPer = 0.08, genetypebarPer = 0.03,
              fontsize.x = 20, fontsize.y = 16,
              fontsize.legend1 = 10,
              geneTypeLegendOn = F
              )

getGoiDotplot(resDir = 'misc/visualization_results', rds = rds1,
              newAnnotation = F,
              goiFname = 'misc/marker_genes_CC.xlsx',
              expCondCheck = 'expCond2', expCondCheckFname = 'dotplot_donor',
              cellcluster = c('T/NK1', 'T/NK2', 'T/NK3'),
              geneTypeOrder = c('S phase', 'G2/M phase'),
              # expCondReorderLevels = paste(rep(names(table(Idents(rds1))), each=3), rep(c('I', 'A', 'F'), time = length(names(table(Idents(rds1))))), sep = '_'),
              expCondReorderLevels = paste(rep(rev(c('T/NK1', 'T/NK2', 'T/NK3')), each=3), rep(rev(c('I', 'A', 'F')), time = 3), sep = '_'),
              # expCond = c('I', 'A'),
              # expCondReorderLevels = paste(rep(rev(c('T/NK1', 'T/NK2', 'T/NK3')), each=2), rep(rev(c('I', 'A')), time = 3), sep = '_'),
              dotPlotFnamePrefix = 'cc_TNKonly',
              # dotPlotFnamePrefix = 'cc_TNKonly_IAonly',
              dotPlotMinExpCutoff = 0, dotPlotMaxExpCutoff = 4,
              dotPlotWidth = 26, dotPlotHeight = 6,
              legendPer = 0.08, genetypebarPer = 0.03,
              fontsize.x = 20, fontsize.y = 16,
              fontangle.x = 90, fontangle.y = 0,
              fontsize.legend1 = 10,
              geneTypeLegendOn = T, fontsize.legend2 = 20
)

getGoiDotplot(resDir = 'misc/visualization_results', rds = rds1,
              newAnnotation = F,
              goiFname = 'misc/marker_genes_CC.xlsx',
              expCondCheck = 'expCond2', expCondCheckFname = 'dotplot_tissue',
              geneTypeOrder = c('S phase', 'G2/M phase'),
              cellcluster = c('T/NK1', 'T/NK2', 'T/NK3'),
              expCond = c('I', 'A'),
              expCondReorderLevels = paste(rep(rev(c('T/NK1', 'T/NK2', 'T/NK3')), each=2), rep(rev(c('I', 'A')), time = 3), sep = '_'),
              dotPlotFnamePrefix = 'cc_TNKonly_IAonly',
              dotPlotMinExpCutoff = 0, dotPlotMaxExpCutoff = 4,
              dotPlotWidth = 26, dotPlotHeight = 6,
              legendPer = 0.08, genetypebarPer = 0.03,
              fontsize.x = 20, fontsize.y = 16,
              fontangle.x = 90, fontangle.y = 0,
              fontsize.legend1 = 10,
              geneTypeLegendOn = T, fontsize.legend2 = 20)

## -------------------------------------------------------------------------------------------------------------------- ##
##### 5.4.3 heatmap visualization
getGoiHeatmap(heatmap.view = 'idents', resDir = 'misc/visualization_results', rds = rds1,
              newAnnotation = F, geneNames = c('CD44', 'ITGA6', 'EPCAM', 'KRT5', 'KRT7', 'TUBB4A', 'TUBB4B', 'PAX8'),
              expCondCheck = 'expCond2',
              expCondCheckFname = 'heatmap_tissue',
              plotFnamePrefix = 'geneset1_SEs_cells',
              display = 'cell')

getGoiHeatmap(heatmap.view = 'idents', resDir = 'misc/visualization_results', rds = rds1,
              newAnnotation = F, geneNames = c('CD44', 'ITGA6', 'EPCAM', 'KRT5', 'KRT7', 'TUBB4A', 'TUBB4B', 'PAX8'),
              expCondCheck = 'expCond1',
              expCondCheckFname = 'heatmap_tissue',
              plotFnamePrefix = 'geneset1_SEs_cells',
              display = 'expCond.merge', debug = T)

getGoiHeatmap(heatmap.view = 'expCond.idents', resDir = 'misc/visualization_results', rds = rds1,
              newAnnotation = F, geneNames = c('CD44', 'ITGA6', 'EPCAM', 'KRT5', 'KRT7', 'TUBB4A', 'TUBB4B', 'PAX8'),
              expCondCheck = 'expCond2',
              expCondCheckFname = 'heatmap_tissue',
              plotFnamePrefix = 'geneset1_expCond_idents_merge_SECEonly',
              expCondReorderLevels = c('I CE', 'A CE', 'F CE', 'I SE', 'A SE', 'F SE'),
              display = 'expCond.merge', draw.lines = F,
              cellcluster = c('CE', 'SE'),
              plotWidth = 6, plotHeight = 4, fontangle.top = 90, debug = T)

getGoiHeatmap(heatmap.view = 'expCond.idents', resDir = 'misc/visualization_results', rds = rds1,
              newAnnotation = F, geneNames = c('CD44', 'ITGA6', 'EPCAM', 'KRT5', 'KRT7', 'TUBB4A', 'TUBB4B', 'PAX8'),
              expCondCheck = 'expCond2',
              expCondCheckFname = 'heatmap_tissue',
              plotFnamePrefix = 'geneset1_expCond_idents_merge_SECEonly_IAonly',
              expCondReorderLevels = c('I CE', 'A CE', 'I SE', 'A SE'),
              display = 'sample.merge', draw.lines = F,
              cellcluster = c('CE', 'SE'),
              expCond = c('I', 'A'),
              fontsize.top = 2, fontsize.y = 12, barHeight.top = 0.05,
              plotWidth = 5, plotHeight = 3, fontangle.top = 90, debug = T)

getGoiHeatmap(heatmap.view = 'expCond', resDir = 'misc/visualization_results', rds = rds1,
              newAnnotation = F, geneNames = c('CD44', 'ITGA6', 'EPCAM', 'KRT5', 'KRT7', 'TUBB4A', 'TUBB4B', 'PAX8'),
              expCondCheck = 'expCond2',
              expCondCheckFname = 'heatmap_tissue',
              plotFnamePrefix = 'geneset1_expCond_idents_merge_SECEmerge2',
              expCondReorderLevels = c('I', 'A', 'F'),
              display = 'expCond.merge', draw.lines = F,
              cellcluster = c('CE', 'SE'),
              plotWidth = 6, plotHeight = 4, fontangle.top = 90, debug = T)

getGoiHeatmap(heatmap.view = 'idents', resDir = 'misc/visualization_results', rds = rds1,
              newAnnotation = F, geneNames = c('CD44', 'ITGA6', 'EPCAM', 'KRT5', 'KRT7', 'TUBB4A', 'TUBB4B', 'PAX8'),
              expCondCheck = 'expCond2',
              expCondCheckFname = 'heatmap_tissue',
              plotFnamePrefix = 'geneset1_sampleMerge',
              display = 'sample.merge')

getGoiHeatmap(heatmap.view = 'idents', resDir = 'misc/visualization_results', rds = rds1,
              newAnnotation = F, geneNames = c('CD44', 'ITGA6', 'EPCAM', 'KRT5', 'KRT7', 'TUBB4A', 'TUBB4B', 'PAX8'),
              expCondCheck = 'expCond2',
              expCondCheckFname = 'heatmap_tissue',
              cellcluster = c('SE', 'CE', 'ST1', 'ST5'),
              expCondReorderLevels = c('SE', 'CE', 'ST1', 'ST5'),
              plotFnamePrefix = 'geneset1_sampleMerge2',
              display = 'sample.merge', debug = T)

## -------------------------------------------------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------------------------------------------------- ##
##### 5.5.1 DE analysis on cells
de.res1 <- getClusterExpCondDe(resDir = 'misc/', rds = rds1,
                               newAnnotation = F, cellcluster = 'MP', deMethod = 'wilcox',
                               expCondCheck = 'expCond2', expCondCheckFname = 'IFcomp_wilcox', compGroup = 'I/F')

de.res2 <- getClusterExpCondDe(resDir = 'misc/', rds = rds1,
                               newAnnotation = F, deMethod = 'MAST',
                               outputExcel = as.logical(F),
                               expCondCheck = 'expCond2', expCondCheckFname = 'IFcomp_allCells_MAST',
                               compGroup = 'I/F')

de.res3 <- getClusterExpCondDe(resDir = 'misc/', rds = rds1,
                               cellcluster = 'MP',
                               newAnnotation = F, deMethod = 'MAST',
                               outputExcel = as.logical(F),
                               expCondCheck = 'expCond2', expCondCheckFname = 'IFcomp_MAST_MP',
                               compGroup = 'I/F')

source('/Users/yanli/Desktop/scRIC_development/scRICA_fnPseudotime_exploration/twee.R')
twee(path = './misc/results_wOrgClusterAnnotation_DEGs')

head(de.res1$clusterDeMarkers$MP)
de.res1$clusterDeResSummary
head(de.res1$clusterTopDeMarkers$MP)
##### 5.5.1 DE analysis on pseudo-bulk counts
meta1 <- data.frame('bulksamp' = names(table(rds1$expCond)),
                    'gender' = c('F', 'F', 'F', 'M', 'M', 'M', 'M', 'M', 'F', 'F', 'M', 'M', 'M', 'M', 'M', 'M', 'F', 'F' ),
                    'bmi' = c(23.5, 23.5, 23.5, 18, 18, 28.4, 28.4, 28.4, 32, 32, 20, 20, 20, 30, 30, 30, 21, 21))
de.res3 <- getClusterExpCondDe(resDir = 'misc/', rds = rds1,
                               newAnnotation = F, cellcluster = 'MP',
                               deMethod = 'DESeq2.bulk',
                               deseq2bulk.metaCovariateInput = meta1,
                               expCondCheck = 'expCond2', expCondCheckFname = 'IFcomp_DEseq2bulk',
                               compGroup = 'I/F')
## -------------------------------------------------------------------------------------------------------------------- ##
##### 5.5.2 Over expressed genes
oe.res1 <- getExpCondClusterMarkers(resDir = 'misc/', rds = rds1,
                                    deMethod = 'MAST',
                                    cellcluster = c('T/NK1', 'T/NK2', 'T/NK3', 'MP', 'B/P'),
                                    newAnnotation = F,
                                    expCondCheck = 'expCond2',
                                    expCondCheckFname = 'tissue_IR_overExpressedGenes')

twee(path = './misc/')

head(oe.res1$fullRes$A)
head(oe.res1$sigUp$A)
table(oe.res1$sigUp$A$cluster)
## -------------------------------------------------------------------------------------------------------------------- ##
##### 5.5.3 HIPPO cells sub-clustering analysis
## hippo clustering
getHippoRes(resDir = 'misc/', rds = rds1, newAnnotation = F,
            expCondCheck = 'expCond2', expCond = 'F',
            hippoResNamePrefix = 'Fonly_TNK23',
            cellcluster = c('T/NK2', 'T/NK3'), noClusters = 5)

getHippoRes(resDir = 'misc/', rds = rds1, newAnnotation = F,
            hippoResNamePrefix = 'Allcells_TNK',
            cellcluster = c('T/NK1', 'T/NK2', 'T/NK3'), noClusters = 4)
#
# getHippoRes(resDir = 'misc/', rds = rds1, newAnnotation = F,
#             hippoResNamePrefix = 'Allcells_TNK23',
#             cellcluster = c('T/NK2', 'T/NK3'), noClusters = 5)
#
# getHippoRes(resDir = 'misc/', rds = rds1, newAnnotation = F,
#             hippoResNamePrefix = 'Allcells_TNK23_222',
#             cellcluster = c('T/NK2', 'T/NK3'), noClusters = 5)
#
# getHippoRes(resDir = 'misc/', rds = rds1, newAnnotation = F,
#             hippoResNamePrefix = 'Allcells_TNK2',
#             cellcluster = c('T/NK2'), noClusters = 3)

twee(path = './misc/')
## update cell clustering
rds.hippo2 <- updateHippoIdents(resDir = 'misc/', rds = rds1, resSaveFname = 'hippo_tnk.Rds',
                               hippoResList = list('TNK_hippo' = 'misc/hippo_results/Allcells_TNK_k4/lightHIPPOres_Allcells_TNK_k4.Rdata'),
                               hippoResK = c(4),
                               newAnnotation = F )
table(Idents(rds.hippo2))
## -------------------------------------------------------------------------------------------------------------------- ##
##### 5.5.4 pseudo-time trajectory analysis with slingshot
?getExpCondClusterPseudotime
psdt.res1 <- getExpCondClusterPseudotime(resDir = 'misc/visualization_results', rds = rds1,
                                         expCondCheck = 'expCond2', expCondCheckFname = 'tissueSep_SE_Aonly',expCond = c('A'),
                                         cellcluster = c('SE'),
                                         newAnnotation = F)

psdt.res2 <- getExpCondClusterPseudotime(resDir = 'misc/visualization_results', rds = rds1,
                                         expCondCheck = 'expCond2', expCondCheckFname = 'tissueSep_en1pv1',
                                         cellcluster = c('EN1', 'P/V1'),
                                         newAnnotation = F)

twee(path = './misc/visualization_results/results_wOrgClusterAnnotation_pseudoTime/')
## -------------------------------------------------------------------------------------------------------------------- ##


