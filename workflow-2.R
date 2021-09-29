
#set working directory
setwd("C:/Users/dedePC/Desktop/patway_analysis_IBD_data")

#load required libraries
library(dplyr)
library(RColorBrewer)
library(org.Hs.eg.db)
library(tidyverse)
library(EnhancedVolcano)
library(VennDiagram)
library(rWikiPathways)
library(data.table)
library(RCy3)

# Differential expression data visualization
# We will use a publicly available dataset, which identified two different 
# subtypes of IBD - CD and UC. DEG analysis was performed before then we will use the results file from prevous workflow
# First, let's import the data and use a volcano plot to visualize the result 
# of the differential gene expression analysis result, and use a Venn diagram 
# to study how many differentially expressed genes are shared between the subtypes.

#CD read, the datasets may change depends on our observation (ileum or rectum)
dataset.CD <- read.delim("data/table_CD_Rectum_vs_nonIBD_Rectum.tab")
#UC read
dataset.UC <- read.delim("data/table_UC_Rectum_vs_nonIBD_Rectum.tab")

#get entrez ID of each gene symbols for each disease type
hs <- org.Hs.eg.db
entrezID <- AnnotationDbi::select(hs, 
            keys = dataset.CD$X,
            columns = c("ENTREZID", "SYMBOL"),
            keytype = "SYMBOL")

#filter out double gene symbols
entrezID<- entrezID %>% distinct(entrezID$SYMBOL, .keep_all = TRUE)

#get back up data before starting
dataset.b  <- dataset.CD
dataset.b2 <- dataset.UC

# add entrezIDs to each dataset
dataset.CD <- cbind(entrezID$ENTREZID,dataset.CD)
dataset.UC <- cbind(entrezID$ENTREZID,dataset.UC)

#change column names
colnames(dataset.CD)[1] <- "ENTREZ.ID"
colnames(dataset.CD)[2] <- "SYMBOL"
colnames(dataset.UC)[1] <- "ENTREZ.ID"
colnames(dataset.UC)[2] <- "SYMBOL"

# filter genes without Entrez Gene identifier
dataset.CD <- dataset.CD %>% tidyr::drop_na(ENTREZ.ID)
dataset.UC <- dataset.UC %>% tidyr::drop_na(ENTREZ.ID)
#remove some unused columns
dataset.CD <- subset( dataset.CD, select = -c(3,5,6,7,9 ) )
dataset.UC <- subset( dataset.UC, select = -c(3,5,6,7,9 ) )

#output folder should be created beforehand run the code
png('output/volcanoplot_CD_rectum.png')
EnhancedVolcano(dataset.CD , title = "Chrons Disease on rectum", lab = dataset.CD$SYMBOL, 
                labSize = 3, x = 'log2FoldChange', y = 'pvalue', pCutoff = 0.05, FCcutoff = 1)
dev.off()

png('output/volcanoplot_UC_rectum.png')
EnhancedVolcano(dataset.UC, title = "Ulcerative Colitis Disease on rectum", lab = dataset.UC$SYMBOL, 
                labSize = 3, x = 'log2FoldChange', y = 'pvalue', pCutoff = 0.05, FCcutoff = 1)
dev.off()

#list of all deg from two disease types 
deg.CD  <-unique(dataset.CD[!is.na(dataset.CD$pvalue) & dataset.CD$pvalue < 0.05 & abs(dataset.CD$log2FoldChange) > 1,c(1,2)])
CD.up   <-unique(dataset.CD[!is.na(dataset.CD$pvalue) & dataset.CD$pvalue < 0.05 & dataset.CD$log2FoldChange > 1,c(1,2)])
CD.down <-unique(dataset.CD[!is.na(dataset.CD$pvalue) & dataset.CD$pvalue < 0.05 & dataset.CD$log2FoldChange < -1,c(1,2)])

deg.UC  <-unique(dataset.UC[!is.na(dataset.UC$pvalue) & dataset.UC$pvalue < 0.05 & abs(dataset.UC$log2FoldChange) > 1,c(1,2)])
UC.up   <-unique(dataset.UC[!is.na(dataset.UC$pvalue) & dataset.UC$pvalue < 0.05 & dataset.UC$log2FoldChange > 1,c(1,2)])
UC.down <-unique(dataset.UC[!is.na(dataset.UC$pvalue) & dataset.UC$pvalue < 0.05 & dataset.UC$log2FoldChange < -1,c(1,2)])


venn.diagram(x = list(CD.up$ENTREZ.ID, CD.down$ENTREZ.ID, UC.up$ENTREZ.ID, UC.down$ENTREZ.ID),
             category.names = c("CD up", "CD down" ,"UC up", "UC down"),
             filename = 'output/venn_genes.png',
             output=FALSE,
             col=c("#440154ff","#440154ff", '#21908dff','#21908dff'),
             fill = c(alpha("#440154ff",0.3),alpha("#440154ff",0.3), alpha('#21908dff',0.3),alpha('#21908dff',0.3)),
             cex = 1.5,
)

## Pathway enrichment analysis
#We will perform pathway enrichment with the gene sets of all pathway models in WikiPathways (human only).
gmt <- rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens",format = "gmt")
wp2gene   <- readPathwayGMT(gmt)
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
bkgd.genes <- unique(dataset.CD[,c(1,2)])# genes are commonn for both disease type CD and UC

# The clusterProfiler R-package is used to perform overrepresentation analysis (ORA). 
# The function can be easily replaced to use other enrichment methods (GSEA / rSEA / etc). 
# We will run the analysis separately for CD and UC subtype.

##################################FOR CD DISEASE##############################################
ewp.CD <- clusterProfiler::enricher(
  deg.CD$ENTREZ.ID,
  universe = bkgd.genes$ENTREZ.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.02,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)
ewp.CD.res <- as.data.frame(ewp.CD) 

# number of genes measured in pathways
length(ewp.CD@universe)
# number of DEG in pathways
length(deg.CD$ENTREZ.ID[deg.CD$ENTREZ.ID %in% unique(wp2gene$gene)])
num.pathways.CD <- dim(ewp.CD.res)[1]

# export enrichment result
png('output/CD_barplot.png', width = 1200, height=1000)
ggplot(ewp.CD[1:num.pathways.CD], aes(x=reorder(Description, -pvalue), y=Count)) +
  geom_bar(stat ="identity", fill="#BA8CD7") +
  coord_flip() +
  labs(x="", y="CD DEG gene count", fill="") +
  theme(axis.text=element_text(size=25)) + 
  theme(legend.position="none")
dev.off()
write.table(ewp.CD.res, file="output/CD_enrich_res.txt", sep="\t", quote=FALSE, row.names = FALSE)

# > Interpretation
# - **Q5**: How many pathways are altered in the CD subtype and how do they link to 
#IBD  (expected or unexpected)?


############################FOR UC DISEASE#################################
ewp.UC <- clusterProfiler::enricher(
  deg.UC$ENTREZ.ID,
  universe = bkgd.genes$ENTREZ.ID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.02,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)
ewp.UC.res <- as.data.frame(ewp.UC) 

# number of genes measured in pathways
length(ewp.UC@universe)
# number of DEG in pathways
length(deg.UC$ENTREZ.ID[deg.UC$ENTREZ.ID %in% unique(wp2gene$gene)])
num.pathways.UC <- dim(ewp.UC.res)[1]

# export enrichment result
png('output/UC_barplot.png', width = 1200, height=1000)
ggplot(ewp.UC[1:num.pathways.UC], aes(x=reorder(Description, -pvalue), y=Count)) +
  geom_bar(stat ="identity", fill="#BA8CD7") +
  coord_flip() +
  labs(x="", y="UC DEG gene count", fill="") +
  theme(axis.text=element_text(size=25)) + 
  theme(legend.position="none")
dev.off()
write.table(ewp.UC.res, file="output/UC_enrich_res.txt", sep="\t", quote=FALSE, row.names = FALSE)

venn.diagram(x = list(ewp.CD.res$ID, ewp.UC.res$ID),
             category.names = c("CD" , "UC"),
             filename = 'output/venn_pathways.png',
             output=TRUE,
             col=c("#440154ff", '#21908dff'),
             fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
             cex = 1.5,
)

#to find common altered pathways in both subtypes intersection of two results
common.pathways <- intersect(ewp.CD.res$ID, ewp.UC.res$ID) 
res <- ewp.CD.res[common.pathways,]
write.table(res, file='output/common_altered-pathways.txt', sep = "\t", quote = FALSE, row.names = FALSE)

#############Pathway visualization##############
# The pathways can then be visualized with the gene expression data as shown with the 
# "Overview of proinflammatory and profibrotic mediators" (WP5095) pathway from WikiPathways. 
# The pathway was altered in both subtypes. 
# interest and visualize the data on that pathway. 

#before starting we should merge data table of CD and UC disease 
data.IBD <- merge(dataset.CD, dataset.UC, by = "ENTREZ.ID")
colnames(data.IBD)[3] <- "logFC_CD"
colnames(data.IBD)[6] <- "logFC_UC"
colnames(data.IBD)[4] <- "pvalue_CD"
colnames(data.IBD)[7] <- "pvalue_UC"

RCy3::cytoscapePing()
RCy3::installApp(c("wikipathways","CyTargetLinker"))

RCy3::commandsRun('wikipathways import-as-pathway id=WP5095') 
toggleGraphicsDetails()
RCy3::mapTableColumn("Ensembl", "Homo sapiens", "Ensembl", "Entrez Gene")
loadTableData(data.IBD, data.key.column = "ENTREZ.ID", table.key.column = "Entrez Gene")

RCy3::installApp("enhancedGraphics")
RCy3::copyVisualStyle("WikiPathways", "my_style_heatmap")

RCy3::setNodeCustomHeatMapChart(c("logFC_CD","logFC_UC"), slot = 2, style.name = "my_style_heatmap", 
                                colors = c("#CC3300","#FFFFFF","#6699FF","#CCCCCC"))

RCy3::setVisualStyle("my_style_heatmap")

#Saving output
png.file <- file.path("output/PathwayVisualization.png")
exportImage(png.file,'PNG', zoom = 500)
cys.file <- file.path(getwd(),"output/PathwayVisualization.cys")
saveSession(cys.file)
#comment following line if you want to manipulate the visualization in Cytoscape
RCy3::closeSession(save.before.closing = F)


##########################Pathway overlap visualization#################################
# There is often crosstalk and overlap between pathways enriched in gene expression analyses. 
# The following step visualizes the overlap between the enriched pathways in a pathway-gene network. 
# The genes not present in any pathway are included in the visualization but can be filtered 
#in a follow-up step if preferred. 

pwy <- unique(ewp.CD.res[,c(1,2)])#enriched pathways in CD 
colnames(pwy) <- c("id","label")
pwy$type <- 'pathway'
edges <- wpid2gene[wpid2gene$wpid %in% pwy$id,]
colnames(edges) <- c("source", "target")
genes <- unique(deg.CD)
colnames(genes) <- c("id","label")
genes <- transform(genes, id = as.character(id))
genes$type <- 'gene'
edges <- unique(edges[edges$target %in% genes$id,])
genes <- genes[genes$id %in% edges$target,]
nodes <- dplyr::bind_rows(genes, pwy)
rownames(nodes) <- NULL
createNetworkFromDataFrames(nodes=nodes,edges=edges,title="Pathway-Gene-Associations", collection="PathwayGeneCrosstalk")
loadTableData(data.IBD, data.key.column = "ENTREZ.ID", table.key.column = "id")

# Visual style
RCy3::copyVisualStyle("default","wp.vis")
RCy3::setNodeLabelMapping("label", style.name="wp.vis")
RCy3::lockNodeDimensions(TRUE, style.name="wp.vis")
RCy3::setNodeShapeMapping('type', c('gene','pathway'), c("ellipse","hexagon"), style.name="wp.vis")
RCy3::setNodeSizeMapping('type', c('gene','pathway'), c(40,25), mapping.type = "d", style.name = "wp.vis")
data.values<-c(-1,0,1) 
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("logFC_CD", data.values, node.colors, default.color = "#99FF99", style.name = "wp.vis")
RCy3::setVisualStyle("wp.vis")
RCy3::toggleGraphicsDetails()

# Saving output
svg.file <- file.path(getwd(), "output/PathwayCrosstalk.svg")
exportImage(svg.file,'SVG')
png.file <- file.path(getwd(), "output/PathwayCrosstalk.png")
exportImage(png.file,'PNG', zoom = 500)
cys.file <- file.path(getwd(), "output/PathwayCrosstalk.cys")
saveSession(cys.file) 
#comment following line if you want to manipulate the visualization in Cytoscape
#RCy3::closeSession(save.before.closing = F)

##############################Drug target information#####################################
# Next, we will add information about known drug-target interactions for the genes 
# in the affected pathways using information from DrugBank using the CyTargetLinker app.
# We will show this for the UC subtype. 

RCy3::cytoscapePing()
installApp('CyTargetLinker') 

pwy <- unique(ewp.CD.res[,c(1,2)])
colnames(pwy) <- c("id","label")
pwy$type <- 'pathway'
edges <- wpid2gene[wpid2gene$wpid %in% pwy$id,]
colnames(edges) <- c("source", "target")
genes <- unique(deg.CD)
colnames(genes) <- c("id","label")
genes <- transform(genes, id = as.character(id))
genes$type <- 'gene'
edges <- unique(edges[edges$target %in% genes$id,])
genes <- genes[genes$id %in% edges$target,]
nodes <- dplyr::bind_rows(genes, pwy)
rownames(nodes) <- NULL
createNetworkFromDataFrames(nodes=nodes,edges=edges,title="Pathway-Gene-Associations", collection="PathwayGeneCrosstalk")
loadTableData(data.IBD, data.key.column = "ENTREZ.ID", table.key.column = "id")

# Visual style
RCy3::copyVisualStyle("default","wp.vis")
RCy3::setNodeLabelMapping("label", style.name="wp.vis")
RCy3::lockNodeDimensions(TRUE, style.name="wp.vis")
RCy3::setNodeShapeMapping('type', c('gene','pathway'), c("ellipse","hexagon"), style.name="wp.vis")
RCy3::setNodeSizeMapping('type', c('gene','pathway'), c(40,25), mapping.type = "d", style.name = "wp.vis")
data.values<-c(-1,0,1) 
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("logFC_UC", data.values, node.colors, default.color = "#99FF99", style.name = "wp.vis")
RCy3::setVisualStyle("wp.vis")
RCy3::toggleGraphicsDetails()

drugbank <- file.path(getwd(), "data/drugbank-5.1.0.xgmml")

# run CyTargetLinker
commandsRun(paste0('cytargetlinker extend idAttribute="id" linkSetFiles="', drugbank, '"'))
commandsRun('cytargetlinker applyLayout network="current"')
RCy3::setVisualStyle("wp.vis")

#let's change the visualization of the drugs in the network using the ByPass option
selected <- RCy3::selectNodes(nodes="drug", by.col = "CTL.Type")
RCy3::setNodeShapeBypass(node.names = selected$nodes, new.shapes = "Triangle")
RCy3::setNodeColorBypass(node.names = selected$nodes, "#FFFFCE")
RCy3::setNodeBorderColorBypass(node.names = selected$nodes, "#000000")
RCy3::setNodeBorderWidthBypass(node.names = selected$nodes, 4)
RCy3::clearSelection()
RCy3::toggleGraphicsDetails()

png.file <- file.path(getwd(), "output/drug_target.png")
exportImage(png.file,'PNG', zoom = 500)






