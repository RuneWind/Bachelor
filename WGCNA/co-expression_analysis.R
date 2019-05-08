print("********** STARTING WGCNA ***********")
library("WGCNA")
library("dynamicTreeCut")

enableWGCNAThreads(nThreads = 16)
getwd()
setwd("/Users/Rune/Documents/Skole/Uni/S6/Bachelor/WGCNA")
print("\\\\\\ Loading data //////")
datExpr = read.delim(file="datExpr.txt", header=TRUE, sep="\t") # Filtered with goodSamplesGenes
print(paste("Working on", as.character(nrow(datExpr)), "samples and", as.character(ncol(datExpr)), "genes"))


SubGeneNames = colnames(datExpr)
print("\\\\\\ Data loaded //////")

# Choosing soft-threshold to fit a scale-free topology to the network
print("\\\\\\ Preparing thresholds //////")
powers = c(c(1:10), seq(from=12, to=20, by=2))
sft = pickSoftThreshold(datExpr, dataIsExpr=TRUE, powerVector=powers, corFnc=cor, corOptions=list(use="p"), networkType="unsigned")

# Plot thresholds
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

print("\\\\\\ Plotting thresholds //////")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n", main=paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
par(mfrow=c(1,1))
print("\\\\\\ Threshold calculations and plooting done //////")

softPower = 3 # The chosen softpower
# Calculate adjacency matrix
print("\\\\\\ Caltulating adjacency matrix //////")
adj = adjacency(datExpr, type="unsigned", power=softPower)
print("\\\\\\ Adjacency matrix done! //////")

#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
print("\\\\\\ Building TOM //////")
TOM = TOMsimilarityFromExpr(datExpr, networkType="unsigned", TOMType="unsigned", power=softPower)
colnames(TOM) = rownames(TOM) = SubGeneNames
dissTOM = 1-TOM
print("\\\\\\ TOM + dissTOM finished //////")

#hierarchical clustering of the genes based on the TOM dissimilarity measure
print("\\\\\\ Hierarchial clustering //////")
geneTree = hclust(as.dist(dissTOM),method="average")
print("\\\\\\ Clustering done //////")

#plot the resulting clustering tree (dendrogram)
print("\\\\\\ Plotting dendrogram //////")
plot(geneTree, xlab="", sub="", cex=0.3)

# Module identification using dynamic tree cut
print("\\\\\\ Using cutTreeDynamic //////")
dynamicMods = cutreeDynamic(dendro=geneTree,  method="tree", minClusterSize=30);
#dynamicMods = cutreeDynamic(dendro=geneTree, distM=dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
print(paste(100 - (table(dynamicColors)[["grey"]]/ncol(datExpr))*100, "% (", table(dynamicColors != "grey")[["TRUE"]], ") of genes assigned to a module", sep=""))

print("\\\\\\ Plotting DendroAndColors //////")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05, main="Gene dendrogram and module colors")

#Visualize the Tom plot
print("\\\\\\ TOM plot //////")
plotTOM = dissTOM^7 # Brings out module structure
png(filename="full_TOM.png",type="cairo", height=960, width=960)
TOMplot(plotTOM,  dendro=geneTree, as.character(dynamicColors))
dev.off()
print("\\\\\\ Done plotting //////")

#discard the unassigned genes, and focus on the rest
print("\\\\\\ DendroAndColors for modules only //////")
restGenes= (dynamicColors != "grey")
diss1=1-TOMsimilarityFromExpr(datExpr[,restGenes], power = softPower)

colnames(diss1) =rownames(diss1) =SubGeneNames[restGenes]
hier1 = hclust(as.dist(diss1), method="average")
png(filename="dendroAndColors.png",type="cairo", height=960, width=960)
plotDendroAndColors(hier1, dynamicColors[restGenes], "Dynamic Tree Cut", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05, main="Gene dendrogram and module colors")
dev.off()
print("\\\\\\ Plotting done! //////")

#set the diagonal of the dissimilarity to NA 
diag(diss1) = NA;

#Visualize the Tom plot of modules only
print("\\\\\\ TOM plot of modules only //////")
plotDiss1 = diss1^4 # Brings out module structure
png(filename="module_TOM.png", type="cairo", height=960, width=960)
TOMplot(plotDiss1, dendro=hier1, as.character(dynamicColors[restGenes]))
dev.off()
print("\\\\\\ Done plotting! //////")

# Extract modules
print("\\\\\\ Writing table for genes with in module... //////")
module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  print(color)
  module=SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}
print("\\\\\\ Done writing tables! //////")

# Expression patterns of genes
print("\\\\\\ Heatmap of expression patterns //////")
module.order <- unlist(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))
m = log10(datExpr[,module.order])

png(filename="expressionHeatmap.png",type="cairo", height=960, width=960)
heatmap(t(m), Rowv=NA, Colv=NA, col=gray.colors(100), labRow=NA, scale="row", RowSideColors=dynamicColors[module.order])
dev.off()
print("\\\\\\ Heatmap done! //////")

# Quantify module similarity by eigengene correlation
print("\\\\\\ Eigengene similarity //////")
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
plotEigengeneNetworks(MEs, "", marDendro=c(0,3,1,5), marHeatmap=c(2,2,1,0), plotHeatmaps=FALSE)
print("\\\\\\ Eigengene similarity plot finished //////")

print("### WGCNA FINISHED ###")