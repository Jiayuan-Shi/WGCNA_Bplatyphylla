library(WGCNA)
library(dplyr)
library(data.table)

### Read_expression_matrix
x<-fread("expression_matrix.txt")
tx<-as.data.frame(t(x[,-1]))
colnames(tx)<-x$ID

### Read_phyno_matrix
y<-read.table("phyno.txt",row.names=1,header=T,comment.char = "",check.names=F)

### goodSamplesGenes
gsg<-goodSamplesGenes(tx, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(tx)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(tx)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  tx = tx[gsg$goodSamples, gsg$goodGenes]
}

###filter
n<-nrow(tx)
meanFPKM<-0.5
#constant_for_select
tx[n+1,]<-apply(tx[c(1:nrow(tx)),],2,mean)
tx<-tx[1:n,tx[n+1,] > meanFPKM]

### for meanFpkm in row n+1 and it must be above what you set--select meanFpkm>opt$meanFpkm(by rp)
fpkmSamples = rownames(tx)
traitSamples =rownames(y)
traitRows = match(fpkmSamples, traitSamples)
datTraits = y[traitRows,]

### Choose a set of soft thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(tx, powerVector = powers)

## Plot the results:
pdf(file="threase_hold.pdf",width=24,height=18)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.65,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=1,col="red")
dev.off()

###ajacency_matrix
myadjacency <- adjacency(tx, power = 14)
TOM = TOMsimilarity(myadjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");

###cluster__tree
pdf(file="1_Gene clustering on TOM-based dissimilarity.pdf",width=24,height=18)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
pdf(file="2_Dynamic Tree Cut.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
MEList = moduleEigengenes(tx, colors = dynamicColors)
MEs = MEList$eigengenes%>%orderMEs()

###Rename to moduleColors
moduleColors = dynamicColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
nGenes = ncol(tx)
nSamples = nrow(tx)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

###draw relation plot
pdf(file="3_Module-trait relationships.pdf",width=10,height=10)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

###calculate_MM&GS
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(tx, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
#names of those trait
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(tx, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

###for_export_module&phynotype_relation_plot
for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      #sizeGrWindow(7, 7)
      pdf(file=paste("4_", trait, "_", module,"_Module membership vs gene significance.pdf",sep=""),width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}
names(tx)
probes = names(tx)

##exportMMGS
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "10_GS_and_MM.xls",sep="\t",row.names=F)

##hubgene
module = "blue"
column = match(module, modNames)
moduleGenes = moduleColors==module
lightyellow_module<-as.data.frame(dimnames(data.frame(tx))[[2]][moduleGenes])
names(lightyellow_module)="genename"
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
lightyellow_MMGS<-as.data.frame(cbind(MM,GS))
rownames(lightyellow_MMGS)=lightyellow_module$genename
hub_b<-abs(lightyellow_MMGS$MM)>0.8&abs(lightyellow_MMGS$GS)>0.2

mymodule<-subset(lightyellow_MMGS, abs(lightyellow_MMGS$MM)>0.8&abs(lightyellow_MMGS$GS)>0.2)
write.table(mymodule, "hubgene_MMGS_blue_0802.txt",quote = F,col.names = F,sep = "\t")