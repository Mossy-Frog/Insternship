library(WGCNA)


DfComplete<-read.table("C:/Users/maxence.lechemia/Documents/Multi-omics/normalized counts by deseq2 with proteo without outliers.txt")



DfComplete<-t(DfComplete)
gsg <- goodSamplesGenes(datExpr = DfComplete, verbose = 3)


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(DfComplete)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(DfComplete)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  DfComplete = DfComplete[gsg$goodSamples, gsg$goodGenes]
}






sampleTree = hclust(dist(DfComplete), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)



powers = c(c(1:10), seq(from = 12, to=50, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(DfComplete, powerVector = powers, verbose = 5)



#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")


plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

DfComplete<-as.data.frame(DfComplete)
for ( i in 1:length(DfComplete))
{
  for (j in 1:length(DfComplete[,i]))
  {
    DfComplete[j,i]<-as.numeric(DfComplete[j,i])
  }
}

#Plateau reached at around power 50, but power is between 1 and 30 so 30
# DfComplete[] <- lapply(DfComplete, as.numeric)
# DfComplete<-data.matrix(DfComplete)
net = blockwiseModules(DfComplete, power = 16,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
# open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "MultiOmicWGCNAtranscriproteo.RData")




#############test 2.b part wgcna#########

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(DfComplete, powerVector = powers, verbose = 5)
#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.65,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



softPower = 9;
adjacency = adjacency(DfComplete, power = softPower)

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


# Calculate eigengenes
MEList = moduleEigengenes(DfComplete, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
#sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.1
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(DfComplete, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs


#sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "C:/Users/maxence.lechemia/Documents/scripts formalisés/multiomicGTPwithoutoutliers.RData")





####################Part 2 C of WGCNA###########################
DfComplete<-as.data.frame(DfComplete)
for ( i in 1:length(DfComplete))
{
  for (j in 1:length(DfComplete[,i]))
  {
    DfComplete[j,i]<-as.numeric(DfComplete[j,i])
  }
}


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(DfComplete, powerVector = powers, verbose = 5)
# Plot the results:
#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

bwnet = blockwiseModules(DfComplete, maxBlockSize = 2000,
                         power = 16, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "femaleMouseTOM-blockwise",
                         verbose = 3)

# Load the results of single-block analysis
# load(file = "FemaleLiver-02-networkConstruction-auto.RData");
# Relabel blockwise modules
bwLabels = matchLabels(bwnet$colors, moduleLabels);
# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)

table(bwLabels)

# open a graphics window
#sizeGrWindow(6,6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#sizeGrWindow(12,9)
plotDendroAndColors(geneTree,
                    cbind(moduleColors, bwModuleColors),
                    c("Single block", "2 blocks"),
                    main = "Single block gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

singleBlockMEs = moduleEigengenes(DfComplete, moduleColors)$eigengenes;
blockwiseMEs = moduleEigengenes(DfComplete, bwModuleColors)$eigengenes;

# single2blockwise = match(names(singleBlockMEs), names(blockwiseMEs))
# signif(diag(cor(blockwiseMEs[, single2blockwise], singleBlockMEs)), 3)
# 
# 
# 




#######################################Part 3 WGCNA#####################



# Load network data saved in the second part.
lnames = load(file = "C:/Users/maxence.lechemia/Documents/scripts formalisés/multiomicGTPwithoutoutliers.RData");
lnames

# Define numbers of genes and samples
nGenes = ncol(DfComplete);
nSamples = nrow(DfComplete);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(DfComplete, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

datTraits <-as.data.frame(MEs)
datTraits <-datTraits[,-c(2)]
print(names(datTraits))
MetaData<-read.csv2("C:/Users/maxence.lechemia/Documents/Multi-omics/antibio resistance per strain.csv")
datTraits[,colnames(MetaData)]<-""

for (i in 1:length(rownames(datTraits)))
{
  for (j in 1:length(MetaData$name))
  {
    if (grepl(MetaData$name[j],rownames(datTraits)[i]))
    {
      datTraits[i,c(colnames(MetaData))]<-MetaData[j,]
    }
  }
}

# datTraits[c(1:6,19:30),c(colnames(MetaData))]<-MetaData[1,]
# datTraits[c(7:12),c(colnames(MetaData))]<-MetaData[4,]
# datTraits[c(13:18),c(colnames(MetaData))]<-MetaData[3,]
# datTraits[c(43:48),c(colnames(MetaData))]<-MetaData[5,]
# datTraits[c(31:36),c(colnames(MetaData))]<-MetaData[2,]
# datTraits[c(37:42),c(colnames(MetaData))]<-MetaData[6,]
# datTraits[c(49:60),c(colnames(MetaData))]<-MetaData[2,]
# datTraits$name[which(datTraits$name=="A2")]<-1
# datTraits$name[which(datTraits$name=="B1")]<-2
# datTraits$name[which(datTraits$name=="A2-AST2")]<-1
# datTraits$name[which(datTraits$name=="A2-AST17")]<-1
# datTraits$name[which(datTraits$name=="B1-AST2")]<-2
# datTraits$name[which(datTraits$name=="B1-AST17")]<-2
# datTraits[c(1:3,7:12),c(colnames(MetaData))]<-MetaData[1,]
# datTraits[c(4:6,13:18),c(colnames(MetaData))]<-MetaData[2,]
# datTraits$name[which(datTraits$name=="A2")]<-1
# datTraits$name[which(datTraits$name=="B1")]<-2
datTraits<-datTraits[,-c(1:(length(MEs0))-1)]
datTraits[c("repiq","treatment","type")]<-""
datTraits$treatment[c(1:68)]<-1
datTraits$treatment[c(7:12,31:36)]<-2
datTraits$treatment[c(13:18,37:42)]<-3
datTraits$repiq[c(1:68)]<-1
datTraits$repiq[c(19:24,43:48,57,60,62)]<-2
datTraits$type[c(1:54)]<-1
datTraits$type[c(55:62)]<-2
datTraits$type[c(63:68)]<-3
# datTraits$name[c(61:78)]<-3
# datTraits$name[c(79:84)]<-4
datTraits$name<- as.numeric(factor(datTraits$name))

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
#sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
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

# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$CAZ);
names(weight) = "CAZ"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(DfComplete, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)); 
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(DfComplete, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")



geneModuleMembership<-as.data.frame(geneModuleMembership)
for ( i in 1:length(geneModuleMembership))
{
  for (j in 1:length(geneModuleMembership[,i]))
  {
    geneModuleMembership[j,i]<-as.numeric(geneModuleMembership[j,i])
  }
}

module = "pink"
column = match(module, modNames);
moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
par(mfrow = c(1,1));
a<-verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),ylim = c(0.7,1),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for GENTAMYCIN",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
a<-text(abs(geneModuleMembership[moduleGenes, column]),
     abs(geneTraitSignificance[moduleGenes, 1]),names(DfComplete)[moduleColors=="pink"])
print(a)



# Recalculate module eigengenes
MEs = moduleEigengenes(DfComplete, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$name);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
#sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)



#sizeGrWindow(1000,1000);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)


# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)









output<-"C:/Users/maxence.lechemia/Documents/Multi-omics/ClusterisationGeneGTPwithoutoutliers/"



if (!dir.exists(output)){
  dir.create(output)
} else {
  print("Dir already exists!")
}  



ClusterName<-names(table(moduleColors))
for ( i in 1:length(ClusterName))
{
  Listgene<-names(DfComplete)[moduleColors==ClusterName[i]]
  write.csv(Listgene,file = paste0(output,ClusterName[i],".csv"),row.names = F)
}



write.csv(geneModuleMembership,"C:/Users/maxence.lechemia/Documents/Multi-omics/genemoduleMembershipGTPwithoutoutliers.csv")











ClusterName<-names(table(moduleColors))



for ( i in 1:length(ClusterName))
{
  module<-ClusterName[i]
  pdf(file = paste0(output,module,".pdf"))
  for ( j in 1:length(datTraits))
  {
    weight<-as.data.frame(datTraits[,j])
    names(weight)<-colnames(datTraits)[j]
    modNames = substring(names(MEs), 3)
    geneModuleMembership = as.data.frame(cor(DfComplete, MEs, use = "p"));
    MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)); 
    names(geneModuleMembership) = paste("MM", modNames, sep="");
    names(MMPvalue) = paste("p.MM", modNames, sep="");
    geneTraitSignificance = as.data.frame(cor(DfComplete, weight, use = "p"));
    GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
    names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
    names(GSPvalue) = paste("p.GS.", names(weight), sep="")
    geneModuleMembership<-as.data.frame(geneModuleMembership)
    for ( k in 1:length(geneModuleMembership))
    {
      for (l in 1:length(geneModuleMembership[,k]))
      {
        geneModuleMembership[l,k]<-as.numeric(geneModuleMembership[l,k])
      }
    }
    column = match(module, modNames);
    moduleGenes = moduleColors==module;
    #sizeGrWindow(7, 7);
    par(mfrow = c(1,1));
    a<-verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),ylim = c(0.7,1),
                          abs(geneTraitSignificance[moduleGenes, 1]),
                          xlab = paste("Module Membership in", module, "module"),
                          ylab = paste("Gene significance for ",colnames(datTraits)[j]),
                          main = paste("Module membership vs. gene significance\n"),
                          cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    text(abs(geneModuleMembership[moduleGenes, column]),
            abs(geneTraitSignificance[moduleGenes, 1]),names(DfComplete)[moduleColors=="pink"])
  }
  dev.off()
}

