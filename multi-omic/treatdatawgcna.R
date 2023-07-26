
file<-"C:/Users/maxence.lechemia/Documents/Multi-omics/ClusterisationGeneGTPwithoutoutliers//"
NotreDbAnnotated<-read.csv("C:/Users/maxence.lechemia/Documents/Multi-omics/AnnotatedDbWithNewKeyWords.csv")

ListKeyWords<-colnames(NotreDbAnnotated)[19:61]


DfClusterGene<-data.frame(ListKeyWords)
DfClusterGene<-t(DfClusterGene)
colnames(DfClusterGene)<-DfClusterGene[1,]

ListClusters<-list.files(file)
LenCluster<-c()
# DfClusterGene<-rbind(DfClusterGene,DfClusterGene[1,])
for (i in 1:length(ListClusters))
{
  ListGene<-read.csv(paste0(file,ListClusters[i]))
  LenCluster<-c(LenCluster,length(ListGene$x))
  newdf<-NotreDbAnnotated[which(NotreDbAnnotated$go_locus_tag %in% c(ListGene$x)),c(19:61)]
  for (j in 1:length(newdf))
  {
    newdf[,j]<-as.integer(as.logical(newdf[,j]))
  }
  DfClusterGene<-rbind(DfClusterGene,mapply(sum,newdf))
}

DfClusterGene<-DfClusterGene[-1,]


rownames(DfClusterGene)<-ListClusters

for (i in 1:length(rownames(DfClusterGene)))
{
  rownamestrsplit<-strsplit(rownames(DfClusterGene)[i],".cs")
  rownames(DfClusterGene)[i]<-rownamestrsplit[[1]]
}

DfClusterGeneProp<-as.data.frame(DfClusterGene)

for (i in 1:length(DfClusterGeneProp[,1]))
{
  print(i)
  for ( j in 1:length(DfClusterGeneProp[1,]))
  {
    print(j)
    DfClusterGeneProp[i,j]<-(as.numeric(DfClusterGeneProp[i,j])/as.numeric(LenCluster[i]))*100
  }
}

DfClusterGenePropZtransfo<-as.data.frame(DfClusterGene)
for ( i in 1:length(DfClusterGenePropZtransfo[,1]))
{
  DfClusterGenePropZtransfo[i,]<-scale(as.numeric(DfClusterGenePropZtransfo[i,]))
}


DfClusterGenePropminmax<-as.data.frame(DfClusterGene)
for ( i in 1:length(DfClusterGenePropminmax[,1]))
{
  DfClusterGenePropminmax[i,]<-(as.numeric(DfClusterGenePropminmax[i,])-min(as.numeric(DfClusterGenePropminmax[i,])))/(max(as.numeric(DfClusterGenePropminmax[i,]))-min(as.numeric(DfClusterGenePropminmax[i,])))
}

matrix <- as.numeric(as.matrix(DfClusterGenePropminmax))
dim(matrix) <- dim(DfClusterGenePropminmax)
colnames(matrix)<-colnames(DfClusterGenePropminmax)

library("FactoMineR")
library("factoextra")
res.pca <- PCA(matrix, graph = FALSE)

fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)


library(pheatmap)
pheatmap(matrix[,-3],labels_row = rownames(DfClusterGene))



















##############Weighted keyword correlation###################

geneModuleMembership<-read.csv("C:/Users/maxence.lechemia/Documents/Multi-omics/genemoduleMembershipGTPwithoutoutliers.csv")


rownames(geneModuleMembership)<-geneModuleMembership$X
geneModuleMembership<-geneModuleMembership[-1]


colnames(geneModuleMembership) <- gsub("^.{0,2}", "", colnames(geneModuleMembership))


geneModuleMembership<-abs(geneModuleMembership)

ListKeyWords<-colnames(NotreDbAnnotated)[19:61]


WeightedDfClusterGene<-data.frame(ListKeyWords)
WeightedDfClusterGene<-t(WeightedDfClusterGene)
colnames(WeightedDfClusterGene)<-WeightedDfClusterGene[1,]

ListClusters<-list.files(file)
LenCluster<-c()
# WeightedDfClusterGene<-rbind(WeightedDfClusterGene,WeightedDfClusterGene[1,])
for (i in 1:length(ListClusters))
{
  ListGene<-read.csv(paste0(file,ListClusters[i]))
  tempstrsplit<-strsplit(ListClusters[i],".cs")
  ListClusters[i]<-tempstrsplit[[1]]
  LenCluster<-c(LenCluster,length(ListGene$x))
  newdf<-NotreDbAnnotated[,c(19:61)]
  for (j in 1:length(newdf))
  {
    newdf[,j]<-as.integer(as.logical(newdf[,j]))
  }
  WeightedDf <- data.frame(mapply(`*`,newdf,geneModuleMembership[ListClusters[i]]))
  WeightedDfClusterGene<-rbind(WeightedDfClusterGene,mapply(sum,WeightedDf))
}

WeightedDfClusterGene<-WeightedDfClusterGene[-1,]


rownames(WeightedDfClusterGene)<-ListClusters

for (i in 1:length(rownames(WeightedDfClusterGene)))
{
  rownamestrsplit<-strsplit(rownames(WeightedDfClusterGene)[i],".cs")
  rownames(WeightedDfClusterGene)[i]<-rownamestrsplit[[1]]
}


WeightedDfClusterGenePropMinMax<-as.data.frame(WeightedDfClusterGene)
for ( i in 1:length(WeightedDfClusterGenePropMinMax[,1]))
{
  WeightedDfClusterGenePropMinMax[i,]<-(as.numeric(WeightedDfClusterGenePropMinMax[i,])-min(as.numeric(WeightedDfClusterGenePropMinMax[i,])))/(max(as.numeric(WeightedDfClusterGenePropMinMax[i,]))-min(as.numeric(WeightedDfClusterGenePropMinMax[i,])))
}






matrix <- as.numeric(as.matrix(WeightedDfClusterGene))
dim(matrix) <- dim(WeightedDfClusterGene)
colnames(matrix)<-colnames(WeightedDfClusterGene)


res.pca <- PCA(matrix, graph = FALSE)

fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)



pheatmap(matrix[,-3],labels_row = rownames(DfClusterGene))