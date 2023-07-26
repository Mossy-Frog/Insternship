library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(stringr)
library(vsn)


print_explo_plot<-function(vsd,
                           dds_filtered,
                           UserInputColNames,
                           UserInputMostImportantAttribute,
                           SepRowHeatmap,
                           SepColHeatmap,
                           output)
{
  MostImportantAttribute<-which(UserInputColNames==UserInputMostImportantAttribute)
  pdf(file= paste0(output,"ExploratoireAuto.pdf"),onefile = TRUE )
  meanSdPlot(assay(vsd),ranks=F)
  # meanSdPlot(assay(rld),ranks=F)
  
  dds<-estimateSizeFactors(dds_filtered)
  
  
  sampleDists_VSD <- dist(t(assay(vsd)))
  # sampleDists_rlog <- dist(t(assay(rld)))
  
  
  cat_df<-data.frame("Remove"=c(rep("",length(vsd@colData@listData[[1]]))))
  
  for (i in 1:(length(vsd@colData@listData)-1))
  {
    cat_df[names(vsd@colData@listData)[i]]<-vsd@colData@listData[[i]]
  }
  
  
  cat_df<-subset(cat_df,select=UserInputColNames)
  row.names(cat_df) = row.names(colData(vsd))
  
  
  if(SepRowHeatmap==0)
  {
    SepRowHeatmap<-3
  }
  
  if(SepColHeatmap==0)
  {
    SepColHeatmap<-8
  }
  
  
  sampleDistMatrix_vsd <- as.matrix( sampleDists_VSD )
  colnames(sampleDistMatrix_vsd) <- colnames(vsd)
  rownames(sampleDistMatrix_vsd) <- NULL
  
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix_vsd,
           clustering_distance_rows = sampleDists_VSD,
           clustering_distance_cols = sampleDists_VSD,
           col = colors, #les noms sont bien au bon endroit
           main = "Heatmap echantillons~souche, VST",
           cutree_cols = SepColHeatmap, cutree_rows= SepRowHeatmap,
           annotation_col = cat_df)
  
  
  
  poisd <- PoissonDistance(t(counts(dds_filtered)))
  samplePoisDistMatrix <- as.matrix( poisd$dd )
  
  colnames(samplePoisDistMatrix) <- colnames(dds_filtered)
  rownames(samplePoisDistMatrix) <- NULL
  pheatmap(samplePoisDistMatrix,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,
           col = colors, main="Heatmap echantillons~souche, dist de Poisson",
           cutree_cols = SepColHeatmap, cutree_rows= SepRowHeatmap,
           annotation_col = cat_df)
  if (MostImportantAttribute==0)
  {
    pcaData<-plotPCA(vsd, intgroup = as.character(names(vsd@colData@listData)[1:(length(vsd@colData@listData)-1)]))+
      ggtitle("PCA, transformation VST",returnData=T)  
  }
  
  if (MostImportantAttribute!=0)
  {
    temp_plot<-plotPCA(vsd, intgroup = as.character(names(vsd@colData@listData)[MostImportantAttribute]))+
      ggtitle("PCA, transformation VST") 
    pcaData <- plotPCA(vsd, intgroup = as.character(names(vsd@colData@listData)[MostImportantAttribute]),returnData=T)
    print(temp_plot)
  }
  
  # ggplot(pcaData, aes(x = PC1, y = PC2)) +
  #   geom_point(size =3) +
  #   coord_fixed() +  ggtitle("PCA with VST data")
  
  gpca <- glmpca(counts(dds), L=2)
  gpca.dat <- gpca$factors
  for (i in 1:(length(dds@colData@listData)-1))
  {
    gpca.dat[names(dds@colData@listData)[i]]<-dds@colData@listData[i]
  }
  
  for (i in 1:(length(dds@colData@listData)-1))
  {
    if (MostImportantAttribute!=i)
    {
      temp_plot<-ggplot(gpca.dat, aes_string(x = "dim1", y = "dim2", color =names(dds@colData@listData)[MostImportantAttribute], shape =names(dds@colData@listData)[i])
      ) + geom_point(size =3) + coord_fixed() + 
        ggtitle("glmpca - Generalized PCA")  
      print(temp_plot)
    }
  }
  
  
  mds <- as.data.frame(colData(vsd))  %>%
    cbind(cmdscale(sampleDistMatrix_vsd))
  
  
  for (i in 1:(length(dds@colData@listData)-1))
  {
    if (MostImportantAttribute!=i)
    {
      temp_plot<-ggplot(mds, aes_string(x = "`1`", y = "`2`", color =names(dds@colData@listData)[MostImportantAttribute], shape =names(dds@colData@listData)[i])) +
        geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
      print(temp_plot)
    }
  }
  
  mdsPois <- as.data.frame(colData(dds)) %>%
    cbind(cmdscale(samplePoisDistMatrix))
  
  for (i in 1:(length(dds@colData@listData)-1))
  {
    if (MostImportantAttribute!=i)
    {
      temp_plot<-ggplot(mdsPois, aes_string(x = "`1`", y = "`2`", color =names(dds@colData@listData)[MostImportantAttribute], shape =names(dds@colData@listData)[i])) +
        geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")
      print(temp_plot)
    }
  }
  dev.off()
}

