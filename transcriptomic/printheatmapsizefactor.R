library(DESeq2)
library(magrittr)
library(apeglm)
library(data.table)
library(readr)
library(ggplot2)
library(DESeq2)
library(genefilter)
library(pheatmap)
library(gridExtra)
source("C:/Users/maxence.lechemia/Documents/scripts formalisés/trimtablebyrow.R")
source("C:/Users/maxence.lechemia/Documents/scripts formalisés/rowlocustogene.R")

print_heatmap_sizefactor<-function(df,
                                   vsd,
                                   PathAnnotated,
                                   listcol,
                                   SelectGene=1,
                                   TableSelection=c(),
                                   NumberSelect=100,
                                   ColTarget=c(),
                                   Flatten=TRUE,
                                   UserInputMostImportantAttribute,
                                   ColDf,
                                   RowDf,
                                   BoolClusterCol=T,
                                   BoolClusterRow=T)
{
  topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE),5000)
  topVarGenes<-topVarGenes[!duplicated(topVarGenes)]
  print(topVarGenes)
  mat  <- assay(vsd)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  print(mat)
  NotreDbAnnotee<-read.csv(PathAnnotated)
  
  if (SelectGene==1)#Takes the first *NumberSelect* genes in the table and print it and does not take Table Selection
  {
    mat<-mat[c(1:NumberSelect),]
    print("PING AAAAAAAAAAAAAAAAAAAAAAAAAAAAAH")
  }
  else if(SelectGene==2)#Takes the gene in the column ColTarget from the PATH OF TABLE TableSelection
  {
    TableSelection<-read_table(TableSelection)
    mat<-trim_table_by_row(Trimed =  mat,Trimer = TableSelection[ColTarget])
  }
  else if(SelectGene==3)#Takes the gene from the list TableSelection
  {
    mat<-trim_table_by_row(Trimed =  mat,Trimer = TableSelection)
  }
  print(mat)
  #Permet de ne mettre que les genes qui sont annotés au moins une fois sur notre dbannotée
  
  if (SelectGene>1)
  {
    for (i in row.names(mat))
    {
      count<-0
      for (j in listcol)
      {
        if (length(which(NotreDbAnnotee$go_locus_tag==i))!=0)
        {
          if (NotreDbAnnotee[which(NotreDbAnnotee$go_locus_tag==i),j]==TRUE)
          {
            count<-count+1
            print("ping")
          }
          else
          {
          }
        }
        else
        {
        }
      }
      if(count<1)
      {
        mat<-mat[!(row.names(mat) %in% c(i)),]
      }
    }
  }
  anno_row <- as.data.frame(NotreDbAnnotee[which(NotreDbAnnotee$go_locus_tag%in% row.names(mat)),listcol])
  for (i in 1:length(anno_row))
  {
    anno_row[,i]<-as.numeric(anno_row[,i])
    print(i)
  }
  rownames(anno_row) = row.names(mat)[which(row.names(mat)%in% NotreDbAnnotee$go_locus_tag)]
  for (i in row.names(mat))
  {
    if (length(which(NotreDbAnnotee$go_locus_tag==i))!=0)
    {
      
      if (NotreDbAnnotee[which(NotreDbAnnotee$go_locus_tag==i),"gene"]!="")
      {
        row.names(mat)[which(row.names(mat)==i)]<-NotreDbAnnotee[which(NotreDbAnnotee$go_locus_tag==i),"gene"]
      }
      else
      {
        
      }
    }
    else
    {
      
    }
  }
  nthroot = function(x,n) {
    (abs(x)^(1/n))*sign(x)
  }
  
  print(mat)
  print(row.names(mat))
  # mat  <- mat - rowMeans(mat)
  
  for (i in 1:length(df[,1]))
  {
    print(i)
    for (j in 1:length(mat[,i]))
    {
      mat[j,i]<-nthroot(mat[j,i],4)
    }
  }
  anno <- as.data.frame(colData(vsd)[, c(UserInputMostImportantAttribute,"sizeFactor")])
  for (i in 1:length(ColDf))
  {
    anno[ColDf[i]]<-df[ColDf[i]] 
  }
  
  anno_row <- row_locus_to_gene(df = anno_row,NotreDbAnnotee = NotreDbAnnotee)
  print("mat c'est ça:")
  print(mat)
  print(anno)
  print(anno_row)
  jpeg(filename = "C:/Users/maxence.lechemia/Documents/plot.jpeg")
  a<-pheatmap(mat,
           annotation_col = anno,
           annotation_row = anno_row,
           fontsize=10,
           cluster_cols = BoolClusterCol,
           cluster_rows = BoolClusterRow)
  print(a)
  dev.off()
  return(a)
  
}
