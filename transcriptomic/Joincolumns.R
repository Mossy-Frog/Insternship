Go_df<- read.csv("C:/Users/maxence.lechemia/Documents/Transcriptomique/parsed-gff-with-GOterms.csv")


# test<-c("V4","V5")
# 
# print(Go_df[test[1]])

join_column<- function(df,ListCol,NameCol,overwrite,erase)#Si overwrite, donne la priorité aux données par ordre décroissant dans la ListCol
#Si Erase, enlève les colonnes données dans listcol
{
  df[NameCol]<-c("")
  for (i in ListCol)
  {
    for (j in 1:length(df[[i]]))
    {
      if (overwrite)
      {
        if ((is.na(df[j,i])==FALSE)&(identical(df[j,i],"")==FALSE))
        {
          df[j,NameCol]<-df[j,i]
          print(j) 
        }
      }
      else
      {
        if ((is.na(df[j,i])==FALSE)&(identical(df[j,i],"")==FALSE))
        {
          df[j,NameCol]<-paste0(toString(df[j,NameCol]),",",toString(df[j,i]))
          print(j) 
        }
      }
    }
  }
  if (erase)
  {
    df <- df[,!names(df) %in% ListCol]
  }
  return(df)
}




NewGo<- join_column(df = Go_df,ListCol = c("locus_tag.x.x","locus_tag.y.x","old_locus_tag"),NameCol = "go_locus_tag",overwrite = TRUE,erase = T)

NewGo<- join_column(df = NewGo,ListCol = c("locus_tag.x.y","locus_tag.y.y"),NameCol = "new_locus_tag",overwrite = TRUE,erase = T)

NewGo<- join_column(df = NewGo,ListCol = c("product.x","product.y"),NameCol = "product",overwrite = F,erase = T)

NewGo<- join_column(df = NewGo,ListCol = c("gene.x","gene.y"),NameCol = "gene",overwrite = T,erase = T)

NewGo<- join_column(df = NewGo,ListCol = c("gene_synonym.x","gene_synonym.y"),NameCol = "gene_synonym",overwrite = FALSE,erase = T)

NewGo<- join_column(df = NewGo,ListCol = c("protein_id.x"),NameCol = "protein_id_EMBL",overwrite = T,erase = T)

NewGo<- join_column(df = NewGo,ListCol = c("protein_id.y"),NameCol = "protein_id_Refseq",overwrite = T,erase = T)

NewGo<- join_column(df = NewGo,ListCol = c("Note.x","Note.y"),NameCol = "Note",overwrite = FALSE,erase = TRUE)


write.csv(NewGo, "C:/Users/maxence.lechemia/Documents/Transcriptomique/parsed-gff-with-GOterms-joined-cols.csv", row.names=FALSE)
