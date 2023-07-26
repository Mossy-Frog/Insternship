#Convert the rownames of the dataframe containing locus tags, into gene names according to the annotation dataset
#Auteur: Maxence LECHEMIA
row_locus_to_gene <- function(df,NotreDbAnnotee)
{
  failsafe<-c("")
  for (i in rownames(df))
  {
    if (length(which(NotreDbAnnotee$go_locus_tag==i))!=0)
    {
      if (NotreDbAnnotee[which(NotreDbAnnotee$go_locus_tag==i),"gene"]!="" & (NotreDbAnnotee[which(NotreDbAnnotee$go_locus_tag==i),"gene"] %in% failsafe)==F)
      {
        print(NotreDbAnnotee[which(NotreDbAnnotee$go_locus_tag==i),"gene"])
        row.names(df)[which(row.names(df)==i)]<-NotreDbAnnotee[which(NotreDbAnnotee$go_locus_tag==i),"gene"]
        failsafe<-append(failsafe,NotreDbAnnotee[which(NotreDbAnnotee$go_locus_tag==i),"gene"])
      }
    }
  }
  return(df)
}