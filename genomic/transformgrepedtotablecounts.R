


chemin<-"C:/Users/maxence.lechemia/Documents/Génomique/output usable/"
list_files <- list.files(path =chemin)

NotreDbAnnotated<-read.csv("C:/Users/maxence.lechemia/Documents/Transcriptomique/GFF parsed CSV/parsed_gff_gene_K96243_2015_refseq.csv")

output_dir<-"C:/Users/maxence.lechemia/Documents/Génomique/Genomiccounts"

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir already exists!")
}


for (i in list_files)
{
  SubChemin<-list.files(path =paste0(chemin,i)) #Grep of non synonimous (genomique grep) is first, insertdelet second and stopcodon third
  RawGenomicTableNonSynonimous<-read.table(paste0(chemin,i,"/",SubChemin[1]))
  RawGenomicTableInsertDelet<-read.table(paste0(chemin,i,"/",SubChemin[2]))
  GenomicTableNonSynonimous<-data.frame("LocusTag"=NotreDbAnnotated$locus_tag," "=0)
  GenomicTableInsertDelet<-data.frame("LocusTag"=NotreDbAnnotated$locus_tag," "=0)
  print(i)
  for (j in 1:length(NotreDbAnnotated$locus_tag))
  {
    for (k in 1:length(RawGenomicTableNonSynonimous$V1))
    {
      if ( RawGenomicTableNonSynonimous$V1[k] == NotreDbAnnotated$locus_tag[j] |
           RawGenomicTableNonSynonimous$V1[k] == NotreDbAnnotated$gene[j] |
           grepl(RawGenomicTableNonSynonimous$V1[k],NotreDbAnnotated$gene_synonym[j]))
      {
        GenomicTableNonSynonimous$X.[j]<-GenomicTableNonSynonimous$X.[j]+1
        count<-1
      }
    }
    for (k in 1:length(RawGenomicTableInsertDelet$V1))
    {
      if ( RawGenomicTableInsertDelet$V1[k] == NotreDbAnnotated$locus_tag[j] |
           RawGenomicTableInsertDelet$V1[k] == NotreDbAnnotated$gene[j] |
           grepl(RawGenomicTableInsertDelet$V1[k],NotreDbAnnotated$gene_synonym[j]))
      {
        GenomicTableInsertDelet$X.[j]<-GenomicTableInsertDelet$X.[j]+1
        count<-1
      }
    }
  }
  write.table(GenomicTableInsertDelet, 
              file=paste0(output_dir,"/",i,"-InserDelet.txt")
              , append = FALSE,
              sep = "  ",
              dec = ".",
              row.names = F,
              col.names =F)
  write.table(GenomicTableNonSynonimous, 
              file=paste0(output_dir,"/",i,"-NonSynonimous.txt")
              , append = FALSE,
              sep = "  ",
              dec = ".",
              row.names = F,
              col.names =F)
}



RawGenomicTable<-read.table("C:/Users/maxence.lechemia/Documents/Génomique/Genomics4Maxence/11_K96243vsCori/output/genomiquegrep-C.txt")

GenomicTable<-data.frame("LocusTag"=NotreDbAnnotated$locus_tag," "=0)

for (i in 1:length(RawGenomicTable$V1))
{
  for (j in 1:length(NotreDbAnnotated$locus_tag))
  {
    count<-0
    if ( RawGenomicTable$V1[i] == NotreDbAnnotated$locus_tag[j] |
         RawGenomicTable$V1[i] == NotreDbAnnotated$gene[j] |
         grepl(RawGenomicTable$V1[i],NotreDbAnnotated$gene_synonym[j]))
    {
      GenomicTable$X.[j]<-GenomicTable$X.[j]+1
      count<-1
      break 
    }
  }
  if(count==0)
  {
    print(RawGenomicTable$V1[i])
  }
}
length(RawGenomicTable$V1)
sum(GenomicTable$X.)

colnames(GenomicTable)<-NA

write.csv(GenomicTable,"C:/Users/maxence.lechemia/Documents/Génomique/counts mutated genomic/C-ori.csv")
