library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(stringr)
library(vsn)



#####Use of DESEQ2 ot normalize the couts of different types. 



source("C:/Users/maxence.lechemia/Documents/scripts formalisés/createdfdeseq.R")

chemin<-"C:/Users/maxence.lechemia/Documents/Multi-omics/ordered formated fused genomic and all transcripto outliers removed/"
df<-create_df_deseq(chemin = chemin,
                    UserInput = list(list("A2","B1","C","PDUC"),list("AST2","AST17"),list("R5","R10"),list("seq2"),list("Transcripto","Proteomic","Génomique")),
                    UserInputColNames = list("strain","treatment","repiq","sequencage","omictype"),
                    UserInputIfEmpty =list("B1","Notreatment","norepiq","seq1","Transcripto"),
                    UserInputRelevel= list("B1","Notreatment","norepiq","seq1","Transcripto")) #Creation of the dataframe for the deseq2

FormulaDES<-paste0("~","strain" )#formula

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = df,
                                       directory = chemin,
                                       design = as.formula(FormulaDES))#create the table of counts for deseq2

dds <- estimateSizeFactors(ddsHTSeq)#formatage
sizeFactors(dds)#formatage

normalized_counts <- counts(dds, normalized=TRUE)#normalization

write.table(normalized_counts, file="C:/Users/maxence.lechemia/Documents/Multi-omics/normalized counts by deseq2 with proteo without outliers.txt", sep="\t", quote=F, col.names=NA)
