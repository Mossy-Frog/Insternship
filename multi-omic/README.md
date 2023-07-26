# Multi Omic

This file contains the scripts used for the multi-omic analysis with the WGCNA package.

The numerous attempts of WGCNA analysis with different types of data are present in the scripts "multiomic wgcna".

The script "transform non conform to conform" is used to format the table of counts, in order that each biological layers are in the exact same format.

The script "Normalizecount" is used to test a Z transformation to normalize each table of counts

The script "DESeq2normalization" presents a way to normalize each table of counts using the DESeq2 package. It wields better results than the Z transformation.

The "treatdatawgcna" is an attempt to further analyse clusters after the WGCNA analysis.
