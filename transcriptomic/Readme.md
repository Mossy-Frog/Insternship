# Transcriptomic


This file contains every scripts used for the GFF annotation and the transcriptomic automation.




```mermaid
graph LR;
    GFFparser-->Joindf;
    Joindf-->ScriptgetAPIGOterms;
    ScriptgetAPIGOterms-->Joincolumns;
    Joincolumns-->Cleancolumns;
    Cleancolumns-->MergeRows;
    MergeRows-->GotermTreatment;
```
### Workflow of the creation of a GFF annotated with both GO-terms and their associated keywords






```mermaid
flowchart LR
autoexplo --> autowritecsvstat --> autoheatmaplogfoldchange -->autoheatmapsizefactor

```
### Workflow of the transcriptomic analysis automation.







```mermaid
flowchart LR
subgraph autoexplo
  direction TB
  createdfdeseq --> createmainvsd
  createmainvsd --> printexploplot
  
end
```
### Functions used by the autoexplo script, outputing a pdf providing information on the dataset using DESeq2.






```mermaid
flowchart LR
subgraph autowritecsvstat
  direction TB
  createpossibilitytable --> createdfdeseq
  createdfdeseq --> creatematrixlevel
  creatematrixlevel --> releveling
  releveling --> rowtosubset
  rowtosubset --> Statisticdeseq
end
```
### Functions used by the autowritecsvstat script, outputing numerous csv files, containing comparisons between each type of sample or conditions from the dataset.







```mermaid
flowchart LR
subgraph autoheatmaplogfoldchange
  direction RL
  Annotation ~~~ rowtosubset
end
```
### Function and data used by the autoheatmaplogfoldchange, outputing a JPEG containing an annotated heatmap corresponding to an inputted comparison CSV file.






```mermaid
flowchart LR
subgraph autoheatmapsizefactor
  direction TB
  createdfdeseq --> createmainvsd
  Annotation ~~~ createmainvsd
  createmainvsd -->printheatmapsizefactor

end
```
### Functions and data used by the autoheatmapsizefactor, outputing a JPEG file containing an annotated heatmap corresponding to the counts analysis by DESeq2

