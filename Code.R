source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("org.Hs.eg.db")

# Read Data from CSV file
df <- read.csv("sample_rnaseq_data.csv", header = TRUE, row.names = "gene_id")
head(df)

subsetDF1 <-  subset(df, rownames(df) %in% c("ENSG00000000003", "ENSG00000000457", "ENSG00000000460",
                                             "ENSG00000001167", "ENSG00000002822", "ENSG00000003096", 
                                             "ENSG00000003509", "ENSG00000004866", "ENSG00000004975",
                                             "ENSG00000005189", "ENSG00000005801", "ENSG00000006016"))
                     
head(subsetDF1)
rownames(subsetDF1)

library(edgeR)
library(org.Hs.eg.db)
subsetDF1$Symbol <- mapIds(org.Hs.eg.db, group=group, rownames(subsetDF1),
                               keytype="ENSEMBL", column="SYMBOL")

head(subsetDF1)

rownames(subsetDF1) <- subsetDF1$Symbol # Use Gene Symbol as rownames 
subsetDF1 <- subsetDF1[ , !(names(subsetDF1) %in% c("Symbol"))] #remove Symbol column

logsubsetDF1 <- log(subsetDF1)
logsubsetDF1

logsubsetDF1 <- logsubsetDF1[is.finite(rowSums(logsubsetDF1)),] # Remove -Inf / Na value
logsubsetDF1

library(gplots)

png(filename = "heatmapDemo.png", # output format and filename of the figure.
     bg = "transparent", # background color
     width = 700, # figure width
     height = 700, # figure height
     pointsize = 18 # base font size
     )

col.pan <- colorpanel(100, "blue", "white", "red") # define colour panel
heatmap.2(as.matrix(logsubsetDF1), 
          main="R Heatmap demostration",  # Title of the chart
          col=col.pan, # color panel defined above
          Rowv=TRUE, # determines if and how the row dendrogram should be reordered.
          Colv = FALSE, # determines if and how the column dendrogram should be reordered. 
          cexCol=1.4, # font size for the labels on column
          trace="none", # character string indicating whether a solid "trace" line should be drawn across 'row's or down 'column's, 'both' or 'none'.
          dendrogram="row", # character string indicating whether to draw 'none', 'row', 'column' or 'both' dendrograms.
          density.info="histogram", # character string indicating whether to superimpose a 'histogram', a 'density' plot, or no plot ('none') on the color-key.
          margin=c(10,9)) # numeric vector of length 2 containing the margins for column and row names, respectively.

dev.off() # shutdown the figure output device when it is done