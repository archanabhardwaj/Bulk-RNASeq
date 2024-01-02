library(scatterplot3d)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(data.table)
library(ComplexHeatmap)


#### Figure 1B####

##Code to generate figure 1B

rlog <- read.table(file="rlog_gene_symbol", header = TRUE)
data_n <- data.frame(rlog) 
Hc <- rep("Hc", times=10)
LC <- rep("LC", times=8)
auCC <- rep("auCC", times=4)
disease <- factor(c(Hc, LC, auCC))
my_design <- as.data.frame(disease)
row.names(my_design) <- patients
pca <- prcomp(t(data_n))
percentage <- round(((pca$sdev^2) / (sum(pca$sdev^2))) * 100, 2)
pca_data <- data.frame(pca$x, SampleType=my_design$disease)
pca_data <- read.table(file="pca_Data", header = TRUE)
colors <- c("lightslateblue", "forestgreen", "darkorange2")

pdf(file = "Figure1B.pdf", width=8, height=8)
s3d <- scatterplot3d(pca_data[,1:3], pch = 1, color=colors, tick.marks = TRUE, label.tick.marks = TRUE, grid=TRUE, box=TRUE,
)

legend(s3d$xyz.convert(7.5, 3.5, 4), legend = levels(pca_data$SampleType),
      col =  c("lightslateblue", "forestgreen", "darkorange2"), pch = 1)
         
dev.off()




##Code to generate figure 1D
data <- read.table(file="DEG_Ana/result/disease_Hc_vs_LC_G1", header = TRUE)
data <- subset(data, padj < 0.05)
d_down <-  subset(data, log2FoldChange < -1.5) 
d_up <-  subset(data, log2FoldChange > 1.5) 
dt <- list(DOWN= rownames(d_down), UP= rownames(d_up))

# conversion of gene symbol to entrezid for KEGG analysis ## 
ids_down <- bitr(d_down$Gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
ids_up <- bitr(d_up$Gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
dt <- list(DOWN= ids_down$ENTREZID, UP= ids_up$ENTREZID)
#ck <- compareCluster(geneCluster = dt, fun = enrichKEGG)
# Using gene symbol for GO analysis ## 
ck_go <- compareCluster(geneCluster = dt, fun = enrichGO,OrgDb    = org.Hs.eg.db,
               ont      = "BP",
             #  level    = 3,
               readable = TRUE)

pdf(file = "Figure1D.pdf",
    width = 7, # The width of the plot in inches
    height = 8) # The height of the plot in inches


dotplot(ck_go, showCategory=10,font.size=10)

dev.off()



##Code to generate figure 1E

res_sorted_f1 <- read.table(file="disease_Hc_vs_LC_G1", header = TRUE) # DEG_Ana/result/disease_Hc_vs_LC_G1
res_sel_f1 <- subset(res_sorted_f1, padj < 0.05)
selectedRowsf <- res_sel_f1[1:50,]
res_sel_50  <- selectedRowsf

meta <- read.table(file="metadata_all", header = TRUE)
meta_fil <- subset(meta, groups %in% c("Hc","LC1","LC2"))
## read normalized data # 
rlog <- read.table(file="rlog_gene_symbol", header = TRUE)
mat  <- rlog[rownames(res_sel_50) ,rownames(meta_fil) ] 
mat_new_n   <- mat - rowMeans(mat)
row.names(mat_new_n) <- res_sel_50$SYMBOL
mat_new_n <- mat_new
anf <- data.frame(meta_fil[, c("groups")] )


colours <- list(
  'group' = c("Hc" = "forestgreen", 'LC1' = 'darkorange2', 'LC2' = 'Sienna'))
  
# create column annotation 
  
colAnn <- HeatmapAnnotation(df = meta_fil,
  which = 'col',
  col = colours,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'))


pdf(file="Figure1E.pdf",  width = 7, height = 10)

Heatmap(as.matrix(mat_new_n), name = "Score",top_annotation = colAnn,  

    show_row_names = TRUE,border = TRUE,cluster_columns = FALSE,  cluster_rows = TRUE,     show_column_names = TRUE,  column_split =meta_fil$Condition
  )
     
dev.off()
     
     



##Code to generate figure 1F and 1G

process_data <- function(file_path) {
  res_sorted <- read.table(file = file_path, header = TRUE)
  res_sel <- subset(res_sorted, padj < 0.05)
d_down <-  subset(res_sel, log2FoldChange < -1.5) 
d_up <-  subset(res_sel, log2FoldChange > 1.5) 
dt <- list(DOWN= rownames(d_down), UP= rownames(d_up))
  return(dt)
}
list_deg_g2 <- process_data("disease_HC_vs_LC_G2")
list_deg_g3 <- process_data("disease_LC1_vs_LC2")
# Using gene symbol for GO analysis ## 
ck_go <- compareCluster(geneCluster = list_deg_g2 , fun = enrichGO,OrgDb    = org.Hs.eg.db,
               ont      = "BP",
             #  level    = 3,
               readable = TRUE)
pdf(file = "Figure1F.pdf",
    width = 7, # The width of the plot in inches
    height = 8) # The height of the plot in inches
dotplot(ck_go, showCategory=10,font.size=10)

dev.off()

ck_go <- compareCluster(geneCluster = list_deg_g3 , fun = enrichGO,OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               readable = TRUE)        
pdf(file = "Figure1G.pdf",
    width = 7, # The width of the plot in inches
    height = 8) # The height of the plot in inches
dotplot(ck_go, showCategory=10,font.size=10)

dev.off()


##Code to generate figure 1H
res_sorted_f1 <- read.table(file="disease_Hc_vs_LC_G1", header = TRUE) # DEG_Ana/result/disease_Hc_vs_LC_G1
res_sel_f1 <- subset(res_sorted_f1, padj < 0.05)
selectedRowsf2 <- res_sel_f1[1:50,]

res_sorted_f2 <- read.table(file="disease_Hc_vs_LC_G2", header = TRUE) # DEG_Ana/result/disease_Hc_vs_LC_G1
res_sel_f2 <- subset(res_sorted_f2, padj < 0.05)
selectedRowsf2 <- res_sel_f1[1:50,]

gene_list <- unique(c(rownames(selectedRowsf1), rownames(selectedRowsf2))
meta <- read.table(file="metadata_all", header = TRUE) # read normalized data from all samples 
meta_fil <- subset(meta, groups %in% c("Hc","LC1","LC2"))
## read normalized data # 
rlog <- read.table(file="rlog_gene_symbol", header = TRUE)
mat  <- rlog[gene_list ,rownames(meta_fil) ] 
mat_new_n   <- mat - rowMeans(mat)

colours <- list(
  'group' = c("Hc" = "forestgreen", 'LC1' = 'darkorange2', 'LC2' = 'Sienna'))
  
# create column annotation 
colAnn <- HeatmapAnnotation(df = meta_fil,
  which = 'col',
  col = colours,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'))


pdf(file="Figure1H.pdf",  width = 7, height = 10)
Heatmap(as.matrix(mat_new_n), name = "Score",top_annotation = colAnn,  
 show_row_names = TRUE,border = TRUE,cluster_columns = FALSE, cluster_rows = TRUE,    show_column_names = TRUE,  column_split =meta_fil$Condition
  )
     
dev.off()
     
