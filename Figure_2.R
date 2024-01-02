
library(data.table)
library(ComplexHeatmap)
library(ggpval)
library(ggplot2)
library(ggsignif)
library(ggpubr)
 library(tidyverse)
 library(data.table)
 library(reshape2)
library(VennDetail)
library(clusterProfiler)
library(org.Hs.eg.db)



##Code to generate figure 2A

# Function to process deg list  #
process_data <- function(file_path) {
  res_sorted <- read.table(file = file_path, header = TRUE)
  res_sel <- subset(res_sorted, padj < 0.05)
  selectedRowsf <- res_sel 
  return(selectedRowsf)
}

# Process data
res_sel_50_g1 <- process_data("disease_Hc_vs_LC_G1") # DEG_Ana/result/disease_Hc_vs_LC_G1
res_sel_50_g2 <- process_data("disease_Hc_vs_LC_G2") # DEG_Ana/result/disease_Hc_vs_LC_G2
res_sel_50  <- rbind(res_sel_50_g1 ,res_sel_50_g2)
ig_genes <- read.table(file="IG_genes.txt", header = TRUE)
res_sel_50 <- res_sel_50[!res_sel_50$SYMBOL %in% ig_genes$genes,]
meta_fil <- subset(read.table(file = "metadata_lc_groups", header = TRUE), groups %in% c("Hc", "LC1", "LC2", "CC"))
rlog <- read.table(file="rlog_norm", header = TRUE)   
mat  <- rlog[rownames(res_sel_50) ,rownames(meta_fil) ] 
mat_new_n   <- mat - rowMeans(mat)
row.names(mat_new_n) <- res_sel_50$SYMBOL
levels(meta_fil) <- c("Hc", "LC",  "LC2", "CC")
colours <- list(
  'group' = c("Hc" = "forestgreen", 'LC1' = 'darkorange2', 'LC2' = 'Sienna', 'CC' = 'purple'))
colAnn <- HeatmapAnnotation(df = anno_new1,
  which = 'col',
  col = colours,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'))

pdf(file="Figure2A.pdf",  width = 7, height = 10)
Heatmap(as.matrix(mat_new_n), name = "Score",top_annotation = colAnn,  
show_row_names = TRUE,border = TRUE,cluster_columns = FALSE,     show_column_names = TRUE,  column_split =meta_fil$Condition,
)
dev.off()


##Code to generate figure 2B 
process_data_down <- function(file_path) {
res_sorted <- read.table(file = file_path, header = TRUE)
res_sel <- subset(res_sorted, padj < 0.05)
d_down <-  subset(res_sel, log2FoldChange < -1.5) 
dt <- rownames(d_down)
  return(dt)
}

process_data_up <- function(file_path) {
res_sorted <- read.table(file = file_path, header = TRUE)
res_sel <- subset(res_sorted, padj < 0.05)
d_up <-  subset(res_sel, log2FoldChange > 1.5) 
dt <- rownames(d_up)
return(dt)
}

list_deg_g2_up <- process_data_up("disease_CC_vs_LC_G1") # DEG_Ana/result/disease_Hc_vs_LC_G1
list_deg_g3_up <- process_data_up("disease_CC_vs_LC_G2")# DEG_Ana/result/disease_Hc_vs_LC_G2
list_deg_g2_down <- process_data_down("disease_CC_vs_LC_G1") # DEG_Ana/result/disease_Hc_vs_LC_G1
list_deg_g3_down <- process_data_down("disease_CC_vs_LC_G2")# DEG_Ana/result/disease_Hc_vs_LC_G2
# Using gene symbol for GO analysis ## 
dt <- list("CD_cs_LC1_Up"=list_deg_g2_up,
"CD_cs_LC1_Down"=list_deg_g2_down,
"CD_cs_LC2_Up"=list_deg_g3_up,
"CD_cs_LC2_Down"=list_deg_g3_down)
ck_go <- compareCluster(geneCluster = dt , fun = enrichGO,OrgDb    = org.Hs.eg.db,
               ont      = "BP",
             #  level    = 3,
               readable = TRUE)
pdf(file = "Figure2B.pdf",
    width = 7, # The width of the plot in inches
    height = 8) # The height of the plot in inches
dotplot(ck_go, showCategory=5,font.size=10)

dev.off()





##Code to generate figure 2C

# Process data
res_sel_50_g1 <- process_data("disease_Hc_vs_LC_G1") # DEG_Ana/result/disease_Hc_vs_LC_G1
res_sel_50_g2 <- process_data("disease_Hc_vs_LC_G2") # DEG_Ana/result/disease_Hc_vs_LC_G2
res_sel_50  <- rbind(res_sel_50_g1 ,res_sel_50_g2)
sel_genes <-read.table(file="cxl_gene_list", header = TRUE) 
res_sel_50 <- res_sel_50[!res_sel_50$SYMBOL %in% sel_genes$Genes,]
meta_fil <- subset(read.table(file = "metadata_all", header = TRUE), groups %in% c("Hc", "LC1", "LC2", "CC"))
rlog <- read.table(file="rlog_gene_symbol", header = TRUE)
mat  <- rlog[rownames(res_sel_50) ,rownames(meta_fil) ] 
mat_new_n   <- mat - rowMeans(mat)
row.names(mat_new_n) <- res_sel_50$SYMBOL
levels(meta_fil) <- c("Hc", "LC",  "LC2", "CC")
colours <- list(
  'group' = c("Hc" = "forestgreen", 'LC1' = 'darkorange2', 'LC2' = 'Sienna', 'auCC' = 'purple'))
colAnn <- HeatmapAnnotation(df = anno_new1,
  which = 'col',
  col = colours,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'))
pdf(file="Figure2C.pdf",  width = 7, height = 10)
Heatmap(as.matrix(mat_new_n), name = "Score",top_annotation = colAnn,  
show_row_names = TRUE,border = TRUE,cluster_columns = FALSE,     show_column_names = TRUE,  column_split =mat_new_n$Condition,
)
dev.off()





##Code to generate figure 2D
res_sel_50_g1 <- process_data("disease_Hc_vs_LC_G1") # DEG_Ana/result/disease_Hc_vs_LC_G1
res_sel_50_g2 <- process_data("disease_Hc_vs_LC_G2") # DEG_Ana/result/disease_Hc_vs_LC_G2
selectedRowsf1 <- res_sel_50_g1[1:25,]
selectedRowsf2 <- res_sel_50_g2[1:25,]
res_sel_50  <- rbind(selectedRowsf1 ,selectedRowsf2)
meta <- read.table(file="metadata_allgroups", header = TRUE)
meta_fil <- subset(meta, groups %in% c("Hc","LC1","LC2"))
## read normalized data # 
rlog <- read.table(file="rlog_gene_symbol", header = TRUE)
mat  <- rlog[rownames(res_sel_50) ,rownames(meta_fil) ] 
mat_new_n   <- mat - rowMeans(mat)
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


pdf(file="Figure2D.pdf",  width = 7, height = 10)
Heatmap(as.matrix(mat_new_n), name = "Score",top_annotation = colAnn,  
    show_row_names = TRUE,border = TRUE,cluster_columns = FALSE,  cluster_rows = TRUE,     show_column_names = TRUE,  column_split =meta_fil$Condition ) 
dev.off()
     
     
##Code to generate figure 2E
data <- read.table(file="deconvolution_result.txt", header = TRUE, sep = "\t", row.names = 1)
data <- data[,1:22]
data_g <- melt(data)
data_g$Groupf <- as.factor(data_g$cell_types) 
p <- ggplot(data_g, aes(x=Groupf, y=value, color = Groupf)) +
geom_bar(position="dodge", stat="identity")  + 
theme_bw(base_size = 8)
pdf(file="Figure_2E.pdf",  width = 7, height = 10)
p + theme_classic()
dev.off()
 

##Code to generate figure 2f
c_palette <- c("forestgreen", "darkorange2", "Sienna", "magenta4")
data_g <- melt(data)
data_g$Groupf <- as.factor(data_g$Group) 
data_g$Group <- factor(data_g$Groupf , levels=c("HC", "LC1", "LC2", "CC"))
data_g_m1 <-  subset(data_g, variable %in% c("Macrophages_M1","T_cells_CD4_memory_resting", "B_cells_memory","T_cells_CD4_naive"))
stat.test <- compare_means(value ~ Groupf,  data =data_g,  group.by = "variable",ref.group = "HC",
               method = "wilcox.test", p.adjust.method=  "BH")
               
pdf(file="Figure_2F.pdf",  width = 4, height = 4)
p1 <- ggplot(data_g_m1, aes(x=Groupf, y=value, color = Groupf)) +
geom_boxplot(outlier.shape = NA) +
#geom_jitter(alpha=0.4, width=0.2, aes(color=Groupf)) +
geom_jitter(width=0.15)+
scale_colour_manual(values= c_palette) +
facet_wrap(~ variable) +
theme_bw(base_size = 8) 
p1
dev.off()


##Code to generate figure 2G

##################### intersection upset plot  ########################################
cc_lc1  <- read.table(file="disease_Hc_vs_LC_G1", header = TRUE)  # DEG_Ana/result/disease_Hc_vs_LC_G1
cc_lc1<- subset(cc_lc1, padj < 0.05)
cc_lc1 <- filter(cc_lc1, log2FoldChange < 0 | log2FoldChange  > 0 ,padj  < 0.05)
cc_lc2  <- read.table(file="disease_Hc_vs_LC_G2", header = TRUE) 
cc_lc2 <- subset(cc_lc2, padj < 0.05)
cc_lc2 <- filter(cc_lc2, log2FoldChange < 0 | log2FoldChange  > 0,padj  < 0.05)
cc_lc3  <- read.table(file="disease_Hc_vs_CC", header = TRUE) 
cc_lc3 <- subset(cc_lc3, padj < 0.05)
cc_lc3 <- filter(cc_lc3, log2FoldChange < 0 | log2FoldChange  > 0,padj  < 0.05)

ven <- venndetail(list("HC_vs_LC1" =cc_lc1$SYMBOL,
"HC_vs_LC2" =cc_lc2$SYMBOL, 
                    "HC_vs_CC" =cc_lc3$Gene)
 )
pdf(file="Figure_2G.pdf",  width = 7, height = 6)  
plot(ven, type = "upset")   
dev.off()


##Code to generate figure 2H

data_matrix <- result(ven, wide = TRUE)
# extract for unique sets
du2 <- subset(data_matrix ,  SharedSets < 2) 
g1 <- subset(du2, HC_vs_LC1 ==1 )
# extract for shared among two sets 
du3 <- subset(data_matrix,  SharedSets < 3) 
g2 <- subset(du3, HC_vs_LC1 ==1 & HC_vs_CC ==1 )
# extract for shared among two sets 
du3 <- subset(data_matrix,  SharedSets == 3) 
g3 <- d3
gene_list <- list("LC1_vs_Hc" = g1$Detail ,"shared_LC1_and_CC" = g2$Detail, "shared_among_all" = g3$Detail)
formula_res <- compareCluster(gene_list, fun="enrichGO",OrgDb = org.Hs.eg.db,keyType="SYMBOL", ont = "BP")
pdf(file="Figure_2H.pdf",  width = 8, height = 8)
dotplot(formula_res, font.size=8, showCategory = 5) #+ facet_grid(~category)

dev.off()



##Code to generate figure 2I

# extract binary matrix ### 
##################### GO based on unique and common  ########################################
data_matrix <- result(ven, wide = TRUE)
# extract for unique sets
du2 <- subset(data_matrix ,  SharedSets < 2) 
g1 <- subset(du2, HC_vs_LC2 ==1 )
cc_lc2f = cc_lc2[!duplicated(cc_lc2$SYMBOL),]
g1_deg  <- subset( cc_lc2f, SYMBOL %in% g1$Detail)
    
# filter for fold change ### 
Hc_vs_LC_G2_down_unique <- subset(g1_deg, log2FoldChange < -1.5 )
Hc_vs_LC_G2_up_unique   <- subset(g1_deg, log2FoldChange > 1.5 )
 
# extract for shared among three sets 
du3 <- subset(data_matrix,  SharedSets < 3) 
g2 <- subset(du3, HC_vs_LC2 ==1 & HC_vs_CC ==1 )
cc_lc3f = cc_lc3[!duplicated(cc_lc3$Gene),]
g2_deg  <- subset( cc_lc3f, Gene %in% g2$Detail)

# filter for fold change ###    
Hc_vs_auCC_down_shared <- subset(g2_deg, log2FoldChange < -1.5 )
Hc_vs_auCC_up_shared   <- subset(g2_deg, log2FoldChange > 1.5 )
  
g3 <- subset(du2, HC_vs_LC2 ==1 & HC_vs_CC ==1 )
cc_lc2f = cc_lc2[!duplicated(cc_lc2$SYMBOL),]
g3_deg  <- subset(cc_lc2f, SYMBOL %in% g3$Detail)
     
# filter for fold change ###    
Hc_vs_LC_G2_down_shared <- subset(g3_deg, log2FoldChange < -1.5 )
Hc_vs_LC_G2_up_shared <- subset(g3_deg, log2FoldChange > 1.5 )

gene_list <- list("Hc_LC_G2_d_un" = Hc_vs_LC_G2_down_unique$SYMBOL ,"Hc_LC_G2_up_un" = Hc_vs_LC_G2_up_unique$SYMBOL,
 "Hc_auCC_d_sh" = Hc_vs_auCC_down_shared$Gene,   "Hc_auCC_up_sh" = Hc_vs_auCC_up_shared$Gene,
"Hc_LC_G2_d_sh" = Hc_vs_LC_G2_down_shared$SYMBOL,  "Hc_LC_G2_up_sh" = Hc_vs_LC_G2_up_shared$SYMBOL)
     
formula_res <- compareCluster(g12, fun="enrichGO",OrgDb = org.Hs.eg.db,keyType="SYMBOL", ont = "BP")
pdf(file="Figure_2I.pdf",  width = 8, height = 8)
dotplot(formula_res, font.size=8, showCategory = 5) #+ facet_grid(~category)

dev.off()



##Code to generate figure 2J
##################### correlation ########################################
data <- read.table(file="rlog_norm", header = TRUE) 
sel_genes <- read.table(file="GeneList_Imm_channels_correlations.txt", header = TRUE)  # selective gene list # 
sel_rows <- subset(data,data$genesymbol %in% sel_genes$Gene)
rownames(sel_rows) <- sel_rows$genesymbol
sel_rows  <- data.frame(t(sel_rows))
# get hclust from lower type ### 
pdf(file="Figure_2J.pdf",  width = 8, height = 14)
corrplot(cor(sel_rows, method = 'square', order = 'hclust', type = 'lower', diag = FALSE,addrect = 4)
 dev.off()
     
     
     


