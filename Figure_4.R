
#Figure 4 H #### 
library(pheatmap)
library(ggplot2)
library(VennDetail)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ComplexHeatmap)
library(stringr)

   
##Code to generate figure 4A

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


# process data
process_data <- function(file_path) {
  res_sorted <- read.table(file = file_path, header = TRUE)
  res_sel <- subset(res_sorted, padj < 0.05)
  selectedRowsf <- res_sel 
  return(selectedRowsf)
}



g1_down <- process_data_down("disease_CD_vs_LC") # DEG_Ana/result/disease_CD_vs_Hc
g2_down <- process_data_down("disease_UC_vs_LC") # DEG_Ana/result/disease_UC_vs_Hc
g1_up <- process_data_up("disease_CD_vs_LC") # DEG_Ana/result/disease_CD_vs_Hc
g2_up <- process_data_up("disease_UC_vs_LC") # DEG_Ana/result/disease_UC_vs_Hc
g1_downf <- bitr(g1_down$symbol, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
g2_downf  <- bitr(g2_down$symbol, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
g1_upf  <- bitr(g1_up$symbol, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
g2_upf <- bitr(g2_up$symbol, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

dt <- list("CD_DOWN"= g1_downf$ENTREZID,"CD_UP"= g1_upf$ENTREZID, "UC_DOWN"= g2_downf$ENTREZID,"UC_UP"= g2_upf$ENTREZID      )
ck_go <- compareCluster(geneCluster = dt, fun = enrichGO,OrgDb    = org.Hs.eg.db,
               ont      = "BP",
             #  level    = 3,
               readable = TRUE)
       
        
pdf(file = "Figure_4A.pdf",
    width = 7, # The width of the plot in inches
    height = 8) # The height of the plot in inches
dotplot(ck_go, showCategory=10,font.size=6,by='GeneRatio')
dev.off()

    
   
   
   
##Code to generate figure 4B####
   
## read deg list ### 
cd1 <- read.table(file="disease_CD_vs_LC2.txt", header = TRUE)
cd1f <- filter(cd1, log2FoldChange < -5 | log2FoldChange  > 5,padj  < 0.05)
cd1ff <- cd1f[cd1f$symbol %like% "^MIR",] 
cd2 <- read.table(file="disease_UC_vs_LC2.txt", header = TRUE)
cd2f <- filter(cd2, log2FoldChange < -5 | log2FoldChange  > 5,padj  < 0.05)
cd2ff <- cd2f[cd2f$symbol %like% "^MIR",] 
selectived_mi <- rbind(cd2ff, cd1ff)

## read rlog data # 
d <- read.table(file="rlog_gene_symbol", header = TRUE)
d$id <- rownames(d)
mi_rows <- subset(d,d$gene_name %in% unique(selectived_mi$symbol))

meta <- read.table(file="metadata_all", header = TRUE) # read normalized data from all samples 
meta$sample <- rownames(meta)
df <- subset(meta, meta$disease %in% c("LC1", "LC2", "CD", "UC"))
mi_data <- mi_rows[,c(df$sample,"gene_name")]
mat_n  <- mi_data - rowMeans(mi_data)
#custom order of labels ## 
order <- c("LC2056" ,"LC2053" ,"LC4002" ,"LC4004", "LC41" ,  "LC55"  , "LC2051" ,"LC2063" ,"C1675" , "C1492"  ,"D0742" , "D0743", "C1661",  "C1666" , "C2485" , "C1658"  )
meta  <- meta[order,]
mat_new_n <- mat_n[, order]
colours <- list(
  'disease' = c('LC1' = 'darkorange2', 'LC2' = 'orange', 'UC' = 'purple', 'CD' = 'pink' ))
meta$disease <- factor(meta$disease, levels = c("LC1", "LC2", "CD", "UC"))
colAnn <- HeatmapAnnotation(df = meta,
  which = 'col',
  col = colours,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'))
pdf(file="Fig_4B.pdf",  width = 7, height = 10)
Heatmap(as.matrix(mat_new_n ), name = "Score",top_annotation = colAnn,  
show_row_names = TRUE,border = TRUE,cluster_columns = FALSE,cluster_rows = TRUE

)
     
dev.off()
     
     
     
     
    
##Code to generate figure 4C

######mirna for pathway analysis ####### 
 
d <- read.table(file="miRNA_gene_list.txt", sep = "\t", header = TRUE) # custom gene list 
ent <- bitr(d$symbol, fromType='SYMBOL', toType='ENTREZID',OrgDb = org.Hs.eg.db)
re <- enrichWP(ent$ENTREZID, "Homo sapiens")

pdf(file = "Figure_2C.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 6) # The height of the plot in inches

dotplot(re)
dev.off()
write.csv(file="wiki_database_annotation_from_clusterProfiler.csv", data.frame(re))

     
     
   
##Code to generate figure 4D
# combined up and down for intersection ### 
res_sel_50_g1 <- process_data("disease_CD_vs_Hc") # DEG_Ana/result/disease_CD_vs_Hc
res_sel_50_g2 <- process_data("disease_UC_vs_Hc") # DEG_Ana/result/disease_UC_vs_Hc
res_sel_50_g3 <- process_data("disease_LC2_vs_Hc") # DEG_Ana/result/disease_LC2_vs_Hc
gene_list <- list("CD_vs_Hc" =rownames(res_sel_50_g1), "UC_vs_Hc" =  rownames(res_sel_50_g2), "LC2_vs_Hc" = rownames(res_sel_50_g3))
pdf(file="Fig_4D.pdf",  width = 7, height = 10)
plot(ven, type = "upset")
dev.off()
     
     
   
##Code to generate figure 4E
## process for both up and down seperately 
 
g1_down <- process_data_down("disease_CD_vs_Hc") # DEG_Ana/result/disease_CD_vs_Hc
g2_down <- process_data_down("disease_UC_vs_Hc") # DEG_Ana/result/disease_UC_vs_Hc
g3_down <- process_data_down("disease_LC2_vs_Hc") # DEG_Ana/result/disease_LC2_vs_Hc
g1_up <- process_data_up("disease_CD_vs_Hc") # DEG_Ana/result/disease_CD_vs_Hc
g2_up <- process_data_up("disease_UC_vs_Hc") # DEG_Ana/result/disease_UC_vs_Hc
g3_up <- process_data_up("disease_LC2_vs_Hc") # DEG_Ana/result/disease_LC2_vs_Hc
 

gene_list <- list("CD_vs_Hc_down" = g1_down, "CD_vs_Hc_up" =  g1_up, 
"UC_vs_Hc_down" = g2_down, "UC_vs_Hc_up" =  g2_up, 
"LC2_vs_Hc_down" = g3_down, "LC2_vs_Hc_up" =  g3_up)

ck_go <- compareCluster(geneCluster = gene_list , fun = enrichGO,OrgDb    = org.Hs.eg.db,
               ont      = "BP",
             #  level    = 3,
               readable = TRUE)
          
pdf(file = "Figure_4E.pdf",
    width = 7, # The width of the plot in inches
    height = 8) # The height of the plot in inches
dotplot(ck_go, showCategory=10,font.size=10)
dev.off()

 
 
 
#####Code to generate figure 4F
# read data from ven interesection among LC2_vs_UC, LC2_UC_CD### 
g1 =read.table(file="DEGs/results/Intersection/gene_list_LC2_UC_shared_down.txt", header = TRUE)
g2 =read.table(file="DEGs/results/Intersection/gene_list_LC2_UC_shared_up.txt", header = TRUE)
g3 =read.table(file="DEGs/results/Intersection/gene_list_LC2_UC_CD_shared_down.txt", header = TRUE)
g4 =read.table(file="DEGs/results/Intersection/gene_list_LC2_UC_CD_shared_up.txt", header = TRUE)
g5 =read.table(file="DEGs/results/Intersection/gene_list_LC2_CD_down.txt", header = TRUE)
g6 =read.table(file="DEGs/results/Intersection/gene_list_LC2_CD_up.txt", header = TRUE)
all_g <- list("Down_LC2_UC" = g1$gene, 
"Up_LC2_UC" = g2$gene, 
"Down_LC2_UC_CD" = g3$gene, 
"Up_LC2_UC_CD" = g4$gene, 
"Down_LC2_CD" = g5$gene, 
"Down_LC2_CD" = g6$gene)

formula_res <- compareCluster( all_g, fun="enrichGO",OrgDb = org.Hs.eg.db,keyType="SYMBOL", ont = "BP")
pdf(file="Figure_4F.pdf",  width = 6, height =9)
dotplot(formula_res, font.size=8, showCategory = 5) #+ facet_grid(~category)
dev.off()


     
#####Code to generate figure 4G
dg <- read.table(file="rlog_gene_symbol", header = TRUE)
genes <- read.table(file="selective.txt", header = TRUE) 
df <- subset(dg, dg$gene_name %in% genes$gene)
data =  df[!duplicated(df$gene_name),]
meta <- read.table(file="metadata_all", header = TRUE) # read normalized data from all samples 
meta$sample <- rownames(meta )
df <- subset(meta , meta$disease %in% c("Hc", "LC1", "LC2", "CC", "CD", "UC"))
dataf <- data[,df$sample]
rownames(dataf) <-  data$gene_name
mat_n  <- dataf - rowMeans(dataf)

order  <- c("Hc8K11" , "Hc18K21", "HcV2" ,   "HcV3"  ,  "HcV4"  ,  "HcV5"  ,  "HcV6"   , "HcV7"   , "HcV8"  ,  "HcV9" , "LC2056" , "LC2053" , "LC4002" , "LC4004"  , "LC41"  ,  "LC55" ,   "LC2051" , "LC2063" , "auCC7" ,  "auCC13"  ,"auCC15" , "auCC18" , "C1675" ,  "C1492"  , "D0742"  ,"D0743" ,  "C1661"  , "C1666"  , "C2485"  , "C1658"  )
    
meta_order <- data.frame(meta[order,])
mat_new_n <- mat_n[,order]   
colours <- list(
  'disease' = c("Hc" = "forestgreen", "LC1" = "darkorange2", "LC2" = "Sienna","CC"= "lightslateblue", "CD"= "magenta4", "UC" =  "grey"))
meta_order$disease <- factor(meta_order$disease, levels = c("Hc", "LC1", "LC2", "CC", "CD", "UC"))
colAnn <- HeatmapAnnotation(df = meta_order,
  which = 'col',
  col = colours,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'))
pdf(file="Fig_4g.pdf",  width = 7, height = 10)
Heatmap(as.matrix(mat_new_n , name = "Score",top_annotation = colAnn,show_row_names = TRUE,border = TRUE,cluster_columns = FALSE,cluster_rows = TRUE)
dev.off()
     
   
   
   
   
   
   
   
   
   
   
   






























