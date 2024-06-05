#######Overlap of cell-specific SAS and cancer DEGs#########
##Epithelial cells as an example
###Epithelial cell-specific down-regulation of SAS and cancer DEGs overlap
load("/boot3/cjl/Aging/data/UCSC/GTEX_RNAseq_data/细胞类型识别/组织统计/DEGs_list.Rdata")
tumor_tissue <- as.data.frame(read.csv("/boot3/cjl/Aging/data/UCSC/GTEX_RNAseq_data/细胞类型识别/组织统计/tumor_tissue.csv",header = T,sep = ","))
DEGs_tolal <- data.frame(module_names=c(1),
                         down_age_gene=c(1),
                         tumor_names=c(1),
                         tumor_Up_DEGs=c(1),
                         down_age_tumor_up_DEGs=c(1),
                         down_age_tumor_up_DEGs_P=c(1)
)
#aging_down_cancer_up
intersect_gene <- list()
for (x in tumor_tissue$tumor) {
  Down_up_go.obj <- testGeneOverlap(newGeneOverlap(EpiCell_Down_More2,DEGs_list[[paste0(x,"_DEGs_up")]],genome.size = tumor_tissue[which(tumor_tissue$tumor == x),4]))
  intersect_gene <- c(intersect_gene,list(getIntersection(Down_up_go.obj)))
  intersection_num <- length(getIntersection(Down_up_go.obj))
  Pvalue <- getPval(Down_up_go.obj)
  DEGs_tolal <- rbind(DEGs_tolal,c("Epithelial_down",length(EpiCell_Down_More2),x,length(DEGs_list[[paste0(x,"_DEGs_up")]]),intersection_num,Pvalue))
}
names(intersect_gene) <- tumor_tissue$tumor

save(intersect_gene,file="/boot3/cjl/Aging/data/UCSC/GTEX_RNAseq_data/细胞类型识别/组织统计/Epithelial_cell/Epithelial_cancer_intersect.rda")

a <- as.data.frame(table(unlist(intersect_gene)))
#aging_down_cancer_down
DEGs_down_tolal <- data.frame(
  tumor_down_DEGs=c(1),
  down_age_tumor_down_DEGs=c(1),
  down_age_tumor_down_DEGs_P=c(1)
)
for (x in tumor_tissue$tumor) {
  Down_up_go.obj <- testGeneOverlap(newGeneOverlap(EpiCell_Down_More2,DEGs_list[[paste0(x,"_DEGs_down")]],genome.size = tumor_tissue[which(tumor_tissue$tumor == x),4]))
  intersection_num <- length(getIntersection(Down_up_go.obj))
  Pvalue <- getPval(Down_up_go.obj)
  DEGs_down_tolal <- rbind(DEGs_down_tolal,c(length(DEGs_list[[paste0(x,"_DEGs_down")]]),intersection_num,Pvalue))
}
DEGs_tolal <- cbind(DEGs_tolal,DEGs_down_tolal)
DEGs_tolal <- DEGs_tolal[-1,]

paixu <- as.data.frame(read.csv("/boot3/cjl/Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/paixu.csv",header = T,sep = ","))
DEGs_tolal <- DEGs_tolal[match(paixu$tumor,DEGs_tolal$tumor_names),]

data <- DEGs_tolal[,c(4:9)]
data$tumor_Up_DEGs <- paste0(DEGs_tolal$tumor_names," (",data$tumor_Up_DEGs,")")
data$tumor_down_DEGs <- paste0(DEGs_tolal$tumor_names," (",data$tumor_down_DEGs,")")
####Fig.4C------
##aging_down_cancer_up_plot
library(pheatmap)
library(ggthemes)
Down_up_data <- data[,c(1,2,3)]
Down_up_data$down_age_tumor_up_DEGs <- as.numeric(Down_up_data$down_age_tumor_up_DEGs)
bk <- c(seq(0,4500,by=100))
p1 <- pheatmap(as.matrix(Down_up_data[,2]),
               scale = "none",
               color = c(colorRampPalette(colors = c("white","#B31312"))(length(bk))),
               cluster_rows = F,
               cluster_cols = F,
               fontsize = 12,
               cellwidth = 28,#小格子的宽度
               cellheight = 20, #小格子的长度
               legend_breaks=seq(0),
               display_numbers = TRUE,
               number_format = "%.0f",
               legend = FALSE,
               fontsize_col = 10,
               labels_row = Down_up_data$tumor_Up_DEGs
)

Down_up_data$down_age_tumor_up_DEGs_P <- as.numeric(Down_up_data$down_age_tumor_up_DEGs_P)
Down_up_data$down_age_tumor_up_DEGs_P <- -1*log10(as.numeric(Down_up_data$down_age_tumor_up_DEGs_P))
Down_up_data$tp <- ifelse(Down_up_data$down_age_tumor_up_DEGs_P < -log10(0.05),"",
                          ifelse(Down_up_data$down_age_tumor_up_DEGs_P >= -log10(0.05) & Down_up_data$down_age_tumor_up_DEGs_P<=-log10(0.01),"*",
                                 ifelse(Down_up_data$down_age_tumor_up_DEGs_P >= -log10(0.01) & Down_up_data$down_age_tumor_up_DEGs_P<=-log10(0.001),"**","***")))
Down_up_data$tumor_Up_DEGs <- factor(Down_up_data$tumor_Up_DEGs,levels =rev(Down_up_data$tumor_Up_DEGs))
ggplot(Down_up_data,aes(x=tumor_Up_DEGs,y=down_age_tumor_up_DEGs_P))+
  geom_bar(stat="identity",position="dodge",fill="#B31312")+
  coord_flip()+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme(panel.background = element_blank(),
        axis.title.x=element_blank())
#aging_down_cancer_down_plot
Down_down_data <- data[,c(4,5,6)]
Down_down_data$down_age_tumor_down_DEGs <- as.numeric(Down_down_data$down_age_tumor_down_DEGs)
bk <- c(seq(0,4500,by=100))
pheatmap(as.matrix(Down_down_data[,2]),
         scale = "none",
         color = c(colorRampPalette(colors = c("white","#3887BE"))(length(bk))),
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 12,
         cellwidth = 28,#小格子的宽度
         cellheight = 20, #小格子的长度
         legend_breaks=seq(0),
         display_numbers = TRUE,
         number_format = "%.0f",
         legend = FALSE,
         fontsize_col = 10,
         labels_row = Down_down_data$tumor_down_DEGs
)
Down_down_data$down_age_tumor_down_DEGs_P <- -1*log10(Down_down_data$down_age_tumor_down_DEGs_P)
Down_down_data$tp <- ifelse(Down_down_data$down_age_tumor_down_DEGs_P < -log10(0.05),"",
                            ifelse(Down_down_data$down_age_tumor_down_DEGs_P >= -log10(0.05) & Down_down_data$down_age_tumor_down_DEGs_P<=-log10(0.01),"*",
                                   ifelse(Down_down_data$down_age_tumor_down_DEGs_P >= -log10(0.01) & Down_down_data$down_age_tumor_down_DEGs_P<=-log10(0.001),"**","***")))
Down_down_data$tumor_down_DEGs <- factor(Down_down_data$tumor_down_DEGs,levels =rev(Down_down_data$tumor_down_DEGs))
ggplot(Down_down_data,aes(x=tumor_down_DEGs,y=down_age_tumor_down_DEGs_P))+
  geom_bar(stat="identity",position="dodge",fill="#3887BE")+
  coord_flip()+
  geom_text(aes(label = tp), vjust = 1,color="red",angle=90)+
  theme(panel.background = element_blank(),
        axis.title.x=element_blank())

####---Fig.4E-----
Epithelial_cancer <- c("UCEC","UCS","BRCA","LUAD","LUSC","COAD","READ","OV","PAAD","TGCT","THCA")
Epithelial_cancer_intergene <- intersect_gene[Epithelial_cancer]
Epithelial_cancer_order <- Epithelial_cancer_intergene[order(sapply(Epithelial_cancer_intergene, length),decreasing = T)]
pubudata <- data.frame(matrix(NA,nrow=382,ncol = 11))
rownames(pubudata) <- data$Var1
colnames(pubudata) <- names(Epithelial_cancer_order)
for(x in 1:11){
  pubudata[Epithelial_cancer_order[[x]],x] <- "Up"
}
pubudata[is.na(pubudata)] = ""
col = c("Up" = "#B31312")
library()
alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),   
  Up = alter_graphic("rect", fill = col["Up"])
)
col_order <- colnames(pubudata)
oncoPrint(pubudata, col = col,
          alter_fun = alter_fun,show_row_names = F,column_order = col_order,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                             foo1 = 1:11,
                                             col = list(foo1 = c("1"="#B31312",
                                                                 "2"="#B92625",
                                                                 "3"="#BF3938",
                                                                 "4"="#C54C4B",
                                                                 "5"="#CB5F5E",
                                                                 "6"="#D17271",
                                                                 "7"="#D88585",
                                                                 "8"="#DE9898",
                                                                 "9"="#E4ABAB",
                                                                 "10"="#EABEBE",
                                                                 "11"="#F0D1D1")),
                                             gp = gpar(col = "black")
          )
)

####Fig.4F------
library(pheatmap)
gene <- c("SDC1","MEX3A","PHLDA2")
tumor_tissue <- as.data.frame(read.csv("./Aging/data/UCSC/TCGA_RNAseq_data/DEseq2循环/tumor_tissue.csv",header = T,sep = ","))
all_gene_log2fc <- data.frame(matrix(NA, nrow = 1, ncol = 24))
colnames(all_gene_log2fc) <- tumor_tissue$tumor
all_gene_P <- data.frame(matrix(NA, nrow = 1, ncol = 24))
colnames(all_gene_P) <- tumor_tissue$tumor
for(j in 1:length(gene)){
  all_gene_log2fc[j,] <- 0
  all_gene_P[j,] <- 0
  gene <- gene[j]
  for (i in 1:length(tumor_tissue$tumor)) {
    DEseq2data <- as.data.frame(read.table(paste0("./Aging/data/UCSC/TCGA_RNAseq_data/",tumor_tissue[i,1],"/DEGs_all_info.txt"),header=T,row.names=1,check.names=F))
    all_gene_log2fc[j,i] <- DEseq2data[gene,]$log2FoldChange
    all_gene_P[j,i] <- ifelse(abs(all_gene_log2fc[j,i])<= 1.5,0,DEseq2data[gene,]$pvalue)
    
  }
}
rownames(all_gene_log2fc) <- gene
all_gene_log2fc[is.na(all_gene_log2fc)] <- 0
rownames(all_gene_P) <- gene
all_gene_P[is.na(all_gene_P)] <- 0
display_number <- all_gene_P
if (!is.null(all_gene_P)){
  sssmt <- all_gene_P>0&all_gene_P<0.001
  all_gene_P[sssmt] <-'***'
  ssmt <- all_gene_P >0.001& all_gene_P< 0.01
  all_gene_P[ssmt] <-'**'
  smt <- all_gene_P >0.01& all_gene_P <0.05
  all_gene_P[smt] <- '*'
  all_gene_P[!ssmt&!smt&!sssmt]<- ''
} else {
  all_gene_P <- F
}
bk <- c(seq(-9,-1.5,by=0.01),seq(-1.4,1.4,by=0.01),seq(1.5,9,by=0.01))
pheatmap::pheatmap(all_gene_log2fc,
                   scale= 'none',
                   color= c(colorRampPalette(colors = c("#203972","#a8ccda"))(length(bk)*2.5/6),colorRampPalette(colors = c("#EEF5FF","#EEF5FF"))(length(bk)/6),colorRampPalette(colors = c("#F78CA2","#872341"))(length(bk)*2.5/6)),
                   legend_breaks=seq(-9,9,1.5),
                   breaks=bk,
                   display_numbers = all_gene_P,number_color = "white")


