library(ggplot2)
######Sharing of SAS from different cells between tissues######
#--Fig.S5B-
####Take endothelial cells for example
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/组织统计/Snakey_data.rda")
Endothelial_cell <- Sankey_data[which(Sankey_data$cell_name == "Endothelial cell"),]
#Up
Endothelial_cell_Up <- Endothelial_cell[which(Endothelial_cell$moduleType== "Up"),]
EndCell_Up <- all_tissue_module_list[Endothelial_cell$pathway[which(Endothelial_cell$moduleType== "Up")]]
Endothelial_cell_all_tissue <- list()
for (onetissue in unique(Endothelial_cell_Up$tissue)) {
  Endothelial_cell_one <- Endothelial_cell_Up[which(Endothelial_cell_Up$tissue == onetissue),]
  Endothelial_cell_one_tissue <- list(unique(unlist(EndCell_Up[Endothelial_cell_one$pathway])))
  names(Endothelial_cell_one_tissue) <- onetissue
  Endothelial_cell_all_tissue <- c(Endothelial_cell_all_tissue,Endothelial_cell_one_tissue)
}
EndCell_Up_total <- as.data.frame(table(as.data.frame(unlist(Endothelial_cell_all_tissue))))
EndCell_Up_More2 <- as.character(EndCell_Up_total[which(EndCell_Up_total$Freq>=2),]$Var1)
EndCell_Up_total_freq <- as.data.frame(table(EndCell_Up_total$Freq))
colnames(EndCell_Up_total_freq) <- c("pinshu","cishu")
#Down
Endothelial_cell_Down <- Endothelial_cell[which(Endothelial_cell$moduleType== "Down"),]
EndCell_Down <- all_tissue_module_list[Endothelial_cell$pathway[which(Endothelial_cell$moduleType== "Down")]]
Endothelial_cell_all_tissue <- list()
for (onetissue in unique(Endothelial_cell_Down$tissue)) {
  Endothelial_cell_one <- Endothelial_cell_Down[which(Endothelial_cell_Down$tissue == onetissue),]
  Endothelial_cell_one_tissue <- list(unique(unlist(EndCell_Down[Endothelial_cell_one$pathway])))
  names(Endothelial_cell_one_tissue) <- onetissue
  Endothelial_cell_all_tissue <- c(Endothelial_cell_all_tissue,Endothelial_cell_one_tissue)
}

EndCell_Down_total <- as.data.frame(table(as.data.frame(unlist(Endothelial_cell_all_tissue))))
EndCell_Down_More2 <- as.character(EndCell_Down_total[which(EndCell_Down_total$Freq>=2),]$Var1)
EndCell_Down_total_freq <- as.data.frame(table(EndCell_Down_total$Freq))
colnames(EndCell_Down_total_freq) <- c("pinshu","cishu")
#两组放在一起画
two_group_data <- cbind(EndCell_Up_total_freq,EndCell_Down_total_freq[,2])
colnames(two_group_data) <- c("tissue_num","Up","Down")
df1 <- reshape2::melt(two_group_data,id.vars = 'tissue_num')
df1$gene <- c(434,152,25,5,274,47,9,0)
ggplot(df1, aes(
  x = factor(tissue_num,levels = unique(tissue_num)),
  y = ifelse(variable == "Up", value, -value),
  fill = variable)) +
  geom_col(width = .85,alpha=0.5)+
  geom_col(
    aes(x = factor(tissue_num,levels = unique(tissue_num)),
        y = ifelse(variable == "Up", gene, -gene),
        fill = variable), width = .5 
  )+
  scale_fill_manual(values = c("#EE3432","#1B64A4"))+   
  geom_text(                                                  
    aes(label=value,                                        
    ),
    size=3)+
  ylab("Number of gene")+xlab("Number of tissues with sharing genes")+
  scale_y_continuous(                                       
    labels = abs,                            
    expand = expansion(mult = c(0.1, 0.1)))+
  theme_bw()+
  theme(panel.grid=element_blank (),legend.position = "none")

#---Fig.2E-----
setwd("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Tissuesstatistics/CellTypeEnrich")
#Up and Down regulated
for (cell in c("Endothelial","Fibroblast","Tcell")) {
  Cell_up <- read.csv(paste0(cell,"_Up_Enrich.csv"),header = T)
  Cell_down <- read.csv(paste0(cell,"_Down_Enrich.csv"),header = T)
  Cell_data <- rbind(Cell_up[,c(4,5,6)],Cell_down[,c(4,5,6)])
  Cell_data$up_down <- c(rep("Up",length(rownames(Cell_up))),rep("Down",length(rownames(Cell_down))))
  Cell_data <- Cell_data %>% arrange(up_down,desc(LogP))
  Cell_data$LogP <- ifelse(Cell_data$up_down == "Down",Cell_data$LogP,-Cell_data$LogP)
  Cell_data$Description <- factor(Cell_data$Description,Cell_data$Description)
  x_max = max(abs(Cell_data$LogP))
  color_list <- c("#B0CAE6","#E5B3B3")
  
  p <- ggplot(Endothelial,aes(x =LogP, y = Description, fill = up_down)) + 
    geom_col() + 
    scale_x_continuous(limits = c(-x_max,x_max)) +
    scale_fill_manual(values = color_list) +
    geom_text(data = Endothelial[which(Endothelial$up_down == "Up"),],
              aes(x = -0.5, y = Description, label = Description,color = up_down),
              size = 4,
              hjust = 1) + 
    geom_text(data = Endothelial[which(Endothelial$up_down == "Down"),],
              aes(x = 0.5, y = Description, label = Description,color = up_down),
              size = 4,
              hjust = 0)  + 
    scale_colour_manual(values = color_list) +
    labs(x="log(Pvalue)", y=NULL, title="enriched pathway") +
    annotate("text", x = 25, y = 15, label = "Up", size = 10, fontface = "bold", color="#E5B3B3") +
    annotate("text", x = -25, y =5, label = "Down", size =10, fontface = "bold", color="#B0CAE6") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 25,hjust=0.5), 
      axis.text.x = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = 'none'
    )
  pdf(paste0("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Tissuesstatistics/CellTypeEnrich",cell,"_enrichplot.pdf"),width=8,height = 6)
  print(p)
  dev.off()
}
#Down regulated
for (cell in c("Epithelial","MastCell","Macrophage")) {
  Cell_data <- read.csv(paste0(cell,"_Down_Enrich.csv"),header = T)
  Cell_data$up_down <- c(rep("Down",length(rownames(Cell_up))))
  Cell_data <- Cell_data %>% arrange(up_down,desc(LogP))
  Cell_data$LogP <- ifelse(Cell_data$up_down == "Down",Cell_data$LogP,-Cell_data$LogP)
  Cell_data$Description <- factor(Cell_data$Description,Cell_data$Description)
  x_max = max(abs(Cell_data$LogP))
  color_list <- c("#B0CAE6","#E5B3B3")
  
  p <- ggplot(Endothelial,aes(x =LogP, y = Description, fill = up_down)) + 
    geom_col() + 
    scale_x_continuous(limits = c(-x_max,x_max)) +
    scale_fill_manual(values = color_list) +
    geom_text(data = Endothelial[which(Endothelial$up_down == "Up"),],
              aes(x = -0.5, y = Description, label = Description,color = up_down),
              size = 4,
              hjust = 1) + 
    geom_text(data = Endothelial[which(Endothelial$up_down == "Down"),],
              aes(x = 0.5, y = Description, label = Description,color = up_down),
              size = 4,
              hjust = 0)  + 
    scale_colour_manual(values = color_list) +
    labs(x="log(Pvalue)", y=NULL, title="enriched pathway") +
    annotate("text", x = 25, y = 15, label = "Up", size = 10, fontface = "bold", color="#E5B3B3") +
    annotate("text", x = -25, y =5, label = "Down", size =10, fontface = "bold", color="#B0CAE6") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 25,hjust=0.5), 
      axis.text.x = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = 'none'
    )
  pdf(paste0("./Aging/data/UCSC/GTEX_RNAseq_data/CellType_ident/Tissuesstatistics/CellTypeEnrich",cell,"_enrichplot.pdf"),width=8,height = 6)
  print(p)
  dev.off()
}