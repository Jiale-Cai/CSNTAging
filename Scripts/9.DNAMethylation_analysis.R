########Age-differential methylation analysis of epithelial cells UTDA genes in normal tissues######
library(broom)
library(biomaRt)
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/细胞类型识别/组织统计/Epithelial_cell/Epithelial_cancer_intersect.rda")
Epithelial_cancer <- intersect_gene[c("BRCA","LUAD","LUSC","COAD","READ","OV","TGCT")]
Normal_tissue_Epithelial <- data.frame(cancaer=c("BRCA","LUAD","LUSC","COAD","READ","OV","TGCT"),
                                       tissue=c("Breast","Lung","Lung","Colon","Colon","Ovary","Testis"))
gene <- read.table("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/gencode_v42_annotation.txt",sep = "\t",header = T)
#get epithelial cell UTDA genes methylation matrices
cancer <- Normal_tissue_Epithelial$cancaer[2]
for (cancer in Normal_tissue_Epithelial$cancaer) {
  tissue <- Normal_tissue_Epithelial[Normal_tissue_Epithelial$cancaer == cancer,2]
  Methy_df <- as.data.frame(fread(paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/",tissue,"_Methy.txt"),sep = "\t"))
  res_pmatch <- read.table(paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UTDA_cg_site/",cancer,"res_pmatch.txt"),sep = "\t",header = T)
  rownames(Methy_df) <- Methy_df$cg
  Methy_df <- Methy_df[,-1]
  Methy_df_gn <- Methy_df[match(res_pmatch$cg_id,rownames(Methy_df),nomatch = 0),]
  gene_cg <- res_pmatch[match(rownames(Methy_df),res_pmatch$cg_id,nomatch = 0),]
  tissue_Methy_beta <- c()
  onegene <- unique(gene_cg$gene)[1]
  for (onegene in unique(gene_cg$gene)) {
    cg_str <- gene_cg[which(gene_cg$gene == onegene),"cg_id"]
    methy_data <- Methy_df_gn[cg_str,]
    methy_res <- colSums(methy_data)/length(cg_str)
    tissue_Methy_beta <- rbind(tissue_Methy_beta,methy_res)
    print(onegene)
  }
  rownames(tissue_Methy_beta) <- unique(gene_cg$gene)
  Epithelial_down_gene_methy_df <- cbind(rownames(tissue_Methy_beta),tissue_Methy_beta)
  colnames(Epithelial_down_gene_methy_df)[1] <- "gene_name"
  write.table(Epithelial_down_gene_methy_df,paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/Epithelial_",cancer,"_MethyBeta.txt"),sep = "\t",quote = F,row.names = F)
}

#####Identifying genes that increase with age in normal tissues####
methy_clinic <- read.table("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/GTEX_Methy_Clin_data.txt",sep = "\t",header = T)
for (cancer in Normal_tissue_Epithelial$cancaer) {
  Epithelial_Methy_df <- as.data.frame(fread(paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/Epithelial_",cancer,"_MethyBeta.txt"),sep = "\t")) 
  rownames(Epithelial_Methy_df) <- Epithelial_Methy_df$V1
  Epithelial_Methy_df <- Epithelial_Methy_df[,-1]
  clinic_methy_epi <- methy_clinic[match(colnames(Epithelial_Methy_df),methy_clinic$Sample_title,nomatch = 0),]
    results <- data.frame(gene=NA,estimate=NA,std.error=NA,statistic=NA,p.value=NA)
  for (gene in rownames(Epithelial_Methy_df)) {
    Epithelial_Methy_tmp <- Epithelial_Methy_df[rownames(Epithelial_Methy_df) == gene,]
    Epithelial_Methy_tmp <- as.data.frame(t(Epithelial_Methy_tmp))
    colnames(Epithelial_Methy_tmp) <- "gene"
    df_tmp <- merge(Epithelial_Methy_tmp, clinic_methy_epi[,c(1,2,4)], by.x = "row.names", by.y = "Sample_title")
    colnames(df_tmp) <- c("patient", colnames(df_tmp)[2:ncol(df_tmp)])
    lm_fit <- lm(formula = as.formula(gene ~ age + Sex), data=df_tmp)  # fit linear model
    summary(lm_fit)
    result <- tidy(lm_fit)
    result_df <- as.data.frame(result[result$term == "age",])
    result_df$term <- gene
    colnames(result_df) <- c("gene", "estimate", "std.error", "statistic", "p.value")
    results <- rbind(results,result_df)
  }
  results <- results[-1,]
  results$q.value <- p.adjust(results$p.value, method = "BH")
  results$Sig <- ifelse(results$q.value < 0.05, TRUE, FALSE)
  write.table(results,paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/",cancer,"_results.txt"),sep = "\t",quote = F,col.names = T,row.names = F)
}

#statistical results
Age_Methy_DEGs <- data.frame(gene_nums=NA,Up_nums=NA,Down_nums=NA,Nosig=NA)
for (cancer in Normal_tissue_Epithelial$cancaer) {
  df <- read.table(paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/",cancer,"_results.txt"),sep = "\t",header = T)
  df$type <- NA
  for(i in 1:nrow(df)){
    if(df$q.value[i] < 0.05 & df$estimate[i] < 0){
      df$type[i] <- "down"
    } else if(df$q.value[i] < 0.05 & df$estimate[i] > 0){
      df$type[i] <- "up"
    } else {
      df$type[i] <- "Nosig"
    }
  }
  
  df_total <- data.frame(gene_nums=length(rownames(df)),
                         Up_nums=length(rownames(df)[which(df$type == "up")]),
                         Down_nums=length(rownames(df)[which(df$type == "down")]),
                         Nosig=length(rownames(df)[which(df$type == "Nosig")]))
  rownames(df_total) <- cancer
  Age_Methy_DEGs <- rbind(Age_Methy_DEGs,df_total)
}
###Fig.5A-----
Age_Methy_DEGs <- Age_Methy_DEGs[-1,]
write.table(Age_Methy_DEGs,paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/Age_Methy_DEGs.csv"),sep = ",",quote = F,row.names = T,col.names = T)
Age_Methy_DEGs <- read.table(paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/Age_Methy_DEGs.csv"),sep = ",",header = T)
data <- reshape2::melt(Age_Methy_DEGs[,-2])
data$cancer <- factor(data$cancer,levels = rev(Age_Methy_DEGs$cancer))
data$variable <- factor(data$variable,levels = c("Up_nums","Down_nums","Nosig"))
paired=c("#DD5746","#4793AF","#FFEFEF")
ggplot(data,aes(cancer,value,fill=variable))+
  labs(x='tissue',cex=10,y='case')+theme_test(base_size = 20)+
  scale_fill_manual(values = paired)+
  geom_bar(stat="identity",position="fill",width = 0.8,size = 0.25)+
  theme(axis.text.x=element_text(angle=90, vjust=0.5))




########Tumor-differential methylation analysis of epithelial cells UTDA genes in cancers######
load(file="./Aging/data/UCSC/GTEX_RNAseq_data/细胞类型识别/组织统计/Epithelial_cell/Epithelial_cancer_intersect.rda")
Epithelial_sig <- c("UCEC","BRCA","LUAD","LUSC","COAD","READ","PAAD","THCA")
Epithelial_cancer <- intersect_gene[Epithelial_sig]
for (cancer in Epithelial_sig) {
  TCGA_methy <- as.data.frame(fread(paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UTDA_methy/",cancer,"methybeta.txt"),sep = "\t"))
  rownames(TCGA_methy) <- TCGA_methy$gene_name
  TCGA_methy <- TCGA_methy[,-1]
  type <- ifelse(substr(colnames(TCGA_methy),14,14) == "0","Tumor","Normal")
  sample_type <- data.frame(sample_id=colnames(TCGA_methy),
                            type= type)
  sample_type$type <- factor(sample_type$type,levels = c("Tumor","Normal"))
  Methy_data=as.data.frame(t(TCGA_methy))
  sample_id <- rownames(Methy_data)
  Methy_data <- as.data.frame(lapply(Methy_data,as.numeric))
  rownames(Methy_data) <- sample_id
  identical(rownames(Methy_data),sample_type$sample_id)
  
  Methy_data <- cbind(Methy_data,sample_type$type)
  colnames(Methy_data)<- gsub("\\.", "-",colnames(Methy_data))
  genes <- rownames(TCGA_methy)
  total<-data.frame(gene_name=genes,p.value=NA,detalR=NA)
  #Rank sum test to obtain p-value
  for (gene in genes) {
    print(gene)
    res<- wilcox.test(unlist(Methy_data[,gene])~type,data = Methy_data)###rank-sum test
    total[which(total$gene_name == gene),2]<- res$p.value
  }
  total$q.value <- p.adjust(total$p.value, method = "BH")
  cancer_methy <- TCGA_methy[,sample_type$sample_id[which(sample_type$type == "Tumor")]]
  normal_methy <- TCGA_methy[,sample_type$sample_id[which(sample_type$type == "Normal")]]
  #deltaR
  for (gene in genes) {
    print(gene)
    Onecancer_Methy <- cancer_methy[gene,]
    Onenormal_Methy <- normal_methy[gene,]
    resDetalR <- rowMeans(Onecancer_Methy)-rowMeans(Onenormal_Methy)
    total[which(total$gene_name == gene),3] <- resDetalR
  }
  total$Sig <- ifelse(total$q.value < 0.05, TRUE, FALSE)
  total$type <- NA
  for(i in 1:nrow(total)){
    if(total$q.value[i] < 0.05 & total$detalR[i] < 0){
      total$type[i] <- "down"
    } else if(total$q.value[i] < 0.05 & total$detalR[i] > 0){
      total$type[i] <- "up"
    } else {
      total$type[i] <- "Nosig"
    }
  }
  total$type <- "Nosige"
  total$p.value <- 0
  total$detalR <- 0
  write.table(total,paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UTDA_methy/Dif_analysis/",cancer,"_diffRes.txt"),sep = "\t",quote = F,row.names = F,col.names = T)
}

#statistical results
All_diff_Methy <- data.frame(gene_nums=NA,Up_nums=NA,Down_nums=NA,Nosig=NA)
cancer <- tumor_normal[1]
for(cancer in tumor_normal){
  one_cancer_total <-  read.table(paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UTDA_methy/Dif_analysis/",cancer,"_diffRes.txt"),sep = "\t",header = T)
  diff_methy_total <- data.frame(gene_nums=length(rownames(one_cancer_total)),
                                 Up_nums=length(rownames(one_cancer_total)[which(one_cancer_total$type == "up")]),
                                 Down_nums=length(rownames(one_cancer_total)[which(one_cancer_total$type == "down")]),
                                 Nosig=length(rownames(one_cancer_total)[which(one_cancer_total$type == "Nosig")]))
  rownames(diff_methy_total) <- cancer
  All_diff_Methy <- rbind(All_diff_Methy,diff_methy_total)
}

All_diff_Methy <- All_diff_Methy[-1,]

####Fig.5B------
data <- reshape2::melt(All_diff_Methy[,-2])
data$cancer <- factor(data$cancer,levels = rev(All_diff_Methy$cancer))
data$variable <- factor(data$variable,levels = c("Down_nums","Up_nums","Nosig"))
paired=c("#4793AF","#DD5746","#FFEFEF")
ggplot(data,aes(cancer,value,fill=variable))+
  labs(x='tissue',cex=10,y='case')+theme_test(base_size = 20)+
  scale_fill_manual(values = paired)+
  geom_bar(stat="identity",position="fill",color ="black",width = 0.8,size = 0.25)+
  theme(axis.text.x=element_text(angle=90, vjust=0.5))



#####Age DMGs and cancer DMGs overlap####
library(GeneOverlap)
load("./Aging/data/UCSC/GTEX_RNAseq_data/aging_gene_list_new.Rdata")
all_inter_list <- list()
for (cancer in c("LUAD","LUSC","COAD","READ")) {
  #Normal
  df <- read.table(paste0("./Aging/data/UCSC/DNA Methy/GTEx_Methylation/Normal_tisse_Methy_all/UTDA_epithelial_methy/",cancer,"_results.txt"),sep = "\t",header = T)
  df$type <- NA
  for(i in 1:nrow(df)){
    if(df$q.value[i] < 0.05 & df$estimate[i] < 0){
      df$type[i] <- "down"
    } else if(df$q.value[i] < 0.05 & df$estimate[i] > 0){
      df$type[i] <- "up"
    } else {
      df$type[i] <- "Nosig"
    }
  }
  #cancer
  one_cancer_total <-  read.table(paste0("./Aging/data/UCSC/DNA Methy/TCGA_Methylation/UTDA_methy/Dif_analysis/",cancer,"_diffRes.txt"),sep = "\t",header = T)
  UpInAge_methy <- df$gene[which(df$type=="up")]
  DownIncancer_methy <- one_cancer_total$gene_name[which(one_cancer_total$type=="down")]
  genes <- intersect(UpInAge_methy,DownIncancer_methy)
  inter_list <- list(genes)
  names(inter_list) <- cancer
  all_inter_list <- c(all_inter_list,inter_list)
}










