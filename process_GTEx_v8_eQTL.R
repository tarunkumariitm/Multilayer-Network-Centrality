######
# This script builds multitissue gene-gene coexpression network byusing GTEx data.
# The 




options(warn=-1)
#install.packages("R.utils")
#BiocManager::install("EnsDb.Hsapiens.v86")#BiocManager::install("EnsDb.Hsapiens.v79")
library(R.utils)
library(EnsDb.Hsapiens.v86)
setwd("/data/users/cs15d017/backup_from_home/v8_experiments/")
cov_path = ("/data/users/cs15d017/backup_from_home/v8_experiments/code/GTEx_Analysis_v8_eQTL_covariates/")

#tissue_list = c("Adipose_Subcutaneous", "Adrenal_Gland", )

##########################################################


remove_extra_ensembl_ids <- function(tissue_df, ensembl_gene_list, HGNC_gene_mapping)
{
  for (i in length(ensembl_gene_list):1)
  {
    temp = HGNC_gene_mapping[which(HGNC_gene_mapping == ensembl_gene_list[i]), 2]
    if (length(temp)<1)
    {
      tissue_df <- tissue_df[-c(i),]
      print(i)
    }
  }
  tissue_df
}

##########################################################

add_column <- function(tissue_incomplete, column_name)
{
  num_genes = dim(tissue_incomplete)[1]
  random_vector = runif(n = num_genes, min = 0.0001, max = 0.00011)
  df = as.data.frame(tissue_incomplete)
  names(random_vector) <- column_name
  #tissue_incomplete$column_name = random_vector
  tissue_complete = cbind(df, random_vector)
  names(tissue_complete)[names(tissue_complete) == "random_vector"] <- column_name
  data.matrix(tissue_complete)
}


#### ---------------------------------------- ###################
# Read ggene expression data for a particular tssiue from GTEx data  

get_tissue <- function(tissue_name)
{
  if(tissue_name=="pituitary")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Pituitary", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Pituitary", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Pituitary",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)
  }else if(tissue_name == "adrenal")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Adrenal_Gland", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Adrenal_Gland", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Adrenal_Gland",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)
  }else if(tissue_name == "adipose")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Adipose_Visceral_Omentum", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Adipose_Visceral_Omentum", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Adipose_Visceral_Omentum",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)
  }else if(tissue_name == "muscle")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Muscle_Skeletal", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Muscle_Skeletal", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Muscle_Skeletal",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)    
  }else if(tissue_name == "heart")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Heart_Atrial_Appendage", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Heart_Atrial_Appendage", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Heart_Atrial_Appendage",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)    
  }else if(tissue_name == "kidney")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Kidney_Cortex", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Kidney_Cortex", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Kidney_Cortex",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)    
  }else if(tissue_name == "liver")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Liver", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Liver", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Liver",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)    
  }else if(tissue_name == "thyroid")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Thyroid", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Thyroid", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Thyroid",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)    
  }else if(tissue_name == "hypothalamus")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Brain_Hypothalamus", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Brain_Hypothalamus", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Brain_Hypothalamus",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)    
    
  }else if(tissue_name == "ovary")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Ovary", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Ovary", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Ovary",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)    
  }else if(tissue_name == "intestine")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Small_Intestine_Terminal_Ileum", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Small_Intestine_Terminal_Ileum", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Small_Intestine_Terminal_Ileum",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)    
  }else if(tissue_name == "stomach")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Stomach", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Stomach", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Stomach",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)    
  }else if(tissue_name == "pancreas")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Pancreas", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Pancreas", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Pancreas",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)    
  }else if(tissue_name == "uterus")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Uterus", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Uterus", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Uterus",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)    
  }else if(tissue_name == "breast")
  {
    file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Breast_Mammary_Tissue", ".v8.normalized_expression.bed.gz", sep = "")
    data = try(gunzip(file, remove=FALSE), silent = TRUE)
    decompressed_file = paste("/data/users/cs15d017/backup_from_home/v8_experiments/data/GTEx_Analysis_v8_eQTL_expression_matrices/", "Breast_Mammary_Tissue", ".v8.normalized_expression.bed", sep = "")
    covariates = read.delim(paste0(cov_path,"Breast_Mammary_Tissue",".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)    
  }else
  {
    print("Tissue not found")
  }
  
  #----------- Processing Tissue data --------#

  #----------- Adjust for covariates ---------#
  
  tissue_df <- as.data.frame(read.csv(decompressed_file,header = T, sep="\t",stringsAsFactors=FALSE))
  ensembl_gene_list = strtrim(tissue_df$gene_id, 15)
  HGNC_gene_mapping = ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl_gene_list, keytype = "GENEID", columns = c("SYMBOL"))
  
  tissue_df = remove_extra_ensembl_ids(tissue_df, ensembl_gene_list, HGNC_gene_mapping)
  headers = ensembldb::select(EnsDb.Hsapiens.v86, keys= strtrim(tissue_df$gene_id, 15), keytype = "GENEID", columns = c("SYMBOL"))[,2]
  
  
  gene_tpm = t(tissue_df[,5:dim(tissue_df)[2]])
  colnames(gene_tpm) = headers
  
  gene_tpm = t(gene_tpm)
  colnames(gene_tpm) = gsub("\\.", "-", colnames(gene_tpm))
  
  #covariates = read.delim(paste0(cov_path,tissue_list[2],".v8.covariates.txt"), sep = '\t', header = T, stringsAsFactors = F, check.names = F)
  gene_tpm <- gene_tpm[,colnames(covariates)[-1]]
  
  old_row <- covariates$ID
  old_col <- colnames(covariates)
  cov2 <- data.frame( t(covariates[,-1]))
  colnames(cov2) <- c(old_row)
  
  i <- 1
  residuals <- apply(gene_tpm, 1, function(y) {
    # y <- unname(unlist(genematrix[i, -c(1)]))
    if ((i%%1000) == 0)
    {
    cat(i, '\n')
    }
    i <<- i+1
    temp <- data.frame(y=y, cov2)
    r <- residuals(lm(y~., temp))
    # cat(str(r))
    return(r)
  })

  residuals
}


############################################################



print("preparing adipose tissue")
adipose_vis <- get_tissue("adipose")
print("preparing liver tissue")
liver <- get_tissue("liver")
print("preparing muscle tissue")
muscle <- get_tissue("muscle")
print("preparing intestine tissue")
intestine <- get_tissue("intestine")
print("preparing thyroid tissue")
thyroid <- get_tissue("thyroid")
print("preparing ovary tissue")
ovary <- get_tissue("ovary")
print("preparing breast tissue")
breast <- get_tissue("breast")
print("preparing hypothalamus tissue")
hypothalamus <- get_tissue("hypothalamus")
hypothalamus = add_column(hypothalamus, "LEP")
print("preparing kidney tissue")
kidney <- get_tissue("kidney")
print("preparing pancreas tissue")
pancreas <- get_tissue("pancreas")
print("preparing pituitary tissue")
pituitary <- get_tissue("pituitary")
print("preparing adrenal tissue")
adrenal <- get_tissue("adrenal")
print("preparing uterus tissue")
uterus <- get_tissue("uterus")
print("preparing stomach tissue")
stomach <- get_tissue("stomach")
print("preparing heart tissue")
heart <- get_tissue("heart")

#-----------------------------------------------#

##### Obtain common rows and column IDs corresponding to each hormone (tissue pair)

acth = Reduce(intersect, list(rownames(pituitary),rownames(adrenal)))#  
adiponectin = Reduce(intersect, list(rownames(adipose_vis),rownames(muscle)))#
adrenaline = Reduce(intersect, list(rownames(adrenal),rownames(heart)))#
adrenocorticotropic = Reduce(intersect, list(rownames(pituitary),rownames(adrenal)))#
aldosterone = Reduce(intersect, list(rownames(adrenal),rownames(kidney)))#
angiotensin = Reduce(intersect, list(rownames(liver),rownames(pituitary)))#
anp = Reduce(intersect, list(rownames(heart),rownames(kidney)))#
calcitonin = Reduce(intersect, list(rownames(thyroid),rownames(kidney)))#
corticotropin = Reduce(intersect, list(rownames(hypothalamus),rownames(pituitary)))#
cortisol = Reduce(intersect, list(rownames(adrenal),rownames(liver)))#
crh = Reduce(intersect, list(rownames(hypothalamus),rownames(pituitary)))#
estradiol = Reduce(intersect, list(rownames(ovary),rownames(pituitary)))#
follitropin = Reduce(intersect, list(rownames(pituitary),rownames(ovary)))#
gastrin = Reduce(intersect, list(rownames(intestine),rownames(stomach)))#
ghr = Reduce(intersect, list(rownames(hypothalamus),rownames(pituitary)))#
ghrelin = Reduce(intersect, list(rownames(stomach),rownames(hypothalamus)))#
glp1 = Reduce(intersect, list(rownames(intestine),rownames(pancreas)))#
glucagon = Reduce(intersect, list(rownames(pancreas),rownames(liver)))#
hcg = Reduce(intersect, list(rownames(pituitary),rownames(ovary)))#
insulin = Reduce(intersect, list(rownames(pancreas),rownames(muscle)))#
leptin = Reduce(intersect, list(rownames(adipose_vis),rownames(hypothalamus)))#
luteinizing = Reduce(intersect, list(rownames(pituitary),rownames(ovary)))#
norepinephrine = Reduce(intersect, list(rownames(adrenal),rownames(intestine)))#
oxytocin = Reduce(intersect, list(rownames(hypothalamus),rownames(uterus)))#
progesterone = Reduce(intersect, list(rownames(ovary),rownames(uterus)))#
prolactin = Reduce(intersect, list(rownames(pituitary),rownames(breast)))#
relaxin = Reduce(intersect, list(rownames(ovary),rownames(uterus)))#
somatostatin = Reduce(intersect, list(rownames(hypothalamus),rownames(pituitary)))#
somatotrophin = Reduce(intersect, list(rownames(pituitary),rownames(liver)))#
thyrotropin = Reduce(intersect, list(rownames(pituitary),rownames(thyroid)))#
thyroxin = Reduce(intersect, list(rownames(thyroid),rownames(liver)))#
trh = Reduce(intersect, list(rownames(hypothalamus),rownames(pituitary)))#
tsh = Reduce(intersect, list(rownames(pituitary),rownames(thyroid)))#
vasopressin = Reduce(intersect, list(rownames(hypothalamus),rownames(kidney)))#
vitamind = Reduce(intersect, list(rownames(kidney),rownames(intestine)))#


acth_cols = Reduce(intersect, list(colnames(pituitary),colnames(adrenal)))#  
adiponectin_cols = Reduce(intersect, list(colnames(adipose_vis),colnames(muscle)))#
adrenaline_cols = Reduce(intersect, list(colnames(adrenal),colnames(heart)))#
adrenocorticotropic_cols = Reduce(intersect, list(colnames(pituitary),colnames(adrenal)))#
aldosterone_cols = Reduce(intersect, list(colnames(adrenal),colnames(kidney)))#
angiotensin_cols = Reduce(intersect, list(colnames(liver),colnames(pituitary)))#
anp_cols = Reduce(intersect, list(colnames(heart),colnames(kidney)))#
calcitonin_cols = Reduce(intersect, list(colnames(thyroid),colnames(kidney)))#
corticotropin_cols = Reduce(intersect, list(colnames(hypothalamus),colnames(pituitary)))#
cortisol_cols = Reduce(intersect, list(colnames(adrenal),colnames(liver)))#
crh_cols = Reduce(intersect, list(colnames(hypothalamus),colnames(pituitary)))#
estradiol_cols = Reduce(intersect, list(colnames(ovary),colnames(pituitary)))#
follitropin_cols = Reduce(intersect, list(colnames(pituitary),colnames(ovary)))#
gastrin_cols = Reduce(intersect, list(colnames(intestine),colnames(stomach)))#
ghr_cols = Reduce(intersect, list(colnames(hypothalamus),colnames(pituitary)))#
ghrelin_cols = Reduce(intersect, list(colnames(stomach),colnames(hypothalamus)))#
glp1_cols = Reduce(intersect, list(colnames(intestine),colnames(pancreas)))#
glucagon_cols = Reduce(intersect, list(colnames(pancreas),colnames(liver)))#
hcg_cols = Reduce(intersect, list(colnames(pituitary),colnames(ovary)))#
insulin_cols = Reduce(intersect, list(colnames(pancreas),colnames(muscle)))#
leptin_cols = Reduce(intersect, list(colnames(adipose_vis),colnames(hypothalamus)))#
luteinizing_cols = Reduce(intersect, list(colnames(pituitary),colnames(ovary)))#
norepinephrine_cols = Reduce(intersect, list(colnames(adrenal),colnames(intestine)))#
oxytocin_cols = Reduce(intersect, list(colnames(hypothalamus),colnames(uterus)))#
progesterone_cols = Reduce(intersect, list(colnames(ovary),colnames(uterus)))#
prolactin_cols = Reduce(intersect, list(colnames(pituitary),colnames(breast)))#
relaxin_cols = Reduce(intersect, list(colnames(ovary),colnames(uterus)))#
somatostatin_cols = Reduce(intersect, list(colnames(hypothalamus),colnames(pituitary)))#
somatotrophin_cols = Reduce(intersect, list(colnames(pituitary),colnames(liver)))#
thyrotropin_cols = Reduce(intersect, list(colnames(pituitary),colnames(thyroid)))#
thyroxin_cols = Reduce(intersect, list(colnames(thyroid),colnames(liver)))#
trh_cols = Reduce(intersect, list(colnames(hypothalamus),colnames(pituitary)))#
tsh_cols = Reduce(intersect, list(colnames(pituitary),colnames(thyroid)))#
vasopressin_cols = Reduce(intersect, list(colnames(hypothalamus),colnames(kidney)))#
vitamind_cols = Reduce(intersect, list(colnames(kidney),colnames(intestine)))#




###### Module to adjust p-values
# In a supraadjacency matrix, p-values are adjusted for each block independently

adjust_p <- function(n, hormone, prob_matrix)
{
  L <- 2
  N <- L * n
  prob_matrix_adj <- matrix(0,N,N)
  current_prob = matrix(0,n,n)
  print("part 1 done")
  if(TRUE)
  {
    for (i in 1:L)
    {
      for (j in 1:L)
      {
        print("dimension")
        print(dim(current_prob))
        print(((i-1)*n)+1)
        print(((i)*n))
        print(((j-1)*n)+1)
        print(((j)*n))
        current_prob = prob_matrix[(((i-1)*n)+1):(((i)*n)),(((j-1)*n)+1):(((j)*n))]
        prob_matrix_adj[(((i-1)*n)+1):(((i)*n)),(((j-1)*n)+1):(((j)*n))] <- matrix(p.adjust(as.vector(current_prob), method='fdr'),ncol=n)
        print("i")
        print(i)
      }
    }
  }
  print("writing probability matrix")
  #write.table(prob_matrix_adj,file=paste("./data/", hormone, "_eQTL_spearman_prob_matrix_adj_removed_CF_lt_R_latest_15k.csv", sep=""), append = FALSE, sep = ",",row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))
  prob_matrix_adj
}


count_cf <- function(tissue)
{
  library(sva)
  library("WGCNA", quietly = T)
  drops <- c("tiss_dat...2.")
  tissue = tissue[, !(names(tissue) %in% drops)]
  tissue_variances <- apply(X=tissue, MARGIN=2, FUN=var)
  tissue_sorted <- sort(tissue_variances, decreasing=TRUE, index.return=TRUE)$ix[1:10000]
  tissue = tissue[,tissue_sorted]
  mod=matrix(1,nrow=dim(tissue)[1],ncol=1)
  colnames(mod)="Intercept"
  nsv=num.sv(t(tissue),mod, method = "be") ## num.sv requires data matrix with features(genes) in the rows and samples in the column
  print(paste("Number of PCs estimated to be removed:", nsv))
}


#-----------------------#
# Finds correlation matrix and writes it in the memory.

write_corr <- function(hormone)
{
  library(WGCNA)
  library(Hmisc)
  
  if (hormone == "acth")
  {
    source_tissue = pituitary[acth,]
    target_tissue = adrenal[acth,]
    source_tissue= source_tissue[,acth_cols]
    target_tissue= target_tissue[,acth_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
    
  }
  if (hormone == "adiponectin")
  {
    source_tissue = adipose_vis[adiponectin,]
    target_tissue = muscle[adiponectin,]
    source_tissue= source_tissue[,adiponectin_cols]
    target_tissue= target_tissue[,adiponectin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  
  if (hormone == "adrenaline")
  {
    source_tissue = adrenal[adrenaline,]
    target_tissue = heart[adrenaline,]
    source_tissue= source_tissue[,adrenaline_cols]
    target_tissue= target_tissue[,adrenaline_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "adrenocorticotropic")
  {
    source_tissue = pituitary[adrenocorticotropic,]
    target_tissue = adrenal[adrenocorticotropic,]
    source_tissue= source_tissue[,adrenocorticotropic_cols]
    target_tissue= target_tissue[,adrenocorticotropic_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "aldosterone")
  {
    source_tissue = adrenal[aldosterone,]
    target_tissue = kidney[aldosterone,]
    source_tissue= source_tissue[,aldosterone_cols]
    target_tissue= target_tissue[,aldosterone_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "angiotensin")
  {
    source_tissue = liver[angiotensin,]
    target_tissue = pituitary[angiotensin,]
    source_tissue= source_tissue[,angiotensin_cols]
    target_tissue= target_tissue[,angiotensin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "anp")
  {
    source_tissue = heart[anp,]
    target_tissue = kidney[anp,]
    source_tissue= source_tissue[,anp_cols]
    target_tissue= target_tissue[,anp_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "calcitonin")
  {
    source_tissue = thyroid[calcitonin,]
    target_tissue = kidney[calcitonin,]
    source_tissue= source_tissue[,calcitonin_cols]
    target_tissue= target_tissue[,calcitonin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "corticotropin")
  {
    source_tissue = hypothalamus[corticotropin,]
    target_tissue = pituitary[corticotropin,]
    source_tissue= source_tissue[,corticotropin_cols]
    target_tissue= target_tissue[,corticotropin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "cortisol")
  {
    source_tissue = adrenal[cortisol,]
    target_tissue = liver[cortisol,]
    source_tissue= source_tissue[,cortisol_cols]
    target_tissue= target_tissue[,cortisol_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "crh")
  {
    source_tissue = hypothalamus[crh,]
    target_tissue = pituitary[crh,]
    source_tissue= source_tissue[,crh_cols]
    target_tissue= target_tissue[,crh_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "estradiol")
  {
    source_tissue = ovary[estradiol,]
    target_tissue = pituitary[estradiol,]
    source_tissue= source_tissue[,estradiol_cols]
    target_tissue= target_tissue[,estradiol_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "follitropin")
  {
    source_tissue = pituitary[follitropin,]
    target_tissue = ovary[follitropin,]
    source_tissue= source_tissue[,follitropin_cols]
    target_tissue= target_tissue[,follitropin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "gastrin")
  {
    source_tissue = intestine[gastrin,]
    target_tissue = stomach[gastrin,]
    source_tissue= source_tissue[,gastrin_cols]
    target_tissue= target_tissue[,gastrin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "ghr")
  {
    source_tissue = hypothalamus[ghr,]
    target_tissue = pituitary[ghr,]
    source_tissue= source_tissue[,ghr_cols]
    target_tissue= target_tissue[,ghr_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "ghrelin")
  {
    source_tissue = stomach[ghrelin,]
    target_tissue = hypothalamus[ghrelin,]
    source_tissue= source_tissue[,ghrelin_cols]
    target_tissue= target_tissue[,ghrelin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "glp1")
  {
    source_tissue = intestine[glp1,]
    target_tissue = pancreas[glp1,]
    source_tissue= source_tissue[,glp1_cols]
    target_tissue= target_tissue[,glp1_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "glucagon")
  {
    source_tissue = pancreas[glucagon,]
    target_tissue = liver[glucagon,]
    source_tissue= source_tissue[,glucagon_cols]
    target_tissue= target_tissue[,glucagon_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  
  if (hormone == "hcg")
  {
    source_tissue = pituitary[hcg,]
    target_tissue = ovary[hcg,]
    source_tissue= source_tissue[,hcg_cols]
    target_tissue= target_tissue[,hcg_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "insulin")
  {
    source_tissue= pancreas[insulin,]
    target_tissue= muscle[insulin,]
    source_tissue= source_tissue[,insulin_cols]
    target_tissue= target_tissue[,insulin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "leptin")
  {
    source_tissue = adipose_vis[leptin,]
    target_tissue = hypothalamus[leptin,]
    source_tissue= source_tissue[,leptin_cols]
    target_tissue= target_tissue[,leptin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  
  if (hormone == "luteinizing")
  {
    source_tissue = pituitary[luteinizing,]
    target_tissue = ovary[luteinizing,]
    source_tissue= source_tissue[,luteinizing_cols]
    target_tissue= target_tissue[,luteinizing_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "norepinephrine")
  {
    source_tissue = adrenal[norepinephrine,]
    target_tissue = intestine[norepinephrine,]
    source_tissue= source_tissue[,norepinephrine_cols]
    target_tissue= target_tissue[,norepinephrine_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "oxytocin")
  {
    source_tissue = hypothalamus[oxytocin,]
    target_tissue = uterus[oxytocin,]
    source_tissue= source_tissue[,oxytocin_cols]
    target_tissue= target_tissue[,oxytocin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "progesterone")
  {
    source_tissue = ovary[progesterone,]
    target_tissue = uterus[progesterone,]
    source_tissue= source_tissue[,progesterone_cols]
    target_tissue= target_tissue[,progesterone_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "prolactin")
  {
    source_tissue = pituitary[prolactin,]
    target_tissue = breast[prolactin,]
    source_tissue= source_tissue[,prolactin_cols]
    target_tissue= target_tissue[,prolactin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "relaxin")
  {
    source_tissue = ovary[relaxin,]
    target_tissue = uterus[relaxin,]
    source_tissue= source_tissue[,relaxin_cols]
    target_tissue= target_tissue[,relaxin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "somatostatin")
  {
    source_tissue = hypothalamus[somatostatin,]
    target_tissue = pituitary[somatostatin,]
    source_tissue= source_tissue[,somatostatin_cols]
    target_tissue= target_tissue[,somatostatin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "somatotrophin")
  {
    source_tissue = pituitary[somatotrophin,]
    target_tissue = liver[somatotrophin,]
    source_tissue= source_tissue[,somatotrophin_cols]
    target_tissue= target_tissue[,somatotrophin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "thyrotropin")
  {
    source_tissue = pituitary[thyrotropin,]
    target_tissue = thyroid[thyrotropin,]
    source_tissue= source_tissue[,thyrotropin_cols]
    target_tissue= target_tissue[,thyrotropin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "thyroxin")
  {
    source_tissue = thyroid[thyroxin,]
    target_tissue = liver[thyroxin,]
    source_tissue= source_tissue[,thyroxin_cols]
    target_tissue= target_tissue[,thyroxin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "trh")
  {
    source_tissue = hypothalamus[trh,]
    target_tissue = pituitary[trh,]
    source_tissue= source_tissue[,trh_cols]
    target_tissue= target_tissue[,trh_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "tsh")
  {
    source_tissue = pituitary[tsh,]
    target_tissue = thyroid[tsh,]
    source_tissue= source_tissue[,tsh_cols]
    target_tissue= target_tissue[,tsh_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "vasopressin")
  {
    source_tissue = hypothalamus[vasopressin,]
    target_tissue = kidney[vasopressin,]
    source_tissue= source_tissue[,vasopressin_cols]
    target_tissue= target_tissue[,vasopressin_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  if (hormone == "vitamind")
  {
    source_tissue = kidney[vitamind,]
    target_tissue = intestine[vitamind,]
    source_tissue= source_tissue[,vitamind_cols]
    target_tissue= target_tissue[,vitamind_cols]
    print(dim(source_tissue))
    print(dim(target_tissue))
  }
  
  print(hormone)
  snap_genes = c(read.csv(paste("/data/users/cs15d017/gene_mapping/unique_genes_",hormone,".csv", sep = ""), sep=",", header = FALSE, stringsAsFactors = FALSE))
  
  #drops <- c("tiss_dat...2.")
  #source_tissue = source_tissue[ , !(names(source_tissue) %in% drops)]
  #target_tissue = target_tissue[ , !(names(target_tissue) %in% drops)]
  snap_genes_inGTEx = intersect(snap_genes, colnames(source_tissue))
  snap_genes_inGTEx_ind = match(snap_genes_inGTEx, colnames(source_tissue))
  #source_tissue = source_tissue[, colSums(source_tissue != 0) > 0]
  source_variances <- apply(X=source_tissue, MARGIN=2, FUN=var)
  source_sorted <- sort(source_variances, decreasing=TRUE, index.return=TRUE)$ix[1:10000]
  #target_tissue = target_tissue[, colSums(target_tissue != 0) > 0]
  target_variances <- apply(X=target_tissue, MARGIN=2, FUN=var)
  target_sorted <- sort(target_variances, decreasing=TRUE, index.return=TRUE)$ix[1:10000]
  #snap_genes_inGTEx = intersect(snap_genes, colnames(source_tissue))
  common_columns <- intersect(source_sorted, target_sorted)
  union_columns <- union(union(source_sorted, target_sorted),snap_genes_inGTEx_ind)
  source_tissue = source_tissue[,union_columns]
  target_tissue = target_tissue[,union_columns]
  print("source and target tissues found")
  write.csv(colnames(source_tissue), file = paste("./data/", hormone, "_eQTL_gene_list_with_snap_v8_lt_latest_15k.csv", sep=""), row.names = FALSE)
  
  df = cbind(source_tissue, target_tissue)
  
  spearman_corraltion = rcorr(as.matrix(df), type = "spearman")
  corr_matrix <- spearman_corraltion$r
  prob_matrix = spearman_corraltion$P
  print("spearman correlation calculated")
  print(dim(prob_matrix))
  n = length(colnames(source_tissue))
  corr_matrix[is.na(corr_matrix)] <- 0
  prob_matrix[is.na(prob_matrix)] <- 0.99
  prob_matrix_adj = adjust_p(n, hormone, prob_matrix)
  #prob_matrix_adj = prob_matrix
  
  write.table(corr_matrix, file=paste("./data/", hormone, "_eQTL_spearman_corr_matrix_with_snap_removed_CF_lt_R_latest_15k.csv", sep=""), append = FALSE, sep = ",",row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))
  
  print("finding s_sec scores")
  source_genes_full = apply(unique(read.csv(paste("./data/", hormone, "_source_genes_full_latest.csv", sep=""), header = FALSE, stringsAsFactors = FALSE)), 2, toupper)
  target_genes_full = apply(unique(read.csv(paste("./data/", hormone, "_target_genes_full_latest.csv", sep=""), header = FALSE, stringsAsFactors = FALSE)), 2, toupper)
  current_hormone_genes = colnames(source_tissue)
  print(hormone)
  print("length of current_hormone_genes - genes per tissue")
  length(current_hormone_genes)
  print("No. of samples")
  print(dim(source_tissue)[1])
  source_genes = intersect(source_genes_full, current_hormone_genes)
  write.csv(source_genes, file=paste("./data/", hormone, "_eQTL_source_genes_with_snap_lt_human_latest_15k.csv", sep=""), row.names = FALSE)
  print("final source genes")
  print(length(source_genes))
  target_genes = intersect(target_genes_full, current_hormone_genes)
  write.csv(target_genes, file=paste("./data/", hormone, "_eQTL_target_genes_with_snap_lt_human_latest_15k.csv", sep=""), row.names = FALSE)
  print("final target genes")
  print(length(target_genes))
  prob_matrix[prob_matrix == 0] <- 0.00001
  s_sec_indices = match(target_genes, current_hormone_genes)
  s_sec_matrix = prob_matrix[1:n, n+s_sec_indices]
  corr_s_sec = corr_matrix[1:n, n+s_sec_indices]
  prob_s_sec = prob_matrix[1:n, n+s_sec_indices]
  if (dim(as.matrix(s_sec_matrix))[2]<2){
    s_sec = -log(s_sec_matrix)
  }else{
    s_sec = rowSums(-log(s_sec_matrix))
  }
  
  source_s_sec_indices = match(source_genes, current_hormone_genes)
  source_s_sec_matrix = prob_matrix[n:(2*n), source_s_sec_indices]
  if (dim(as.matrix(source_s_sec_matrix))[2]<2){
    source_s_sec = -log(source_s_sec_matrix)
  }else{
    source_s_sec = rowSums(-log(source_s_sec_matrix))
  }
  write.csv(s_sec, file=paste("./output/", hormone, "_eQTL_s_sec_spearman_with_snap_removed_CF_lt_v8_latest_15k.csv", sep=""))
  write.csv(source_s_sec, file=paste("./output/", hormone, "_eQTL_source_s_sec_spearman_with_snap_removed_CF_lt_v8_latest_15k.csv", sep=""))
}



#-------- run the following block with the list of hormones for which you want to obtain the supraadjacency matrix


if(TRUE){
  ptm <- proc.time()
  #hormones = c('acth', 'adiponectin', 'adrenaline', 'adrenocorticotropic', 'aldosterone', 'angiotensin', 'anp', 'calcitonin', 'cortisol', 'crh', 'estradiol', 'follitropin', 'gastrin', 'ghr', 'ghrelin', 'glp1', 'glucagon', 'hcg', 'insulin', 'leptin', 'luteinizing', 'norepinephrine', 'oxytocin', 'progesterone', 'prolactin', 'relaxin', 'somatostatin', 'somatotrophin', 'thyrotropin', 'thyroxin', 'trh', 'tsh', 'vasopressin', 'vitamind')
  hormones = c('adrenaline', 'aldosterone', 'angiotensin', 'cortisol', 'estradiol', 'glucagon', 'insulin', 'leptin', 'norepinephrine', 'progesterone', 'somatotrophin', 'thyroxin', 'vitamind')
  library(foreach)
  library(doParallel)
  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(7) #not to overload your computer: max should be cores - 1
  registerDoParallel(cl)
  
  print("starting parallel code")
  foreach(i=1:length(hormones)) %dopar% {
    write_corr(hormones[i]) #calling a function
    #do other things if you want
  }
  ##stop cluster
  stopCluster(cl)
  print("code complete")
  print(proc.time() - ptm)
}

#print("starting write_corr function")
#write_corr("insulin")
print("Done, done. Yayyy")