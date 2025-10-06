library(PARSE)
library(Seurat)
library(SingleCellExperiment)

load(file = "islet/raw/raw_islet.Rdata")
Z = read.table("islet/raw/GSE81608_series_matrix.txt",header = TRUE)
Z0 = t(Z)[-1,c(10:15)]
phenotypic = as.data.frame(Z0)
rownames(phenotypic) = 1:1600
colnames(phenotypic) = c("donor id","condition","age","ethnicity","gender","cell subtype")
phenotypic[,1] = sapply(strsplit(Z0[,1],": "),function(x){x[2]})
phenotypic[,2] = sapply(strsplit(Z0[,2],": "),function(x){x[2]})
phenotypic[,3] = as.numeric(sapply(strsplit(Z0[,3],": "),function(x){x[2]}))
phenotypic[,4] = sapply(strsplit(Z0[,4],": "),function(x){x[2]})
phenotypic[,5] = sapply(strsplit(Z0[,5],": "),function(x){x[2]})
phenotypic[,6] = sapply(strsplit(Z0[,6],": "),function(x){x[2]})

load("islet/raw/islet_gene.rda")
mu_name = islet_gene
gene_info = read.csv("islet/raw/homologene.data",sep='',header = FALSE)
gene_select = as.numeric(mu_name)
gene_pos = rep(NA,length(gene_select))
gene_no = as.numeric(gene_info[,3])


for(i in 1:length(gene_select)){
  gene_pos[i] = which.min(abs(gene_no-gene_select[i]))
}


human_islets = phenotypic[,-1]
true.label = factor(human_islets$`cell subtype`)
zero_count <- colSums( y == 0 )
n = nrow(y)

##########################
# Islet Data
#########################
data = log(y+1)
data = data[,(n - zero_count) >= 0.01 * n]
dim(data)
saveRDS(data,"islet/islet_data.rds")

##########################
# Normalized Islet Data
#########################
seurat_obj <- CreateSeuratObject(counts = t(y[,(n - zero_count) >= 0.01 * n]))

# Normalize using defaults (LogNormalize, scale.factor=1e4)
seurat_obj <- NormalizeData(seurat_obj)

# Optional: inspect normalized values
data_norm <- t(GetAssayData(seurat_obj, layer = "data"))
dim(data_norm)

saveRDS(data_norm,"islet/islet_data_norm.rds")

##########################
# Gene Name
#########################
gene_name = gene_info[gene_pos,4]
gene_name = gene_name[(n - zero_count) >= 0.01 * n]
saveRDS(gene_name, "islet/gene_name.rds")


