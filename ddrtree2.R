library(monocle)
library(RColorBrewer)
library(ggplot2)
library(Seurat)
library(igraph)


# load data table
rawdata = read.csv("C:\\Users\\WF\\Desktop\\ddrtree\\18046-data.csv", header = TRUE)
rownames(rawdata) = rawdata$person_id_new
data10k = read.csv("C:\\Users\\WF\\Desktop\\ddrtree\\new_data1.csv", header = TRUE)
data10k$gender = ifelse(rawdata[data10k$person_id_new,"gender"]=="男",1,0)
data10k$age = rawdata[data10k$person_id_new,"age"]
data10k = data10k[!is.na(data10k$gender) & !is.na(data10k$age),]
# features = c('person_id_new','sbp','dbp','HbA1c','total_cholesterol','triglycerides','Creatinine','ALT','HDLC')
# feature_table = t(data10k[features])
# 
# # remove >5 s.d.
# means = rowMeans(feature_table)
# std = apply(feature_table, 1, sd, na.rm=TRUE)
# qc = apply(feature_table, 2, function(v) all(abs(v-means)<=5*std))
# feature_table = feature_table[,qc]
# 
# # rank normalization
# ranked_matrix = t(apply(feature_table, 1, rank))
# ranked_matrix = feature_table
# ranked_matrix = as.matrix(ranked_matrix)
# 
# for (i in 1:8){
#   model = lm(feature_table[i,] ~ data10k$gender[qc] + data10k$age[qc])
#   feature_table[i,] = model$residuals
# }


features = c('person_id_new', 'sbp', 'dbp', 'HbA1c', 'total_cholesterol', 'triglycerides', 'Creatinine', 'ALT', 'HDLC')
feature_table = data10k[features]

# 提取 person_id_new 并从后续计算中排除
person_id_new = feature_table$person_id_new
feature_table = feature_table[,-1]

# 转置以进行后续计算
feature_table = t(feature_table)

# 移除标准差大于5的观测值
means = rowMeans(feature_table, na.rm=TRUE)
std = apply(feature_table, 1, sd, na.rm=TRUE)
qc = apply(feature_table, 2, function(v) all(abs(v-means)<=5*std))
feature_table = feature_table[, qc]

# 排名归一化
ranked_matrix = t(apply(feature_table, 1, rank))

# 使用线性模型调整特征
for (i in 1:nrow(feature_table)) {
  model = lm(ranked_matrix[i, ] ~ data10k$gender[qc] + data10k$age[qc])
  feature_table[i, ] = model$residuals
}

# 将 person_id_new 添加回处理后的数据
processed_data = cbind(person_id_new[qc], t(feature_table))

# 保存为 CSV 文件
write.csv(processed_data, file = "processed_data.csv", row.names = FALSE)





# meta data
seurat_obj <- CreateSeuratObject(counts = feature_table)
pd <- new('AnnotatedDataFrame', data = seurat_obj@meta.data)
fData <- data.frame(gene_short_name = features, row.names = features)
fd <- new('AnnotatedDataFrame', data = fData)

# reduce dimension
cds <- newCellDataSet(as(ranked_matrix, "sparseMatrix"), phenoData = pd,featureData = fd)
cds <- estimateSizeFactors(cds)
cds$Size_Factor <- 1
cds <- reduceDimension(cds, max_components = 2, method='DDRTree')
# 提取降维后的坐标
reduced_dims <- reducedDimS(cds)
x_coords <- reduced_dims[1,]
y_coords <- reduced_dims[2,]
print(x_coords)
print(y_coords)
# 打印前几个坐标点（例如前5个）
# 创建一个包含坐标的新数据框
coordinates_df <- data.frame(X = reduced_dims[1,], Y = reduced_dims[2,])

# 将数据框写入 CSV 文件
write.csv(coordinates_df, file = "reduced_dimensions.csv", row.names = FALSE)








# viz
plot_cell_trajectory(cds, show_tree = FALSE, color_by = rank(data10k$gender[qc])) + labs(color = "hdlc") + 
  scale_color_gradient(low = "green", high = "#8B008B")#,limits = c(0,4))

# pairwise Spearman plot

pheatmap::pheatmap(abs(cor(data10k[features], method = "spearman")), display_numbers = TRUE)

