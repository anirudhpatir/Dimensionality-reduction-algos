# Initialize --------------------------------------------------------------

# Library
library(NMF)
library(jackstraw)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(gridExtra)

# Data
data(esGolub)
dim(esGolub)

# PCA ---------------------------------------------------------------------

# remove genes with 0 SD
rm = which(apply(exprs(esGolub),1,function(x)sd(x)==0))
esGolub = esGolub[-rm,]

# PCA
pca = prcomp(t(exprs(esGolub)), scale = T)

# Variance explained
par(mfrow = c(1,2))
plot(summary(pca)$importance[2,], type = "o")
plot(summary(pca)$importance[3,], type = "o")

# relevant PCs
pca_sig = permutationPA(exprs(esGolub))

# Cos2 - importance of loadings
pca_cos2 = pca$rotation^2
pca_contribution = apply(pca_cos2,2,function(x) x*100/sum(x))
plot(sort(pca_contribution[,2], decreasing = T), type = "l")
abline(h = 1/ncol(exprs(esGolub)), col = "red")
sum(pca_contribution[,2]>1/ncol(exprs(esGolub)))

# Sginifcanty contributig features
thresh=1/ncol(exprs(esGolub))
thresh=0.14
pca_features = list()
for(i in 1:pca_sig$r){
  tmp = which(pca_contribution[,i]>thresh)
  pca_features[[i]] = rownames(pca_contribution)[tmp]
}


# Visualizing -------------------------------------------------------------

# Annotation
annoc = data.frame(type=esGolub$ALL.AML,cell=esGolub$Cell)
rownames(annoc) = colnames(esGolub)
annor = data.frame(factor(melt(pca_features)[,-1]))
rownames(annor) = unlist(pca_features)
colnames(annor) = "pc"

# Plot important PCs
pheatmap(t(pca$x[,1:pca_sig$r]),
         cluster_rows = T,cluster_cols = F,
         annotation_col = annoc)
ggplot(data.frame(pca$x), aes(PC3,PC5,col=anno$type,shape=anno$cell))+
  geom_point()+theme_light()


# Plot all important features
pheatmap(exprs(esGolub)[rownames(annor),],
         cluster_rows = F,cluster_cols = F,
         annotation_col = annoc, 
         annotation_row = annor,
         scale = "row")

# Correlation between important features, specific to PC 'i'
i = 1
pc_features = pca_features[[i]]
pc_data = exprs(esGolub)[pc_features,]
pc_features_cor = cor(t(pc_data))
pc_features_cor[pc_features_cor<=0]=0
feature_clus = pheatmap(pc_features_cor)

## Incase positive and negatively weighted features
feature_clus = cutree(feature_clus$tree_col,k=2)
annor_pc = data.frame(clus=factor(feature_clus))

# Plot features, specific to PC 'i'
pheatmap(pc_data[order(annor_pc$clus),],
         cluster_rows = F,cluster_cols = F,
         annotation_col = annoc,annotation_row = annor_pc,
         scale = "row")
pc_feature_pos = melt(t(pc_data[rownames(annor_pc),]))
pc_feature_pos$Var1 = as.numeric(as.factor(pc_feature_pos$Var1))
a1=ggplot(pc_feature_pos, aes(Var1,value,fill=Var2))+
  geom_area(show.legend = F)+theme_light()

## Incase positive and negatively weighted features
pc_feature_pos = melt(t(pc_data[rownames(annor_pc)[which(annor_pc$clus==1)],]))
pc_feature_pos$Var1 = as.numeric(as.factor(pc_feature_pos$Var1))
a1=ggplot(pc_feature_pos, aes(Var1,value,fill=Var2))+
  geom_area(show.legend = F)+theme_light()+labs(x = "Index",title="Negative weighted features")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
pc_feature_neg = melt(t(pc_data[rownames(annor_pc)[which(annor_pc$clus==2)],]))
pc_feature_neg$Var1 = as.numeric(as.factor(pc_feature_neg$Var1))
a2=ggplot(pc_feature_neg, aes(Var1,value,fill=Var2))+
  geom_area(show.legend = F)+theme_light()+labs(x = "Index",title="Positive weighted features")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
grid.arrange(a1,a2)
