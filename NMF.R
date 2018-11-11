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

# NMF ---------------------------------------------------------------------

# remove genes with 0 SD
rm = which(apply(exprs(esGolub),1,function(x)sd(x)==0))
esGolub = esGolub[-rm,]
dt = as.matrix(exprs(esGolub))

# NMF
nmf_dt = nmf(dt, rank = 3)

# matrices
w = basis(nmf_dt)
h = coef(nmf_dt)

# Important features
nmf_f = extractFeatures(nmf_dt)
  
# Visualizing -------------------------------------------------------------

# Annotation
annoc = data.frame(type=esGolub$ALL.AML,cell=esGolub$Cell)
rownames(annoc) = colnames(esGolub)
annor = data.frame(factor(melt(nmf_f)[,-1]))
rownames(annor) = rownames(exprs(esGolub))[unlist(nmf_f)]
colnames(annor) = "NMF"

# Plot important NMFs
pheatmap(h,
         cluster_rows = T,cluster_cols = F,
         annotation_col = annoc)

# Plot all important features
pheatmap(exprs(esGolub)[rownames(annor),],
         cluster_rows = F,cluster_cols = F,
         annotation_col = annoc, 
         annotation_row = annor,
         scale = "row")

# Correlation between important features, specific to PC 'i'
i = 1
nmf_features = nmf_f[[i]]
nmf_data = exprs(esGolub)[nmf_features,]
nmf_features_cor = cor(t(nmf_data))
nmf_features_cor[nmf_features_cor<=0]=0
pheatmap(nmf_features_cor)

# Plot features, specific to PC 'i'
pheatmap(nmf_data,
         cluster_rows = T,cluster_cols = F,
         annotation_col = annoc,
         scale = "row")
nmf_feature_pos = melt(t(nmf_data))
nmf_feature_pos$Var1 = as.numeric(as.factor(nmf_feature_pos$Var1))
ggplot(nmf_feature_pos, aes(Var1,value,fill=Var2))+
  geom_area(show.legend = F)+theme_light()
