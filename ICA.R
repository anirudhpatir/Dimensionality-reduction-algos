# Initialize --------------------------------------------------------------

# Library
library(NMF)
library(ica)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(gridExtra)

# Data
data(esGolub)
dim(esGolub)

# ICA ---------------------------------------------------------------------

# remove genes with 0 SD
rm = which(apply(exprs(esGolub),1,function(x)sd(x)==0))
esGolub = esGolub[-rm,]
dt = as.matrix(exprs(esGolub))

# iCA
nc = 3
ica = icafast(dt,nc=nc)

# significant feature for each component
thresh=0.001
ica_s = ica$S^2
ica_s = ica_s/sum(ica_s)
plot(sort(ica_s),type="l")
rownames(ica_s) = rownames(dt)
ica_features = list()
for(i in 1:nc){
  tmp = which(ica_s[,i]>thresh)
  ica_features[[i]] = rownames(ica_s)[tmp]
}

# Visualizing -------------------------------------------------------------

# Annotation
annoc = data.frame(type=esGolub$ALL.AML,cell=esGolub$Cell)
rownames(annoc) = colnames(esGolub)
annor = data.frame(factor(melt(ica_features)[,-1]))
rownames(annor) = unlist(ica_features)
colnames(annor) = "ic"

# Plot ICAs
colnames(ica$W) = colnames(esGolub)
pheatmap(ica$S, cluster_cols = F,cluster_rows = T)
pheatmap(ica$W, cluster_cols = F,cluster_rows = T,
         annotation_col = annoc)

# Plot all important features
pheatmap(exprs(esGolub)[unique(unlist(ica_features)),],
         cluster_rows = T,cluster_cols = F,
         annotation_col = annoc,
         #annotation_row = annor,
         scale = "row")

# Correlation between important features, specific to PC 'i'
i = 1
ic_features = ica_features[[i]]
ic_data = exprs(esGolub)[ic_features,]
ic_features_cor = cor(t(ic_data))
ic_features_cor[ic_features_cor<=0]=0
feature_clus = pheatmap(ic_features_cor)

## Incase positive and negatively weighted features
feature_clus = cutree(feature_clus$tree_col,k=2)
annor_ic = data.frame(clus=factor(feature_clus))

# Plot features, specific to PC 'i'
pheatmap(ic_data,
         cluster_rows = F,cluster_cols = F,
         annotation_col = annoc,
         scale = "row")
pc_feature = melt(t(ic_data))
pc_feature$Var1 = as.numeric(as.factor(pc_feature$Var1))
ggplot(pc_feature, aes(Var1,value,fill=Var2))+
  geom_area(show.legend = F)+theme_light()

## Incase positive and negatively weighted features
ic_feature_pos = melt(t(ic_data[rownames(annor_ic)[which(annor_ic$clus==1)],]))
ic_feature_pos$Var1 = as.numeric(as.factor(ic_feature_pos$Var1))
a1=ggplot(ic_feature_pos, aes(Var1,value,fill=Var2))+
  geom_area(show.legend = F)+theme_light()+labs(x = "Index",title="Negative weighted features")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
ic_feature_neg = melt(t(ic_data[rownames(annor_ic)[which(annor_ic$clus==2)],]))
ic_feature_neg$Var1 = as.numeric(as.factor(ic_feature_neg$Var1))
a2=ggplot(ic_feature_neg, aes(Var1,value,fill=Var2))+
  geom_area(show.legend = F)+theme_light()+labs(x = "Index",title="Positive weighted features")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
grid.arrange(a1,a2)


