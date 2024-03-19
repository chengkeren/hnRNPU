##该代码主要对两个及两个以上单细胞样本进行合并,并提供了定义分组找差异基因的方法
##载入seurat包
library(dplyr)
library(Seurat)
setwd("G:/20240113_SQY_RIPseq/scRNA-seq")
##读入pbmc数据
load('germ_cell_cluster_new.Rds')


pbmc.data <- Read10X(data.dir = "ctrl_filtered_feature_bc_matrix/")
colnames(pbmc.data) <- paste('ctrl', colnames(pbmc.data), sep = '_')
##创建Seurat对象与数据过滤
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "ctrl", min.cells = 3, min.features = 200)

##读入pbmc1数据
pbmc1.data <- Read10X(data.dir = "cko_filtered_feature_bc_matrix/")
colnames(pbmc1.data) <- paste('cko', colnames(pbmc1.data), sep = '_')
##创建Seurat对象与数据过滤
pbmc1 <- CreateSeuratObject(counts = pbmc1.data, project = "cko", min.cells = 3, min.features = 200)

merged<-merge(pbmc,pbmc1)

####如果3个样进行合并
##读入pbmc2数据
##pbmc2.data <- Read10X(data.dir = "D:/BC9_TUMOR1/")
##创建Seurat对象与数据过滤
##pbmc2 <- CreateSeuratObject(counts = pbmc2.data, project = "TUMOR3", min.cells = 3, min.features = 200)
##merged<-merge(pbmc,c(pbmc1,pbmc2))

##计算每个细胞的线粒体基因转录本数的百分比（%）,使用[[ ]] 操作符存放到metadata中
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^mt-")

##展示基因及线粒体百分比
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(merged, features = "percent.mt",y.max=20)#y.max=20通过调整y.max更好的展示线粒体数量分布


plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

##过滤细胞：保留gene数大于200小于5000的细胞；目的是去掉空GEMs和1个GEMs包含2个以上细胞的数据；而保留线粒体基因的转录本数低于10%的细胞,为了过滤掉死细胞等低质量的细胞数据。
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 15)

##表达量数据标准化,LogNormalize的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000 )
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
#merged <- NormalizeData(merged) 或者用默认的

##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

##提取表达量变变化最高的10个基因；
top10 <- head(VariableFeatures(merged), 10)
top10

plot1 <- VariableFeaturePlot(merged)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

##PCA分析：
#PCA分析数据准备阶段，需要使用ScaleData()进行数据归一化。
##对所有基因进行归一化的方法如下：
#all.genes <- rownames(merged)
#merged <- ScaleData(merged, features = all.genes)

##为了加快速度，用默认参数，选取标准化高变基因（2000个）,速度更快。
merged <- ScaleData(merged)

##如果要消线粒体的效应
#merged <- ScaleData(merged, vars.to.regress = "percent.mt")

##线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
merged <- RunPCA(merged, features = VariableFeatures(object = merged))

##检查PCA分群结果, 这里只展示前5个PC,每个PC只显示5个基因；
print(merged[["pca"]], dims = 1:5, nfeatures = 5)

##展示主成分基因分值
VizDimLoadings(merged, dims = 1:2, reduction = "pca")

##绘制pca散点图
DimPlot(merged, reduction = "pca")

##画第1个或15个主成分的热图；
DimHeatmap(merged, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(merged, dims = 1:15, cells = 500, balanced = TRUE)

##确定数据集的分群个数 
#方法1：Jackstraw置换检验算法；重复取样（原数据的1%）,重跑PCA,鉴定p-value较小的PC；计算‘null distribution’(即零假设成立时)时的基因scores。
#merged <- JackStraw(merged, num.replicate = 100)
#merged <- ScoreJackStraw(merged, dims = 1:20)
#JackStrawPlot(merged, dims = 1:15)

#方法2：肘部图（碎石图）,基于每个主成分对方差解释率的排名。
ElbowPlot(merged, ndims = 50)

##分群个数这里选择10,建议尝试选择多个主成分个数做下游分析,对整体影响不大；在选择此参数时,建议选择偏高的数字（为了获取更多的稀有分群,“宁滥勿缺”）；有些亚群很罕见,如果没有先验知识,很难将这种大小的数据集与背景噪声区分开来。

##非线性降维（UMAP/tSNE)基于PCA空间中的欧氏距离计算nearest neighbor graph,优化任意两个细胞间的距离权重（输入上一步得到的PC维数）。
merged <- FindNeighbors(merged, dims = 1:50)

##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
merged <- FindClusters(merged, resolution = 1.3)

##使用Idents（）函数可查看不同细胞的分群；
head(Idents(merged), 5)

##Seurat提供了几种非线性降维的方法进行数据可视化（在低维空间把相似的细胞聚在一起）,比如UMAP和t-SNE,运行UMAP需要先安装'umap-learn'包,这里不做介绍。
merged <- RunTSNE(merged, dims = 1:50)

##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；

DimPlot(merged,reduction = "tsne",label = TRUE,pt.size = 0.5)
DimPlot(merged,reduction = "tsne",label = TRUE,group.by="orig.ident",pt.size = 0.5)

gene<-c("Dazl", "Ddx4","Etv4","Gfra1","Id4","Stra8","Utf1","Zbtb16","Uchl1","Etv5", "Nanos2", "Dmrt1", "Kit", "Rhox13", "Sall4","Lhx1","Ret","Pou3f1","Oct4","Cdh1", "Lin28a","Sohlh1","Sohlh2", "Nanos3","Piwil4","Dnmt3l","Etv4","Foxo1","Rhox13","Amh", "Wt1","Sox9","Clu","Rhox8","Ctsl","Cst12","Amhr2","Gstm6","Cldn11","Inhba","Inhbb","Gata4","Trf","Testin","N-cadherin","ZO-1","Hsd3b1","Cyp11a1","Star", "Cyp17a1","Acta2","Myh11","Tagln","Fbxl22",
        "Cd68","Cd52","Lyz2", "Adgre1","Tk1","Csf1r","F13a1","Mrc1","Pecam1","Cd34","Esam","Flt1", "Ly6c1","Igf1","Cnmd","Igkc","Ighm")
gene2<-c("Upk1b", "Upk3b", "Lgals7", "Aldh1a2", "Mt2", "Ildr2", "Krt18", "Gpm6a", "Plxna4", "Krt7,Nr5a2", "Bex4", "Lect1", "Fst", "Hsd17b1", "Slc18a2", "Inha", "Serpine2", "Ivns1abp", "Fam13a", "Col1a2", "Col3a1", "Col1a1", "Ogn", "Bgn", "Tcf21", "Dcn", "Pdgfra", "Lum", "Mgp", "Lyz2", "Laptm5", "H2.Aa", "Cd74", "Fcer1g", "Ctss", "H2.Eb1", "C1qb", "C1qc", "Cd52", "Kdr", "Mmrn2", "Cdh5", "Egfl7", "Pecam1", "Cldn5", "Esam", "Flt1", "Cd93", "Ctla2a", "Gm15698", "Gdf9", "H1foo", "Padi6", "Ooep", "Oosp1", "Rfpl4", "Tcl1", "Khdc1b", "Nlrp14")
DoHeatmap(merged, features = gene2) + NoLegend()
save(merged,file = "merge.RData")
sub_pbmc1<-subset(merged, idents = c(5,9,17,21))
DimPlot(sub_pbmc1,reduction = "tsne",label = TRUE,pt.size = 0.5)

##每个cluster按照marker规定细胞类型
new.cluster.ids <- c("Stroma","unknown","PTM","unknown","Sertoli Cells","Sertoli Cells","Germ Cells","Macrophage","PTM","Germ Cells","Germ Cells","PTM","Stroma","Leydig cells","Endothelial cells","Macrophage","Macrophage","unknown","unknown","unknown","unknown","unknown","unknown")
names(new.cluster.ids) <- levels(merged)
merged_cluster <- RenameIdents(merged, new.cluster.ids)
DimPlot(merged_cluster,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(merged_cluster,reduction = "tsne",label = TRUE,group.by="orig.ident",pt.size = 1.5)
DoHeatmap(merged_cluster, features = gene) + NoLegend()
FeaturePlot(merged_cluster,features = "Hnrnpu",cols = c("gray", "red"))
VlnPlot(merged_cluster, pt.size = 0,features ="Hnrnpu")





#提取germ cell可能群体

germ_cell<-subset(merged_cluster, idents = c("Germ Cells"))
sub_pbmc<-germ_cell
DimPlot(sub_pbmc, reduction = "tsne",label = TRUE, pt.size = 1.5)

##########对取出的cluster重新进行分群
sub_pbmc <- FindVariableFeatures(sub_pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sub_pbmc)
sub_pbmc <- ScaleData(sub_pbmc, features = all.genes)

sub_pbmc <- RunPCA(sub_pbmc, features = VariableFeatures(object = sub_pbmc))

sub_pbmc <- FindNeighbors(sub_pbmc, dims = 1:50)

sub_pbmc <- FindClusters(sub_pbmc, resolution = 0.5)

sub_pbmc <- RunTSNE(sub_pbmc, dims = 1:50)
DimPlot(sub_pbmc,reduction = "tsne",label = TRUE,group.by="orig.ident",pt.size = 0.5)
DimPlot(sub_pbmc, reduction = "tsne",label = TRUE, pt.size = 1.5)
sub_pbmc.markers <- FindAllMarkers(sub_pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sub_pbmc.markers,file = "subgermcell_markers.csv")
gene<-c("Dnmt3a","Dnmt3l","Kdm1b","Piwil4", "Rac3","Spocd1","Utf1","Cdkn1c","Dapk1","Ly6h","Palld","Ccnd2","Etv4","Mmp9",
        "Ttc28","Gfra1","Id4","Psmd7","Sfrp2", "Etv5", "Lhx1","Ret","Nanos3","Neurog3","Rarg","Rragd", "Ddit4", "Egr4","Dmrtb1","Kit","Rhox13","Stra8", "Dnmt3b")
DoHeatmap(sub_pbmc, features = gene) + NoLegend()

new.cluster.ids <- c("proSG","RNA metabolism related","RNA metabolism related","Diff_SPG","SSC","proSG","unknown")
names(new.cluster.ids) <- levels(sub_pbmc)
germ_cell_cluster<- RenameIdents(sub_pbmc, new.cluster.ids)
DimPlot(germ_cell_cluster,reduction = "tsne",label = FALSE,pt.size = 1.5)
DimPlot(germ_cell_cluster,reduction = "tsne",label = TRUE,group.by="orig.ident",pt.size = 1.5)
germ_genes<-c("Cdkn1c","Ccnd2","Ttc28","Zfp36l1","Gfra1","Psmd7","Sfrp2", "Etv5", "Lhx1","Ret","Nanos3","Rragd", "Ddit4","Dmrtb1","Kit","Rhox13","Stra8", "Dnmt3b")
DoHeatmap(germ_cell_cluster, features = germ_genes) + NoLegend()

table(Idents(germ_cell_cluster))
table(germ_cell_cluster$orig.ident)
prop.table(table(germ_cell_cluster$orig.ident))
cell.prop<-as.data.frame(prop.table(table(Idents(germ_cell_cluster), germ_cell_cluster$orig.ident)))
colnames(cell.prop)<-c("cluster","group","proportion")
write.csv(cell.prop,file = "germcell_prop.csv")
DotPlot(germ_cell_cluster, features = c("Cdkn1c","Ccnd2","Ttc28","Zfp36l1","Gfra1","Psmd7","Sfrp2", "Etv5", "Lhx1","Ret","Nanos3","Rragd", "Ddit4","Dmrtb1","Kit","Rhox13","Stra8", "Dnmt3b"),cols = c("blue", "red"))
save(germ_cell_cluster,file = "sub_germ_cluster.Rdata")
levels(germ_cell_cluster)<-c("proSG","SSC","Diff_SPG","RNA metabolism related","unknown")
FeaturePlot(germ_cell_cluster, features = c("Cdkn1c","Ccnd2","Ttc28","Zfp36l1","Gfra1","Psmd7","Sfrp2", "Etv5", "Lhx1","Ret","Nanos3","Rragd", "Ddit4","Dmrtb1","Kit","Rhox13","Stra8", "Dnmt3b"),cols = c("gray", "red"))
VlnPlot(germ_cell_cluster, pt.size = 0,features = c("Cdkn1c","Ccnd2","Ttc28","Zfp36l1","Gfra1","Psmd7","Sfrp2", "Etv5", "Lhx1","Ret","Nanos3","Rragd", "Ddit4","Dmrtb1","Kit","Rhox13","Stra8", "Dnmt3b") )









