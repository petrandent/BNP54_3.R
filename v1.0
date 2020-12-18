library(gplots)
library(RColorBrewer)
library(gProfileR)
library(randomForest)
library(caret)
library(NbClust)

#file import
filename="GeneExpressionDataset_normalized.tsv"
df=read.csv(filename,header=T,sep="\t")

#conditions and genes
condition=c(rep("WT",10),rep("TG",10),rep("TherA",10),rep("TherB",10),rep("TherC",10),rep("TherD",10))
genes=df[,1]

#function that takes the rownumber of the gene in df and returns anova results
anova_fun=function(x){
gene_exps=df[x,2:ncol(df)]
df_gene=data.frame(Expression=as.numeric(gene_exps),Condition=condition)
fit=aov(Expression~Condition,data = df_gene)
results=TukeyHSD(fit)[[1]][c("WT-TG","WT-TherA","WT-TherB","WT-TherC","WT-TherD"),c("diff","p adj")]
return(results)
}

#this list contains diff and p_adj values for each WT-treatment combination for each gene.the result is a list of genes with the above mentioned values.
print("If you have a 2GB RAM this will take 10mins...")
result_list=lapply(1:nrow(df),anova_fun)
names(result_list)=genes

#diff and p adj values for all genes per treatment
wt_tg=lapply(result_list,function(x) x[1,])
wt_thera=lapply(result_list,function(x) x[2,])
wt_therb=lapply(result_list,function(x) x[3,])
wt_therc=lapply(result_list,function(x) x[4,])
wt_therd=lapply(result_list,function(x) x[5,])


#function that calculates names of the DEGS for each treatment
diff_exp_genes=function(anova_list){ 
  de=character(0)
  for (i in 1:length(anova_list)){
  if (abs(anova_list[[i]][1])>=1 & anova_list[[i]][2]<=0.05) de=c(de,names(anova_list)[i])
  }
  return(de)
}

#application of diff_exp_genes functions to each treatment
de_wt_tg=diff_exp_genes(wt_tg)
de_wt_thera=diff_exp_genes(wt_thera)
de_wt_therb=diff_exp_genes(wt_therb)
de_wt_therc=diff_exp_genes(wt_therc)
de_wt_therd=diff_exp_genes(wt_therd)

#this contains the union of DEG for all 5 treatments.
DEG_union=union(de_wt_tg,union(de_wt_thera,union(de_wt_therb,union(de_wt_therc,de_wt_therd))))

#removal of redundant variables
rm(de_wt_tg,de_wt_thera,de_wt_therb,de_wt_therc,de_wt_therd,wt_tg,wt_thera,wt_therb,wt_therc,wt_therd,filename)

#this will create a txt file with the answer to 1st question.
write.csv(DEG_union,file="DEG_genes.txt",row.names = F)

#list that contains anova results for DE genes only
result_list2=result_list[DEG_union]

#change of sign as WT is considered the control group and treatments the test groups, in this context a + sign indicates overexpression and a negative underexpression.
diff_list2=lapply(result_list2,function(x) -x[,1]) 
long_vector=unlist(diff_list2)
matrix1=matrix(long_vector,ncol=5,byrow = T)
colnames(matrix1)=c("TG-WT","TherA-WT","TherB-WT","TherC-WT","TherD-WT" )
rownames(matrix1)=names(diff_list2)

jpeg("Heatmap.jpg")
heatmap.2(matrix1,density.info = "none",key.title="log2FC", col=rev(brewer.pal(9,"Spectral")),trace="none",main="DE genes",labRow = F,cexCol = 1.2,ylab="genes",margins = c(5,2),adjCol = c(0.5,0.1),lwid = c(1,3),ColSideColors = c("Red","Green","Yellow","Blue","Orange"),srtCol = 0,offsetCol = c(0.1,0.1,0.1,0.1,0.1))
dev.off()


#K-means
num_clust=NbClust(matrix1,method="ward.D2")
#most optimal number of groups is 2,followed by 3, followed by 5.From recommended number of groups in the question the chosen one is 5.
keymeans=kmeans(matrix1,centers=5,algorithm="Lloyd",iter.max=30)
cluster1_genes=names(which(keymeans$cluster==1))
cluster2_genes=names(which(keymeans$cluster==2))
cluster3_genes=names(which(keymeans$cluster==3))
cluster4_genes=names(which(keymeans$cluster==4))
cluster5_genes=names(which(keymeans$cluster==5))

#creating lists of the genes in each cluster in the clusters.txt file.
write.table(x=c("cluster 1",cluster1_genes,"\n\n"),file="clusters.txt",sep="\t",row.names = F)
write.table(x=c("cluster 2",cluster2_genes,"\n\n"),file="clusters.txt",sep="\t",row.names = F,append = T)
write.table(x=c("cluster 3",cluster3_genes,"\n\n"),file="clusters.txt",sep="\t",row.names = F,append = T)
write.table(x=c("cluster 4",cluster4_genes,"\n\n"),file="clusters.txt",sep="\t",row.names = F,append = T)
write.table(x=c("cluster 5",cluster5_genes,"\n\n"),file="clusters.txt",sep="\t",row.names = F,append = T)


#functional gene Analysis
funcenr1<-gprofiler(query=as.character(cluster1_genes), organism= "mmusculus", src_filter=c("GO:BP", "KEGG"))

funcenr2<-gprofiler(query=as.character(cluster2_genes), organism= "mmusculus", src_filter=c("GO:BP", "KEGG"))

funcenr3<-gprofiler(query=as.character(cluster3_genes), organism= "mmusculus", src_filter=c("GO:BP", "KEGG"))

funcenr4<-gprofiler(query=as.character(cluster4_genes), organism= "mmusculus", src_filter=c("GO:BP", "KEGG"))

funcenr5<-gprofiler(query=as.character(cluster5_genes), organism= "mmusculus", src_filter=c("GO:BP", "KEGG"))

#filtering of gene analysis according to the criteria specified by the question.
filter_funcenr1=funcenr1[funcenr1$term.size<=200 &
                         funcenr1$p.value<=0.01,c(12,3)]

filter_funcenr2=funcenr2[funcenr2$term.size<=200 &
                           funcenr2$p.value<=0.01,c(12,3)]

filter_funcenr3=funcenr3[funcenr3$term.size<=200 &
                           funcenr3$p.value<=0.01,c(12,3)]

filter_funcenr4=funcenr4[funcenr4$term.size<=200 &
                           funcenr4$p.value<=0.01,c(12,3)]

filter_funcenr5=funcenr5[funcenr5$term.size<=200 &
                           funcenr5$p.value<=0.01,c(12,3)]


#sorting of functional gene analysis results by p-value
final1=filter_funcenr1[order(filter_funcenr1$p.value,decreasing = F),]
final2=filter_funcenr2[order(filter_funcenr2$p.value,decreasing = F),]
final3=filter_funcenr3[order(filter_funcenr3$p.value,decreasing = F),]
final4=filter_funcenr4[order(filter_funcenr4$p.value,decreasing = F),]
final5=filter_funcenr5[order(filter_funcenr5$p.value,decreasing = F),]

#creating enrichment.txt file containing results
write.csv(x=c("cluster 1",final1),file="enrichment.txt",row.names = F)
write.table(x=c("cluster 2",final2),file="enrichment.txt",sep="\t",row.names = F,append = T)
write.table(x=c("cluster 3",final3),file="enrichment.txt",sep="\t",row.names = F,append = T)
write.table(x=c("cluster 4",final4),file="enrichment.txt",sep="\t",row.names = F,append = T)
write.table(x=c("cluster 5",final5),file="enrichment.txt",sep="\t",row.names = F,append = T)


#transforming df to a form that can be used with RF. It has DEGs + state columns.
DEG_df=df[df$Gene %in% DEG_union,]
DEG_df2=t(data.matrix(DEG_df[,2:61]))
colnames(DEG_df2)=make.names(DEG_df$Gene)
DEG_df3=as.data.frame(DEG_df2)
DEG_df3$state=as.factor(condition)

#RF
RFmodel=randomForest(state~.,data = DEG_df3,ntree=1000,mtry=3)
RF_CM=RFmodel$confusion #confusion matrix
write.table(RF_CM,"RFtable.txt") #creating RFtable.txt file

#applying varImp and sorting genes by "Overall" DESC
gene_ginidecrease=varImp(RFmodel,useModel="rf")
gene_ginidecrease$genes=rownames(gene_ginidecrease)
gene_sorted=gene_ginidecrease[order(gene_ginidecrease$Overall,decreasing = T),]

#50 most important genes
gene50=gene_sorted[1:50,2]
write(gene50,"fifty_genes.txt")
