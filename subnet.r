#2021-1-27
#利用random walker 识别网络中的子网：


#识别子网
backgroundNet <- read.table("E:\\课题\\眼科\\AMD\\backgroundNet.txt", sep="\t", header=F, as.is=T)
DEgenes <- read.table("E:\\课题\\眼科\\AMD\\DEgenes.txt", sep="\t", header=F, as.is=T)
DEmiRs <- read.table("E:\\课题\\眼科\\AMD\\DEmiRs.txt", sep="\t", header=F, as.is=T)
DE_gmi <- union(unlist(DEgenes), unlist(DEmiRs))
#差异节点及其直接邻居
#subnet <- backgroundNet[union(which(backgroundNet[,1] %in% as.character(unlist(DEgenes))), which(backgroundNet[,2] %in% as.character(unlist(DEmiRs)))),]

#backgroundNet的邻接矩阵
backgroundNet_nodes <- union(backgroundNet[,1], backgroundNet[,2])
backgroundNet_adj <- matrix(0, nrow=length(backgroundNet_nodes), ncol=length(backgroundNet_nodes))
rownames(backgroundNet_adj) <- backgroundNet_nodes
colnames(backgroundNet_adj) <- backgroundNet_nodes
#无向图
for (i in 1:dim(backgroundNet)[1])
	{
	row_num <- match(backgroundNet[i,1], rownames(backgroundNet_adj))
	col_num <- match(backgroundNet[i,2], colnames(backgroundNet_adj))
	backgroundNet_adj[row_num, col_num] <- 1
	backgroundNet_adj[col_num, row_num] <- 1
	print(i/dim(backgroundNet)[1])
	}

#随机游走
library(RWOAG)
rw_seeds <- intersect(DE_gmi, rownames(backgroundNet_adj))
visProbs <- randomWalk(backgroundNet_adj, rw_seeds)
genelist <- names(visProbs)[order(visProbs, decreasing=T)]

#取randomwalk的genelist的topn基因构成subnet，直到首次达到连通图
n <- 62
l <- matrix(0, nrow=dim(backgroundNet)[1])
genelist_n <- genelist[1:n]
for (i in 1:dim(backgroundNet)[1]){
	l[i,1] <- length(intersect(c(as.character(backgroundNet[i,1]),as.character(backgroundNet[i,2])), genelist_n))
	}
subnet <- backgroundNet[which(l==2),]
subnet_nodes <- union(subnet[,1], subnet[,2])
length(subnet_nodes)
#邻接矩阵
subnet_adj <- matrix(0, nrow=length(subnet_nodes), ncol=length(subnet_nodes))
rownames(subnet_adj) <- subnet_nodes
colnames(subnet_adj) <- subnet_nodes
#subnet 表示互作关系对
#无向图
for (i in 1:dim(subnet)[1])
	{
	row_num <- match(subnet[i,1], rownames(subnet_adj))
	col_num <- match(subnet[i,2], colnames(subnet_adj))
	subnet_adj[row_num, col_num] <- 1
	subnet_adj[col_num, row_num] <- 1
	print(i/dim(subnet)[1])
	}
#利用dfs判断subnet是否是连通图
library(graph)
library(RBGL)
x <- list()
length(x) <- nrow(subnet_adj)
names(x) <- rownames(subnet_adj)
for (i in 1:nrow(subnet_adj)) {x[[i]] <- list(edges=which(subnet_adj[i, ]==1))}
subnet_g1 <- graphNEL(nodes=rownames(subnet_adj), edgeL=x, edgemode="undirected")
isConnected(subnet_g1)
#DFS只能做connect graph
#DFS(subnet_g1, genelist_n[1])
#n为62时首次达到连通图
write.table(subnet, "E:\\课题\\眼科\\AMD\\subnet.txt", sep="\t", row.names=F, col.names=F, quote=F)

#subnet中的gene做功能富集分析
library("org.Hs.eg.db")
library("clusterProfiler")
subnet_gene<-read.table("E:\\课题\\眼科\\AMD\\fig1\\subnet_gene.txt",header=F,sep="\t")
head(subnet_gene)
subnet_gene_id<-bitr(subnet_gene[,1], fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=org.Hs.eg.db)[,2]
ego <- enrichGO(gene = subnet_gene_id, OrgDb = 'org.Hs.eg.db', keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, readable=T)
ego
emapplot(ego)
head(ego)
write.table(ego, "E:\\课题\\眼科\\AMD\\fig1\\ego.txt", sep="\t", quote=F)
barplot(ego, showCategory=20, title="Enriched GO terms")
dotplot(ego, showCategory=20, title="Enriched GO terms")

#识别子网中入度为0和出度为0的节点
subnet_allnodes <- union(subnet[,1], subnet[,2])
#write.table(subnet_allnodes, "E:\\课题\\眼科\\AMD\\subnet_allnodes.txt", sep="\t", row.names=F, col.names=F, quote=F)
subnet_outdegree0 <- as.data.frame(subnet_allnodes[which(!(subnet_allnodes %in% unique(subnet[,1])))])
subnet_indegree0 <- as.data.frame(subnet_allnodes[which(!(subnet_allnodes %in% unique(subnet[,2])))])
write.table(subnet_outdegree0, "E:\\课题\\眼科\\AMD\\subnet_outdegree0.txt", row.names=F, col.names=F)
write.table(subnet_indegree0, "E:\\课题\\眼科\\AMD\\subnet_indegree0.txt", row.names=F, col.names=F)
subnet_seeds <- as.data.frame(union(as.character(unlist(subnet_outdegree0)), as.character(unlist(subnet_indegree0))))
write.table(subnet_seeds, "E:\\课题\\眼科\\AMD\\subnet_seeds.txt", quote=F, row.names=F, col.names=F)

#读入路径
path <- read.csv("E:\\课题\\眼科\\AMD\\fig2\\path.txt", sep="\t", header=F, fill=T)
#AMD_genes <- read.table("E:\\课题\\眼科\\AMD\\AMD genes.txt", header=F, as.is=T)
#AMD_miRs <- read.table("E:\\课题\\眼科\\AMD\\AMD miRNAs.txt", header=F, as.is=T)
#AMD_gmi <- union(unlist(AMD_genes), unlist(AMD_miRs))
path_length <- apply(path, 1, function(x){return(sum(!(x=="")))}) #路径的长度
#write.table(path_length, "E:\\课题\\眼科\\AMD\\fig2\\path_length.txt", sep="\t", row.names=F, col.names=F)
mcn <- read.table("E:\\课题\\眼科\\AMD\\manually curated nodes.txt", sep="\t", header=F)
path_length_known <- apply(path, 1, function(x){return(length(intersect(x, as.character(mcn[,1]))))}) #路径中含有的已知AMD基因和miRNA的数目
#write.table(path_length_known, "E:\\课题\\眼科\\AMD\\fig2\\path_length_known.txt", sep="\t", row.names=F, col.names=F)
path_length_known_percent <- path_length_known/path_length #路径中已知AMD基因和miRNA的百分比
#write.table(path_length_known_percent, "E:\\课题\\眼科\\AMD\\fig2\\path_length_known_percent.txt", sep="\t", row.names=F, col.names=F)
path_length_de <- apply(path, 1, function(x){return(length(intersect(x, as.character(DE_gmi))))}) #路径中含有的差异基因/miRNA的数目
#write.table(path_length_de, "E:\\课题\\眼科\\AMD\\fig2\\path_length_de.txt", sep="\t", row.names=F, col.names=F)


#The fraction of top nodes ranked by the degree in backgroundNet and subnet
d1<-read.table("E:\\My project\\眼科\\AMD\\fig1\\backgroundNet_ranked by degree.txt", sep="\t", header=F, as.is=T) #d1为按degree降序排列的backgroundNet中的节点是否是差异的基因/miRNA（0/1）
d2<-read.table("E:\\My project\\眼科\\AMD\\fig1\\subnet_ranked by degree.txt", sep="\t", header=F, as.is=T) #d2为按degree降序排列的subnet中的节点是否是差异的基因/miRNA（0/1）

for (i in 1:nrow(d2)){
	d2[i,2]<-sum(d2[1:i,1])/i
	}
plot(d2[,2])
plot(d2[round(nrow(d2)*seq(0.01,1,0.01)),2], type="l") #将节点分成100等份，即subnet中top1%，top2%，...，top100%的节点中差异基因/miRNA节点所占的百分比

for (i in 1:nrow(d1)){
	d1[i,2]<-sum(d1[1:i,1])/i
	}
lines(d1[,2])
lines(d1[round(nrow(d1)*seq(0.01,1,0.01)),2]) #将节点分成100等份，即backgroundNet网络中top1%，top2%，...，top100%的节点中差异基因/miRNA节点所占的百分比

#随机1000次，保持节点数目与子网节点数目相同，形成随机子网，然后还是统计degree，按degree降序排列，计算top1%，top2%，...，top100%的节点中差异基因/miRNA节点所占的百分比。随机1000次的结果取平均，画成一条线。
subnet_random_nodes <- sample(backgroundNet_nodes)[1:length(subnet_allnodes)]
subnet_random <- backgroundNet[intersect(which(backgroundNet[,1] %in% subnet_random_nodes), which(backgroundNet[,2] %in% subnet_random_nodes)),]
dim(subnet_random)
#temp <- table(c(subnet_random[,1],subnet_random[,2]))
#subnet_random_nodes_dd <- names(temp)[order(temp,decreasing=T)] #按degree decreasing排列
#为了评价随机找的子网的网络连通性差（都是散点，基本连接不成网络，所以计算网络的聚类系数，注意将孤立节点的聚类系数设为0）
#计算聚类系数，clustering coefficient，CCv=n/Ck2=2n/k(k-1)





#检验powerlaw分布
x<-log(1:11,10)
y<-log(table(as.numeric(table(c(subnet[,1],subnet[,2])))),10)
plot(x,y,type="p",xlab="degree",ylab="number of nodes")
z=lm(y~x)
lines(x,fitted(z))




