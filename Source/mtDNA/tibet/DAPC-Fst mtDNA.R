library(adegenet)
library(dartR)
library(hierfstat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(raster)

theme_white <- function (base_size = 11, base_family = "") {
  theme_bw() %+replace% 
    theme(
      panel.grid.major  = element_line(color = "white"),
      panel.background = element_rect(color="white", fill = NA),
      panel.border = element_rect(color = "black",fill=NA),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text.y = element_blank(),
      axis.text.x = element_text(color = "black",size=10),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 16)
    )
}

#DAPC ANALYSIS
dist_group=matrix(,16,29)
gene_list=c("ATP6","ATP8", "CO1","CO2","CO3","CYB","ND1","ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "12S_rRNA", "16S_rRNA", "D-loop")
for(i in 1:length(gene_list)){
  data=fasta2genlight(file=paste("gene",i,".fasta", sep = ""), parallel = FALSE)
  samples=read.table("pop.txt")
  o=match(data$ind.names, samples$V1)
  pop=factor(samples$V2[o]) 
  data$pop=pop 
  dapc1=dapc(data, n.pca=100, n.da=2) 
  
  #centroids coordinates
  grp.coord=data.frame(dapc1$grp.coord)
  grp.coord$POP=row.names(grp.coord)
  
  #euclidean distance
  dist_group[i,1]=gene_list[i]
  dist_group[i,2]=round(pointDistance(grp.coord[grp.coord$POP=="0-500", 1:2], grp.coord[grp.coord$POP=="501-1000", 1:2], lonlat = F),3)
  dist_group[i,3]=round(pointDistance(grp.coord[grp.coord$POP=="0-500", 1:2], grp.coord[grp.coord$POP=="1001-1500", 1:2], lonlat = F),3)
  dist_group[i,4]=round(pointDistance(grp.coord[grp.coord$POP=="0-500", 1:2], grp.coord[grp.coord$POP=="1501-2000", 1:2], lonlat = F),3)
  dist_group[i,5]=round(pointDistance(grp.coord[grp.coord$POP=="0-500", 1:2], grp.coord[grp.coord$POP=="2001-2500", 1:2], lonlat = F),3)
  dist_group[i,6]=round(pointDistance(grp.coord[grp.coord$POP=="0-500", 1:2], grp.coord[grp.coord$POP=="2501-3000", 1:2], lonlat = F),3)
  dist_group[i,7]=round(pointDistance(grp.coord[grp.coord$POP=="0-500", 1:2], grp.coord[grp.coord$POP=="3001-4000", 1:2], lonlat = F),3)
  dist_group[i,8]=round(pointDistance(grp.coord[grp.coord$POP=="0-500", 1:2], grp.coord[grp.coord$POP=="4001", 1:2], lonlat = F),3)
  dist_group[i,9]=round(pointDistance(grp.coord[grp.coord$POP=="501-1000", 1:2], grp.coord[grp.coord$POP=="1001-1500", 1:2], lonlat = F),3)
  dist_group[i,10]=round(pointDistance(grp.coord[grp.coord$POP=="501-1000", 1:2], grp.coord[grp.coord$POP=="1501-2000", 1:2], lonlat = F),3)
  dist_group[i,11]=round(pointDistance(grp.coord[grp.coord$POP=="501-1000", 1:2], grp.coord[grp.coord$POP=="2001-2500", 1:2], lonlat = F),3)
  dist_group[i,12]=round(pointDistance(grp.coord[grp.coord$POP=="501-1000", 1:2], grp.coord[grp.coord$POP=="2501-3000", 1:2], lonlat = F),3)
  dist_group[i,13]=round(pointDistance(grp.coord[grp.coord$POP=="501-1000", 1:2], grp.coord[grp.coord$POP=="3001-4000", 1:2], lonlat = F),3)
  dist_group[i,14]=round(pointDistance(grp.coord[grp.coord$POP=="501-1000", 1:2], grp.coord[grp.coord$POP=="4001", 1:2], lonlat = F),3)
  dist_group[i,15]=round(pointDistance(grp.coord[grp.coord$POP=="1001-1500", 1:2], grp.coord[grp.coord$POP=="1501-2000", 1:2], lonlat = F),3)
  dist_group[i,16]=round(pointDistance(grp.coord[grp.coord$POP=="1001-1500", 1:2], grp.coord[grp.coord$POP=="2001-2500", 1:2], lonlat = F),3)
  dist_group[i,17]=round(pointDistance(grp.coord[grp.coord$POP=="1001-1500", 1:2], grp.coord[grp.coord$POP=="2501-3000", 1:2], lonlat = F),3)
  dist_group[i,18]=round(pointDistance(grp.coord[grp.coord$POP=="1001-1500", 1:2], grp.coord[grp.coord$POP=="3001-4000", 1:2], lonlat = F),3)
  dist_group[i,19]=round(pointDistance(grp.coord[grp.coord$POP=="1001-1500", 1:2], grp.coord[grp.coord$POP=="4001", 1:2], lonlat = F),3)
  dist_group[i,20]=round(pointDistance(grp.coord[grp.coord$POP=="1501-2000", 1:2], grp.coord[grp.coord$POP=="2001-2500", 1:2], lonlat = F),3)
  dist_group[i,21]=round(pointDistance(grp.coord[grp.coord$POP=="1501-2000", 1:2], grp.coord[grp.coord$POP=="2501-3000", 1:2], lonlat = F),3)
  dist_group[i,22]=round(pointDistance(grp.coord[grp.coord$POP=="1501-2000", 1:2], grp.coord[grp.coord$POP=="3001-4000", 1:2], lonlat = F),3)
  dist_group[i,23]=round(pointDistance(grp.coord[grp.coord$POP=="1501-2000", 1:2], grp.coord[grp.coord$POP=="4001", 1:2], lonlat = F),3)
  dist_group[i,24]=round(pointDistance(grp.coord[grp.coord$POP=="2001-2500", 1:2], grp.coord[grp.coord$POP=="2501-3000", 1:2], lonlat = F),3)
  dist_group[i,25]=round(pointDistance(grp.coord[grp.coord$POP=="2001-2500", 1:2], grp.coord[grp.coord$POP=="3001-4000", 1:2], lonlat = F),3)
  dist_group[i,26]=round(pointDistance(grp.coord[grp.coord$POP=="2001-2500", 1:2], grp.coord[grp.coord$POP=="4001", 1:2], lonlat = F),3)
  dist_group[i,27]=round(pointDistance(grp.coord[grp.coord$POP=="2501-3000", 1:2], grp.coord[grp.coord$POP=="3001-4000", 1:2], lonlat = F),3)
  dist_group[i,28]=round(pointDistance(grp.coord[grp.coord$POP=="2501-3000", 1:2], grp.coord[grp.coord$POP=="4001", 1:2], lonlat = F),3)
  dist_group[i,29]=round(pointDistance(grp.coord[grp.coord$POP=="3001-4000", 1:2], grp.coord[grp.coord$POP=="4001", 1:2], lonlat = F),3)
  
  #DAPC plot
  pdf(file = paste("DAPC_",gene_list[i],".pdf", sep = ""))
  scatter(dapc1, scree.da=F, bg="white", pch=20,  
          cell=1, cstar=0, solid=.8, cex=1.5, clab=1, leg=TRUE, sub=paste("DAPC on ", gene_list[i], sep = ""), possub = "topleft")
  dev.off()
}

dist_group=data.frame(dist_group)
colnames(dist_group)[1:29]=c("Gene",
                             "0-500 vs 501-1000","0-500 vs 1001-1500","0-500 vs 1501-2000","0-500 vs 2001-2500","0-500 vs 2501-3000","0-500 vs 3001-4000","0-500 vs 4001",
                             "501-1000 vs 1001-1500","501-1000 vs 1501-2000","501-1000 vs 2001-2500","501-1000 vs 2501-3000","501-1000 vs 3001-4000","501-1000 vs 4001",
                             "1001-1500 vs 1501-2000","1001-1500 vs 2001-2500","1001-1500 vs 2501-3000","1001-1500 vs 3001-4000","1001-1500 vs 4001",
                             "1501-2000 vs 2001-2500","1501-2000 vs 2501-3000","1501-2000 vs 3001-4000","1501-2000 vs 4001",
                             "2001-2500 vs 2501-3000","2001-2500 vs 3001-4000","2001-2500 vs 4001",
                             "2501-3000 vs 3001-4000","2501-3000 vs 4001",
                             "3001-4000 vs 4001")

write.table(dist_group, "DPCA distance between centroids.txt", col.names = T, sep="\t", row.names = F, quote = F)


rm(list = setdiff(ls(), c("gene_list", "theme_white")))

#FST ANALYSIS
multi.page=vector(mode="list",length = length(gene_list) * 2)

for(i in 1:length(gene_list)){
  data=fasta2genlight(file=paste("gene",i,".fasta", sep = ""), parallel = FALSE)
  samples=read.table("pop.txt")
  o=match(data$ind.names, samples$V1)
  pop=factor(samples$V2[o]) 
  data$pop=pop
  df=gl2gi(data) #convert from genlight to genind object
  fst=pairwise.neifst(genind2hierfstat(df))
  fst=data.frame(fst)
  fst$POP=row.names(fst)
  fst[,1:8]=round(fst[,1:8],3)
  
  group_0_500_vs_501_1000=matrix(,200,2)
  group_0_500_vs_1001_1500=matrix(,200,2)
  group_0_500_vs_1501_2000=matrix(,200,2)
  group_0_500_vs_2001_2500=matrix(,200,2)
  group_0_500_vs_2501_3000=matrix(,200,2)
  group_0_500_vs_3001_4000=matrix(,200,2)
  group_0_500_vs_4001=matrix(,200,2)
  group_501_1000_vs_1001_1500=matrix(,200,2)
  group_501_1000_vs_1501_2000=matrix(,200,2)
  group_501_1000_vs_2001_2500=matrix(,200,2)
  group_501_1000_vs_2501_3000=matrix(,200,2)
  group_501_1000_vs_3001_4000=matrix(,200,2)
  group_501_1000_vs_4001=matrix(,200,2)
  group_1001_1500_vs_1501_2000=matrix(,200,2)
  group_1001_1500_vs_2001_2500=matrix(,200,2)
  group_1001_1500_vs_2501_3000=matrix(,200,2)
  group_1001_1500_vs_3001_4000=matrix(,200,2)
  group_1001_1500_vs_4001=matrix(,200,2)
  group_1501_2000_vs_2001_2500=matrix(,200,2)
  group_1501_2000_vs_2501_3000=matrix(,200,2)
  group_1501_2000_vs_3001_4000=matrix(,200,2)
  group_1501_2000_vs_4001=matrix(,200,2)
  group_2001_2500_vs_2501_3000=matrix(,200,2)
  group_2001_2500_vs_3001_4000=matrix(,200,2)
  group_2001_2500_vs_4001=matrix(,200,2)
  group_2501_3000_vs_3001_4000=matrix(,200,2)
  group_2501_3000_vs_4001=matrix(,200,2)
  group_3001_4000_vs_4001=matrix(,200,2)
  
  for(j in 1:nrow(group_0_500_vs_501_1000)){
    pop(df) = sample(pop(df))
    fst_toss=pairwise.neifst(genind2hierfstat(df))
    fst_toss=data.frame(fst_toss)
    fst_toss$POP=row.names(fst_toss)
    group_0_500_vs_501_1000[j,1]=j
    group_0_500_vs_1001_1500[j,1]=j
    group_0_500_vs_1501_2000[j,1]=j
    group_0_500_vs_2001_2500[j,1]=j
    group_0_500_vs_2501_3000[j,1]=j
    group_0_500_vs_3001_4000[j,1]=j
    group_0_500_vs_4001[j,1]=j
    group_501_1000_vs_1001_1500[j,1]=j
    group_501_1000_vs_1501_2000[j,1]=j
    group_501_1000_vs_2001_2500[j,1]=j
    group_501_1000_vs_2501_3000[j,1]=j
    group_501_1000_vs_3001_4000[j,1]=j
    group_501_1000_vs_4001[j,1]=j
    group_1001_1500_vs_1501_2000[j,1]=j
    group_1001_1500_vs_2001_2500[j,1]=j
    group_1001_1500_vs_2501_3000[j,1]=j
    group_1001_1500_vs_3001_4000[j,1]=j
    group_1001_1500_vs_4001[j,1]=j
    group_1501_2000_vs_2001_2500[j,1]=j
    group_1501_2000_vs_2501_3000[j,1]=j
    group_1501_2000_vs_3001_4000[j,1]=j
    group_1501_2000_vs_4001[j,1]=j
    group_2001_2500_vs_2501_3000[j,1]=j
    group_2001_2500_vs_3001_4000[j,1]=j
    group_2001_2500_vs_4001[j,1]=j
    group_2501_3000_vs_3001_4000[j,1]=j
    group_2501_3000_vs_4001[j,1]=j
    group_3001_4000_vs_4001[j,1]=j
    
    group_0_500_vs_501_1000[j,2]=fst_toss[fst_toss$POP=="0-500",names(fst_toss)=="X501.1000"]
    group_0_500_vs_1001_1500[j,2]=fst_toss[fst_toss$POP=="0-500",names(fst_toss)=="X1001.1500"]
    group_0_500_vs_1501_2000[j,2]=fst_toss[fst_toss$POP=="0-500",names(fst_toss)=="X1501.2000"]
    group_0_500_vs_2001_2500[j,2]=fst_toss[fst_toss$POP=="0-500",names(fst_toss)=="X2001.2500"]
    group_0_500_vs_2501_3000[j,2]=fst_toss[fst_toss$POP=="0-500",names(fst_toss)=="X2501.3000"]
    group_0_500_vs_3001_4000[j,2]=fst_toss[fst_toss$POP=="0-500",names(fst_toss)=="X3001.4000"]
    group_0_500_vs_4001[j,2]=fst_toss[fst_toss$POP=="0-500",names(fst_toss)=="X4001"]
    group_501_1000_vs_1001_1500[j,2]=fst_toss[fst_toss$POP=="501-1000",names(fst_toss)=="X1001.1500"]
    group_501_1000_vs_1501_2000[j,2]=fst_toss[fst_toss$POP=="501-1000",names(fst_toss)=="X1501.2000"]
    group_501_1000_vs_2001_2500[j,2]=fst_toss[fst_toss$POP=="501-1000",names(fst_toss)=="X2001.2500"]
    group_501_1000_vs_2501_3000[j,2]=fst_toss[fst_toss$POP=="501-1000",names(fst_toss)=="X2501.3000"]
    group_501_1000_vs_3001_4000[j,2]=fst_toss[fst_toss$POP=="501-1000",names(fst_toss)=="X3001.4000"]
    group_501_1000_vs_4001[j,2]=fst_toss[fst_toss$POP=="501-1000",names(fst_toss)=="X4001"]
    group_1001_1500_vs_1501_2000[j,2]=fst_toss[fst_toss$POP=="1001-1500",names(fst_toss)=="X1501.2000"]
    group_1001_1500_vs_2001_2500[j,2]=fst_toss[fst_toss$POP=="1001-1500",names(fst_toss)=="X2001.2500"]
    group_1001_1500_vs_2501_3000[j,2]=fst_toss[fst_toss$POP=="1001-1500",names(fst_toss)=="X2501.3000"]
    group_1001_1500_vs_3001_4000[j,2]=fst_toss[fst_toss$POP=="1001-1500",names(fst_toss)=="X3001.4000"]
    group_1001_1500_vs_4001[j,2]=fst_toss[fst_toss$POP=="1001-1500",names(fst_toss)=="X4001"]
    group_1501_2000_vs_2001_2500[j,2]=fst_toss[fst_toss$POP=="1501-2000",names(fst_toss)=="X2001.2500"]
    group_1501_2000_vs_2501_3000[j,2]=fst_toss[fst_toss$POP=="1501-2000",names(fst_toss)=="X2501.3000"]
    group_1501_2000_vs_3001_4000[j,2]=fst_toss[fst_toss$POP=="1501-2000",names(fst_toss)=="X3001.4000"]
    group_1501_2000_vs_4001[j,2]=fst_toss[fst_toss$POP=="1501-2000",names(fst_toss)=="X4001"]
    group_2001_2500_vs_2501_3000[j,2]=fst_toss[fst_toss$POP=="2001-2500",names(fst_toss)=="X2501.3000"]
    group_2001_2500_vs_3001_4000[j,2]=fst_toss[fst_toss$POP=="2001-2500",names(fst_toss)=="X3001.4000"]
    group_2001_2500_vs_4001[j,2]=fst_toss[fst_toss$POP=="2001-2500",names(fst_toss)=="X4001"]
    group_2501_3000_vs_3001_4000[j,2]=fst_toss[fst_toss$POP=="2501-3000",names(fst_toss)=="X3001.4000"]
    group_2501_3000_vs_4001[j,2]=fst_toss[fst_toss$POP=="2501-3000",names(fst_toss)=="X4001"]
    group_3001_4000_vs_4001[j,2]=fst_toss[fst_toss$POP=="3001-4000",names(fst_toss)=="X4001"]

    if(j%%10==0){print(j)}
  }
  
  group_0_500_vs_501_1000=data.frame(group_0_500_vs_501_1000)
  colnames(group_0_500_vs_501_1000)[1:2]=c("Toss","Fst")
  group_0_500_vs_1001_1500=data.frame(group_0_500_vs_1001_1500)
  colnames(group_0_500_vs_1001_1500)[1:2]=c("Toss","Fst")
  group_0_500_vs_1501_2000=data.frame(group_0_500_vs_1501_2000)
  colnames(group_0_500_vs_1501_2000)[1:2]=c("Toss","Fst")
  group_0_500_vs_2001_2500=data.frame(group_0_500_vs_2001_2500)
  colnames(group_0_500_vs_2001_2500)[1:2]=c("Toss","Fst")
  group_0_500_vs_2501_3000=data.frame(group_0_500_vs_2501_3000)
  colnames(group_0_500_vs_2501_3000)[1:2]=c("Toss","Fst")
  group_0_500_vs_3001_4000=data.frame(group_0_500_vs_3001_4000)
  colnames(group_0_500_vs_3001_4000)[1:2]=c("Toss","Fst")
  group_0_500_vs_4001=data.frame(group_0_500_vs_4001)
  colnames(group_0_500_vs_4001)[1:2]=c("Toss","Fst")
  group_501_1000_vs_1001_1500=data.frame(group_501_1000_vs_1001_1500)
  colnames(group_501_1000_vs_1001_1500)[1:2]=c("Toss","Fst")
  group_501_1000_vs_1501_2000=data.frame(group_501_1000_vs_1501_2000)
  colnames(group_501_1000_vs_1501_2000)[1:2]=c("Toss","Fst")
  group_501_1000_vs_2001_2500=data.frame(group_501_1000_vs_2001_2500)
  colnames(group_501_1000_vs_2001_2500)[1:2]=c("Toss","Fst")
  group_501_1000_vs_2501_3000=data.frame(group_501_1000_vs_2501_3000)
  colnames(group_501_1000_vs_2501_3000)[1:2]=c("Toss","Fst")
  group_501_1000_vs_3001_4000=data.frame(group_501_1000_vs_3001_4000)
  colnames(group_501_1000_vs_3001_4000)[1:2]=c("Toss","Fst")
  group_501_1000_vs_4001=data.frame(group_501_1000_vs_4001)
  colnames(group_501_1000_vs_4001)[1:2]=c("Toss","Fst")
  group_1001_1500_vs_1501_2000=data.frame(group_1001_1500_vs_1501_2000)
  colnames(group_1001_1500_vs_1501_2000)[1:2]=c("Toss","Fst")
  group_1001_1500_vs_2001_2500=data.frame(group_1001_1500_vs_2001_2500)
  colnames(group_1001_1500_vs_2001_2500)[1:2]=c("Toss","Fst")
  group_1001_1500_vs_2501_3000=data.frame(group_1001_1500_vs_2501_3000)
  colnames(group_1001_1500_vs_2501_3000)[1:2]=c("Toss","Fst")
  group_1001_1500_vs_3001_4000=data.frame(group_1001_1500_vs_3001_4000)
  colnames(group_1001_1500_vs_3001_4000)[1:2]=c("Toss","Fst")
  group_1001_1500_vs_4001=data.frame(group_1001_1500_vs_4001)
  colnames(group_1001_1500_vs_4001)[1:2]=c("Toss","Fst")
  group_1501_2000_vs_2001_2500=data.frame(group_1501_2000_vs_2001_2500)
  colnames(group_1501_2000_vs_2001_2500)[1:2]=c("Toss","Fst")
  group_1501_2000_vs_2501_3000=data.frame(group_1501_2000_vs_2501_3000)
  colnames(group_1501_2000_vs_2501_3000)[1:2]=c("Toss","Fst")
  group_1501_2000_vs_3001_4000=data.frame(group_1501_2000_vs_3001_4000)
  colnames(group_1501_2000_vs_3001_4000)[1:2]=c("Toss","Fst")
  group_1501_2000_vs_4001=data.frame(group_1501_2000_vs_4001)
  colnames(group_1501_2000_vs_4001)[1:2]=c("Toss","Fst")
  group_2001_2500_vs_2501_3000=data.frame(group_2001_2500_vs_2501_3000)
  colnames(group_2001_2500_vs_2501_3000)[1:2]=c("Toss","Fst")
  group_2001_2500_vs_3001_4000=data.frame(group_2001_2500_vs_3001_4000)
  colnames(group_2001_2500_vs_3001_4000)[1:2]=c("Toss","Fst")
  group_2001_2500_vs_4001=data.frame(group_2001_2500_vs_4001)
  colnames(group_2001_2500_vs_4001)[1:2]=c("Toss","Fst")
  group_2501_3000_vs_3001_4000=data.frame(group_2501_3000_vs_3001_4000)
  colnames(group_2501_3000_vs_3001_4000)[1:2]=c("Toss","Fst")
  group_2501_3000_vs_4001=data.frame(group_2501_3000_vs_4001)
  colnames(group_2501_3000_vs_4001)[1:2]=c("Toss","Fst")
  group_3001_4000_vs_4001=data.frame(group_3001_4000_vs_4001)
  colnames(group_3001_4000_vs_4001)[1:2]=c("Toss","Fst")
  
  p.value=matrix(,28,2)
  p.value[1,1]="0-500 vs 501-1000"
  p.value[2,1]="0-500 vs 1001-1500"
  p.value[3,1]="0-500 vs 1501-2000"
  p.value[4,1]="0-500 vs 2001-2500"
  p.value[5,1]="0-500 vs 2501-3000"
  p.value[6,1]="0-500 vs 3001-4000"
  p.value[7,1]="0-500 vs 4001"
  p.value[8,1]="501-1000 vs 1001-1500"
  p.value[9,1]="501-1000 vs 1501-2000"
  p.value[10,1]="501-1000 vs 2001-2500"
  p.value[11,1]="501-1000 vs 2501-3000"
  p.value[12,1]="501-1000 vs 3001-4000"
  p.value[13,1]="501-1000 vs 4001"
  p.value[14,1]="1001-1500 vs 1501-2000"
  p.value[15,1]="1001-1500 vs 2001-2500"
  p.value[16,1]="1001-1500 vs 2501-3000"
  p.value[17,1]="1001-1500 vs 3001-4000"
  p.value[18,1]="1001-1500 vs 4001"
  p.value[19,1]="1501-2000 vs 2001-2500"
  p.value[20,1]="1501-2000 vs 2501-3000"
  p.value[21,1]="1501-2000 vs 3001-4000"
  p.value[22,1]="1501-2000 vs 4001"
  p.value[23,1]="2001-2500 vs 2501-3000"
  p.value[24,1]="2001-2500 vs 3001-4000"
  p.value[25,1]="2001-2500 vs 4001"
  p.value[26,1]="2501-3000 vs 3001-4000"
  p.value[27,1]="2501-3000 vs 4001"
  p.value[28,1]="3001-4000 vs 4001"

  p.value[1,2]=1-ecdf(group_0_500_vs_501_1000$Fst)(fst[fst$POP=="0-500", names(fst)=="X501.1000"])
  p.value[2,2]=1-ecdf(group_0_500_vs_1001_1500$Fst)(fst[fst$POP=="0-500", names(fst)=="X1001.1500"])
  p.value[3,2]=1-ecdf(group_0_500_vs_1501_2000$Fst)(fst[fst$POP=="0-500", names(fst)=="X1501.2000"])
  p.value[4,2]=1-ecdf(group_0_500_vs_2001_2500$Fst)(fst[fst$POP=="0-500", names(fst)=="X2001.2500"])
  p.value[5,2]=1-ecdf(group_0_500_vs_2501_3000$Fst)(fst[fst$POP=="0-500", names(fst)=="X2501.3000"])
  p.value[6,2]=1-ecdf(group_0_500_vs_3001_4000$Fst)(fst[fst$POP=="0-500", names(fst)=="X3001.4000"])
  p.value[7,2]=1-ecdf(group_0_500_vs_4001$Fst)(fst[fst$POP=="0-500", names(fst)=="X4001"])
  p.value[8,2]=1-ecdf(group_501_1000_vs_1001_1500$Fst)(fst[fst$POP=="501-1000", names(fst)=="X1001.1500"])
  p.value[9,2]=1-ecdf(group_501_1000_vs_1501_2000$Fst)(fst[fst$POP=="501-1000", names(fst)=="X1501.2000"])
  p.value[10,2]=1-ecdf(group_501_1000_vs_2001_2500$Fst)(fst[fst$POP=="501-1000", names(fst)=="X2001.2500"])
  p.value[11,2]=1-ecdf(group_501_1000_vs_2501_3000$Fst)(fst[fst$POP=="501-1000", names(fst)=="X2501.3000"])
  p.value[12,2]=1-ecdf(group_501_1000_vs_3001_4000$Fst)(fst[fst$POP=="501-1000", names(fst)=="X3001.4000"])
  p.value[13,2]=1-ecdf(group_501_1000_vs_4001$Fst)(fst[fst$POP=="501-1000", names(fst)=="X4001"])
  p.value[14,2]=1-ecdf(group_1001_1500_vs_1501_2000$Fst)(fst[fst$POP=="1001-1500", names(fst)=="X1501.2000"])
  p.value[15,2]=1-ecdf(group_1001_1500_vs_2001_2500$Fst)(fst[fst$POP=="1001-1500", names(fst)=="X2001.2500"])
  p.value[16,2]=1-ecdf(group_1001_1500_vs_2501_3000$Fst)(fst[fst$POP=="1001-1500", names(fst)=="X2501.3000"])
  p.value[17,2]=1-ecdf(group_1001_1500_vs_3001_4000$Fst)(fst[fst$POP=="1001-1500", names(fst)=="X3001.4000"])
  p.value[18,2]=1-ecdf(group_1001_1500_vs_4001$Fst)(fst[fst$POP=="1001-1500", names(fst)=="X4001"])
  p.value[19,2]=1-ecdf(group_1501_2000_vs_2001_2500$Fst)(fst[fst$POP=="1501-2000", names(fst)=="X2001.2500"])
  p.value[20,2]=1-ecdf(group_1501_2000_vs_2501_3000$Fst)(fst[fst$POP=="1501-2000", names(fst)=="X2501.3000"])
  p.value[21,2]=1-ecdf(group_1501_2000_vs_3001_4000$Fst)(fst[fst$POP=="1501-2000", names(fst)=="X3001.4000"])
  p.value[22,2]=1-ecdf(group_1501_2000_vs_4001$Fst)(fst[fst$POP=="1501-2000", names(fst)=="X4001"])
  p.value[23,2]=1-ecdf(group_2001_2500_vs_2501_3000$Fst)(fst[fst$POP=="2001-2500", names(fst)=="X2501.3000"])
  p.value[24,2]=1-ecdf(group_2001_2500_vs_3001_4000$Fst)(fst[fst$POP=="2001-2500", names(fst)=="X3001.4000"])
  p.value[25,2]=1-ecdf(group_2001_2500_vs_4001$Fst)(fst[fst$POP=="2001-2500", names(fst)=="X4001"])
  p.value[26,2]=1-ecdf(group_2501_3000_vs_3001_4000$Fst)(fst[fst$POP=="2501-3000", names(fst)=="X3001.4000"])
  p.value[27,2]=1-ecdf(group_2501_3000_vs_4001$Fst)(fst[fst$POP=="2501-3000", names(fst)=="X4001"])
  p.value[28,2]=1-ecdf(group_3001_4000_vs_4001$Fst)(fst[fst$POP=="3001-4000", names(fst)=="X4001"])
  
  p.value=data.frame(p.value)
  colnames(p.value)[1]="Analysis"
  colnames(p.value)[2]="P-value"
  p.value$`P-value`=as.numeric(as.character(p.value$`P-value`))
  p.value[,2]=round(p.value[,2],3)
  
  #plot Fst results
  p1=ggdensity(group_0_500_vs_501_1000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="0-500", names(fst)=="501-1000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="0-500 vs 501-1000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p1
  
  p2=ggdensity(group_0_500_vs_1001_1500, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="0-500", names(fst)=="1001-1500"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="0-500 vs 1001-1500", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p2
  
  p3=ggdensity(group_0_500_vs_1501_2000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="0-500", names(fst)=="1501-2000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="0-500 vs 1501-2000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p3
  
  p4=ggdensity(group_0_500_vs_2001_2500, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="0-500", names(fst)=="2001-2500"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="0-500 vs 2001-2500", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p4
  
  p5=ggdensity(group_0_500_vs_2501_3000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="0-500", names(fst)=="2501-3000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="0-500 vs 2501-3000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p5
  
  p6=ggdensity(group_0_500_vs_3001_4000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="0-500", names(fst)=="3001-4000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="0-500 vs 3001-4000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p6
  
  p7=ggdensity(group_0_500_vs_4001, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="0-500", names(fst)=="4001"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="0-500 vs 4001", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p7
  
  p8=ggdensity(group_501_1000_vs_1001_1500, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="501-1000", names(fst)=="1001-1500"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="501-1000 vs 1001-1500", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p8
  
  p9=ggdensity(group_501_1000_vs_1501_2000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="501-1000", names(fst)=="1501-2000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="501-1000 vs 1501-2000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p9
  
  p10=ggdensity(group_501_1000_vs_2001_2500, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="501-1000", names(fst)=="2001-2500"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="501-1000 vs 2001-2500", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p10
  
  p11=ggdensity(group_501_1000_vs_2501_3000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="501-1000", names(fst)=="2501-3000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="501-1000 vs 2501-3000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p11
  
  p12=ggdensity(group_501_1000_vs_3001_4000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="501-1000", names(fst)=="3001-4000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="501-1000 vs 3001-4000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p12
  
  p13=ggdensity(group_501_1000_vs_4001, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="501-1000", names(fst)=="4001"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="501-1000 vs 4001", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p13
  
  p14=ggdensity(group_1001_1500_vs_1501_2000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="1001-1500", names(fst)=="1501-2000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="1001-1500 vs 1501-2000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p14
  
  p15=ggdensity(group_1001_1500_vs_2001_2500, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="1001-1500", names(fst)=="2001-2500"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="1001-1500 vs 2001-2500", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p15
  
  p16=ggdensity(group_1001_1500_vs_2501_3000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="1001-1500", names(fst)=="2501-3000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="1001-1500 vs 2501-3000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p16
  
  p17=ggdensity(group_1001_1500_vs_3001_4000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="1001-1500", names(fst)=="3001-4000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="1001-1500 vs 3001-4000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p17
  
  p18=ggdensity(group_1001_1500_vs_4001, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="1001-1500", names(fst)=="4001"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="1001-1500 vs 4001", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p18
  
  p19=ggdensity(group_1501_2000_vs_2001_2500, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="1501-2000", names(fst)=="2001-2500"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="1501-2000 vs 2001-2500", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p19
  
  p20=ggdensity(group_1501_2000_vs_2501_3000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="1501-2000", names(fst)=="2501-3000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="1501-2000 vs 2501-3000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p20
  
  p21=ggdensity(group_1501_2000_vs_3001_4000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="1501-2000", names(fst)=="3001-4000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="1501-2000 vs 3001-4000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p21
  
  p22=ggdensity(group_1501_2000_vs_4001, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="1501-2000", names(fst)=="4001"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="1501-2000 vs 4001", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p22
  
  p23=ggdensity(group_2001_2500_vs_2501_3000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="2001-2500", names(fst)=="2501-3000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="2001-2500 vs 2501-3000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p23
  
  p24=ggdensity(group_2001_2500_vs_3001_4000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="2001-2500", names(fst)=="3001-4000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="2001-2500 vs 3001-4000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p24
  
  p25=ggdensity(group_2001_2500_vs_4001, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="2001-2500", names(fst)=="4001"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="2001-2500 vs 4001", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p25
  
  p26=ggdensity(group_2501_3000_vs_3001_4000, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="2501-3000", names(fst)=="3001-4000"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="2501-3000 vs 3001-4000", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p26
  
  p27=ggdensity(group_2501_3000_vs_4001, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="2501-3000", names(fst)=="4001"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="2501-3000 vs 4001", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p27
  
  p28=ggdensity(group_3001_4000_vs_4001, x="Fst", y="..density..", color = c("blue"),fill=c("blue") ,rug = T)+
    geom_vline(xintercept = fst[fst$POP=="3001-4000", names(fst)=="4001"], color = "blue", size=0.5, linetype="dashed")+
    labs(title="3001-4000 vs 4001", x ="Fst", y = "Density")+
    geom_density(adjust=2)+
    theme_white()
  p28
  
  text=paste("                                                                                        ",gene_list[i], sep= "")
  text.gene <- ggparagraph(text = text, face = "italic", size = 20, color = "black")
  
  stable.gene2 = ggtexttable(p.value, rows = NULL, theme = ttheme(base_size = 5))
  
  gga1=ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, 
                 p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, 
                 ncol = 7, nrow = 4)
  gga4=ggarrange(text.gene, gga1, ncol = 1, nrow = 2, heights = c(0.05, 1))
  multi.page[[2*i-1]]=gga4
  gga4=ggarrange(text.gene, stable.gene2, ncol = 1, nrow = 2, heights = c(0.05, 1))
  multi.page[[2*i]]=gga4
}
ggexport(multi.page, filename = "fst mitochondrial DNA.pdf", heigth=30, width = 15)
