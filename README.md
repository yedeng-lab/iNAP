# iNAP
Integrated Network Analysis Pipeline (https://inap.denglab.org.cn/)

**Notice: The iNAP website is updated and please make sure you are visiting our [new iNAP website](https://inap.denglab.org.cn/).**

Please see more details in the [iMeta paper](http://www.imeta.science/index.html) and [Chinese Version](https://mp.weixin.qq.com/s/C8w1_EzO9v0-hkOebpE3JQ): 

**Kai Feng, Xi Peng, Zheng Zhang, Songsong Gu, Qing He, Wenli Shen, Zhujun Wang, Danrui Wang, Qiulong Hu, Yan Li, Shang Wang, Ye Deng. iNAP: an integrated Network Analysis Pipeline for microbiome studies. iMeta. 2022,1:e13. https://doi.org/10.1002/imt2.13**

**iNAP procedure file**: https://inap.denglab.org.cn/static/iNAP-Denglab-Jan.2022.pdf

**iNAP tutorial video**: in [Youtube](https://youtu.be/lCb-Nsp5bwM) or [Bilibili](https://www.bilibili.com/video/BV1a3411p72v/).

**Network visualization video**: in [Youtube](https://youtu.be/jAYexCTZYlI) or [Bilibili](https://www.bilibili.com/video/BV1LT4y1i7Sy/).

This is a descriptive file to run some network calculations for iNAP, including SparCC, SPIEC-EASI and eLSA/LA. The codes here are mainly dumped from the original workflow of each method and modified minorly for iNAP.

## SparCC
**Details for SparCC can be found at [dlegor/SparCC](https://github.com/dlegor/SparCC), which provides a SparCC script for python3 version. The installation and relevant introduction is available as well.**  

**Here we only introduce brief script to obtain SparCC results after the Step of Majority selection in iNAP.**  
> Plant_6_microbe_8.txt is the downloaded file from iNAP after majority selection. The default values can be found in iNAP.
```shell
SparCC.py Plant_6_microbe_8.txt  -i 20 -x 10 -t 0.1 -c SparCC_correlation.txt 
mkdir Plant_6_microbe_8
MakeBootstraps.py Plant_6_microbe_8.txt -n 100 -t permutation_#.txt -p Plant_6_microbe_8/ 
for i in `seq 0 99`; 
do 
  SparCC.py Plant_6_microbe_8/permutation_$i.txt -i 20 -x 10 -t 0.1 -c Plant_6_microbe_8/perm_cor_$i.txt  
done
PseudoPvals.py SparCC_correlation.txt Plant_6_microbe_8/perm_cor_#.txt 100 -o SparCC_pseudo_p_value.txt -t two_sided 
```
- 20: Number of inference iteration to average over
- 10: Number of exclusion iterations to remove strongly correlated pairs
- 0.1: Correlation strengh exclusion threshold
- 100: Number of shuffled times
- two_sided: one_sided or two_sided


## SPIEC-EASI
**Details for SPIEC-EASI can be found at [SPIEC-EASI](https://github.com/zdk123/SpiecEasi), which provides detailed installation and relevant introduction.**  

**Here we only introduce brief script to obtain SPIEC-EASI results of interdomain interactions in iNAP.**
> use 6 for plant abundance data and 8 for microbial data as majority selection criterion.
```Rscript
library(SpiecEasi)
majority_otu <- 8
majority_plant <- 6
otu_table <- read.table("otu.txt",header = T, row.names = 1,sep = "\t")
plant_table <- read.table("plant.txt",header = T, row.names = 1,sep = "\t")
samp.name <- intersect(colnames(otu_table),colnames(plant_table))
otu_table[is.na(otu_table)] = 0
otu_table2<-otu_table[,samp.name] 
otu_table2[otu_table2>0] <- 1
otu_table_MV <- otu_table[rowSums(otu_table2)>= majority_otu,]
plant_table[is.na(plant_table)] = 0
plant_table2<-plant_table[,samp.name] 
plant_table2[plant_table2>0] <- 1
plant_table_MV <- plant_table[rowSums(plant_table2)>= majority_plant,]
rename_otu<-function(x,y){
  result<-paste(x,"_",y,sep = "")
}
rownames(otu_table_MV)<-sapply(rownames(otu_table_MV),rename_otu,y="M")
rownames(plant_table_MV)<-sapply(rownames(plant_table_MV),rename_otu,y="P")
sp.res <- multi.spiec.easi(datalist=list(t(otu_table_MV),t(plant_table_MV)),method="glasso",sel.criterion = "stars",verbose=FALSE,pulsar.select = TRUE,lambda.min.ratio=0.1,nlambda=50,pulsar.params = list(rep.num=50,thresh=0.05,ncores=8),cov.output=TRUE) # change ncores=1 if met errors in windows
opticov <- getOptiCov(sp.res)
colnames(opticov) = rownames(opticov) <- colnames(sp.res[["est"]][["data"]])
write.table(as.matrix(opticov),file = "SpiecEasi_result.txt",sep = "\t",quote=FALSE,col.names=NA)
```


## eLSA/LA
**Details for eLSA/LA can be found at [elsa](https://bitbucket.org/charade/elsa/wiki/Home), which provides detailed installation and relevant introduction.**  

**Here we only introduce brief script to obtain LSA results after the Step of Majority selection in iNAP.**
```bash
lsa_compute Plant_6_microbe_8.txt LSA_calculation_result.txt -s 4 -r 1 -d 3 -p perm -x 100 -b 0 -f none -t simple -n robustZ -q scipy 
```
- **4**: number of spots
- **1**: number of replicates for each time spot
- **3**: maximum time delay
- **perm**: Use permutation; **theo** Theoretical approximation; **mix** Use theoretical approximation for pre-screening and then use permutation;
- **100**: Permutation number 100 or precision=0.01/permutation for p-value estimation
- **0**: Number of bootstraps for 95% confidence interval estimation
- **none**: fill up with zeros; **zero** fill with zero order splines; **linear** fill up with linear splines;**quadratic** fill up with quadratic spline;**cubic** fill up with cubic spline; **slinear** fill up with slinear; **nearest** fill up with nearest neighbor
- **simple**: simple averaging; **SD** standard deviation weighted averaging; **Med** simple Median; **MAD** median absolute deviation weighted median;
- **robustZ**: percentileZ normalization + robust estimates (with perm, mix and theo, and must use this for theo and mix, default); **percentile** percentile normalization, including zeros (only with perm); **percentileZ** percentile normalization + Z-normalization; **pnz** percentile normalization, excluding zeros (only with perm); **none** none;
- **scipy**: Qvalue calculation method

## Pearson/Spearman correlations
**Here we provided one alternative correlation calculation scripts for Pearson/Spearman methods using 'WGCNA' package in R program. **

**This script would get similar correlation matrix as in iNAP pipeline with 'keep blank' for missing values, and the P values were adjusted using 'fdr' methods, which depends on you. The P values in iNAP were not adjusted, but you can using similar methods as shown here to adjust.**
```Rscript
library(WGCNA)
otu <- read.table("Plant_6_microbe_8.txt",header = T,row.names = 1,sep="\t")
corres <- corAndPvalue(t(otu),use="pairwise.complete.obs",method = c("spearman"))
corres.r <- corres$cor
corres.p <- corres$p
corres.p.adj <- p.adjust(corres.p[upper.tri(corres.p,diag = F)],method = "fdr")
corres.p.adj.matrix <- matrix(0,nrow=nrow(corres.p),ncol=ncol(corres.p))
corres.p.adj.matrix[lower.tri(corres.p,diag = F)] <- corres.p.adj
corres.p.adj.matrix<-as.matrix(as.dist(corres.p.adj.matrix))
colnames(corres.p.adj.matrix) = rownames(corres.p.adj.matrix) <- rownames(corres.p)
write.table(corres.r,"Correlation_matrix_r_values.txt",sep="\t",col.names = NA)
write.table(corres.p.adj.matrix,"Correlation_matrix_p_adj_values.txt",sep="\t",col.names = NA)
```
- **spearman**: using "pearson" or "spearman"
- **fdr**: using "fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY" or other adjustment methods

