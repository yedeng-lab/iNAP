# iNAP
Integrated Network Analysis Pipeline (http://mem.rcees.ac.cn:8081/)

This is a descriptive file to run some network calculations for iNAP, including SparCC, SPIEC-EASI and eLSA/LA. The codes here are mainly dumped from the original workflow of each method and modified minorly for iNAP.

## SparCC
**Details for SparCC can be found at [https://github.com/dlegor/SparCC](dlegor/SparCC), which provides a SparCC script for python3 version. The installation and relevant introduction is available as well.**  

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
**Details for SPIEC-EASI can be found at [https://github.com/zdk123/SpiecEasi](SPIEC-EASI), which provides detailed installation and relevant introduction.**  

**Here we only introduce brief script to obtain SPIEC-EASI results of interdomain interactions in iNAP.**
> use 6 for plant abundance data and 8 for microbial data as majority selection criterion.
```Rscript
majority_otu <- 8
majority_plant <- 6
otu_table <- read.table("otu.txt",header = T, row.names = 1,sep = "\t")
plant_table <- read.table(inputPlant,header = T, row.names = 1,sep = "\t")
samp.name <- interaction(colnames(otu_table),colnames(plant_table))
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
sp.res <- multi.spiec.easi(datalist=list(t(otu_table_MV),t(plant_table_MV)),method="glasso",sel.criterion = "stars",verbose=FALSE,pulsar.select = TRUE,lambda.min.ratio=0.1,nlambda=50,pulsar.params = list(rep.num=50,thresh=0.05,ncores=8),cov.output=TRUE)
opticov <- getOptiCov(sp.res)
colnames(opticov) = rownames(opticov) <- colnames(sp.res[["est"]][["data"]])
write.table(as.matrix(opticov),file = "SpiecEasi_result.txt",sep = "\t",quote=FALSE,col.names=NA)
```


## eLSA/LA
**Details for eLSA/LA can be found at [https://bitbucket.org/charade/elsa/wiki/Home](elsa), which provides detailed installation and relevant introduction.**  

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


