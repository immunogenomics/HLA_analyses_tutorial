# Tutorial for statistical test for HLA association with complex traits





## HLA association and fine-mapping

This tutorial corresponds to the section "*HLA association and fine-mapping*" in the manuscript.

We use outputs from the HLA imputation tutorial. Other data are in `data_assoc` directory. Most scripts are in `script_assoc` directory.



**Note!**

If you use Minimac3 imputation in the above section, you should first convert the allele name (CHR:POS) to the original name (HLA_A*XX:XX) in the output VCF file.



If you use `SNP2HLA.csh`,  `SNP2HLA.py` or MIS, you do not have to do this procedure.



```bash
refVCF="/data/srlab/ssakaue/share/forSid/Tutorial_1KGonly"
output="/data/srlab/ssg34/HLA_tutorial/data/imputation/MM3/hgdp_all_chr6.hg19.ba.only.GSA.final.SHAPEIT.imputed"

zcat ${refVCF}.vcf.gz | grep -v "#" | awk '{print $2,$3,$4,$5}' > ${refVCF}.converter

python script_assoc/convert_vcf_allele.py ${output}.dose.vcf.gz ${refVCF}.converter data_assoc/converted.info | bgzip -c > data_assoc/converted.vcf.gz
```



`data_assoc/converted.info` provides R2 and AF information embedded in the VCF file.

`data_assoc/converted.vcf.gz` is the output imputed VCF file with corrected variant names.



### Single-marker test

We first convert imputed genotype to dosage txt file by `PLINK2`.

```bash
imputed="data_assoc/converted"

plink2 \
 --vcf ${imputed}.vcf.gz dosage=DS \
 --make-pgen --out ${imputed}

# this will create ${imputed}.{pgen,psam,pvar}

plink2 --pfile ${imputed} \
  --pheno data_assoc/phenotype.txt \
  --covar data_assoc/covariates.txt \
  --pheno-name trait_name \
  --glm omit-ref hide-covar cols=chrom,pos,ref,alt,test,nobs,beta,se,ci,tz,p,a1freqcc,a1freq \
  --ci 0.95 \
  --out single_marker_assoc \
  --covar-variance-standardize

```



`single_marker_assoc.trait_name.glm.logistic.hybrid` is the output association statistics, and you can extract results for HLA alleles by `grep "HLA_"`, HLA amino acids by `grep "AA_"`, and HLA intragenic SNPs by `grep "SNPS_"`.



You can also do this by using custom R script.

```bash
# when you only test for HLA alleles and amino acids (modify as necessry)
cat ${imputed}.pvar | grep -v "#" | grep -E 'HLA_|AA_' | cut -f3 >  test_markers.txt
cat ${imputed}.pvar | grep -v "#" | grep -E 'HLA_|AA_' | awk '{print $3,"T"}' >  test_markers_alleles.txt # this is to make sure to output the dosages of the "presence" of the allele coded as T, but not the "absence" coded as A.

plink2 --vcf ${imputed}.vcf.gz dosage=DS --export A --extract test_markers.txt --export-allele test_markers_alleles.txt --out ${imputed}
```

`${imputed}.raw` is the table of dosages for HLA alleles and amino acids.



(In `R`)

```r
d <- read.table("data_assoc/converted.raw", header=T)
pheno <- read.table("data_assoc/phenotype.txt", header=T)
cov <- read.table("data_assoc/covariates.txt", header=T)

d <- merge(d, pheno, by = "IID")
d <- merge(d, cov, by = "IID")
d$trait_name <- d$trait_name - 1 # cases to be 1 and controls to be 0

# if you want to test for HLA_A*01:01
summary(glm(trait_name ~ HLA_A.01.01_T + sex + PC1 + PC2, data = d, family = binomial))

```



Then, you can see the same statistics found in PLINK output.

```r
Call:
glm(formula = trait_name ~ HLA_A.01.01_T + sex + PC1 + PC2, family = binomial, data = d)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-1.0457  -0.8305  -0.7443   1.4218   1.9433  

Coefficients:
              Estimate Std. Error z value Pr(>|z|)  
(Intercept)   -0.39449    0.23453  -1.682   0.0926 .
HLA_A.01.01_T -0.03740    0.21231  -0.176   0.8602  
sex           -0.38172    0.15030  -2.540   0.0111 *
PC1            0.05466    0.07488   0.730   0.4654  
PC2           -0.17217    0.07520  -2.289   0.0221 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 1068.1  on 906  degrees of freedom
Residual deviance: 1056.0  on 902  degrees of freedom
AIC: 1066

Number of Fisher Scoring iterations: 4
```



### Omnibus test

First, we extract amino acid polymorphisms from the dosage output and creat `*.raw` file by using `PLINK2`.

In this example, we also apply QC to extract any amino acid polymorphisms with *Rsq* > 0.7.

```bash
sed 1d data_assoc/converted.info | grep "^AA_" | awk '{if($5>0.7)print $1}' > QCed_AA_variants.txt
sed 1d data_assoc/converted.info | grep "^AA_" | awk '{if($5>0.7)print $1,"T"}' > QCed_AA_variants_alleles.txt

plink2 --pfile ${imputed} \
  --extract QCed_AA_variants.txt \
  --export-allele QCed_AA_variants_alleles.txt \
  --export A \
  --out ${imputed}.QCed_AA

# output amino acid names without "_T" and by converting "-" with "minus" as a workaround in R
head -n1 ${imputed}.QCed_AA.raw | cut -f7- | tr "\t" "\n" | sed -e "s/_T$//" | sed 's/-/minus/' > ${imputed}.QCed_AA.allele_name.txt


# run omnibus test by using those input files!
python3 script_assoc/run_omnibus_AAtest.py \
  --aaraw ${imputed}.QCed_AA.raw \
  --out test_output \
  --allele ${imputed}.QCed_AA.allele_name.txt \
  --pheno data_assoc/phenotype.txt \
  --cov data_assoc/covariates.txt \
  --phenoname trait_name \
  --covname sex PC1 PC2

```



In this example, `test_output_omnibus_result.txt` is the output statistics with the tested amino acid positions, deviance explained by including the tested amino acid position into the model, and p value from the ANOVA.

```bash
$ head test_output_omnibus_result.txt

ALLELE_NAME	OMNIBUS_DEVIANCE	OMNIBUS_PVALUE
AA_A_minus22_29910338_exon1	1.25834923435173	0.533031574597542
AA_A_minus15_29910359_exon1	3.52122377743922	0.171939623712512
AA_A_minus11_29910371_exon1	2.15789903921791	0.339952451524512
AA_A_minus2_29910398_exon1	4.53914553021627	0.103356328081297
AA_A_9_29910558_exon2	5.54515225056707	0.593743388681586
AA_A_12_29910567_exon2	0.0109563791947949	0.916635505963766
AA_A_17_29910582_exon2	0.0188695939268655	0.890740998058756
AA_A_19_29910588_exon2	0.000823349387019334	0.977108589575259
AA_A_43_29910660_exon2	0.523234132611833	0.769805751833586
```





### Conditional haplotype test

First, we extract two-field allele dosages from the imputed genotype after QC.

```bash
sed 1d data_assoc/converted.info | grep "^HLA_" | cut -f1,5 | grep ":" | awk '{if($2>0.7)print $1}' > QCed_HLA_tf.txt
sed 1d data_assoc/converted.info | grep "^HLA_" | cut -f1,5 | grep ":" | awk '{if($2>0.7)print $1,"T"}' > QCed_HLA_tf_alleles.txt

plink2 --pfile ${imputed} \
  --extract QCed_HLA_tf.txt \
  --export-allele QCed_HLA_tf_alleles.txt \
  --export A \
  --out ${imputed}.QCed_HLA_tf

# output HLA names without "_T", "*", and ":" as a workaround in R
head -n1 ${imputed}.QCed_HLA_tf.raw | cut -f7- | tr "\t" "\n" | sed -e "s/_T$//" | sed 's/*/_/' | sed 's/:/_/'  > ${imputed}.QCed_HLA_tf.allele_name.txt


```



We first perform the same omnibus test for single AA position but using the two-fiedl allele information.

Let's do this for HLA-DRB1 as an example.

```r
library(dplyr)
library(data.table)

# info file summarizes all information on the correspondence between two-field alleles and amino acid residues
info <- readRDS("data_assoc/HLA_DICTIONARY_AA.hg19.imgt3320.AA_tf.in_ref.rds")

HLA="DRB1"
info <- info[info$gene==HLA,]

info$tag <- paste0(info$pos,":",info$AA)

allpos <- sort( unique(info$pos) )

res <- data.frame()
for( pos in allpos ){
   y <- info[ info$pos %in% c(pos), ]
   for( k in 1:length( unique( y$tag )) ){
      ytag <- unique( y$tag )[ k ]
      hap <- ytag
      y_4d <- subset(y, tag == ytag )$hla
      hap_4d <- y_4d
      if( length(hap_4d) > 0 ){
         out <- data.frame(hap, hla = hap_4d, pos )
         res <- rbind(res, out)
   }}}

# res will be used to group two-field alleles based on amino acid residues at each position.


# read imputed two field alleles
dose <- read.table("data_assoc/converted.QCed_HLA_tf.raw", header=T)
dose <- dose[,c(-1,-3:-6)]
allelenames <- as.character(read.table("data_assoc/converted.QCed_HLA_tf.allele_name.txt")$V1)
colnames(dose) <- c("IID",allelenames)

# read phenotype and covariates
pheno <- read.table("data_assoc/phenotype.txt", header=T)
if(max(pheno[,2])==2){pheno[,2]<-pheno[,2] - 1}
cov <- read.table("data_assoc/covariates.txt", header=T)
dat <- merge(dose, pheno, by = "IID")
dat <- merge(dat, cov, by = "IID")

pval_list<-NULL
deviance_list<-NULL

for(thispos in allpos){
  thishaps<-as.character(unique(subset(res,pos==thispos)$hap))
  adopted<-NULL
  for(thishap in thishaps){
    hlas<-as.character(subset(res,hap==thishap)$hla) # extracting all two-field alleles explained by this position-AA residue pairs
    hlas<-hlas[hlas %in% allelenames] # restrict two-field alleles to those in our QCed data
    if(length(hlas)>0){
      dat$thishap <- rowSums(dat[hlas])
      colnames(dat)[ncol(dat)] <- thishap
      adopted<-c(adopted,thishap)
    }
  }
  obj1<-glm(trait_name ~ sex+PC1+PC2,data=dat,family=binomial(link="logit")) # model with covariates
  obj2<-glm(trait_name~as.matrix(dat[adopted])+sex+PC1+PC2,data=dat,family=binomial(link="logit")) # model with this AA position
  Chisqtest <- anova(obj1,obj2,test="Chisq")
  pval<-Chisqtest$`Pr(>Chi)`[2]
  deviance<-Chisqtest$Deviance[2]
  pval_list<-c(pval_list,pval)
  deviance_list<-c(deviance_list,deviance)
}

summary<-data.frame(POSITION=allpos,OMNIBUS_DEVIANCE=deviance_list,OMNIBUS_PVALUE = pval_list)

# this summary is a summary for the first round of conditional haplotype test.
head(summary)

#  POSITION OMNIBUS_DEVIANCE OMNIBUS_PVALUE
#1      -25        0.1598153      0.9232016
#2      -24        0.1573953      0.9243193
#3      -17        2.5568379      0.2784772
#4      -16        0.1598153      0.9232016
#5       -1        0.1670549      0.9198658
#6        4        4.4417129      0.1085161
```



If the strongest association among all the positions was at position 11, we next run similar analyses conditioned on the position 11.

```r
# This script is continued from the above R workspace.

# 
res.prev <- res # transfer res information from the previous round into res.prev

condpos = "11" # this is a position we want to condition on
thishaps<-as.character(unique(subset(res.prev, pos==condpos)$hap))

thishaps

adopted.prev<-NULL
for(thishap in thishaps){
  hlas<-as.character(subset(res,hap==thishap)$hla)
  hlas<-hlas[hlas %in% allelenames]
    if(length(hlas)>0){
    dat$thishap<-rowSums(dat[hlas])
    colnames(dat)[ncol(dat)]<-thishap
    adopted.prev<-c(adopted.prev,thishap)
}
}

adopted.prev
# [1] "11:L" "11:S" "11:V" "11:G" "11:D" "11:P"
# These are the amino acid residues obverved in the data in the previous round at position 11.

# Let's move on to the next round
allpos <- setdiff(allpos, condpos) # all the other positions to analyse

res <- data.frame()
for( pos in allpos ){
   x <- info[ info$pos %in% c(condpos), ] # haplotype information at the position I want to condition on 
   y <- info[ info$pos %in% c(pos), ] # haplotype information at the position I want to analyze
   for( i in 1:length( unique( x$tag )) ){
   for( k in 1:length( unique( y$tag )) ){
      xtag <- unique( x$tag )[ i ]
      ytag <- unique( y$tag )[ k ]
      hap <- paste0(xtag,"_",ytag)
      x_4d <- subset(x, tag == xtag )$hla
      y_4d <- subset(y, tag == ytag )$hla
      hap_4d <- intersect( x_4d, y_4d )
      if( length(hap_4d) > 0 ){
         out <- data.frame(hap, hla = hap_4d, pos )
         res <- rbind(res, out)
    }}}}

head(res)

#         hap            hla pos
#1 11:L_-25:K HLA_DRB1_01_01 -25
#2 11:L_-25:K HLA_DRB1_01_02 -25
#3 11:L_-25:K HLA_DRB1_01_03 -25
#4 11:S_-25:R HLA_DRB1_03_01 -25
#5 11:S_-25:R HLA_DRB1_03_02 -25
#6 11:S_-25:R HLA_DRB1_08_01 -25

# Haplotypes defined by two amino acid positions and their correspondence to the two-field alleles

pval_list<-NULL
deviance_list<-NULL

for(thispos in allpos){
  thishaps<-as.character(unique(subset(res,pos==thispos)$hap))
  adopted<-NULL
  for(thishap in thishaps){
    hlas<-as.character(subset(res,hap==thishap)$hla)
    hlas<-hlas[hlas %in% allelenames]
    if(length(hlas)>0){
      dat$thishap<-rowSums(dat[hlas])
      colnames(dat)[ncol(dat)]<-thishap
      adopted<-c(adopted,thishap)
    }
  }
  obj1<-glm(trait_name ~ as.matrix(dat[adopted.prev])+sex+PC1+PC2,data=dat,family=binomial(link="logit")) # a model including groups defined by the previous round (single position)
  obj2<-glm(trait_name ~ as.matrix(dat[adopted])+sex+PC1+PC2,data=dat,family=binomial(link="logit")) # a model including groups defined by this round (two positions)
  Chisqtest <- anova(obj1,obj2,test="Chisq")
  pval<-Chisqtest$`Pr(>Chi)`[2]
  deviance<-Chisqtest$Deviance[2]
  pval_list<-c(pval_list,pval)
  deviance_list<-c(deviance_list,deviance)
}

summary<-data.frame(POSITION=allpos,OMNIBUS_DEVIANCE=deviance_list,OMNIBUS_PVALUE = pval_list)


# this "summary" is a summary for the second round of conditional haplotype test.
head(summary)

#  POSITION OMNIBUS_DEVIANCE OMNIBUS_PVALUE
#1      -25        2.6065100      0.1064257
#2      -24        2.6065100      0.1064257
#3      -17        0.0000000             NA
#4      -16        2.6065100      0.1064257
#5       -1        0.8824553      0.3475301
#6        4        0.0000000             NA

```



We continue these procedures until we do not get any significant results (`OMNIBUS_PVALUE`).



- Non-additive association test

We use the same data to test non-additive effect.

In this vignette, we generate and use best-guess genotype from the imputation result.

```bash
plink2 --vcf ${imputed}.vcf.gz \
  --extract QCed_HLA_tf.txt \
  --export-allele QCed_HLA_tf_alleles.txt \
  --export A \
  --out ${imputed}.QCed_HLA_tf.best_guess

# Here we do not specify dosage and export the data to a amtrix

# output HLA names without "_T", "*", and ":" as a workaround in R
head -n1 ${imputed}.QCed_HLA_tf.best_guess.raw | cut -f7- | tr "\t" "\n" | sed -e "s/_T$//" | sed 's/*/_/' | sed 's/:/_/'  > ${imputed}.QCed_HLA_tf.best_guess.allele_name.txt

```



Then, in `R`,

```r
dose <- read.table("data_assoc/converted.QCed_HLA_tf.best_guess.raw", header=T)
dose <- dose[,c(-1,-3:-6)]
allelenames <- as.character(read.table("data_assoc/converted.QCed_HLA_tf.best_guess.allele_name.txt")$V1)
colnames(dose) <- c("IID",allelenames)

# read phenotype and covariates
pheno <- read.table("data_assoc/phenotype.txt", header=T)
if(max(pheno[,2])==2){pheno[,2]<-pheno[,2] - 1} # cases to be 1 and controls to be 0

cov <- read.table("data_assoc/covariates.txt", header=T)
dat <- merge(dose, pheno, by = "IID")
dat <- merge(dat, cov, by = "IID")

# if you want to test for HLA_A*01:01's non-additive effect

tested_allele = "HLA_A_01_01"
dat$non_add <- 0
dat[dat[,tested_allele]==1,]$non_add <- 1 # define non-additive term to have 1 if and only if the genotype is heterozygous.

model.0 <- glm(trait_name ~ as.matrix(dat[,tested_allele]) + sex + PC1 + PC2, data = dat, family = binomial)
model.1 <- glm(trait_name ~ as.matrix(dat[,tested_allele]) + non_add + sex + PC1 + PC2, data = dat, family = binomial)

Chisqtest <- anova(model.0,model.1,test="Chisq") # assess the improvement of the model by including non-additive term
model.pval<-Chisqtest$`Pr(>Chi)`[2]
model.deviance<-Chisqtest$Deviance[2]
SUMC<-summary(model.1)$coefficients
main.beta<-SUMC[2,1] # we get coefficients of the additive term
main.se<-SUMC[2,2]
main.p<-SUMC[2,4]
non_add.beta<-SUMC[3,1] # we get coefficients of the non-additive term
non_add.se<-SUMC[3,2]
non_add.p<-SUMC[3,4]
out<-data.frame(allele=tested_allele, anova.deviance=model.deviance, anova.p=model.pval, main.beta=main.beta, main.se=main.se, main.p=main.p, non_add.beta=non_add.beta, non_add.se=non_add.se, non_add.p=non_add.p)

# this "out" summarizes the statistics at this allele (HLA_A*01:01), both additive and non-additive
out

#       allele anova.deviance   anova.p  main.beta   main.se    main.p
#1 HLA_A_01_01      0.1587487 0.6903112 -0.2378009 0.5629893 0.6727405
#  non_add.beta non_add.se non_add.p
#1    0.2325363  0.6026695 0.6996124
```





### Interaction test

Let's test interaction between two HLA alleles of the same gene.

As we explained in the manuscript, it is sometimes important to restrict the analyses to common alleles. The rare * rare allele combination could yield inflated statistics due to the noisy estimate of the effect sizes for both alleles.

If we decide to QC based on MAF > 0.05 for this interaction analysis,

```bash
plink2 --pfile ${imputed} \
  --extract QCed_HLA_tf.txt \
  --export-allele QCed_HLA_tf_alleles.txt \
  --export A \
  --maf 0.05 \
  --out ${imputed}.QCed_HLA_tf.MAF

# output HLA names without "_T", "*", and ":" as a workaround in R
head -n1 ${imputed}.QCed_HLA_tf.MAF.raw | cut -f7- | tr "\t" "\n" | sed -e "s/_T$//" | sed 's/*/_/' | sed 's/:/_/'  > ${imputed}.QCed_HLA_tf.MAF.allele_name.txt


```



And in `R`

```r
dose <- read.table("data_assoc/converted.QCed_HLA_tf.MAF.raw", header=T)
dose <- dose[,c(-1,-3:-6)]
allelenames <- as.character(read.table("data_assoc/converted.QCed_HLA_tf.MAF.allele_name.txt")$V1)
colnames(dose) <- c("IID",allelenames)

# read phenotype and covariates
pheno <- read.table("data_assoc/phenotype.txt", header=T)
if(max(pheno[,2])==2){pheno[,2]<-pheno[,2] - 1} # cases to be 1 and controls to be 0

cov <- read.table("data_assoc/covariates.txt", header=T)
dat <- merge(dose, pheno, by = "IID")
dat <- merge(dat, cov, by = "IID")


# For example, if you want to test for HLA_DRB1*03:01  and HLA_DRB1*15:01

first_allele = "HLA_DRB1_03_01"
seconde_allele = "HLA_DRB1_15_01"

# It is always nice to check the distribution of samples based on these two alleles.

table(round(dat[,c(first_allele, seconde_allele)]))

#              HLA_DRB1_15_01
#HLA_DRB1_03_01   0   1   2
#             0 624  89   4
#             1 156  14   0
#             2  20   0   0

# We confirmed that the number of samples based on allelic combination of HLA-DRB1.

summary(glm(trait_name ~ HLA_DRB1_03_01 + HLA_DRB1_15_01 + HLA_DRB1_03_01*HLA_DRB1_15_01 +  sex + PC1 + PC2, data = dat, family = binomial))


################################
Call:
glm(formula = trait_name ~ HLA_DRB1_03_01 + HLA_DRB1_15_01 + 
    HLA_DRB1_03_01 * HLA_DRB1_15_01 + sex + PC1 + PC2, family = binomial, 
    data = dat)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-1.0787  -0.8305  -0.7422   1.3944   1.9409  

Coefficients:
                              Estimate Std. Error z value Pr(>|z|)  
(Intercept)                   -0.39243    0.23787  -1.650   0.0990 .
HLA_DRB1_03_01                 0.11719    0.16267   0.720   0.4713  
HLA_DRB1_15_01                -0.40575    0.27059  -1.500   0.1337  
sex                           -0.37822    0.15060  -2.511   0.0120 *
PC1                            0.05676    0.07515   0.755   0.4501  
PC2                           -0.17298    0.07519  -2.300   0.0214 *
HLA_DRB1_03_01:HLA_DRB1_15_01  0.33678    0.70759   0.476   0.6341  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 1068.1  on 906  degrees of freedom
Residual deviance: 1052.7  on 900  degrees of freedom
AIC: 1066.7

Number of Fisher Scoring iterations: 4
################################
 
# "HLA_DRB1_03_01:HLA_DRB1_15_01" is showing the effect of interaction between two alleles.
 
model.0 <- glm(trait_name ~ HLA_DRB1_03_01 + HLA_DRB1_15_01 + sex + PC1 + PC2, data = dat, family = binomial) # a model without an interaction term
model.1 <- glm(trait_name ~ HLA_DRB1_03_01 + HLA_DRB1_15_01 + HLA_DRB1_03_01*HLA_DRB1_15_01 + sex + PC1 + PC2, data = dat, family = binomial) # a model with an interaction term
 
anova(model.0,model.1,test="Chisq")
 
################################
Analysis of Deviance Table

Model 1: trait_name ~ HLA_DRB1_03_01 + HLA_DRB1_15_01 + sex + PC1 + PC2
Model 2: trait_name ~ HLA_DRB1_03_01 + HLA_DRB1_15_01 + HLA_DRB1_03_01 * 
    HLA_DRB1_15_01 + sex + PC1 + PC2
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1       901     1052.9                     
2       900     1052.7  1  0.22015   0.6389
################################
```



It is also considered to run permutation tests ( by breaking the relationship between phenotype and genotype ) whether the observed interaction effects could happen by random chance.



### Multi-trait test

```r
library(MVLM)

# Let's use those MAF-QCed HLA alleles in this example.
dose <- read.table("data_assoc/converted.QCed_HLA_tf.MAF.raw", header=T)
dose <- dose[,c(-1,-3:-6)]
allelenames <- as.character(read.table("data_assoc/converted.QCed_HLA_tf.MAF.allele_name.txt")$V1)
colnames(dose) <- c("IID",allelenames)

# read phenotype of a multiple vector (i.e., a matrix)
pheno <- read.table("data_assoc/phenotype_multi.txt", header=T)

head(pheno)

#        IID    trait.A   trait.B     trait.C    trait.D
#1 HGDP00610 0.47607622 0.2454573 0.219182940 0.05928355
#2 HGDP00982 0.20184835 0.4383009 0.007797759 0.35205303
#3 HGDP00001 0.01328931 0.1882389 0.402212096 0.39625968
#4 HGDP01247 0.33558573 0.2070649 0.268822337 0.18852699
#5 HGDP00309 0.31788202 0.2712361 0.067621683 0.34326023
#6 HGDP00786 0.36187225 0.3786774 0.180634545 0.07881575

# This is comprised of multiple traits, A, B, C, and D, and...

all.equal(rowSums(pheno[,2:5]), rep(1,nrow(pheno)))
#[1] TRUE

# The sum of the traits are (almost) equal to 1 in all samples.
# So we assume this phenotype has 3-degrees of freedom, and thus drop trait.D from the analysis in the subsequent models.

# read covariates
cov <- read.table("data_assoc/covariates.txt", header=T)
dat <- merge(dose, pheno, by = "IID")
dat <- merge(dat, cov, by = "IID")

# Let's see the effect of HLA_DRB1*15:01
# We can perform linear regression
model.0 <- lm( cbind(trait.A,trait.B,trait.C) ~ sex + PC1 + PC2, data = dat)
model.1 <- lm( cbind(trait.A,trait.B,trait.C) ~ HLA_DRB1_15_01 + sex + PC1 + PC2, data = dat)

anova(model.0, model.1)

################################
Analysis of Variance Table

Model 1: cbind(trait.A, trait.B, trait.C) ~ sex + PC1 + PC2
Model 2: cbind(trait.A, trait.B, trait.C) ~ HLA_DRB1_15_01 + sex + PC1 + 
    PC2
  Res.Df Df Gen.var.    Pillai approx F num Df den Df Pr(>F)
1    903    0.015147                                        
2    902 -1 0.015150 0.0028114  0.84579      3    900  0.469
################################


# We can also perform MVLM model
mvlm.res.0 <- mvlm( cbind(trait.A,trait.B,trait.C) ~ sex + PC1 + PC2 , data = dat)

mvlm.res.1 <- mvlm( cbind(trait.A,trait.B,trait.C) ~ HLA_DRB1_15_01 + sex + PC1 + PC2, data = dat )
mvlm.res.1$pseudo.rsq["Omnibus Effect",1] 

summary(mvlm.res.1)

################################
                Statistic Numer.DF Pseudo.R.Square p.value    
Omnibus Effect    0.90858        4       4.013e-03 0.52003    
(Intercept)     315.81910        1                 < 1e-20 ***
HLA_DRB1_15_01    0.82089        1       9.064e-04 0.45784    
sex               0.07375        1       8.144e-05 0.96786    
PC1               1.71025        1       1.888e-03 0.16841    
PC2               0.97278        1       1.074e-03 0.38646    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
################################


# Explained variance by HLA_DRB1_15_01 will be..
mvlm.res.1$pseudo.rsq["Omnibus Effect",1]  - mvlm.res.0$pseudo.rsq["Omnibus Effect",1] 

```

