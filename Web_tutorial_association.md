# Tutorial: A statistical genetics guide to identifying HLA alleles driving complex disease





[TOC]

Author: Saori Sakaue (ssakaue@broadinstitute.org)

Lastly updated XX/XX/2022



## Reference panels

Under section "*HLA imputation reference panel*"

Make a table and a link to the ref panels







## Phasing and imputation

Under section "*Tools for genotype phasing and HLA imputation*"

- SNP2HLA bash script (original)
- `SNP2HLA.py` module with adding options for users to specify QC thresholds
- SHAPEIT or EAGLE phasing + minimac3 imputation
- MIS imputation (just show how to prepare input genotype file in VCF)







## HLA association and fine-mapping

Under section "*HLA association and fine-mapping*"



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



- Single-marker test

Convert imputed genotype to dosage txt file by plink2



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



You can also do this by using custom R script

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



- Omnibus test

First, we extract amino acid polymorphisms from the dosage output and creat `*.raw` file by using `PLINK2`.

We also apply QC to extract any amino acid polymorphisms with *Rsq* > 0.7.



```bash
sed 1d data_assoc/converted.info | grep "^AA_" | awk '{print $1,"T"}' > QCed_AA_variants.txt
sed 1d data_assoc/converted.info | grep "^AA_" | awk '{print $1,"T"}' > AA_variants_alleles.txt

plink2 --pfile ${imputed} \
  --extract AA_variants.txt \
  --export-allele AA_variants_alleles.txt \
  --export A \
  --out imp/validation/geno/${id}.v3.AA \



```



- Conditional haplotype test

Notation about QC out individuals



- Non-additive association test
- Interaction test
- Multi-trait test