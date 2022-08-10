#!/usr/bin/env python
import sys
import argparse
import os

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--aaraw', '-a', default=None, type=str,
                    help='Name of the amino acid dosage *.raw file',
                    required=True)
parser.add_argument('--out', '-o', default=None, type=str,
                    help='Prefix of the output file',
                    required=True)
parser.add_argument('--allele', '-l', default=None, type=str,
                    help='Name of the amino acid allele name file',
                    required=True)
parser.add_argument('--pheno', '-p', default=None, type=str,
                    help='Name of the phenotype file',
                    required=True)
parser.add_argument('--cov', '-c', default=None, type=str,
                    help='Name of the covariate file',
                    required=True)
parser.add_argument('--phenoname', '-n', default=None, type=str,
                    help='Name of the phenotype (2nd column of the phenotype file)',
                    required=True)

parser.add_argument('--covname', '-m', default=None, nargs='+',
                    help='Name of the covariates',
                    required=True)

args = parser.parse_args()

cmd = 'echo "library(data.table)" > ' + args.out + "_runcode.R"
os.system(cmd)
cmd = 'echo \'dose<-fread("' + args.aaraw + '", header=T)\' >> ' + args.out + "_runcode.R"
os.system(cmd)
cmd = 'echo \'dose<-dose[,c(-1,-3:-6)]\' >> ' + args.out + "_runcode.R"
os.system(cmd)
cmd = 'echo \'allelenames<-as.character(read.table("' + args.allele + '")$V1)\' >> ' + args.out + "_runcode.R"
os.system(cmd)
cmd = 'echo \'colnames(dose)<-c("IID",allelenames)\' >> ' + args.out + "_runcode.R"
os.system(cmd)
cmd = 'echo \'pheno <- read.table("' + args.pheno + '", header=T)\' >> ' + args.out + "_runcode.R"
os.system(cmd)
cmd = 'echo \'if(max(pheno[,2])==2){pheno[,2]<-pheno[,2] - 1}\' >> ' + args.out + "_runcode.R"
os.system(cmd)
cmd = 'echo \'cov <- read.table("' + args.cov + '", header=T)\' >> ' + args.out + "_runcode.R"
os.system(cmd)
cmd = 'echo \'d <- merge(dose, pheno, by = "IID")\' >> ' + args.out + "_runcode.R"
os.system(cmd)
cmd = 'echo \'d <- merge(d, cov, by = "IID")\' >> ' + args.out + "_runcode.R"
os.system(cmd)


OUT = open(args.out + "_runcode.R", "a")
OUT2 = open(args.out + "_alleles.txt", "w")

print('pval_list<-NULL', file = OUT)
print('deviance_list<-NULL', file = OUT)

AA_allele_dic = {}

with open(args.allele, "r") as f:
	for line in f:
		line = line.rstrip()
		if line.count("_") == 5:
			AA = "_".join(line.split("_")[:5])
			this_allele = line
			AA_allele_dic.setdefault(AA,[]).append(this_allele)
		else:
			AA = line
			this_allele = line
			AA_allele_dic.setdefault(AA,[]).append(this_allele)


for AA in AA_allele_dic:
	num_allele = len(AA_allele_dic[AA])
	print('obj1<-glm(' + args.phenoname + '~' + '+'.join(args.covname) + ', data=d,family=binomial)' ,file = OUT)
	alleles = '+'.join(AA_allele_dic[AA])
	out = 'obj2<-glm(' + args.phenoname + '~' + alleles + '+' + '+'.join(args.covname) + ',data=d,family=binomial(link="logit"))'
	print(out, file = OUT)
	print('Chisqtest <- anova(obj1, obj2, test="Chisq")', file = OUT)
	print('pval <- Chisqtest$`Pr(>Chi)`[2]', file = OUT)
	print('deviance <- Chisqtest$Deviance[2]', file = OUT)
	print('pval_list <- c(pval_list,pval)', file = OUT)
	print('deviance_list <- c(deviance_list,deviance)', file = OUT)
	print(AA, file = OUT2)

OUT.close()
OUT2.close()

cmd = 'echo \'alleles <- as.character(read.table("' + args.out + '_alleles.txt")$V1)\' >> ' + args.out + "_runcode.R"
os.system(cmd)
cmd = 'echo \'summary <- data.frame(ALLELE_NAME = alleles, OMNIBUS_DEVIANCE = deviance_list, OMNIBUS_PVALUE = pval_list)\' >> ' + args.out + "_runcode.R"
os.system(cmd)
cmd = 'echo \'outfile <- "' + args.out + '_omnibus_result.txt"\' >> ' + args.out + "_runcode.R"
os.system(cmd)
cmd = 'echo \'write.table(summary, outfile, sep="\t", quote=F, row.names=F)\' >> ' + args.out + "_runcode.R"
os.system(cmd)

cmd = 'Rscript ' + args.out + "_runcode.R > " + args.out + ".log"
os.system(cmd)


