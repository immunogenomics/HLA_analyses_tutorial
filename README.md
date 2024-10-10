# HLA analyses tutorial
A thorough tutorial on HLA imputation and association, accompanying our manuscript "Tutorial: A statistical genetics guide to identifying HLA alleles driving complex disease"

<div align="center">
<img src="https://raw.githubusercontent.com/immunogenomics/HLA_analyses_tutorial/main/images/for_web_Overview_v4.png" width=85%>
</div>

The tutorial consists of two parts:
- #1. [HLA imputation](https://github.com/immunogenomics/HLA_analyses_tutorial/blob/main/tutorial_HLAQCImputation.ipynb) : We introduce protocols to QC genotype, perform haplotype phasing and HLA imputation. We provide useful scripts and example usage, with example genotype and reference datasets.

- #2. [HLA association and fine-mapping](https://github.com/immunogenomics/HLA_analyses_tutorial/blob/main/tutorial_association.md): We introduce various statistical methods to identify and fine-map disease-associated HLA variations. The HLA imputation results from section #1 will be used. We provide useful scripts with some example phenotype data.

  

### Note for Michigan Imputation Server users

For the amino acid (AA) residue imputation, sometimes the imputed residue is not clear from the post imputation VCF. Please use [this `bim` file](https://github.com/immunogenomics/HLA_analyses_tutorial/blob/main/data/AA_annotated.bim) that summarizes all information about AA residue names (in the 2nd column) that appear on the imputed VCF and the binary allele (residue) difinitions (in the 5th and 6th columns). 

For further naming schemes, please also refer to [this website](https://software.broadinstitute.org/mpg/snp2hla/makereference_manual.html) at the sections of 'Marker Nomenclature'.  

