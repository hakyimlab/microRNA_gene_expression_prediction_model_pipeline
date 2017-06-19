# PredictDB Pipeline

## Introduction 
+ A pipeline used to build PredictDB  

## Prerequisites
1. [Python2.7](http://www.python.org/download/)
2. [R 3.0+](http://www.r-project.org/)
3. [tidyverse](http://tidyverse.org)

## Installation and Setup 

## Run 
+ Open the terminal, and execute:
```
	./run.py \
	--project_name gEUVADIS_LCL_miRNA \
	--molecular_type miRNA \
	--alpha 0.5 \
	--snpset 1KG_snps \
	--expression_path ../input_mirna/expression_phenotypes/miRNA_CTR_Exp.RDS \
	--genotype_path ../input_mirna/genotype/ \
	--gene_annot_path ../input_mirna/gene_annotation/miRBase_miRNA_gene_annotation.RDS \
	--snp_annot_path ../input_mirna/snp_annotation/gEUVADIS.SNP.annotation.RDS 
```

## Output 

