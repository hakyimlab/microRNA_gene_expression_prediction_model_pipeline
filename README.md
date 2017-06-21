# microRNA Gene Expression Prediction Model Pipeline

## Introduction 
+ A pipeline is used to build microRNA gene expression prediction models  

## Requirements 
*[Python2.7](http://www.python.org/download/)
*[R 3.0+](http://www.r-project.org/)
*[data.table](https://github.com/Rdatatable/data.table)
*[glmnet](https://cran.r-project.org/web/packages/glmnet/index.html)
*[qvalue](http://bioconductor.org/packages/release/bioc/html/qvalue.html)
*[dplyr](https://github.com/tidyverse/dplyr)
*[bit64](https://cran.r-project.org/web/packages/bit64/index.html)
*[doMC](https://cran.r-project.org/web/packages/doMC/index.html)

## Installing Pipeline 
```bash
git clone https://github.com/jiamaozheng/microRNA_gene_expression_prediction_model_pipeline
```   

## Supported Input File Formats 
+ [expression file - miRNA_CTR_Exp.RDS](https://s3.amazonaws.com/imlab-jiamaoz/shared/gene_expression_sample.txt)
+ [expression annotation file - miRBase_miRNA_gene_annotation.RDS](https://s3.amazonaws.com/imlab-jiamaoz/shared/gene_annotation_sample.txt)
+ [genotype file - geuvadis_snps.txt ](https://s3.amazonaws.com/imlab-jiamaoz/shared/genotype_sample.txt)
+ [SNP annotation file - gEUVADIS.SNP.annotation.RDS](https://s3.amazonaws.com/imlab-jiamaoz/shared/snp_annotation_sample.txt)


## Command Line Parameters 
  Argument              |  Abbre  | Required | Default  | Description  
  ----------------------| ------- | -------- | -------- | ------------------------
  --project_name	    |  -p     |   Yes    |  None    | Project name (e.g. gEUVADIS_LCL, TCGA or Framingham)
  --molecular_type      |  -m     |   Yes    |  None    | Molecular types (e.g. mRNA, miRNA or shRNA)
  --alpha      	        |  -a     |   Yes    |  None    | Alpha values (e.g. 0.05, 0.5, or 1)
  --snpset    	        |  -s     |   Yes    |  None    | SNP set used for analysis (e.g. 1KG_snps, HapMap_snps) 
  --n_k_folds  	        |  -n     |   No     |  10      | The number of folds for cross-validation (e.g. 10) 
  --fdr_level 	        |  -f     |   No     |  0.05    | FDR used for filtering modelS (e.g. 0.05) 
  --window    	        |  -w     |   No     |  1e6     | The number of bps to +/- transcription start site(TSS)
  --expression_path     |  -e     |   No     |  ''      | User-defined expression file path 
  --genotype_path	    |  -g     |   No     |  ''      | User-defined genotype file path 
  --gene_annot_path	    |  -x     |   No     |  ''      | User-defined gene annotation file path 
  --snp_annot_path	    |  -y     |   No     |  ''      | User-defined snp annotation file path 
  --intermediate_path   |  -i     |   No     |  ''      | User-defined intermediate file path 
  --results_output_path |  -r     |   No     |  ''      | User-defined output file path 


## Running Pipeline 
+ **Example 1** (*basic command, recommended*): 
Navigate to the folder that contains downloaded pipeline `microRNA_gene_expression_prediction_model_pipeline`, create a default input directory by using a [script](https://github.com/jiamaozheng/prediction_model_docker/blob/master/input_directory.sh), navigate to the pipeline folder `microRNA_gene_expression_prediction_model_pipeline`, and then execute the following command with four requried parameters.  
```bash
	./run.py \
	--project_name gEUVADIS_LCL_miRNA \
	--molecular_type miRNA \
	--alpha 0.5 \
	--snpset 1KG_snps \
```
Alternatively, you can fun the following shortcut
```bash
	./run.py -p gEUVADIS_LCL_miRNA -m miRNA -a 0.5 -s 1KG_snps 
```

+ **Example 2** (*basic command + usered-defined input file paths*): 
Navigate to the pipeline folder `microRNA_gene_expression_prediction_model_pipeline`, and then execute the following command with four requried parameters and five user-defined input file option parameters 
```bash
	./run.py \
	--project_name gEUVADIS_LCL_miRNA \
	--molecular_type miRNA \
	--alpha 0.5 \
	--snpset 1KG_snps \

	--expression_path ../input_mirna/expression_phenotypes/miRNA_CTR_Exp.RDS \
	--genotype_path ../input_mirna/genotype/ \
	--gene_annot_path ../input_mirna/gene_annotation/miRBase_miRNA_gene_annotation.RDS \
	--snp_annot_path ../input_mirna/snp_annotation/gEUVADIS.SNP.annotation.RDS \
	--intermediate_path ../input_mirna/intermediate/
```

+ **Example 3** (*basic command + user-defined input file paths + user-provided output file path*): 
Navigate to the pipeline folder `microRNA_gene_expression_prediction_model_pipeline`, and then execute the following command with four requried parameters, five user-defined input file option parameters, and one user-defined output file option parameter
```bash
	./run.py \
	--project_name gEUVADIS_LCL_miRNA \
	--molecular_type miRNA \
	--alpha 0.5 \
	--snpset 1KG_snps \

	--expression_path ../input_mirna/expression_phenotypes/miRNA_CTR_Exp.RDS \
	--genotype_path ../input_mirna/genotype/ \
	--gene_annot_path ../input_mirna/gene_annotation/miRBase_miRNA_gene_annotation.RDS \
	--snp_annot_path ../input_mirna/snp_annotation/gEUVADIS.SNP.annotation.RDS \
	--intermediate_path ../input_mirna/intermediate/ \ 
	--results_output_path ../output_mirna/
```

+ **Example 4** (*basic command + user-defined input file paths + user-defined output file path + default parameters*): 
Navigate to the pipeline folder `microRNA_gene_expression_prediction_model_pipeline`, and then execute the following command with four requried parameters, five user-defined input file option parameters, one user-defined output file option parameter, and three default parameters which can be modified by users if necessary 
```bash
	./run.py \
	--project_name gEUVADIS_LCL_miRNA \
	--molecular_type miRNA \
	--alpha 0.5 \
	--snpset 1KG_snps \

	--expression_path ../input_mirna/expression_phenotypes/miRNA_CTR_Exp.RDS \
	--genotype_path ../input_mirna/genotype/ \
	--gene_annot_path ../input_mirna/gene_annotation/miRBase_miRNA_gene_annotation.RDS \
	--snp_annot_path ../input_mirna/snp_annotation/gEUVADIS.SNP.annotation.RDS 
	--intermediate_path ../input_mirna/intermediate/

	--n_k_folds 10 
	--fdr_level 0.05 
	--window 1e6 
```






