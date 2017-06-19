# Author = 'Jiamao Zheng <jiamaoz@yahoo.com>'

# import libraries and functions 
source('util/imlab_utilities.R')
source('util/generate_models_outputs.R')

# args 
project <- argv[1]
expression_RDS <- argv[2]
geno_file <- argv[3]
gene_annot_RDS <- argv[4]
snp_annot_RDS <- argv[5]
n_k_folds <- as.numeric(argv[6])
alpha <- as.numeric(argv[7])
intermediate_out_dir <- argv[8]
snpset <- argv[9]
window <- as.numeric(argv[10]) 
gene_type <- argv[11]
output_data_dir <- argv[12]
fdr_level <- as.numeric(argv[13]) 
genotype_input_data_dir <- argv[14]
numberof_genotype_files <- ' '

#---------------------------------------------------------------------------
# Step one: to create individual beta, weights, logs and covariances files 
#---------------------------------------------------------------------------
# input: gene expression 
expression = data.frame() 

if (substr(expression_RDS, nchar(expression_RDS)-2, nchar(expression_RDS)) == 'txt'){
	expression = readTextFile(expression_RDS)
} else if (substr(expression_RDS, nchar(expression_RDS)-2, nchar(expression_RDS)) == 'RDS'){
	expression = readRDSFile(expression_RDS)
}
# class(expression) <- 'numeric'
# class(expression) <- 'integer'
n_samples <- nrow(expression) # for sample_info table in the final sqlite db 

# input: gene annotation 
gene_annot = data.frame()

if (substr(gene_annot_RDS, nchar(gene_annot_RDS)-2, nchar(gene_annot_RDS)) == 'txt'){
	gene_annot = readTextFile(gene_annot_RDS)
} else if (substr(gene_annot_RDS, nchar(gene_annot_RDS)-2, nchar(gene_annot_RDS)) == 'RDS'){
	gene_annot = readRDSFile(gene_annot_RDS)
}


# input: snp annotation 
snp_annot = data.frame()

if (substr(snp_annot_RDS, nchar(snp_annot_RDS)-2, nchar(snp_annot_RDS)) == 'txt'){
	snp_annot = readTextFile(snp_annot_RDS)
} else if (substr(snp_annot_RDS, nchar(snp_annot_RDS)-2, nchar(snp_annot_RDS)) == 'RDS'){
	snp_annot = readRDSFile(snp_annot_RDS)
}

# input: genotype 
genotype = data.frame()

if (length(list.files(genotype_input_data_dir)) == 1){
    numberof_genotype_files = 1
} else if (length(list.files(genotype_input_data_dir)) == 22){
    numberof_genotype_files = 22
} else {
    print ('Waring: Only 1 or 22 genotype input files are allowed in the input directory of ../input/genotype/!')
}


if (numberof_genotype_files == 1) {
	genotype = readTextFile(geno_file)
	rownames(genotype) <- genotype$Id 
	genotype <- genotype[, -1]
	genotype <- t(genotype)
}


# create models and covariances 
for (i in 1:22){
	chrom <- i 
	if (numberof_genotype_files == 22) {
		# input: 22 genotype 
		genotype = readTextFile(paste(geno_file, 'chr', chrom, '.txt', sep=''))
		rownames(genotype) <- genotype$Id 
		genotype <- genotype[, -1]
		# genotype <- t(genotype)   
	}

	generate_models(expression, genotype, gene_annot, snp_annot, n_k_folds, alpha, intermediate_out_dir, project, chrom, snpset, window)
}

#--------------------------------------------------------------------------------------
# Step two: to combine all files and generate sqlite db and covariance files   
#--------------------------------------------------------------------------------------
generate_outputs(project, gene_annot, n_samples, n_k_folds, alpha, snpset, fdr_level, gene_type, intermediate_out_dir, output_data_dir)

