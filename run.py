#!/usr/bin/env python

# miRNA usage --------------------------
# ./run.py \
# --project_name gEUVADIS_LCL_miRNA \
# --molecular_type miRNA \
# --alpha 0.5 \
# --snpset 1KG_snps \
# --expression_path ../input_mirna/expression_phenotypes/miRNA_CTR_Exp.RDS \
# --genotype_path ../input_mirna/genotype/ \
# --gene_annot_path ../input_mirna/gene_annotation/miRBase_miRNA_gene_annotation.RDS \
# --snp_annot_path ../input_mirna/snp_annotation/gEUVADIS.SNP.annotation.RDS

# proteomics usage -------------------------
# ./run.py \
# --project_name framingham_proteomics \
# --molecular_type protein_coding \
# --alpha 0.5 \
# --snpset 1KG_snps \
# --expression_path ../input/expression_phenotypes/framingham_immunoassay_final_4751_66.txt \
# --genotype_path ../input/genotype/ \
# --gene_annot_path ../input/gene_annotation/gencode.v26.annotation.parsed.protein_coding.RDS \
# --snp_annot_path ../input/snp_annotation/GTEx.1KB.SNP.annotation.RDS

# import modules and functions 
from util.imlab_utilities import *

__author__ = "Jiamao Zheng <jiamaoz@yahoo.com>"
__version__ = "Revision: 0.01"
__date__ = "Date: 2017-06-08"

class PredictionModel(object):
    """docstring for create_model"""
    def __init__(self):
        # project info 
        self.project = '' 
        self.expression = ''
        self.geno = ''
        self.gene_annot = ''
        self.snp_annot = ''
        self.project_time = ''
        self.hour = ''
        self.minute = '' 
        self.second = ''
        self.project_id = ''

        # input info 
        self.intermediate_out_dir = ''
        self.output_data_dir = ''
        self.logs_data_dir = ''
        self.genotype_input_data_dir = ''
        self.expression_input_data_dir = ''
        self.gene_annotation_input_data_dir = ''
        self.snp_annotation_input_data_dir = ''

        # commond line arguments 
        self.project_name = ''
        self.gene_type = '' 
        self.alpha = ''
        self.snpset = ''
        self.n_k_folds = ''
        self.fdr_level = ''
        self.window = '' 
        # self.numberof_genotype_files = ''

        # optional input file paths 
        self.expression_path = ''
        self.genotype_path = ''
        self.gene_annot_path = ''
        self.snp_annot_path = ''
        self.intermediate_path = ''

        # output 
        self.output_path = ''

        # logs  
        self.msg = ''

    def get_args(self):
        # setup commond line arguments 
        parser = argparse.ArgumentParser()

        # required arguments
        parser.add_argument('-p', '--project_name', required=True, default=None, type=str, help='e.g gEUVADIS, TCGA, or framingham')
        parser.add_argument('-m', '--molecular_type', required=True, default=None, type=str, help='e.g. mRNA, miRNA, or shRNA')
        parser.add_argument('-a', '--alpha', required=True, default=None, type=str, help='alpha values, e.g. 0, 0.05, 0.5, or 1')
        parser.add_argument('-s', '--snpset', required=True, default=None, type=str, help='SNP set used for analysis (e.g. 1KG_snps, HapMap_snps')
        # parser.add_argument('-n', '--numberof_genotype_files', required=True, default=None, type=str, help='The number of genotype files, e.g. 1, 22')

        # not required but has default values. Please include them if you would like to modify the default settings 
        parser.add_argument('-c', '--n_k_folds', required=False, default='10', type=str, help='the number of folds for cross-validation, e.g. 10')
        parser.add_argument('-f', '--fdr_level', required=False, default='0.05', type=str, help='fdr used for model filtering, e.g. 0.05')
        parser.add_argument('-w', '--window', required=False, default='1e6', type=str, help='the number of bps to +/- TSS')

        # not required, and those are optional arguments. Please include them if you would like to provide one or more your own input file paths. If not, please include your input files in script-defined input paths as described in README. 
        parser.add_argument('-e', '--expression_path', required=False, default='', type=str, help='the user-provided path for expression input file')
        parser.add_argument('-g', '--genotype_path', required=False, default='', type=str, help='the user-provided path for genotype input file, use the directory that holds either 1 genotype file or 22 genotype files')
        parser.add_argument('-x', '--gene_annot_path', required=False, default='', type=str, help='the user-provided path for gene annotation input file')
        parser.add_argument('-y', '--snp_annot_path', required=False, default='', type=str, help='the user-provided path for snp annotation input file')
        parser.add_argument('-i', '--intermediate_path', required=False, default='', type=str, help='the user-provided path for holding the intermediate files')

        # not requried, and those are optional arguments. Please include them if you would like to output results to your own path 
        parser.add_argument('-r', '--results_output_path', required=False, default='', type=str, help='the user-provided path for outputing result files')

        # parse the arguments 
        args = parser.parse_args()

        # Change the parameters below to reflect your project requirements
        self.project_name = args.project_name 
        self.gene_type = args.molecular_type 
        self.alpha = args.alpha
        self.snpset = args.snpset
        self.n_k_folds = args.n_k_folds
        self.fdr_level = args.fdr_level
        self.window = args.window 
        # self.numberof_genotype_files = args.numberof_genotype_files
        
        self.expression_path = args.expression_path
        self.genotype_path = args.genotype_path
        self.gene_annot_path = args.gene_annot_path
        self.snp_annot_path = args.snp_annot_path
        self.intermediate_path = args.intermediate_path
        self.output_path = args.results_output_path

    def get_parameters(self):
        # input, output, log, intermediate file paths 
        self.intermediate_out_dir = '../input/intermediate/' # for intermediate data files (individual chromsome files)
        self.expression_input_data_dir = '../input/expression_phenotypes/'
        self.genotype_input_data_dir = '../input/genotype/'
        self.gene_annotation_input_data_dir = '../input/gene_annotation/'
        self.snp_annotation_input_data_dir = '../input/snp_annotation/'
        self.output_data_dir = '../output/'
        self.logs_data_dir = '../logs/'

        #### input file paths 
        ## two options: use user-provided file paths or script default input file paths 
        # gene expression data
        self.expression = get_input_path(self.expression_path, self.expression_input_data_dir, self.expression)

        # genotype data (0-1-2 format)
        self.geno = get_genotype_input_path(self.genotype_path, self.genotype_input_data_dir, self.geno)

        # gene annotation data
        self.gene_annot = get_input_path(self.gene_annot_path, self.gene_annotation_input_data_dir, self.gene_annot)

        # snp annotation data
        self.snp_annot = get_input_path_RDS (self.snp_annot_path, self.snp_annotation_input_data_dir, self.snp_annot)

        # project INFO 
        # self.project_time = datetime.now().strftime('%Y-%m-%d') 
        self.project_time = '2017-06-19'
        self.hour = datetime.now().strftime('%H') 
        self.minute = datetime.now().strftime('%M') 
        self.second = datetime.now().strftime('%S') 
        # self.project_id = str(myuuid.uuid4())
        self.project_id = '5e556ffe-2315-40e8-9d88-11f2e61715cb'
        self.project = self.project_name + '_' + self.gene_type + '_alpha-' + self.alpha + '_snpset-' + self.snpset + '_window-' + self.window + '_cv-folds-' + self.n_k_folds +'_fdr-' + self.fdr_level + '_' + self.project_id + '_' + self.project_time

    def start_logs(self):
        # open logfile and log project info 
        open_log(self.project)

        # project info 
        self.msg = '\n\n' + "---------- PROJECT INFORMATION ----------- " + '\n\n' \
        + "PROJECT NAME: " + self.project_name + '\n\n' \
        + "GENE TYPE: " + self.gene_type + '\n\n' \
        + "ALPHA: " + self.alpha + '\n\n' \
        + "SNPSET: " + self.snpset + '\n\n' \
        + "FDR_LEVEL: " + self.fdr_level + '\n\n' \
        + "WINDOW: " + self.window + '\n\n' \
        + "N_K_FOLDS: " + self.n_k_folds + '\n\n' \
        + "GENE EXPRESSION FILE: " + self.expression + '\n\n' \
        + "GENOTYPE FILE: " + self.geno + '\n\n' \
        + "GENE ANNOTATION FILE: " + self.gene_annot + '\n\n' \
        + "SNP ANNOTATION FILE: " + self.snp_annot + '\n\n' \
        + "PROJECT TIME: " + self.project_time + '_' + self.hour + ':' + self.minute + ':' + self.second + '\n\n' \
        + "PROJECT ID: " + self.project_id + '\n\n' \
        + "PROJECT UNIQUE OUTPUT PREFIX: " + self.project + '\n\n'
        print(self.msg)
        add_log(self.msg)

        # make models 
        self.msg = '\n\n' + "---------- GENERATING MODELS ----------- " + '\n\n' \
        + "MAKING WEIGHTS, BETA, COVARIANCES AND LOGS ..." + '\n\n' 
        print(self.msg)
        add_log(self.msg)

    def end_logs(self):
        # log output information 
        self.msg = '\n\n' + "---------- OUTPUT FILES ----------- " + '\n\n' \
        + "OUTPUT FILES CAN BE FOUDN UNDER THE DIRECTORY OF : " + self.output_data_dir + '\n\n' \
        + "SQLITE DB: " + glob.glob(self.output_data_dir + self.project + '_sqlite.db')[0] + '\n\n' \
        + "COVARIANCES: " + glob.glob(self.output_data_dir + self.project + '_covariances.txt')[0] + '\n\n' \
        + "Done!"
        print(self.msg)
        add_log(self.msg)

        # write and finish logs 
        write_logs()
        finish_log() 
        
    def run(self):
        # start a log file 
        self.start_logs()

        # create a folder for intermediate results
        if len(self.genotype_path) != 0:                               # use user provided input file paths 
            self.intermediate_out_dir = self.intermediate_path + self.project + '/'  
        else: 
            self.intermediate_out_dir = self.intermediate_out_dir + self.project + '/'
        os.system("mkdir " + self.intermediate_out_dir)

        # create a folder for outputs  
        if len(self.output_path) != 0:                               # use user provided input file paths 
            self.output_data_dir = self.input_path + self.project + '/'  
        else: 
            self.output_data_dir = self.output_data_dir + self.project + '/'
        os.system("mkdir " + self.output_data_dir)

        # run prediction models  
        # cmd = 'Rscript %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' %('model.R', self.project, self.expression, self.geno, self.gene_annot, self.snp_annot, self.n_k_folds, self.alpha, self.intermediate_out_dir, self.snpset, self.window, self.gene_type, self.output_data_dir, self.fdr_level, self.genotype_input_data_dir, self.numberof_genotype_files)
        # os.system(cmd)
        subprocess.call(['Rscript', 'model.R', self.project, self.expression, self.geno, self.gene_annot, self.snp_annot, self.n_k_folds, self.alpha, self.intermediate_out_dir, self.snpset, self.window, self.gene_type, self.output_data_dir, self.fdr_level, self.genotype_path])

        # close a log file 
        self.end_logs()


# --------------------------------------
# main functions 
# --------------------------------------
def main():
    #------------------------
    # Instantiate class
    #------------------------
    prediction_models = PredictionModel() 

    #------------------------
    # get arguments 
    #------------------------
    prediction_models.get_args() 

    #-----------------------------
    # setup parameters 
    #-----------------------------
    prediction_models.get_parameters()

    #-----------------------------
    # generate models 
    #-----------------------------
    prediction_models.run()



# initialize the script
if __name__ == '__main__':
    sys.exit(main())

