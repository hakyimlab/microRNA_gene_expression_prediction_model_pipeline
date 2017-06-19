# Author = 'Jiamao Zheng <jiamaoz@yahoo.com>'

# generate model
generate_models <- function(expression, genotype, gene_annot, snp_annot,
    n_k_folds, alpha, out_dir, tis, chrom, snpset, window) {
  # filter gene_annot by chr 
  gene_annot <- subset(gene_annot, gene_annot$chr == chrom)

  # filter snp_annot by chr 
  snp_annot <- subset(snp_annot, snp_annot$chr == chrom)

  rownames(gene_annot) <- gene_annot$gene_id

  # Subset expression data to only include genes with gene_info
  expression <- expression[, intersect(colnames(expression), rownames(gene_annot))]

  exp_samples <- rownames(expression)
  exp_genes <- colnames(expression)
  n_samples <- length(exp_samples)
  n_genes <- length(exp_genes)
  seed <- sample(1:2016, 1)
  log_df <- data.frame(chrom, n_genes, seed, alpha)
  colnames(log_df) <- c('chr', 'n_genes', 'seed_for_cv', 'alpha')
  write.table(log_df, file = out_dir %&% tis %&% '_chr' %&% chrom %&% '_elasticNet_model_log.txt', quote = FALSE, row.names = FALSE, sep = "\t")
  set.seed(seed)
  groupid <- sample(1:n_k_folds, length(exp_samples), replace = TRUE)

  resultsarray <- array(0, c(length(exp_genes), 9))
  dimnames(resultsarray)[[1]] <- exp_genes
  resultscol <- c("gene", "alpha", "cvm", "lambda.iteration", "lambda.min", "n.snps", "R2", "pval", "genename")
  dimnames(resultsarray)[[2]] <- resultscol
  workingbest <- out_dir %&% "working_TW_" %&% tis %&% "_exp_" %&% n_k_folds %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_chr" %&% chrom %&% ".txt"
  write(resultscol, file = workingbest, ncolumns = 9, sep = "\t")

  weightcol <- c("gene","rsid","ref","alt","beta","alpha")
  workingweight <- out_dir %&% "TW_" %&% tis %&% "_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_weights_chr" %&% chrom %&% ".txt"
  write(weightcol, file = workingweight, ncol = 6, sep = "\t")

  covariance_out <- out_dir %&% tis %&% '_chr' %&% chrom %&% '_snpset_' %&% snpset %&% '_alpha_' %&% alpha %&% "_covariances.txt"
  covariance_col <- c('gene','rsid1', 'rsid2', 'covariance')
  write(covariance_col, file = covariance_out, ncolumns = 4, sep = ' ')

  if (length(exp_genes) > 0){
  for (i in 1:length(exp_genes)) {

    cat("Chr", chrom, ": ", i, "/", length(exp_genes), "\n")
    gene <- exp_genes[i]
    print(gene)
    # Reduce genotype data to only include SNPs within specified window of gene.
    geneinfo <- gene_annot[gene,]
    start <- geneinfo$start - window
    end <- geneinfo$end + window
    # Pull cis-SNP info
    cissnps <- subset(snp_annot, snp_annot$pos >= start & snp_annot$pos <= end)

    # Pull cis-SNP genotypes
    cisgenos <- genotype[,intersect(colnames(genotype), cissnps$varID), drop = FALSE]
    # dim(cisgenos)
    # Reduce cisgenos to only include SNPs with at least 1 minor allele in dataset
    cm <- colMeans(cisgenos, na.rm = TRUE)
    minorsnps <- subset(colMeans(cisgenos), cm > 0 & cm < 2)
    minorsnps <- names(minorsnps)
    cisgenos <- cisgenos[,minorsnps, drop = FALSE]
    if (ncol(cisgenos) < 2) {
      # Need 2 or more cis-snps for glmnet.
      bestbetas <- data.frame()
    } else {
      # Pull expression data for gene
      exppheno <- expression[,gene]
      # Scale for fastLmPure to work properly
      exppheno <- scale(exppheno, center = TRUE, scale = TRUE) 
      exppheno[is.na(exppheno)] <- 0
      rownames(exppheno) <- rownames(expression)      
      # Run Cross-Validation
      # parallel = TRUE is slower on tarbell
      bestbetas <- tryCatch(
        { fit <- cv.glmnet(as.matrix(cisgenos),
                          as.vector(exppheno),
                          nfolds = n_k_folds,
                          alpha = alpha,
                          keep = TRUE,
                          foldid = groupid,
                          parallel = TRUE)
          # Pull info from fit to find the best lambda   
          fit.df <- data.frame(fit$cvm, fit$lambda, 1:length(fit$cvm))
          # Needs to be min or max depending on cv measure (MSE min, AUC max, ...)
          best.lam <- fit.df[which.min(fit.df[,1]),]
          cvm.best <- best.lam[,1]
          lambda.best <- best.lam[,2]
          # Position of best lambda in cv.glmnet output
          nrow.best <- best.lam[,3]
          # Get the betas from the best lambda value
          ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best])
          ret[ret == 0.0] <- NA
          # Pull the non-zero betas from model
          as.vector(ret[which(!is.na(ret)),])
        },
        error = function(cond) {
          # Should fire only when all predictors have 0 variance.
          message('Error with gene ' %&% gene %&% ', index ' %&% i)
          message(geterrmessage())
          return(data.frame())
        }
      )
    } 
    if (length(bestbetas) > 0) {
      names(bestbetas) <- rownames(ret)[which(!is.na(ret))]
      # Pull out the predictions at the best lambda value.    
      pred.mat <- fit$fit.preval[,nrow.best]
      res <- summary(lm(exppheno~pred.mat))
      genename <- as.character(gene_annot[gene, 3])
      rsq <- res$r.squared
      pval <- res$coef[2,4]  
      resultsarray[gene,] <- c(gene, alpha, cvm.best, nrow.best, lambda.best, length(bestbetas), rsq, pval, genename)  
      # Output best shrunken betas for PrediXcan
      bestbetalist <- names(bestbetas)
      bestbetainfo <- snp_annot[bestbetalist,]
      betatable <- as.matrix(cbind(bestbetainfo,bestbetas))
      write_covariance(gene, cisgenos, betatable[,"rsid"], betatable[,"varID"], covariance_out)
      # Output "gene", "rsid", "refAllele", "effectAllele", "beta"
      # For future: To change rsid to the chr_pos_ref_alt_build label, change "rsid" below to "varID".
      betafile <- cbind(gene,betatable[,"rsid"],betatable[,"refAllele"],betatable[,"effectAllele"],betatable[,"bestbetas"], alpha)
      # Transposing betafile necessary for correct output from write() function
      write(t(betafile), file = workingweight, ncolumns = 6, append = TRUE, sep = "\t")
      write(resultsarray[gene,], file = workingbest, ncolumns = 9, append = TRUE, sep = "\t")
    } else {
      genename <- as.character(gene_annot[gene,3])
      resultsarray[gene,1] <- gene
      resultsarray[gene,2:8] <- c(alpha,NA,NA,NA,0,NA,NA)
      resultsarray[gene,9] <- genename
    }
  }
  } else {
      cat("Chr", chrom, ": ", i, "/", length(exp_genes), "\n")
  }
  write.table(resultsarray,file=out_dir %&% "TW_" %&% tis %&% "_chr" %&% chrom %&% "_exp_" %&% n_k_folds %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% ".txt",quote=F,row.names=F,sep="\t")
}

# calculate covariance 
write_covariance <- function(gene, cisgenos, model_rsids, model_varIDs, covariance_out) {
  model_geno <- cisgenos[,model_varIDs, drop=FALSE]
  geno_cov <- cov(model_geno)
  cov_df <- data.frame(gene=character(),rsid1=character(),rsid2=character(), covariance=double())
  for (i in 1:length(model_rsids)) {
    for (j in i:length(model_rsids)) {
      cov_df <- rbind(cov_df, data.frame(gene=gene,rsid1=model_rsids[i], rsid2=model_rsids[j], covariance=geno_cov[i,j]))
    }
  }
  write.table(cov_df, file = covariance_out, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
}

# generate outputs
generate_outputs <- function(project, gene_annot, n_samples, n_k_folds, alpha, snpset, fdr_level, gene_type, intermediate_out_dir, output_data_dir) {
  # filter what target you are interested
  gene_annot <- gene_annot %>% filter(gene_type == gene_type)

  # create empty data frame for each set of files (logs, beta, results, covariances) 
  log_files = data.frame()
  beta_files = data.frame()
  results_files = data.frame()
  covariances_files = data.frame()

  # combine all 22 individual files for each set (logs, beta, weights, covariances) 
  out_dir = intermediate_out_dir 
  for (i in 1:22){
    # logs 
    log_files_name = paste(out_dir, project, '_chr', i, '_elasticNet_model_log.txt', sep = '')
    if (file.exists(log_files_name)){
      log_file = readTextFile(log_files_name)
      log_files = rbind(log_files, log_file)
      # file.remove(log_files_name)
    } else {
      message(paste('Warning: log file (', log_files_name,') does not seem to exist!', sep = ''))
    }

    # beta 
    beta_files_name = paste(out_dir, 'TW_', project, '_elasticNet_alpha', alpha, '_', snpset, '_weights_chr', i, '.txt', sep = '')
    if (file.exists(beta_files_name)){
      beta_file = readTextFile(beta_files_name)
      beta_files = rbind(beta_files, beta_file)
      # file.remove(beta_files_name)
    } else {
      message(paste('Warning: beta file (', beta_files_name,') does not seem to exist!', sep = ''))
    }

    # weights 
    results_files_name = paste(out_dir, 'TW_', project, '_chr', i, '_exp_', n_k_folds, '-foldCV_elasticNet_alpha', alpha, '_', snpset, '.txt', sep = '')
    if (file.exists(results_files_name)){
      results_file = readTextFile(results_files_name)
      results_file = na.omit(results_file)  # remove NA 
      results_files = rbind(results_files, results_file)
      # file.remove(results_files_name)
    } else {
      message(paste('Warning: weight file (', results_files_name,') does not seem to exist!', sep = ''))
    }

    # covariances 
    covariances_files_name = paste(out_dir, project, '_chr', i, '_snpset_', snpset, '_alpha_', alpha, '_covariances.txt', sep = '')
    if (file.exists(covariances_files_name)){
      covariances_file = readTextFile(covariances_files_name)
      covariances_files = rbind(covariances_files, covariances_file)
      # file.remove(covariances_files_name)
      } else {
      message(paste('Warning: covariance file (', covariances_files_name,') does not seem to exist!', sep = ''))
    }

      # remove workingbest file 
      workingbest_files_name = paste(out_dir, "working_TW_", project, '_exp_', n_k_folds, '-foldCV_elasticNet_alpha', alpha, '_', snpset, '_chr', i, '.txt', sep = '')
    if (file.exists(workingbest_files_name)){
      # file.remove(workingbest_files_name)
    } else {
      message(paste('Warning: working best file (', workingbest_files_name,') does not seem to exist!', sep = ''))
    }
  }


  #---------------------------------
  # OUTPUT ONE: 1 covariance file
  #---------------------------------
  colnames(covariances_files) = c('GENE', 'RSID1', 'RSID2', 'VALUE')
  write.table(covariances_files, file = paste(output_data_dir, project, '_covariances.txt', sep=''), quote = FALSE, sep = ' ', row.names = FALSE)


  #-------------------------------------------------------------------------------------------------------------
  # OUTPUT TWO: 1 sqlite db that contains four tables including construction, extra, weights and sample info
  #-------------------------------------------------------------------------------------------------------------
  # setup sqlite driver and open connections
  driver <- dbDriver("SQLite")
  new_conn <- dbConnect(drv = driver, dbname = paste(output_data_dir, project, '_sqlite.db', sep=''))

  # construction table 
  log_files <- log_files %>% select(chr, seed_for_cv, n_genes, alpha)
  colnames(log_files) <- c('chr', 'cv_seed', 'n_genes', 'alpha')
  dbWriteTable(new_conn, "construction", log_files)
  dbGetQuery(new_conn, "CREATE INDEX construction_chr ON construction (chr)")

  # extra table
  # calculate fdr, and filter rows.
  extra_df <- results_files %>% filter(gene %in% gene_annot$gene_id)
  extra_filtered <- extra_df
  tryCatch (
      {    # with fdr filtering 
        qobj <- qvalue(extra_df$pval, fdr.level = fdr_level)  
        extra_df$pred.perf.qval <- qobj$qvalues
        extra_df$significant <- qobj$significant
        extra_filtered <- extra_df %>% rename(pred.perf.pval = pval, n.snps.in.model = n.snps, pred.perf.R2 = R2) %>% filter(significant == TRUE) %>% select(-significant)
        extra_filtered = extra_filtered %>% select(gene, n.snps.in.model, pred.perf.R2, pred.perf.pval, genename, pred.perf.qval)
      }, 
      error = function(e) {          # no fdr filtering 
        print(e)
        print('without fdr filtering!')
        extra_filtered <- extra_df %>% rename(pred.perf.pval = pval, n.snps.in.model = n.snps, pred.perf.R2 = R2)
        extra_filtered = extra_filtered %>% select(genename, gene, pred.perf.R2, pred.perf.pval, n.snps.in.model)
      }
    )
  sig_genes <- extra_filtered$gene
  dbWriteTable(new_conn, "extra", extra_filtered)
  dbGetQuery(new_conn, "CREATE INDEX extra_gene ON extra (gene)")

  # weights table
  # drop all rows pertaining to insignificant genes based fdr filtering.
  colnames(beta_files) = c('gene', 'rsid', 'ref_allele', 'eff_allele', 'weight', 'alpha')
  weights_filtered <- beta_files %>% filter(gene %in% sig_genes) %>% select(one_of(c("rsid", "gene", "weight", "ref_allele", "eff_allele")))
  dbWriteTable(new_conn, "weights", weights_filtered)
  dbGetQuery(new_conn, "CREATE INDEX weights_rsid ON weights (rsid)")
  dbGetQuery(new_conn, "CREATE INDEX weights_gene ON weights (gene)")
  dbGetQuery(new_conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")

  # sample info table 
  n_samples <- n_samples
  n_folds_cv <- n_k_folds
  snpset <- snpset
  alpha <- alpha
  metadata_file <- data.frame(n_samples, n_folds_cv, snpset, alpha)
  colnames(metadata_file) <- c("n_samples", "n_folds_cv", "snpset", "alpha")
  dbWriteTable(new_conn, "sample_info", metadata_file)

  # close connections
  dbDisconnect(new_conn)
}


