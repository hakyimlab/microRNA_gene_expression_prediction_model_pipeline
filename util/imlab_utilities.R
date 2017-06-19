# Author = 'Jiamao Zheng <jiamaoz@yahoo.com>'

# load R libraries 
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(glmnet)))
suppressWarnings(suppressMessages(library(qvalue)))
suppressWarnings(suppressMessages(library(RSQLite)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(bit64)))
suppressWarnings(suppressMessages(library(doMC)))
registerDoMC(cores=4)

# arg variables
argv <- commandArgs(trailingOnly = TRUE)

# paste funtion 
"%&%" <- function(a,b) paste(a, b, sep = "")

# a function used to convert data.table to data.frame 
setDF <- function(x) {
        if (!is.data.table(x))
                stop("x must be a data.table")
        setattr(x, "row.names", .set_row_names(nrow(x)))
        setattr(x, "class", "data.frame")
        setattr(x, "sorted", NULL)
        setattr(x, ".internal.selfref", NULL)        
}

# a function used to read RDS file 
readRDSFile <- function(input_filename){
   tryCatch(
		input_data <- readRDS(input_filename), 
		error=function(e){
      message('Error reading with file. Aborting. Detail: ', e)
		}
	)
   return(input_data)
}

# a function used to read txt file using big.table 
readTextFile <- function(input_filename){
     tryCatch(
      {
  			input_data <- fread(input_filename, fill = TRUE)
        setDF(input_data) 
      }, 
			error=function(e){
        message('Error reading with file. Aborting. Detail:  ', e)
			}
	  )
      return(input_data)
}

# a function used to read txt file using big.table 
readTableTextFile <- function(input_filename){
     tryCatch(
      {
        input_data <-read.table(input_filename)
      }, 
      error=function(e){
        message('Error reading with file. Aborting. Detail:  ', e)
      }
    )
      return(input_data)
}