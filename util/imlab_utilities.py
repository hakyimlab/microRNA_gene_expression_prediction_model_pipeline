# All helper functions can be used through the pipeline processing  

__author__ = 'Jiamao Zheng <jiamaoz@yahoo.com>'

# Load modules
from datetime import datetime
import uuid as myuuid 
import subprocess, sys, os, glob, sqlite3, time, argparse, sys 

# log variables 
start_time = time.time() 
logfile = None
current_logs = None
logs_data_dir = '../logs/'

# a function to check whether a file exist or not 
def check_file(file_path):
    if not os.path.exists(file_path): 
      warning = "Please check the existence of file path - " + file_path
      add_log(warning)
      add_log('If not, please create this file directory!' )

#---------------------------------
# functions used to produce logs
#---------------------------------
# create a log file path 
def get_log_path(project):  
    log_path = logs_data_dir + project
    if not os.path.exists(log_path): 
        os.makedirs(log_path)
    return log_path

# create a new log file 
def open_log(project):
    global logfile
    global current_logs

    logfile = open(get_log_path(project) + "/" + project + '.log', "w")
    current_logs = []

# add logs 
def add_log(log_str):
    global current_logs

    current_logs.append(log_str + '\n')

# output log 
def write_logs():
    global logfile
    global current_logs
    global start_time

    runing_time = "\nElapsed Time: " + timeString(time.time() - start_time) + "\n" # calculate how long the program is running
    print(runing_time)
    add_log(runing_time)

    for index in range (len(current_logs)): 
        logfile.write(current_logs[index])  

    current_logs = [] # clear up log 

# close log 
def finish_log():
    global logfile
    logfile.close()

# pretty string for a given number of seconds.
def timeString(seconds):
  tuple = time.gmtime(seconds);
  days = tuple[2] - 1;
  hours = tuple[3];
  mins = tuple[4];
  secs = tuple[5];
  if sum([days,hours,mins,secs]) == 0:
    return "<1s";
  else:
    string = str(days) + "d";
    string += ":" + str(hours) + "h";
    string += ":" + str(mins) + "m";
    string += ":" + str(secs) + "s";
  return string;

# get gene expression, gene annotation and snp annotation input files 
def get_input_path(input_path, input_dir, input_file):
    if len(input_path) != 0:                               # use user provided input file paths 
        input_file = input_path  
    else:                                                            # use script default input file paths 
        if len(glob.glob(input_dir + '*.txt')) == 1: 
            input_file = glob.glob(input_dir + '*.txt')[0]
        else: 
            print(input_file)
            sys.exit('Waring: Only one input file is allowed!')
    return input_file

def get_input_path_RDS(input_path, input_dir, input_file):
    if len(input_path) != 0:                               # use user provided input file paths 
        input_file = input_path  
    else:                                                            # use script default input file paths 
        if len(glob.glob(input_dir + '*.RDS')) == 1: 
            input_file = glob.glob(input_dir + '*.RDS')[0]
        else: 
            print(input_file)
            sys.exit('Waring: Only one input file is allowed!')
    return input_file

# get genotype files 
def get_genotype_input_path(input_path, input_dir, input_file):
    # genotype data (0-1-2 format)
    if len(input_path) != 0:                              # use user provided input file paths 
        if len(glob.glob(input_path + '*.txt')) == 1:
            input_file = glob.glob(input_path + '*.txt')[0]
        elif len(glob.glob(input_path + '*.txt')) == 22:
            input_file = input_path
        else: 
            sys.exit('Waring: Only 1 or 22 genotype input files are allowed!')
    else:                                                         # use script default input file paths 
        if len(glob.glob(input_dir + '*.txt')) == 1:
            input_file = glob.glob(input_dir + '*.txt')[0]
        elif len(glob.glob(input_dir + '*.txt')) == 22:
            input_file = input_dir
        else: 
                sys.exit('Waring: Only 1 or 22 genotype input files are allowed in the input directory of ../input/genotype/!')
    return input_file



