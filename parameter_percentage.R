#percentage_of_parameters: Calculates the parameters (bp, kpbp, hp, kphp, rc, kprc)
#protein_list_input_file: Read a CSV file and extract column 'Accession'

#install.packages('rvest', repos = "http://cran.us.r-project.org")
#install.packages('tidyverse', repos = "http://cran.us.r-project.org")
#install.packages('logging', repos = "http://cran.us.r-project.org")
#install.packages('log4r', repos = "http://cran.us.r-project.org")

#require('rvest')
#require('tidyverse')
#require('logging')
#require('log4r')


library(rvest)
library(tidyverse)
library(logging)
#library(testthat)
library(log4r)


#this function wil log all errors in both console and my_logfile.txt file
log_errors <- function(logFile_path){
  
  my_logfile = logFile_path
  my_console_appender = console_appender(layout = default_log_layout())
  my_file_appender = file_appender(my_logfile, append = TRUE, 
                                   layout = default_log_layout())
  
  my_logger <- log4r::logger(threshold = "INFO", 
                             appenders= list(my_console_appender,my_file_appender)) 
  return (my_logger)
  
}

'This function calculates the 6 parameters for each accession ID (bp, kpbp, hp, kphp, rc, kprc)
  bp-> Betasheet percentage
  kpbp-> known partial betasheet percentage
  hp-> helix percentage
  kphp-> known partial helix percentage
  rc-> random coil percentage
  kprc-> known partial random coil percentage
  Input is accession ID, i.e. P35579
  Output is a vector of 7 items (ID,bp, kpbp, hp, kphp, rc, kprc) of each accession ID'

percentage_of_parameters <- function(ID){
  logFile_path = "./my_logfile.txt"
  #error_message = "this is an error"
  my_logger <- log_errors(logFile_path)
  
  tryCatch(
    {
      #Scrape the following webpage
      url_old="https://www.uniprot.org/uniprot/accession" 
      url_new=gsub("accession", ID, url_old)
      
      #Read the 'secondary structure' table
      sec_structure_table<- url_new %>%
        read_html() %>%
        html_nodes(xpath='//*[@id="secstructure_section"]') %>%
        html_table(fill=TRUE,header=TRUE)
      
      #Read the total length
      seq.length<- url_new %>%
        read_html() %>%
        html_nodes(xpath='//*[@id="secondarystructure"]/svg/text[2]')%>%
        html_text()
      
      #Percentage calculation
      sec_table = as.data.frame(sec_structure_table[[1]][,1])
      #print(sec_table)
      for (keyfeatures in sec_table) {
        keyfeatures = gsub(".*>","",keyfeatures)
      }
      
      sec_table_new = as.data.frame(keyfeatures)
      #keyfeature = as.data.frame(keyfeature)
      
      all_features_length = sec_structure_table[[1]][,5]
      sum_of_all_features_length=sum(all_features_length)
      
      all_features_length_df = as.data.frame(all_features_length)
      
      new_df = cbind(sec_table_new, all_features_length_df)
      
      #keyfeature=gsub(".*>","",sec_structure_table[[1]][,1])
      all_features_length = sec_structure_table[[1]][,5]
      sum_of_all_features_length=sum(all_features_length)
      
      #betaStrand = sum(all_features_length[which(keyfeature=="Beta strandi")])
      betaStrand= sum(new_df$Length[new_df$keyfeatures=="Beta strandi"])
      bp = (betaStrand/ sum_of_all_features_length)*100  
      kpbp = (betaStrand/ (as.integer(seq.length)))*100  
      
      helix= sum(new_df$Length[new_df$keyfeatures=="Helixi"])
      #helix = sum(all_features_length[which(keyfeature=="Helixi")]) 
      hp = (helix/ sum_of_all_features_length)*100    
      kphp = (helix/ (as.integer(seq.length)))*100    
      
      #random_coil = sum(all_features_length[which(keyfeature=="Turni")])
      random_coil= sum(new_df$Length[new_df$keyfeatures=="Turni"])
      rc = (random_coil/ sum_of_all_features_length)*100
      kprc = (random_coil/ (as.integer(seq.length)))*100      
      
      return(c(bp, kpbp, hp, kphp, rc, kprc))
    },
    error=function(error_message) {
      #message(error_message)
      #logerror(error_message)
      log4r::error(my_logger, error_message) 
      #readLines(my_logfile)
      return(NA)
    },
    warnings = function(warning_message){
      #message(warning_message)
      #logerror(warning_message)
      log4r::warnings(my_logger, warning_message)
      #readLines(my_logfile)
      return(NA)
    }
  )
}
percentage_of_parameters("Q99895")

'This function reads a CSV file and 
  extracts the column: Accession
  Input is the directory of a csv file
  Output is a list of accession IDs'
protein_list_input_file <- function(csv_dir){
  
  tryCatch(
    {
      #import csv file
      files = read.csv(file=csv_dir ,header=TRUE)
      #Extract column 'Accession' as a vector
      accession = files %>% pull(Accession) 
      IDs = accession
      ID = as.list(strsplit(IDs, " "))
      return(ID)
    },
    error=function(error_message) {
      #message(error_message)
      #logerror(error_message)
      log4r::error(my_logger, error_message) 
      #readLines(my_logfile)
      return(NA)
    },
    warnings = function(warning_message){
      #message(warning_message)
      #logerror(warning_message)
      log4r::warnings(my_logger, warning_message)
      #readLines(my_logfile)
      return(NA)
    }
  )
}

#protein_list_input_file("./HC1_BB1.csv")

#percentage_of_parameters("P35579")


#this function wraps all other function
#Input is the directory of a csv file
#Output is the dataframe of all parameters along with the accession ID

output_protein_parameters <- function(csv_dir){
  
  logFile_path = "./error_logfile.txt"
  
  my_logger <- log_errors(logFile_path)
  
  #log4r::error(my_logger, error_message)
  #my_logger <- log_errors(logFile_path)
  tryCatch(
    {
      
      ID = protein_list_input_file(csv_dir)
      data = lapply(ID, percentage_of_parameters)
      parameters_in_dataframe = data.frame(ID= sapply( ID, "[", 1), 
                                           bp =as.numeric( sapply( data, "[", 1) ), kpbp =as.numeric( sapply( data, "[", 2) ), 
                                           hp =as.numeric( sapply( data, "[", 3) ), kphp =as.numeric( sapply( data, "[", 4) ),
                                           rc =as.numeric( sapply( data, "[", 5) ), kprc =as.numeric( sapply( data, "[", 6) ))
      #readLines(my_logfile)
      #write.csv(parameters_in_dataframe,"./demo_data.csv", row.names = FALSE)
      return (write.csv(parameters_in_dataframe,"./parameters.csv", row.names = FALSE))
    },
    error=function(error_message) {
      #message(error_message)
      #logerror(error_message)
      log4r::error(my_logger, error_message) 
      #readLines(my_logfile)
      return(NA)
    },
    warnings = function(warning_message){
      #message(warning_message)
      #logerror(warning_message)
      log4r::warnings(my_logger, warning_message)
      #readLines(my_logfile)
      return(NA)
    }
    
  )
  
}
#output_protein_parameters("./HC1_BB1.csv")
#write.csv(output_protein_parameters,"./demo_data.csv", row.names = FALSE)
