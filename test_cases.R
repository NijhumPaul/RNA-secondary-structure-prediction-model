source("parameter_percentage.R")

#this is the right csv directory
output_protein_parameters("./HC1_BB1.csv")

#if the name of the csv file directory is wrong
output_protein_parameters("wrong_csv_file_name")

#if the accession ID does not have any structure table or it is missing any parameters
percentage_of_parameters("Q09666")

#if the accession ID doesn't exist in the website that is extracted
percentage_of_parameters("wrong_accessionID_name")

#if the csv file directory is null
percentage_of_parameters(" ")

#if the csv file directory is null
percentage_of_parameters("")

#if the name of the csv file directory is wrong
protein_list_input_file("wrong_csv_file_name")
