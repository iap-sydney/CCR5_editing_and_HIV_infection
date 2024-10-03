#Steffen Docken
#27-1-23
#This code will run the various functions or the project on modeling how 
#effective various levels of CCR5 KO are in protecting humanized mice from 
#infection

rm(list=ls())
graphics.off()


save_output = TRUE #indicate whether to save output or not

# Load Global Parameters and libraries ------------------------------------

source('setup.R')


# Import and clean data ---------------------------------------------------

source('processing/01_CCR5KO_data_import.R')

# Run analysis ------------------------------------------------------------

source('analysis/01_CCR5_KO_fitting.R')

source('analysis/02_Model_comparison_output.R')

source('analysis/03_survival_bootstrap_plotting.R')

source('analysis/04_VL_plotting.R')

source('analysis/05_Implication_plots.R')
