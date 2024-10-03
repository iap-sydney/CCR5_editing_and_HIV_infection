#Steffen Docken
#13-3-23
# -------------------------------------------------------------------------
#' This file is used to set up everything that's needed across the project
#' It loads libraries, creates functions, sets themes and defaults
#' 
# -------------------------------------------------------------------------


# libraries ---------------------------------------------------------------
library(GGally)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(latex2exp)
library(stats)
library(survival)
library(survminer)
library(lhs)
library(ggnewscale)
library(scales)
library(foreach)
library(doParallel)


# Data location -----------------------------------------------------------

data_file = "Data/2023-02-07_CCR5-KO aggregate data_(with_max_pVLs).csv"

# Global parameters -------------------------------------------------------
Arm_ind_link <- c("100% CCR5-KO" = 1, "75% CCR5-KO" = 2, "50% CCR5-KO" = 3,
                  "25% CCR5-KO" = 4, "0% CCR5-KO" = 5)

num_pred_var = 1 #number of predictor variables to be used
num_model_forms = 9 #number of forms for exponential terms

n_LHS = 1e3 #number of LHS samples
n_pars = 5 #maximum number of parameters considered

LHS_max = 20
LHS_min = -10 #max and min for range of each parameter in LHS sampling

## Generating Latin Hypercube Sample of parameter space to do preliminary search
set.seed(1111) #setting the seed so LHS is reproducible
unit_LHS = randomLHS(n = n_LHS, k = n_pars)
LHS_param_matrix = LHS_min + (LHS_max - LHS_min)*unit_LHS #converting sampling
#from unit interval [0,1] to considered parameter range [-10, 20]
#Each row will correspond to a different sampling and each column a different 
#parameter

n_bootstrap = 1000 #number of bootstrap data sets to be generated
#10 #for debugging #

b_h = 8 # published R_0 estimate for HIV (Ribeiro et al. 2010)

# Color schemes -----------------------------------------------------------
#color palettes:
# The palette with grey:
arm_palette <- c("#018E00", "#10C317", "#7C81F6", "#4293F7", "#1332F5")

names(arm_palette) <- factor(c("Arm=0% CCR5-KO", "Arm=25% CCR5-KO", 
                             "Arm=50% CCR5-KO", "Arm=75% CCR5-KO",
                             "Arm=100% CCR5-KO"))
Surv_arm_palette_surv <- scale_colour_manual(name = "% CCR5-KO",values = arm_palette,
                                      breaks = c("Arm=0% CCR5-KO", "Arm=25% CCR5-KO", 
                                                 "Arm=50% CCR5-KO", "Arm=75% CCR5-KO", 
                                                 "Arm=100% CCR5-KO"),
                                      labels = c("0%", "25%", "50%", "75%", 
                                                 "100%"),
                                      aesthetics = c("color", "fill"),
                                      limits = force)

names(arm_palette) <- factor(c("0% CCR5-KO", "25% CCR5-KO", 
                             "50% CCR5-KO", "75% CCR5-KO",
                             "100% CCR5-KO"))
Surv_arm_palette <- scale_colour_manual(name = "% CCR5-KO",values = arm_palette,
                                           breaks = c("0% CCR5-KO", "25% CCR5-KO", 
                                                      "50% CCR5-KO", "75% CCR5-KO", 
                                                      "100% CCR5-KO"),
                                           labels = c("0% CCR5-KO", "26% CCR5-KO", "54% CCR5-KO",
                                                      "69% CCR5-KO","96% CCR5-KO"),
                                        #labels based on mean percentage
                                           aesthetics = c("color","fill"),
                                           limits = force)

# For Parallelization -----------------------------------------------------

n_cores = detectCores() #number of cores on computer

# Load local functions ----------------------------------------------------
source("R_functions/Prob_infect_functions.R")
source("R_functions/AIC_generic.R")
source("R_functions/param_conf_functions.R")
source("R_functions/surv_KM_bootstrap_functions.R")


# Create output directory structure if it doesn't exist---------------------------------------

if(!dir.exists("output")){
  dir.create("output")
}

output_subdirs = c('Figures', 'Bootstrap_Results')

for (ii in 1:2) {
  if(!dir.exists(paste0("output/", output_subdirs[ii]))){
    dir.create(paste0("output/", output_subdirs[ii]))
  }
}


