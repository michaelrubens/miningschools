This is the replication folder for 
"Management, Productivity, and Technology Choices: Evidence from U.S. Mining Schools". (RAND Journal of Economics)

Author: Michael Rubens (KU Leuven). Contact: michael.rubens@outlook.com

Date written: February 2, 2022

The file 'miningschools_master.do' is the master file that runs all the code.

I/ The folder 'main' contains the code to construct all tables and figures in the main text:

- miningschools_data.dta 	contains the dataset
- miningschools_figures.do 	makes figures 1 & 2
- miningschools_sumstat.do	makes the summary statistics
- miningschools_main.do 	makes tables 1-3 and figure 3
- miningschools_prodfun_bs.do	bootstraps standard errors around the production function coefficients

II/ The folder 'appendix' contains the code to construct all tables and figures in the main text:

- miningschools_appendix.do 	is the main do-file for the appendix tables
- miningschools_acf.do		estimates the production function using Ackerberg, Caves and Frazer (2015)
- miningschools_pf_acf.do	contains the MATA files for ACF production function estimation
- miningschools_acf_bs.do	bootstraps standard errors around the production function coefficients 

The bootstrapping files are currently commented out to run the program quickly. 