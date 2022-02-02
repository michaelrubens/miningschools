/* Management, Productivity, and Technology Choices: Evidence from U.S. Mining Schools

								-	MASTER DO FILE   -

								Michael Rubens (KU Leuven)
===================================================================================================*/

set more off
use ./main/miningschools_data, clear

/*	FIGURES
-------------------------------------------------------------------------------- */	

do ./main/miningschools_figures											// make figures 1,2,4  
	
/*	SUMMARY STATISTICS
-------------------------------------------------------------------------------- */	

do ./main/miningschools_sumstat											// make summary statistics tables
	
/*	RESULTS IN MAIN TEXT
-------------------------------------------------------------------------------- */	

do ./main/miningschools_main												// make tables 1,2,3 and figure 5
exit
	 
/*	RESULTS IN APPENDICES
-------------------------------------------------------------------------------- */	

do ./appendix/miningschools_appendix  											// appendix tables and figures
