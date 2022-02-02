
* Program to get bootstrapped standard errors on cost share estimates

	clear mata
	clear matrix
	local B =  250							// number of bootstrap iterations
	set more off

	use ./main/temp/data_temp, clear	
	drop *bl* *bm*
	set matsize 6000
	save ./main/temp/data_bs, replace
	
	** Matrices to store estimates
	* Production model
	foreach coef in "medbl_cs" "medbm_cs"   {
	matrix _`coef'  = J(`B', 1, 0)
	}

	forvalues b = 1/`B' {
	preserve

	use ./main/temp/data_bs, clear
	set seed `b'
	bsample , cluster(mineid) idcluster(midbis)

	egen medbl_cs = median((w*emp*ndays) / (p*q))	
	 egen medbm_cs = median(qinp/q)	

	foreach var of varlist  medbl_cs medbm_cs {
	egen av`var' = mean(`var')
	mkmat av`var', matrix(av`var')
    matrix _`var'[`b',1] = av`var'[1,1]			// store average markdown for this bootstrap iteration in matrix	
    }	
	restore
	}
	
	keep in 1/`B'

	set matsize  11000
	foreach coef in "medbl_cs" "medbm_cs"  {
	svmat _`coef'  
	}	
	
	sum _medb*
			
	egen bl_cs = sd(_medbl_cs)
	egen bm_cs = sd(_medbm_cs)
	
	estpost su bl_cs bm_cs
	est store A4a_se
	
	
