/* Management, Productivity, and Technology Choices: Evidence from U.S. Mining Schools

							-	RESULTS IN THE MAIN TEXT  -

								Michael Rubens (KU Leuven)
===================================================================================================*/

set more off
*use miningschools_data, clear
 
/* ------------------------------------------------------------------------------
A. Production function (table 2)
-------------------------------------------------------------------------------- */	
	 
	* locomotive dummies
	
	foreach var of varlist locel locst locair {
	gen d`var' = `var'>0
	replace d`var' = . if `var'==.
	}
 
	* take logarithms of variables
 
	gen y = log(q)
	gen l = log(emp   )
	gen inp =  qinp
	gen m = log(inp)
	gen lpowd = log(powd+1)
	gen ke = dlocel 
	gen ks = dlocst
	gen ka = dlocair
	gen x = minsch
	
  	drop if yr<1900	// locomotives unobserved prior to 1900
	
	* calculate variable inputs' cost shares
	
	gen bl_cs = (w*emp*ndays) / (p*q)			// factor revenue share of labor
	gen bm_cs = qinp /q 						// factor revenue share of intermediate inputs
	sum bl_cs bm_cs					
	  
	* Table 1(a): output elasticities of variable inputs
	******************************************************
	
	egen avbl_cs = mean(bl_cs)
	egen avbm_cs = mean(bm_cs)
	egen medbl_cs = median(bl_cs)
	egen medbm_cs = median(bm_cs)
	bys yr: egen medbl_csyr = median(bl_cs)
	bys yr: egen medbm_csyr = median(bm_cs)
	
	label var medbl_cs "Labor"
	label var medbm_cs "Materials"
	twoway connect medbl_csyr medbm_csyr yr, lwidth(thick thick) msize(normalsize normalsize) mfcolor(white white) mcolor(black black) msymbol(square diamond) lcolor(black black) lpattern(solid longdash)  graphregion(color(white)) ytitle("Median revenue share") 
//	graph export .\Paper\fig\fig_revshare.png, replace
		
	drop bl_cs bm_cs
	gen bl_cs = medbl_cs
	gen bm_cs = medbm_cs

	* bootstrapping standard errors around production function coefficients
	
	save ./main/temp/data_temp, replace	
	// do ./main/miningschools_prodfun_bs			// uncomment to bootstrap standard errors
 	use ./main/temp/data_temp, clear

 	estpost su  bl_cs bm_cs 		// store estimates for table 2(a)
	est store A_4a					// 
	
	* Table 1(b): output elasticities of fixed inputs
	******************************************************
	
	preserve
	keep if y~=. & l~=. & m~=. & ke~=. & ks~=. & ka~=. & x~=. & col~=.
	di _N
	restore

	* ols - table 2(b) column (I) 
	
	gen ltfp_cs = y - medbl_cs  *l  - medbm_cs  *m		
	reg ltfp_cs x  ke ka ks yr col		
	local N_pf1 = e(N)

	foreach var of varlist x ke ks ka col yr {
	gen b`var'_cs_1 = _b[`var']
	gen se`var'_cs_1 = _se[`var']
	}
	gen bc_cs_1 = _b[_cons]
	gen bt_cs_1 = byr_cs_1
	
	* fe - table 2(b) column (II)
	
	xtreg ltfp_cs x ke ks ka col yr , fe r		
	
	local N_pf2 = e(N)
	foreach var of varlist x ke ks ka yr col{
	
	gen b`var'_cs_2 = _b[`var']					
	gen se`var'_cs_2 = _se[`var']
	}
	gen bc_cs_2 = _b[_cons]
	gen bt_cs_2 = byr_cs_2
	
	xtset mineid yr
  
	* ar(1) model - table 2(b) column (III)

	gmm (ltfp_cs - {rho}*L.ltfp_cs - ({bke})*(ke-{rho}*L.ke) - {bks}*(ks- {rho}*L.ks)  - {bka}*(ka - {rho}*L.ka)  - {bx}*(x- {rho}*L.x)  - {bc}*(col- {rho}*L.col)  ///
	- {bt}*(yr - {rho}*L.yr)- {c}*(1-{rho}) )    ,  	inst(L.l L.m ke ka ks   x  col  yr     )  
   
	local N_pf3 = e(N)
	mat coef_cs = e(b)
	mat cov_cs = e(V)
	gen bke_cs_3 = coef_cs[1,2]
	gen bks_cs_3 = coef_cs[1,3]
	gen bka_cs_3 = coef_cs[1,4]
	gen bx_cs_3 = coef_cs[1,5]
	gen bcol_cs_3 = coef_cs[1,6]
	gen bt_cs_3 = coef_cs[1,7]
	gen bc_cs_3 = coef_cs[1,8]
	gen seke_cs_3 = sqrt(cov_cs[2,2])
	gen seks_cs_3 = sqrt(cov_cs[3,3])
	gen seka_cs_3 = sqrt(cov_cs[4,4])
	gen sex_cs_3 = sqrt(cov_cs[5,5])
	gen secol_cs_3 = sqrt(cov_cs[6,6])
	gen set_cs_3 = sqrt(cov_cs[7,7])
	gen sec_cs_3 = sqrt(cov_cs[8,8])
 
	save ./main/temp/data_temp, replace	
	
	* r-squared of each regression
	
	forvalues n = 1/2 {
	preserve 
	gen yhat = y - bke_cs_`n'*ke -bks_cs_`n'*ks -bka_cs_`n'*ka -bx_cs_`n'*x -bt_cs_`n'*yr -bc_cs_`n' -bl_cs*l - bm_cs*m - bcol_cs_`n'*col
	reg y yhat
	gen r2 = e(r2)	
	sum r2 
	local r2_pf`n' = round(r2,0.001)
	restore
	}
  
	preserve 
	local n = 3
	gen yhat = y - bke_cs_`n'*ke -bks_cs_`n'*ks -bka_cs_`n'*ka -bx_cs_`n'*x -bt_cs_`n'*yr -bc_cs_`n' -bl_cs*l - bm_cs*m - bcol_cs_`n'*col
	xtset mineid yr
	drop if yhat==. | L.yhat==.
	reg y yhat
	gen r2 = e(r2)	
	sum r2 
	local r2_pf`n' = round(r2,0.001)
	restore
	
	* locomotive coefficients
	
	gen bke = bke_cs_3
	gen bka = bka_cs_3
	gen bks = bks_cs_3
	gen bx = bx_cs_3
	gen bt = bt_cs_3
	gen bcol = bcol_cs_3
 
	forvalues m = 1/3 {
	rename (bx_cs_`m' bke_cs_`m' bka_cs_`m' bks_cs_`m' bcol_cs_`m') (bx_cs bke_cs bka_cs bks_cs bcol_cs)
	estpost su bx_cs bcol_cs bke_cs bka_cs bks_cs 		// store point estimates for table 2b
	est store A_4b_`m'
	drop bx_cs bke_cs bka_cs bks_cs bcol_cs
	}
	
	forvalues m = 1/3 {
	rename (sex_cs_`m' seke_cs_`m' seka_cs_`m' seks_cs_`m' secol_cs_`m') (bx_cs bke_cs bka_cs bks_cs bcol_cs)
//	estpost su bx_cs bcol_cs bke_cs bka_cs bks_cs		// store SE estimates for table 2b
//	est store A_4b_`m'_se
	drop bx_cs bke_cs bka_cs bks_cs bcol_cs
	}
	
	foreach var in "bx_cs" "bke_cs" "bka_cs" "bks_cs" "bcol_cs" {
	gen `var' = .
	}
	label var bx_cs "1(Mining col. grad.)"
	label var bcol_cs "1(Other grad.)"
	label var bke_cs "1(Elec. loc.)"
	label var bka_cs "1(Air loc.)"
	label var bks_cs "1(Steam loc.)"
	label var bl_cs "Labor"
	label var bm_cs "Materials"
	
	  * write table 2(a)
//	esttab  A_4a  A4a_se  using "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\tab\tab4.tex", replace ///
//	mtitle("" "" ) prehead(   &  \multicolumn{2}{c}{ (I)}& \multicolumn{2}{c}{ (II)} & \multicolumn{2}{c}{ (III)} \\\textit{(a) Variable inputs}  & Estimate & S.E. & &&& \\ \hline)  posthead(   \\)    ///
//	cells(mean(fmt(3))) label booktabs nonum noobs collabels(none) gaps f    prefoot(  &&&\\      ) 
	
 
	* write table 2(b)
//	esttab  A_4b_1 A_4b_1_se A_4b_2 A_4b_2_se A_4b_3 A_4b_3_se  using "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\tab\tab4.tex", append ///
//	mtitle("" "" "" "" "" "") prehead(  \hline	 & \multicolumn{6}{c}{ log(TFP)} \\\textit{(b) Fixed inputs}  & Estimate & S.E. & Estimate & S.E.& Estimate & S.E. \\ \hline)  posthead(   \\)    ///
//	cells(mean(fmt(3)))  label booktabs nonum noobs collabels(none) gaps f   prefoot( \\ Model &    \multicolumn{2}{c}{ OLS}& \multicolumn{2}{c}{FE} & \multicolumn{2}{c}{ AR(1)} ///
//	\\Observations &    \multicolumn{2}{c}{`N_pf1'}& \multicolumn{2}{c}{`N_pf2'} & \multicolumn{2}{c}{`N_pf3'}\\ R-squared &    \multicolumn{2}{c}{`r2_pf1'}& \multicolumn{2}{c}{`r2_pf2'} & \multicolumn{2}{c}{`r2_pf3'}  \\)

 /* ------------------------------------------------------------------------------
 B. Interaction effects (table 3a)
-------------------------------------------------------------------------------- */	

	* Table 3(a): Interaction effect
	*********************************

	 gen xke = x*ke
	  
	  gmm (ltfp_cs - {rho}*L.ltfp_cs - ({bke})*(ke-{rho}*L.ke) - {bks}*(ks- {rho}*L.ks)  - {bka}*(ka - {rho}*L.ka)  - {bx}*(x- {rho}*L.x) ///
	 - {bxke}*(xke- {rho}*L.xke)  - {bt}*(yr - {rho}*L.yr)- {bcol}*(col - {rho}*L.col)- {c}*(1-{rho}) )   , ///
		inst(ke  ka   ks   x    yr  xke L.l L.m col )
    
		local N_inter = e(N)
		
		mat inter_coef = e(b)
		mat inter_cov = e(V)
		mat list inter_coef
		
		gen bkxe = inter_coef[1,6]
	 	gen sekxe = sqrt(inter_cov[6,6])

		estpost su bkxe
		est store D_kxe	// store estimates for table 3a
		
		drop bkxe
		rename sekxe bkxe
		estpost su bkxe
		est store D_kxe_se		// store SE estimates for table 3a
		 
/* ------------------------------------------------------------------------------
C. Technology choices (table 1)
-------------------------------------------------------------------------------- */	
	 
	gen const = 1

	drop ltfp*
	gen ltfp = lq - bl_cs*l - bm_cs*m - bke *ke-bks *ks -bka *ka -bx *x -bcol*col -bt  *yr -bc_cs_3 
 	gen tfp=exp(ltfp)

	xtset mineid yr

	label var ke "Electrical"
	label var ks "Steam"
	label var ka "Air"
	
	* Table 1(a): Extensive margin
	*******************************
	
	gen lw = log(w)
	gen lp = log(p)
 
 	* technology choice regressions
	foreach var of varlist ke ks ka {
	xtreg  `var' x  col   i.yr ltfp  l m lw lp , fe r			// wages and prices only vary by year so absorbed by year dummies. xtreg gives wrong SEs but needed to get the within R^2 
	local r2_`var' = round(e(r2_w) ,0.001) 
	areg  `var' x  col   i.yr ltfp  l m lw lp, absorb(mineid) cluster (man_id)		// areg to get correct clustered standard errors
	predict `var'hat1 if x ==1															
	predict `var'hat0 if x ==0															// predicted usage probabilities
	local N_`var' = e(N)
	gen cx`var' = _b[x]
 	gen co`var' = _b[col]
	gen cx`var'_se = _se[x]
 	gen co`var'_se = _se[col]
	}
  	
	estpost su  cxke  coke			// save electrical locomotive estimates for table 1(a) column (I)
	est store C_ke 

	foreach var of varlist ka ks {
	replace cxke = cx`var'
 	replace coke = co`var'
	estpost su  cxke  coke			// same for air and steam locomotives, table 1(a) columns (II)-(III)
	est store C_`var' 
	}

	replace cxke = cxke_se
	replace coke = coke_se
	estpost su  cxke coke
	est store C_ke_se 		// save SEs for electrical locomotives, table 1(a) column (I)

	foreach var of varlist ka ks {
	replace cxke = cx`var'_se
	replace coke = co`var'_se
	estpost su  cxke coke
	est store C_`var'_se		// save SEs for other locomotive types, table 1(a) column (II)-(III)
	}
	
	* average usage of each type
	sum ke ka ks
	foreach var of varlist ke ka ks {
	egen av`var' = mean(`var')
	local av`var' = round(av`var',0.001)
	}		
	
	egen avprobel1 = mean(kehat1)  
	egen avprobel0 = mean(kehat0)  
	
	/*
	* Write table 1a: 
	label var cxke "1(Mining col. grad.) "
	label var coke "1(Other grad.) "
	esttab  C_ke C_ke_se C_ka C_ka_se C_ks C_ks_se using "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\tab\tab1.tex", replace ///
	mtitle("" "" "" "" "" "") prehead(  &  \multicolumn{2}{c}{ (I)} &  \multicolumn{2}{c}{ (II)} &  \multicolumn{2}{c}{ (III)} \\ &  \multicolumn{2}{c}{1(Elec. loc.)}& \multicolumn{2}{c}{1(Air. loc.)} & \multicolumn{2}{c}{1(Steam loc.)} \\  \textit{(a) Extensive margin}& Estimate & S.E. & Estimate & S.E. & Estimate & S.E.\\ \hline) ///
	posthead(   \\)   cells(mean(fmt(3))) label booktabs nonum noobs collabels(none) gaps f   ///
	prefoot(  &&&\\  Average usage & \multicolumn{2}{c}{`avke'}  & \multicolumn{2}{c}{`avka'} & \multicolumn{2}{c}{`avks'} \\  Observations & \multicolumn{2}{c}{  `N_ke'}  & \multicolumn{2}{c}{`N_ka'} & \multicolumn{2}{c}{`N_ks'} ///
	\\  Within R-squared &  \multicolumn{2}{c}{  `r2_ke'}&\multicolumn{2}{c}{`r2_ka'}&\multicolumn{2}{c}{`r2_ks'}\\ &&&\\   ) 
	 */
	 
	* Table 1(b): Intensive margin
	*******************************
	
	foreach var of varlist locel locair locst {		// take logarithms of locomotive counts
	gen l`var' = log(`var')
	}
	rename (llocel llocair llocst)(lke lka lks)
	drop cxk*  cok*
	
	foreach var of varlist lke lks lka {
	xtset mineid yr
	xtreg `var' x col i.yr ltfp  l m, fe r
	local r2_`var' = round(e(r2_w) ,0.001) 
	areg `var' x col i.yr ltfp l m, absorb(mineid) cluster (man_id)					// regress log number of locomotives on mining college graduates and other covaraites
	local N_`var' =  e(N)   
	predict `var'hat, xb															// predicted usage probability
	gen cx`var' = _b[x]
	gen co`var' = _b[col]
	gen cx`var'_se = _se[x]
	gen co`var'_se = _se[col]
	}
  
	estpost su  cxlke colke 
	est store C_lke 

	foreach var of varlist lka lks {
	replace cxlke = cx`var'
	replace colke = co`var'
	estpost su  cxlke  colke
	est store C_`var' 
	}

	replace cxlke = cxlke_se
 	replace colke = colke_se
	estpost su  cxlke  colke
	est store C_lke_se 

	foreach var of varlist lka lks {
	replace cxlke = cx`var'_se
	replace colke = co`var'_se
	estpost su  cxlke  colke
	est store C_`var'_se
	}
	
	label var cxlke "1(Mining col. grad.) "
 	label var colke "1(Other grad.) "
	
	* Write table 1
	/*
	esttab  C_lke C_lke_se C_lka C_lka_se C_lks C_lks_se using "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\tab\tab1.tex", append ///
	mtitle("" "" "" "" "" "") prehead( \hline   &  \multicolumn{2}{c}{log(Elec. loc.)}& \multicolumn{2}{c}{log(Air. loc.)} & \multicolumn{2}{c}{log(Steam loc.)} \\ \textit{(b) Intensive margin}& Estimate & S.E. & Estimate & S.E. & Estimate & S.E.\\ \hline) ///
	posthead(   \\)   cells(mean(fmt(3))) label booktabs nonum noobs collabels(none) gaps f   ///
	prefoot(  &&&\\   Observations & \multicolumn{2}{c}{`N_lke'}  & \multicolumn{2}{c}{`N_lka'} & \multicolumn{2}{c}{`N_lks'} \\  Within R-squared &  \multicolumn{2}{c}{`r2_lke'}&\multicolumn{2}{c}{`r2_lka'}&\multicolumn{2}{c}{`r2_lks'}\\ &&&\\   ) 
	*/
	
		 
/* ---------------------------------
D. Quantify returns to managers
-----------------------------------*/

	gen return = (exp(bke)*avprobel1 + (1-avprobel1))/(exp(bke)*avprobel0 + (1-avprobel0))
	sum return  // factor 1.0298 = increase of 3.0%
	   
	
/* -------------------------------------------------------------------------------------
E. Event study (Figure 5)
-------------------------------------------------------------------------------- */	
 	
	* Locomotive adoption regressions with various lags and leads
	 
 	foreach var of varlist ke ka ks  {
	reg D.`var'  F3D.x F2D.x F1D.x D.x L1D.x L2D.x L3D.x      i.yr   ,  cluster(man_id)
	forvalues n = 1(1)3 {
	local p = `n'+3
	local k = 4-`n'
	gen th_`var'_`k' = _b[F`n'D.x]
	gen se_`var'_`k' = _se[F`n'D.x]
	gen th_`var'_`p' = _b[L`n'D.x]
	gen se_`var'_`p' = _se[L`n'D.x]
	}
	}
	forvalues n = 1(1)6 {
	foreach var in "ke" "ks" "ka" {
	gen ci_`var'_lo_`n' = th_`var'_`n'-se_`var'_`n'*1.96
	gen ci_`var'_hi_`n' = th_`var'_`n'+se_`var'_`n'*1.96
	}
	}		

	foreach var of varlist l {
	reg D.`var'  F3D.x F2D.x F1D.x D.x L1D.x L2D.x L3D.x         ,  cluster(man_id)
	forvalues n = 1(1)3 {
	local p = `n'+3
	local k = 4-`n'
	gen th_`var'_`k' = _b[F`n'D.x]
	gen se_`var'_`k' = _se[F`n'D.x]
	gen th_`var'_`p' = _b[L`n'D.x]
	gen se_`var'_`p' = _se[L`n'D.x]
	}
	}
	forvalues n = 1(1)6 {
	foreach var in "l" {
	gen ci_`var'_lo_`n' = th_`var'_`n'-se_`var'_`n'*1.96
	gen ci_`var'_hi_`n' = th_`var'_`n'+se_`var'_`n'*1.96
	}
	}	
  
	* Plot the event study
	
	preserve 
	keep if _n==1
	keep th* se_* cons
	reshape long th_ke_ th_ks_ th_ka_ th_l_ se_ke_ se_ks_ se_ka_ se_l_, i(cons) j(b)
	gen t = b-4
	
	replace t = t+1 if t>=0	
	foreach var of varlist th*{
	replace `var'=0 if t==0
	}
	foreach var in "ke" "ks" "ka" "l"{
	gen ci_`var'_lo_ = th_`var'_-se_`var'_*1.96
	gen ci_`var'_hi_ = th_`var'_+se_`var'_*1.96
	}
	 
	twoway line th_ke_ t , lwidth(thick) msize(large) msymbol(diamond) lcolor(black)  xlabel(-3(1)3)|| line ci_ke_lo_ ci_ke_hi t , mcolor( red red) lwidth(thick thick)  lcolor( red red)  lpattern( shortdash shortdash)    ///
	xtitle("Time since mining engineer hire") ytitle("Locomotive adoption difference") legend(order(1 "Diff. in elec. adoption rate" 2 "90% CI"))   graphregion(color(white))
*	graph export "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\fig\figure6a.png", replace

	twoway line th_ks_ t , lwidth(thick) msize(large) msymbol(square) lcolor(black) lpattern(solid) mcolor(maroon)  xlabel(-3(1)3)|| line ci_ks_lo_ ci_ks_hi t , mcolor( red red) lwidth(thick thick)  lcolor( red red) lpattern( shortdash shortdash)    	///
	xtitle("Time since mining engineer hire") ytitle("Locomotive adoption difference") legend(order(1 "Diff. in steam adoption rate" 2 "90% CI")) graphregion(color(white))
	* graph export "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\fig\figure6b.png", replace

	twoway line th_ka_ t , lwidth(thick) msize(large) lcolor(black)  lpattern(solid) xlabel(-3(1)3) || line ci_ka_lo_ ci_ka_hi t , mcolor( red red) lwidth(thick thick)  lcolor( red red) lpattern( shortdash shortdash)    ///
	xtitle("Time since mining engineer hire") ytitle("Locomotive adoption difference") legend(order(1 "Diff. in air adoption rate" 2 "90% CI")) graphregion(color(white))
	* graph export "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\fig\figure6c.png", replace

	twoway line th_l_ t , lwidth(thick) msize(large) lcolor(black) lpattern(solid) xlabel(-3(1)3) || line ci_l_lo_ ci_l_hi t , mcolor( red red) lwidth(thick thick)  lcolor( red red) lpattern( shortdash shortdash)    ///
	xtitle("Time since mining engineer hire") ytitle("Labor adoption difference") legend(order(1 "Diff. in air adoption rate" 2 "90% CI")) graphregion(color(white))
	*graph export "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\fig\figure6d.png", replace
	
	restore
	 
	
	
/* -----------------------------------------------------------------------------
F. Information spillovers (Table 3(b) and 3(c)
-------------------------------------------------------------------------------- */	 
	  
	foreach var of varlist ke ks ka {
	bys firmid countyid: egen n`var'fc = sum(`var')		//  number of locomotives at  firms in same county
	bys firmid: egen n`var'f = sum(`var')		// number of locomotives at firm
	gen n`var'fot = n`var'f-n`var'fc				// number of locomotive at firm in other counties
	gen d`var'fot = n`var'fot>0						// firm has locomotives in other counties
	sum d`var'fot
	}
	
	gen xdkefot = x*dkefot
	foreach var of varlist ke ks ka {
	xtreg `var'  x  col  i.yr l m ltfp if d`var'fot ==1, fe r 		// locomotive choice regression if no locomotives in other counties
	local N_`var'ot1 = e(N) 
	local r2_`var'ot1 = round(e(r2_w),0.001)
	areg `var'  x   col  i.yr l m ltfp if d`var'fot ==1, absorb(mineid)  
	gen bx1`var' = _b[x]
	gen sex1`var' = _se[x]
	xtreg `var'  x   col  i.yr l m ltfp if d`var'fot ==0, fe r		// locomotive choice regression if already locomotives in other counties
	local N_`var'ot0 = e(N) 
	local r2_`var'ot0 = round(e(r2_w),0.001)
	areg `var'  x  col  i.yr l m ltfp if d`var'fot ==0, absorb(mineid)  
	gen bx0`var' = _b[x]
	gen sex0`var' = _se[x]
	}
   
	// store the estimates
	
	foreach var of varlist ke ka ks {
	gen bx0 = bx0`var'
	gen bx1 = bx1`var'
	estpost su bx0
	est store E_`var' 
	estpost su  bx1
	est store F_`var' 
	drop bx0 bx1
	rename (  sex0`var' sex1`var') ( bx0 bx1)
	estpost su bx0  
	est store E_`var'_se 
	estpost su  bx1
	est store F_`var'_se 
	drop   bx0 bx1
	}
	gen bx0 = .
	gen bx1 = .
	label var bx0 "1(Mining col. grad.)  "
	label var bx1 "1(Mining col. grad.)   "
	label var bkxe "1(M.C. grad)*1(Elec. loc.)"	

	*Write table 3
	/*
	esttab  D_kxe D_kxe_se   using "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\tab\tab3.tex", replace ///
	mtitle("" "" "" "" "" "") prehead(      &  \multicolumn{2}{c}{(I)}& \multicolumn{2}{c}{(II)} & \multicolumn{2}{c}{(III)} \\ \textit{(a) Different returns} & Estimate & S.E. & &&&\\ \hline) ///
	posthead(   \\)   cells(mean(fmt(3))) label booktabs nonum noobs collabels(none) gaps f   ///
	prefoot(  &&&\\   Observations & \multicolumn{2}{c}{`N_inter'} &&  \\  &&&\\   ) 		
	*/
	
	/*
	esttab  E_ke E_ke_se E_ka E_ka_se E_ks E_ks_se using "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\tab\tab3.tex", append ///
	mtitle("" "" "" "" "" "") prehead( \hline  &  \multicolumn{2}{c}{1(Elec. loc.)}& \multicolumn{2}{c}{1(Air. loc.)} & \multicolumn{2}{c}{1(Steam loc.)} \\  \textit{(b) Loc. not yet used  } & Estimate & S.E. & Estimate & S.E. & Estimate & S.E.\\ \hline) ///
	posthead(   \\)   cells(mean(fmt(3))) label booktabs nonum noobs collabels(none) gaps f   ///
	prefoot(  &&&\\   Observations & \multicolumn{2}{c}{`N_keot0'}  & \multicolumn{2}{c}{`N_kaot0'} & \multicolumn{2}{c}{`N_ksot0'} \\  Within  R-squared &  \multicolumn{2}{c}{`r2_keot0'}&\multicolumn{2}{c}{`r2_kaot0'}&\multicolumn{2}{c}{`r2_ksot0'}\\ &&&\\   ) 
	*/
	/*
	esttab  F_ke F_ke_se F_ka F_ka_se F_ks F_ks_se using "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\tab\tab3.tex", append ///
	mtitle("" "" "" "" "" "") prehead( \hline   &  \multicolumn{2}{c}{1(Elec. loc.)}& \multicolumn{2}{c}{1(Air. loc.)} & \multicolumn{2}{c}{1(Steam loc.)} \\ \textit{(c) Loc. already used  }  & Estimate & S.E. & Estimate & S.E. & Estimate & S.E.\\ \hline) ///
	posthead(   \\)   cells(mean(fmt(3))) label booktabs nonum noobs collabels(none) gaps f   ///
	prefoot(  &&&\\   Observations & \multicolumn{2}{c}{`N_keot1'}  & \multicolumn{2}{c}{`N_kaot1'} & \multicolumn{2}{c}{`N_ksot1'} \\  Within R-squared &  \multicolumn{2}{c}{`r2_keot1'}&\multicolumn{2}{c}{`r2_kaot1'}&\multicolumn{2}{c}{`r2_ksot1'}\\ &&&\\   ) 
	*/
 
