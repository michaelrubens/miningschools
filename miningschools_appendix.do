/* Management, Productivity, and Technology Choices: Evidence from U.S. Mining Schools

							-	RESULTS IN THE APPENDIX  -

								Michael Rubens (KU Leuven)
===================================================================================================*/

	use ./appendix/miningschools_data, clear
	set more off
	
	 
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
	
	
	/* TABLE A4 (I). Auto-regressive approach, Cobb-Douglas
	-------------------------------------------------------*/	

	gmm (y - {rho}*L.y - ({bl})*(l-{rho}*L.l) - ({bm})*(m-{rho}*L.m) - ({bke})*(ke-{rho}*L.ke) - {bks}*(ks- {rho}*L.ks)  - {bka}*(ka - {rho}*L.ka) - {bx}*(x - {rho}*L.x) ///
	- {bcol}*(col - {rho}*L.col) - {bt}*(yr - {rho}*L.yr)- {c}*(1-{rho}) )   ,  	inst( ke ka ks x col L.ke L.ka L.ks L.x  L.m L.l  yr  )
 	local N_ar1 = e(N)

  
	mat coef_2a = e(b)
	mat cov_2a = e(V)
	
	mat list coef_2a

	gen bl_2a = coef_2a[1,2]
	gen bm_2a = coef_2a[1,3]
	gen bke_2a = coef_2a[1,4]
	gen bks_2a = coef_2a[1,5]
	gen bka_2a = coef_2a[1,6]
	gen bx_2a = coef_2a[1,7]
	gen bcol_2a = coef_2a[1,8]
	gen bt_2a = coef_2a[1,9]
	gen bc_2a = coef_2a[1,10]

	mat list cov_2a
	gen sel_2a = sqrt(cov_2a[2,2])
	gen sem_2a = sqrt(cov_2a[3,3])
	gen seke_2a = sqrt(cov_2a[4,4])
	gen seks_2a = sqrt(cov_2a[5,5])
	gen seka_2a = sqrt(cov_2a[6,6])
	gen sex_2a = sqrt(cov_2a[7,7])
	gen secol_2a = sqrt(cov_2a[8,8])
	gen set_2a = sqrt(cov_2a[9,9])
	gen sec_2a = sqrt(cov_2a[10,10])
	 
	* R-squared
	preserve 
	gen tfp = y - bl_2a*l-bm_2a*m - bx_2a*x - bcol_2a*col - bke_2a*ke - bks_2a*ks - bka_2a*ka - bt_2a*yr-bc_2a
	gen yhat = -(- bl_2a*l-bm_2a*m - bx_2a*x  - bke_2a*ke - bcol_2a*col - bks_2a*ks - bka_2a*ka - bt_2a*yr-bc_2a)
	xtset mineid yr
	keep if tfp~=. & L.tfp~=.
	reg y yhat 
	gen r2 = e(r2)	// 0.78
	sum r2 
 	local r2_ar1 = round(r2,0.001)
	gen r2_ar1=substr("`r2_ar1'",1,4)
	local r2_ar1= r2_ar1 
	restore
	 
	 
	/* TABLE A4 (III). Auto-regressive approach, with interaction effects between L, M and locomotives
	----------------------------------------------------------------------------------------------------*/	

	* generate interaction terms
	foreach var of varlist l m x {
	foreach loc of varlist ke ks ka {
	gen `var'_`loc' = `var'*`loc'
	}
	}
	
	xtset mineid yr
	
	gmm (y - {rho}*L.y - ({bl})*(l-{rho}*L.l) - ({bm})*(m- {rho}*L.m) - ({blocel})*(ke-{rho}*L.ke) - {blocst}*(ks- {rho}*L.ks) 	- {blocair}*(ka - {rho}*L.ka) - {bx}*(x - {rho}*L.x) ///
	- ({blke})*(l_ke-{rho}*L.l_ke) - ({blka})*(l_ka-{rho}*L.l_ka)- ({blks})*(l_ks-{rho}*L.l_ks) ///
	- ({bmke})*(m_ke-{rho}*L.m_ke) - ({bmka})*(m_ka-{rho}*L.m_ka)- ({bmks})*(m_ks-{rho}*L.m_ks)- {bt}*(yr - {rho}*L.yr)- {bcol}*(col - {rho}*L.col)- {c}*(1-{rho}) ) ///
	,  	inst( ke ka ks x L.ke  L.ka  L.ks L.x L.m L.l L.l_ke L.l_ka L.l_ks L.m_ke L.m_ka L.m_ks yr col L.col )

	local N_int = e(N)
	mat coef_2b = e(b)
	mat cov_2b = e(V)
	mat list coef_2b

	gen bl_2b = coef_2b[1,2]
	gen bm_2b = coef_2b[1,3]
	gen bke_2b = coef_2b[1,4]
	gen bks_2b = coef_2b[1,5]
	gen bka_2b = coef_2b[1,6]
	gen bx_2b = coef_2b[1,7]
	gen blke_2b = coef_2b[1,8]
	gen blka_2b = coef_2b[1,9]
	gen blks_2b = coef_2b[1,10]
	gen bmke_2b = coef_2b[1,11]
	gen bmka_2b = coef_2b[1,12]
	gen bmks_2b = coef_2b[1,13]
	gen bt_2b = coef_2b[1,14]
	gen bcol_2b = coef_2b[1,15]
	gen bc_2b = coef_2b[1,16]
	  
	mat list cov_2b
	gen sel_2b = sqrt(cov_2b[2,2])
	gen sem_2b = sqrt(cov_2b[3,3])
	gen seke_2b = sqrt(cov_2b[4,4])
	gen seks_2b = sqrt(cov_2b[5,5])
	gen seka_2b = sqrt(cov_2b[6,6])
	gen sex_2b = sqrt(cov_2b[7,7])
	gen selke_2b = sqrt(cov_2b[8,8])
	gen selka_2b = sqrt(cov_2b[9,9])
	gen selks_2b = sqrt(cov_2b[10,10])
	gen semke_2b = sqrt(cov_2b[11,11])
	gen semka_2b = sqrt(cov_2b[12,12])
	gen semks_2b = sqrt(cov_2b[13,13])
	gen set_2b = sqrt(cov_2b[14,14])
	gen secol_2b = sqrt(cov_2b[15,15])
	gen sec_2b = sqrt(cov_2b[16,16])
	  
	* R-squared
	preserve 
	gen tfp = y - bcol_2b*col - bl_2b*l-bm_2b*m - bx_2b*x  - bke_2b*ke - bks_2b*ks - bka_2b*ka - blke_2b*l_ke  - blka_2b*l_ka  - blks_2b*l_ks - bmke_2b*m_ke  - bmka_2b*m_ka  - bmks_2b*m_ks   - bt_2b*yr-bc_2b
	xtset mineid yr
	keep if tfp~=. & L.tfp~=.
	egen avy = mean(y)
	egen tss = sum((y-avy)^2)
	egen rss = sum(tfp^2)
	gen r2 = 1 - rss/tss	// 0.78
	sum r2 
 	local r2_int = round(r2,0.001)
	gen r2_int=substr("`r2_int'",1,4)
	local r2_int= r2_int 
	restore
	 
	** Save production function estimates to table
	
		* table 2a - estimates
	*drop bl bm bke bka bks bt bx blka blke blks bmka bmke bmks
	 foreach var in "l" "m" "ka" "ke" "ks" "x" "t" "col"{
	rename  b`var'_2a b`var'
	}
	estpost su  bx bcol bl bm bke bka bks bt 
	est store A_2a
	drop bx bl bm bke bka bks bt bcol

	* table 2a - SEs
	foreach var in "l" "m" "ka" "ke" "ks" "x" "t" "col"{
	rename  se`var'_2a b`var'
	}
	estpost su  bx bcol bl bm bke bka bks bt 
	est store A_2a_se
	drop bx bl bm bke bka bks bt bcol

	
	* table 2b - estimates
	foreach var in "l" "m" "ka" "ke" "ks" "lka" "lke" "lks" "mka" "mke" "mks" "x" "t" "col"{
	rename  b`var'_2b b`var'
	}
	estpost su  bx bcol bl bm bke bka bks  blka blke blks bmka bmke bmks  bt 
	est store A_2b
	drop bx bl bm bke bka bks bt blka blke blks bmka bmks bmke  bcol

	* table 2b - SEs
	foreach var in "l" "m" "ka" "ke" "ks" "lka" "lke" "lks" "mka" "mke" "mks" "x" "t" "col"{
	rename  se`var'_2b b`var'
	}
	estpost su  bx  bcol bl bm bke bka bks  blka blke blks bmka bmke bmks bt
	est store A_2b_se
  
	
 	/* TABLE A4 (II): ACF(2015) estimator, Cobb-Douglas model 
	----------------------------------------------------------*/	

	save ./appendix/temp_ap/data_temp2, replace	

	do ./appendix/miningschools_acf		// acf estimates
//	do ./appendix/miningschools_acf_bs	// standard errors (bootstrap)	
	 
	use ./appendix/temp_ap/data_temp2, clear
	
	gen const = 1
	merge m:1 const using ./appendix/temp_ap/estimates_acf, nogen // import estimates
	merge m:1 const using ./appendix/temp_ap/estimates_acf_se, nogen	// import SEs
	
	preserve 
	gen tfp = y - bl_acf*l-bm_acf*m - bx_acf*x  - bke_acf*ke - bks_acf*ks - bka_acf*ka  - bt_acf*yr  -bc_acf -bcol_acf*col
	gen yhat = -( bl_acf*l-bm_acf*m - bx_acf*x  - bke_acf*ke - bks_acf*ks - bka_acf*ka  - bt_acf*yr  -bc_acf -bcol_acf*col)
	xtset mineid yr
	drop if tfp==. | L.tfp==.
	reg yhat y 
	gen r2 = e(r2)	 
	sum r2 
	local r2_acf = round(r2,0.001)
	gen r2_acf=substr("`r2_acf'",1,4)
	local r2_acf= r2_acf 
	restore
 
	
	sum *_acf
	 
	gen t = yr
	gen c = 1
	drop bx bl bm bke bka bks bt  bcol
	foreach var of varlist l m ke ks ka x t c  col{
	rename b`var'_acf b`var'
	}	

	estpost su  bx  bcol  bl bm bke bka bks  bt 
	est store A_acf
	 	
 	* R-squared
   	local N_acf = `N_ar1'

	sum bl bm bx bke bka bks bt bc bcol
	
  
	foreach var of varlist l m ke ks ka x t  col {
	rename b`var' b`var'_acf
	rename se`var' b`var'
	}	
	estpost su  bx bcol  bl bm  bke bka bks  bt 
	est store A_acf_se
	
	label var bl "log(Labor)"
	label var bm "log(Materials)"
	label var bke "1(Elec. loc.)"
	label var bka "1(Air loc.)"
	label var bks "1(Steam loc.)"
	label var blke "1(Elec. loc.)*log(Labor)"
	label var blka "1(Air loc.)*log(Labor)"
	label var blks "1(Steam loc.)*log(Labor)"
	label var bmke "1(Elec. loc.)*log(Materials)"
	label var bmka "1(Air loc.)*log(Materials)"
	label var bmks "1(Steam loc.)*log(Materials)"
	label var bx "1(Mining col. grad.)" 
	label var bcol "1(Other grad.)" 
	label var bt "Year"
	
//	esttab  A_2a A_2a_se A_acf A_acf_se  A_2b A_2b_se  using "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\tab\tab2.tex", replace ///
//	mtitle("" "" "" "" "" "" ) prehead( &  \multicolumn{2}{c}{ (I)}& \multicolumn{2}{c}{ (II)} & \multicolumn{2}{c}{ (III)} \\    & \multicolumn{2}{c}{log(Output)}& \multicolumn{2}{c}{log(Output)} & \multicolumn{2}{c}{log(Output)} \\ & Estimate & S.E. & Estimate & S.E.& Estimate & S.E.\\ \hline)  posthead(   \\)    ///
//	cells(mean(fmt(3))) label booktabs nonum noobs collabels(none) gaps f    prefoot(  &&&\\  Model & \multicolumn{2}{c}{Cobb-Douglas}  & \multicolumn{2}{c}{Cobb-Douglas}& \multicolumn{2}{c}{Interaction effects} ///
//	\\  Method & \multicolumn{2}{c}{Autoregressive} & \multicolumn{2}{c}{ACF(2015)} & \multicolumn{2}{c}{Autoregressive}\\ ///
//	Observations & \multicolumn{2}{c}{`N_ar1'} & \multicolumn{2}{c}{`N_acf'}  & \multicolumn{2}{c}{`N_int'}  \\  R-squared &  \multicolumn{2}{c}{`r2_ar1'}&\multicolumn{2}{c}{`r2_acf'}&\multicolumn{2}{c}{`r2_int'}\\ &&&\\   ) 
//	drop bx bl bm bke bka bks bt blka blke blks bmka bmks bmke 
 
	 /* Table A5: Technology choice robustness checks 
	---------------------------------------------------*/	

 /* Common fixed cost component */
 
	use ./appendix/temp_ap/data_temp2, clear
	gen const = 1
	merge m:1 const using ./appendix/temp_ap/estimates_acf, nogen
 	  
	egen bl_cs = median((w*emp*ndays) / (p*q))
	egen bm_cs = median(qinp /q)
 
	gen ltfp_cs = y - bl_cs  *l  - bm_cs  *m		
 
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
 
	gen ltfp = lq - bl_cs*l - bm_cs*m - bke_cs_3 *ke-bks_cs_3 *ks -bka_cs_3 *ka -bx_cs_3 *x -bt_cs_3  *yr -bc_cs_3 -bcol_cs_3 *col
 
	gen tfp=exp(ltfp)

	xtset mineid yr

	label var ke "Electrical"
	label var ks "Steam"
	label var ka "Air"
	
	* Panel a: Extensive margin
	
	xtreg  ke ks ka x col  i.yr ltfp  l m, fe r
	local r2_kebis = round(e(r2_w),0.001)
	gen r2_kebis=substr("`r2_kebis'",1,4)
	local r2_kebis = r2_kebis
	local N_kebis = e(N)
	areg  ke ks ka x col  i.yr ltfp  l m, absorb(mineid) cluster (man_id)
	gen cxke = _b[x]
 	gen coke = _b[col]
	gen cxke_se = _se[x]
 	gen coke_se = _se[col]

	xtreg  ka ks ke x  col i.yr ltfp  l m, fe r
	local r2_kabis = round(e(r2_w),0.001)
	gen r2_kabis=substr("`r2_kabis'",1,4)
	local r2_kabis = r2_kabis
	local N_kabis = e(N)
	areg  ka ks ke x  col i.yr ltfp  l m, absorb(mineid) cluster (man_id)
	gen cxka = _b[x]
 	gen coka = _b[col]
	gen cxka_se = _se[x]
 	gen coka_se = _se[col]

 
	set more off
	xtreg  ks ka ke x  col i.yr ltfp  l m, fe r
	local N_ksbis = e(N)
	local r2_ksbis = round(e(r2_w),0.001)
	gen r2_ksbis=substr("`r2_ksbis'",1,4)
	local r2_ksbis = r2_ksbis

	areg  ks ka ke x col  i.yr ltfp  l m, absorb(mineid) cluster (man_id)
	gen cxks = _b[x]
 	gen coks = _b[col]
	gen cxks_se = _se[x]
 	gen coks_se = _se[col]
 
	estpost su  cxke  coke
	est store Cap_ke 

	foreach var of varlist ka ks {
	replace cxke = cx`var'
	replace coke = co`var'
	estpost su  cxke  coke
	est store Cap_`var' 
	}

	replace cxke = cxke_se
 	replace coke = coke_se
	estpost su  cxke  coke
	est store Cap_ke_se 

	foreach var of varlist ka ks {
	replace cxke = cx`var'_se
 	replace coke = co`var'_se
	estpost su  cxke  coke
	est store Cap_`var'_se
	}
	
	sum ke ka ks
	
	label var cxke "1(Mining col. grad.) "
 	label var coke "1(Other   grad.) "
//	esttab  Cap_ke Cap_ke_se Cap_ka Cap_ka_se Cap_ks Cap_ks_se using "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\tab\tab2_ap.tex", replace ///
//	mtitle("" "" "" "" "" "") prehead(   \textit{(a) Common fixed costs} &  \multicolumn{2}{c}{1(Elec. loc.)}& \multicolumn{2}{c}{1(Air. loc.)} & \multicolumn{2}{c}{1(Steam loc.)} \\ & Estimate & S.E. & Estimate & S.E. & Estimate & S.E.\\ \hline) ///
//	posthead(   \\)   cells(mean(fmt(3))) label booktabs nonum noobs collabels(none) gaps f   ///
//	prefoot(  &&&\\   Observations & \multicolumn{2}{c}{`N_kebis'}  & \multicolumn{2}{c}{`N_kabis'} & \multicolumn{2}{c}{`N_ksbis'} ///
//	\\  Within R-squared &  \multicolumn{2}{c}{`r2_kebis'}&\multicolumn{2}{c}{`r2_kabis'}&\multicolumn{2}{c}{`r2_ksbis'}\\ &&&\\   ) 
  	 
	* Logit model 
	
	encode postoffice, gen(townid)
 
	drop cxke cxka cxks    coke coka coks
	
	set more off
	foreach var of varlist ke ka ks {
	logit   `var' x col  yr l ltfp m  , cluster(man_id)
	local N_`var'log = e(N)
	margins, dydx(x col)  
	mat cx`var' = r(b)
	gen cx`var'=cx`var'[1,1]
	mat secx`var' = r(V)
	gen secx`var' = sqrt(secx`var'[1,1])
	mat co`var' = r(b)
	gen co`var'=co`var'[1,2]
	mat seco`var' = r(V)
	gen seco`var' = sqrt(seco`var'[2,2])
	rename (cx`var' co`var') (cx co)
	estpost su cx co
	est store Cap2_`var'
	drop cx co
	rename (secx`var' seco`var') (cx co)
	estpost su cx co
	est store Cap2_`var'_se
	drop cx co
	}
	
	gen cx = 1	
	gen co = 1
	label var cx "1(Mining col. grad.) "
	label var co "1(Other grad.)"
//	esttab  Cap2_ke Cap2_ke_se Cap2_ka Cap2_ka_se Cap2_ks Cap2_ks_se using "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\tab\tab2_ap.tex", append ///
//	mtitle("" "" "" "" "" "") prehead(  \hline \textit{(b) Logit model} &  \multicolumn{2}{c}{1(Elec. loc.)}& \multicolumn{2}{c}{1(Air. loc.)} & \multicolumn{2}{c}{1(Steam loc.)} \\ & Estimate & S.E. & Estimate & S.E. & Estimate & S.E.\\ \hline) ///
//	posthead(   \\)   cells(mean(fmt(3))) label booktabs nonum noobs collabels(none) gaps f   ///
//	prefoot(  &&&  \\  Observations & \multicolumn{2}{c}{`N_kelog'}  & \multicolumn{2}{c}{`N_kalog'} & \multicolumn{2}{c}{`N_kslog'} ///
//	\\ &&&\\   ) 
			  
	 * Cost dynamics
	 	 
	*drop qcum lqcum
	gen qcum = 0
	xtset mineid yr
	forvalues n = 1/30 {
	replace qcum = qcum+L`n'.q if L`n'.q~= .
	}
	
	gen lqcum=ln(qcum)
	label var lqcum "ln(cumulative Q)"
	gen llp=ln(q/emp)
	label var llp "ln(Q/L)"

	drop cxke cxka cxks    cx  co

	foreach var of varlist  ks ke ka {
	xtreg `var' x col l  m ltfp lqcum  i.yr, fe r
	local r2_`var'cum = round(e(r2_w),0.001)
	local N_`var'cum = e(N)
	areg `var' x col l  m ltfp lqcum i.yr, absorb(mineid) r			
	gen cx`var' = _b[x]
	gen secx`var' = _se[x]
	gen co`var' = _b[col]
	gen seco`var' = _se[col]
	gen cc`var' = _b[lqcum]
	gen secc`var' = _se[lqcum]
	rename (cx`var'   cc`var' co`var') (cx  cc co)
	estpost su cx  cc co
	est store Cap3_`var'
	drop cx  cc co
	rename (secx`var'   secc`var' seco`var')(cx   cc co)
	estpost su cx   cc co
	est store Cap3_`var'_se
	drop cx cc co
	}
 
 
	gen cx = .	
	gen cc = .
	gen co = .


	label var cx "1(Mining col. grad.) "
	label var co "1(Other grad.) "
	label var cc "log(Cum. output) "
//	esttab  Cap3_ke Cap3_ke_se Cap3_ka Cap3_ka_se Cap3_ks Cap3_ks_se using "C:\Users\MichaelRubens\Dropbox\Managerial education\Paper\tab\tab2_ap.tex", append ///
//	mtitle("" "" "" "" "" "") prehead(  \hline \textit{(c) Cost dynamics} &  \multicolumn{2}{c}{1(Elec. loc.)}& \multicolumn{2}{c}{1(Air. loc.)} & \multicolumn{2}{c}{1(Steam loc.)} \\ & Estimate & S.E. & Estimate & S.E. & Estimate & S.E.\\ \hline) ///
//	posthead(   \\)   cells(mean(fmt(3))) label booktabs nonum noobs collabels(none) gaps f   ///
//	prefoot(  &&&  \\  Observations & \multicolumn{2}{c}{`N_kecum'}  & \multicolumn{2}{c}{`N_kacum'} & \multicolumn{2}{c}{`N_kscum'} ///
//	\\  Within R-squared & \multicolumn{2}{c}{`r2_kecum'}  & \multicolumn{2}{c}{`r2_kacum'} & \multicolumn{2}{c}{`r2_kscum'} \\ &&&\\   ) 
	
	 
		
	
