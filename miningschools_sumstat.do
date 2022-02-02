/* Management, Productivity, and Technology Choices: Evidence from U.S. Mining Schools

								-	SUMMARY STATISTICS   -

								Michael Rubens (KU Leuven)
===================================================================================================*/
	
	
	tab minsch				//308 observations for MC graduates
	tab col					//45 observations for OC graduates
	
	preserve
	collapse(sum) minsch col, by(man_id)
	gen dum_minsch = minsch>0
	gen dum_col = col>0
	gen dum_oth = minsch==0 & col==0
	tab dum_minsch			
	tab dum_col				// 17 college graduates, among which 7 MS graduates 
	tab dum_oth
	restore
	
	preserve
	collapse(sum) minsch col, by(mineid )
	gen dum_minsch = minsch>0
	gen dum_col = col>0
	tab dum_minsch			// 61 mines with mc grad, 17 with col grad
	tab dum_col
	restore
	
	preserve
	collapse(sum) minsch col, by(firmid )
	gen dum_minsch = minsch>0
	gen dum_col = col>0
	tab dum_minsch			// 15 firms with mc grad, 12 with col grad
	tab dum_col
	restore
	
* Multi-mine firms

	bys firmid  yr: egen nminef = sum(open)
	gen multimine = nminef>1
	sum multimine 								// 70% of mines owned by multi-mine firms
	
	preserve
	collapse nminef, by(firmid yr)
	gen multimine = nminef>1					
	sum multimine								// 22% of firms own more than one mine
	sum nminef, d								// 2.5 mines per firm
	restore
	
* Number of superintendents by firm and firm-county-year pair

	preserve
	collapse(sum) open, by(firmid  man_id yr)
	gen nman = 1
	collapse(sum) nman, by(firmid  yr)
	sum nman 									// 1.05 superintendent per firm-year, max. 4
	restore
	
	preserve
	keep if multimine == 1
	collapse(sum) open, by(firmid  man_id yr)
	gen nman = 1
	collapse(sum) nman, by(firmid   yr)			// In multi-mine firms, on average 1.21 superintendents
	sum nman 
	restore
 	
* Mine summary statistics - table A1 (a)

	sum q pctship emp powd qinp loc minsch col
	
* Manager summary statistics - table A1 (b)
	
	bys man_id minsch col yr: egen nmineman=sum(open)
	preserve
	collapse age q nmineman, by(man_id minsch col yr)
	sum age q 
	bys minsch col: sum age q nmineman
	restore

	
 * Market structure		

	* State-wide markets

	preserve
	collapse(sum) q , by(firm yr)
	bys yr : egen qyr = sum(q)
	gen ms = q/qyr
	sum ms, d
	sum ms
	sort  ms yr
	restore

	* County-level markets

	preserve
	collapse(sum) q , by(firm countyid yr)
	bys yr countyid : egen qyr = sum(q)
	gen ms = q/qyr
	sum ms, d
	sum ms
	restore

* Simultaneous usage of different locomotive types 

	gen ntype = dlocair+dlocst+dlocair
	replace ntype=. if yr<1900
	gen multitech = ntype>1
	sum multitech

	
	use ./main/miningschools_data, clear
 /*
* how many unique superintendents per year?	
	
	preserve
	collapse mineid , by(sup1 sup2 sup3 yr)
	gen id = _n
	reshape long sup, i(id) j(n)
	drop if sup==""
	quietly bys yr sup:  gen dup = cond(_N==1,0,_n)
	drop if dup>1
	gen nsupyr = 1
	collapse(sum) nsupyr, by(yr)  
	save ./data/tempfiles/nsupyr, replace
	restore
	
	merge m:1 yr using ./data/tempfiles/nsupyr, nogen	

	gen supcostyr = wsup*nsupyr	// annual labor cost of superintendents
	gen lcost = w*emp*ndays	// annual labor cost of production workers
	bys yr: egen lcostyr = sum(lcost)
	
	gen supshare = supcostyr/lcostyr	// 
	preserve
	collapse supshare, by(yr)
	drop if supshare==.
	sum supshare, d		 
	restore
*/
	
