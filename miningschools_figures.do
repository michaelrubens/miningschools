/* Management, Productivity, and Technology Choices: Evidence from U.S. Mining Schools

								-	FIGURES 1 and 2  -

								Michael Rubens (KU Leuven)
===================================================================================================*/

	set more off
	
** Figure 1: Mining manager education


	preserve
	collapse minsch col ngrad, by(yr man_id)
	collapse minsch col ngrad, by(yr)
	label var ngrad "# E.M. graduates"
	label var minsch "E.M. degree"
	label var col "Other college degree"
	twoway bar ngrad yr , yscale(range(0,400) ) lcolor(black) fcolor(white)  || line minsch yr ,xlabel(1898(2)1914) msymbol(square) lcolor(black) ///
	lwidth(thick) msize(large) lpattern(shortdash) yaxis(2)/*legend(size(medlarge))*/ || line col yr, yaxis(2)xlabel(1898(2)1914) lpattern(longdash) lcolor(black)  lwidth(thick)  ytitle("# Mining college graduates ", ///
	  axis(1) margin(large) )    ytitle("Share of mines", axis(2))graphregion(color(white))

  	twoway  connect minsch yr ,xlabel(1898(2)1914) msymbol(square) lcolor(black) mcolor(black) mlwidth(medthick) mfcolor(white) ///
	lwidth(thick) msize(medlarge) lpattern(solid) yaxis(2)/*legend(size(medlarge))*/ ///
	|| connect col yr, yaxis(2)xlabel(1898(2)1914) lpattern(longdash) lcolor(black) msymbol(diamond) msize(medlarge) mcolor(black) mlwidth(medthick) mfcolor(white)  lwidth(thick)   ///
      ytitle("Share of mines", axis(2))graphregion(color(white))
*	graph export .\Paper\fig\figure1.png , replace
	restore	
  
 
	
** Figure 2: Usage rates of mining locomotives

	foreach var of varlist locel locst locair loc {
	gen d`var' = `var'>0
	}

	bys yr: egen qtot = sum(q)
	foreach var of varlist dlocel dlocst dlocair  {
	bys yr: egen pct`var' = mean(`var')
	gen q`var'=q*`var'
	bys yr: egen q`var'tot = sum(q`var')
	gen pctq`var' = q`var'tot/qtot
	}

	foreach var of varlist locel locst locair  {
	bys yr: egen n`var'tot = sum(`var')
	}
	
	label var nloceltot "Electrical"
	label var nlocairtot "Compressed air"
	label var nlocsttot "Steam"
	
	twoway connect nlocairtot  nlocsttot nloceltot yr if yr>1899, lwidth(thick thick thick) mcolor(black black black) mlwidth(medthick medthick medthick) mfcolor(white white white) ///
	lpattern(longdash shortdash solid) lcolor(black black black) msize(medlarge medlarge medlarge) msymbol( diamond triangle  square ) graphregion(color(white)) ytitle("Number of locomotives used") 
*	graph export .\Paper\fig\figure4.png, replace
	 
 
	 
	
