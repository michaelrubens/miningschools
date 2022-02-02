/* Management, Productivity, and Technology Choices: Evidence from U.S. Mining Schools

					-	PRODUCTION FUNCTION ESTIMATION USING ACF (2015): ESTIMATES  -

								Michael Rubens (KU Leuven)
===================================================================================================*/

use ./appendix/temp_ap/data_temp2, clear
save ./appendix/temp_ap/data_ap_temp, replace		// temporary data file

set more off
	 

** First stage estimation 

	xtset mineid yr, year

	* Estimate probability of exit
	gen dx = q~=. & F.q==.
	probit dx  l m ka ke ks x  
	predict px 
	
	
	local I = 3									// 2nd & 3rd order terms for polynomial approximation 
	local J = 3
	local N = 3
	forvalues i=1/`I'{					
		gen l`i' = l^(`i')
		gen m`i' = m^(`i')
	 
			forvalues j=1/`J'{					
			gen l`i'm`j' = l^(`i')*m^(`j')			 
			}
	   }
 
	
	  // ke, ka, ks and x are dummies, so do not enter the polynomials with squared terms:
	   
	  foreach var of varlist l1* l2* l3* m1* m2* m3  {
	  gen `var'ke = `var'*ke
	  gen `var'ka = `var'*ka
	  gen `var'ks = `var'*ks
	  gen `var'x = `var'*x
	  gen `var'keka = `var'*ke*ka
	  gen `var'keks = `var'*ke*ks
	  gen `var'kex = `var'*ke*x
	  gen `var'kaks = `var'*ka*ks
	  gen `var'kax = `var'*ka*x
	  gen `var'ksx = `var'*ks*x
	  gen `var'kekaks = `var'*ke*ka*ks
	  gen `var'kekax = `var'*ke*ka*x
	  gen `var'keksx = `var'*ke*ks*x
	  gen `var'kaksx = `var'*ke*ks*x
	  gen `var'kekaksx = `var'*ke*ka*ks*x
		}	
	  gen keka = ke*ka
	  gen keks = ke*ks
	  gen kex = ke*x
	  gen kaks = ka*ks
	  gen kax = ka*x
	  gen ksx = ks*x
	  gen kekaks = ke*ka*ks
	  gen kekax = ke*ka*x
	  gen keksx = ke*ks*x
	  gen kaksx = ke*ks*x
	  gen kekaksx = ke*ka*ks*x
	   
	  gen xke = kex
	   
	* Estimate phi function
	set matsize 11000
	quietly{
	reg y l1* l2* l3* m1* m2* m3* x col ka* ke* ks*   i.yr  px
	predict double phi if e(sample)==1 
	gen phi_lag=L.phi if e(sample)==1
	}

	
**  Second Stage Estimation 

	gen l_lag = L.l 				// lagged variables
	gen ka_lag = L.ka 
	gen ke_lag = L.ke 
	gen ks_lag = L.ks 
	gen m_lag = L.m 
	gen l_lag2 = l_lag^2
	gen ka_lag2 = ka_lag^2
	gen ke_lag2 = ke_lag^2
	gen ks_lag2 = ks_lag^2
	gen m_lag2 = m_lag^2
	gen l_lagka_lag = l_lag*ka_lag
	gen l_lagke_lag = l_lag*ke_lag
	gen l_lagks_lag = l_lag*ks_lag
	gen lka = l*ka
	gen lke = l*ke
	gen lks = l*ks
	gen l_lagka1 = l_lag*ka
	gen l_lagke1 = l_lag*ke
	gen l_lagks1 = l_lag*ks
	gen col_lag=L.col
	gen x_lag = L.x
	gen yr_lag=L.yr
	gen xke_lag = L.xke

	gen l1m1_lag = l_lag*m_lag
		
	* Averages (for partial derivatives interaction terms)
	
	foreach var of varlist   l m ka ke ks {
	egen av`var' = mean(`var')
	}

	xtset mineid yr
	gen const = 1
	drop if l == .
	drop if ka == . | ke==. | ks==.
	drop if m == .
	drop if x==.
	drop if x_lag==.
	drop if phi == .
	drop if phi_lag == .
	
	di _N

	* GMM estimation
	
	do  ./appendix/miningschools_pf_acf		// run GMM estimation command in mata 

	* Store estimates of this bootstrap iteration in vectors
	
	gen bc_acf = beta_lin[1,1]
	gen bl_acf  = beta_lin[2,1]
	gen bx_acf = beta_lin[3,1]
	gen bka_acf = beta_lin[4,1]
	gen bke_acf = beta_lin[5,1]
	gen bks_acf = beta_lin[6,1]
	gen bt_acf = beta_lin[7,1]
	gen bcol_acf = beta_lin[8,1]
	gen bm_acf = beta_lin[9,1]

*	gen const=1	
	collapse bc_acf bl_acf bka_acf bke_acf bks_acf bm_acf bx_acf  bt_acf bcol_acf  , by(const)		// Save estimates and standard errors in tempfile
 
	save ./appendix/temp_ap/estimates_acf, replace

