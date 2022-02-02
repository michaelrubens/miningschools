/* 
===========================================================
ACF PRODUCTION FUNCTION ESTIMATION: MATA PROGRAM
Code adapted from De Loecker (AEJ Micro-2013) 
===========================================================
*/


/*********** Optimization by Mata****************/



set more off
mata:
mata clear
//OLS estimates used for initial values
	X_main=st_data(., ("const", "l",  "x","ka","ke","ks" , "yr", "col"))
	AVKA=st_data(., ("avka"))
	AVKE=st_data(., ("avke"))
	AVKS=st_data(., ("avks"))
	X_ext=st_data(., ("m"))
	X=(X_main, X_ext)
	Y=st_data(., ("y"))
	beta_init = invsym(X'X)*X'Y
	beta_init = beta_init'
	st_matrix("beta_init", beta_init)
	

	// bringing data from stata; st_view() may be better 
	PHI=st_data(., ("phi")) 
	PHI_lag=st_data(.,("phi_lag"))
	Z_main=st_data(., ("const" ,"l_lag", "x","ka","ke","ks", "yr", "col" ))
	X_main=st_data(., ("const" ,"l",   "x","ka","ke","ks", "yr", "col" ))
	X_ext=st_data(., ("m"))
	X_main_lag=st_data(., ("const","l_lag", "x_lag","ka_lag","ke_lag","ks_lag", "yr_lag" , "col_lag")) 
	X_ext_lag=st_data(., ("m_lag"))
    X=(X_main, X_ext)
	X_lag=(X_main_lag, X_ext_lag)
	Z=(Z_main, X_ext)
	Y=st_data(., ("y"))
	C=st_data(., ("const"))	
	MS=st_data(.,("x"))
	MS_lag=st_data(.,("x_lag"))
	M_lag=st_data(.,("m_lag"))
	
/* De Loecker method - Endogenous Productivity Change*/
// defining GMM_DL_linear function
void GMM_DL_linear(todo, betas, crit, g, H)
{
	external PHI,PHI_lag,Z_main,X_main,X_ext,X_main_lag,X,X_lag,Z,Y,C,MS,MS_lag,M_lag	
		
    // create varaibles to polynomial approximation of g function
	// this time, a set of the variables should include LB, EA, LA
	OMEGA = PHI - X*betas'
	OMEGA_lag = PHI_lag - X_lag*betas'
	OMEGA_lag2 = OMEGA_lag :* OMEGA_lag
	OMEGA_lag3 = OMEGA_lag2 :* OMEGA_lag
	OMEGA_lag_pool = (C, OMEGA_lag,  MS_lag)
	// estimation of g() function 
	g_b = invsym(OMEGA_lag_pool'OMEGA_lag_pool)*OMEGA_lag_pool'OMEGA
	// residual function for GMM
	XI = OMEGA - OMEGA_lag_pool * g_b
	// moment condition is Z'XI
	crit = (Z'XI)'(Z'XI)
} 

S = optimize_init()
optimize_init_evaluator(S, &GMM_DL_linear())
optimize_init_evaluatortype(S, "d0")
optimize_init_technique(S, "nm")
optimize_init_nmsimplexdeltas(S, 0.1)
optimize_init_which(S, "min")
optimize_init_params(S, beta_init)
beta_dl_linear = optimize(S)
j_dl_linear = optimize_result_value(S)
st_matrix("beta_lin", beta_dl_linear')
st_matrix("beta_ols",beta_init)
st_matrix("avka",AVKA)
st_matrix("avke",AVKE)
st_matrix("avks",AVKS)

end

matrix list beta_lin	// 

// output elasticities of variable inputs (partial derivatives)






