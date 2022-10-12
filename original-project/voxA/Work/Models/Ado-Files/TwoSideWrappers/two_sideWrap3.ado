*! two_sideWrap3 v01a  MOW 2013-07-27
program two_sideWrap3, rclass
	version 12.1
	syntax varlist(min=7 max=7 numeric)
	
* Wrapper ado-file for "two_side.ado" (Alan, Honoré, Hu and Leth-Petersen 2011) that:
*
*	(1) calls "two_side"
*	(2) estimates the constant term ("b0")
*	(3) estimates standard deviation of  u_i + e_it ("Sigma")
*	(4) calculates crowd-out-marginal-effects (by calling "margCrowdOutTwoSide3.ado")
*
*	As written, this wrapper is hard-coded for K=3 variables in this specific order: "Ggov Ggov_GovtHigh GovtHigh".

*** Assign variables
*
	tempvar h Ggov Ggov_GovtHigh GovtHigh 	Lower Upper identifier

	clonevar `h'				= `1'
	clonevar `Ggov'				= `2'
	clonevar `Ggov_GovtHigh'	= `3'
	clonevar `GovtHigh'			= `4'
	clonevar `Lower'			= `5'
	clonevar `Upper'			= `6'
	clonevar `identifier'		= `7'

	
*** Get the estimates
*
two_side `0' 	, details vc

	
*** Estimate the constant term
*
	tempname		Xbeta yMean beta0
	
	matrix			b = e(b)
*	matrix list		b
	
	matrix accum	Cov = `Ggov' `Ggov_GovtHigh' `GovtHigh'	, noconstant means(M) deviations
*	matrix list		M
	
	matrix 			Xb = M * b'
*	matrix list		Xb
	
	scalar			`Xbeta' = Xb[1,1]
	
	sum				`h' , meanonly									
		scalar		`yMean' = r(mean)
		
	scalar			`beta0' = `yMean' - `Xbeta'
*	di as result	"Constant term is " `beta0'

	matrix			b0 = `beta0'
	matrix colnames	b0 = _cons
	
	matrix			bWithCon = b , b0					// Append the constant term to the other betas
	
	matrix			b = bWithCon						// Rename the beta vector as simply "b"
	matrix list		b

	
*** Estimate sigma 
*
	tempname	b1 b2 b3 b0 K Nobs Sigma2 Sigma
	
		scalar		`b1' = b[1,1]
		scalar		`b2' = b[1,2]
		scalar		`b3' = b[1,3]
		scalar		`b0' = b[1,4]
		
		scalar		`K' = colsof(b)

	tempvar			xbeta uhat2
	
	quietly gen double 		`xbeta' = `b1'*`Ggov' + `b2'*`Ggov_GovtHigh'  + `b3'*`GovtHigh' + `b0'
	quietly gen double 		`uhat2' = (`h' - `xbeta')^2
	quietly replace			`uhat2' = .				if `h' == `Lower' | `h' == `Upper'
	
	quietly sum `uhat2'
		scalar	`Nobs' = r(N)
		scalar	`Sigma2' = r(sum)/(`Nobs' - `K')
		scalar	`Sigma' = sqrt(`Sigma2')

	di as result "Sigma = " `Sigma'	
	
		
*** Marginal effects - crowd-out
*		
	margCrowdOutTwoSide3 `b1' `b2' `b3' `b0' `Sigma'

	return list
	

	
return scalar coLow		= r(coLow)
return scalar coHigh	= r(coHigh)
return scalar coChange	= r(coChange)
	
end
