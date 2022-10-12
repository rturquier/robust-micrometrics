*! two_sideWrap2hi v01a  MOW 2013-08-31
program two_sideWrap2hi, rclass
	version 12.1
	syntax varlist(min=6 max=6 numeric)
	
* Wrapper ado-file for "two_side.ado" (Alan, Honoré, Hu and Leth-Petersen 2011) that:
*
*	(1) calls "two_side"
*	(2) estimates the constant term ("b0")
*	(3) estimates standard deviation of  u_i + e_it ("Sigma")
*	(4) calculates crowd-out-marginal-effects (by calling "margCrowdOutTwoSide2hi.ado")
*
*	As written, this wrapper is hard-coded for K=2 variables in this specific order: "y Ggov".
*
*	Warning: As written, this wrapper's call to "margCrowdOutTwoSide2hi.ado" is hard-coded specifically for HIGH provision.
*																											^^^^
*																											^^^^
*																											^^^^


*** Assign variables
*
	tempvar h y Ggov 	Lower Upper identifier

	clonevar `h'				= `1'
	clonevar `y'				= `2'
	clonevar `Ggov'				= `3'
	clonevar `Lower'			= `4'
	clonevar `Upper'			= `5'
	clonevar `identifier'		= `6'

	
*** Get the estimates
*
two_side `0' 	, details vc

	
*** Estimate the constant term
*
	tempname		Xbeta yMean beta0
	
	matrix			b = e(b)
*	matrix list		b
	
	matrix accum	Cov = `y' `Ggov'	, noconstant means(M) deviations
*	matrix list		M
	
	matrix 			Xb = M * b'
*	matrix list		Xb
	
	scalar			`Xbeta' = Xb[1,1]
	
	sum				`h' , meanonly									
		scalar		`yMean' = r(mean)
		
	scalar			`beta0' = `yMean' - `Xbeta'
	di as result	"Constant term is " `beta0'

	matrix			b0 = `beta0'
	matrix colnames	b0 = _cons
	
	matrix			bWithCon = b , b0					// Append the constant term to the other betas
	
	matrix			b = bWithCon						// Rename the beta vector as simply "b"
	matrix list		b

	
*** Estimate sigma 
*
	tempname	b1 b2 b0 K Nobs Sigma2 Sigma
	
		scalar		`b1' = b[1,1]
		scalar		`b2' = b[1,2]
		scalar		`b0' = b[1,3]
		
		scalar		`K' = colsof(b)

	tempvar			xbeta uhat2
	
	quietly gen double 		`xbeta' = `b1'*`y' + `b2'*`Ggov'  + `b0'
	quietly gen double 		`uhat2' = (`h' - `xbeta')^2
	quietly replace			`uhat2' = .				if `h' == `Lower' | `h' == `Upper'
	
	quietly sum `uhat2'
		scalar	`Nobs' = r(N)
		scalar	`Sigma2' = r(sum)/(`Nobs' - `K')
		scalar	`Sigma' = sqrt(`Sigma2')

	di as result "Sigma = " `Sigma'	
	
		
*** Marginal effects - crowd-out
*		
	margCrowdOutTwoSide2hi `b1' `b2' `b0' `Sigma'

	return list
	
return scalar incHigh		= r(incHigh)
return scalar coHigh		= r(coHigh)
return scalar q1plusq2High	= r(q1plusq2High)
return scalar q2High		= r(q2High)
	
end
