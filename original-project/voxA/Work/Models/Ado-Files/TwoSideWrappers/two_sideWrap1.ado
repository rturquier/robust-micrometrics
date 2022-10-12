*! two_sideWrap1 v01a  MOW 2013-07-30
program two_sideWrap1, rclass
	version 12.1
	syntax varlist(min=6 max=6 numeric)
	
* Wrapper ado-file for "two_side.ado" (Alan, Honoré, Hu and Leth-Petersen 2011) that:
*
*	(1) calls "two_side"
*	(2) estimates the constant term ("b0")
*	(3) estimates standard deviation of  u_i + e_it ("Sigma")
*	(4) calculates crowd-out-marginal-effects (by calling "margCrowdOutTwoSide1.ado")
*
*	As written, this wrapper is hard-coded for K=1 variable in this specific order: "Ggov".

*** Assign variables
*
	tempvar h Ggov  	Lower Upper identifier Budget

	clonevar `h'				= `1'
	clonevar `Ggov'				= `2'
	clonevar `Lower'			= `3'
	clonevar `Upper'			= `4'
	clonevar `identifier'		= `5'
	clonevar `Budget'			= `6'
	
*** Get the estimates
*
two_side `1' `2' `3' `4' `5' 	, details vc

	
*** Estimate the constant term
*
	tempname		Xbeta yMean beta0
	
	matrix			b = e(b)
*	matrix list		b
	
	matrix accum	Cov = `Ggov', noconstant means(M) deviations
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
	tempname	b1		b0 K Nobs Sigma2 Sigma
	
		scalar		`b1' = b[1,1]
		scalar		`b0' = b[1,2]
		
		scalar		`K' = colsof(b)

	tempvar			xbeta uhat2
	
	quietly gen double 		`xbeta' = `b1'*`Ggov' 			+ `b0'
	quietly gen double 		`uhat2' = (`h' - `xbeta')^2
	quietly replace			`uhat2' = .				if `h' == `Lower' | `h' == `Upper'
	
	quietly sum `uhat2'
		scalar	`Nobs' = r(N)
		scalar	`Sigma2' = r(sum)/(`Nobs' - `K')
		scalar	`Sigma' = sqrt(`Sigma2')

*	di as result "Sigma = " `Sigma'
	
		
*** Marginal effects - crowd-out
*
	tempname G1 Lo1 Up1 G2 Lo2 Up2
	
	sum `Budget'
		local budget1 = r(max)
		local budget2 = r(min)

		
	sum `Ggov'	if `Budget' == `budget1'
		assert r(sd) == 0
		scalar `G1' = r(min)

	sum `Ggov'	if `Budget' == `budget2'
		assert r(sd) == 0
		scalar `G2' = r(min)

		
		
	sum `Lower'	if `Budget' == `budget1'
		assert r(sd) == 0
		scalar `Lo1' = r(min)

	sum `Lower'	if `Budget' == `budget2'
		assert r(sd) == 0
		scalar `Lo2' = r(min)
				
		
		
	sum `Upper'	if `Budget' == `budget1'
		assert r(sd) == 0
		scalar `Up1' = r(min)

	sum `Upper'	if `Budget' == `budget2'
		assert r(sd) == 0
		scalar `Up2' = r(min)
				



	foreach x in G1 Lo1 Up1 G2 Lo2 Up2 {
		di as result "`x' = " ``x''
	}
	
	
		
	margCrowdOutTwoSide1 `b1'  `b0' `Sigma'	`G1' `Lo1' `Up1' `G2' `Lo2' `Up2'

	return list
	

	
return scalar me		= r(me)
	
end
