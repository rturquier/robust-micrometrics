//	file: 
		local Namedo "voxA-mod001C-v01a-RandomEffectsTobit"
		
//	task:	Estimate Random effects Tobit of the Balanced-budget crowd-out model.
//

//	author:
		local aut "\ mow \ 2016-10-04"   

//	Include file to assign directories:
		voxAdir Models
		include "voxA-mod-dirInclude-001A.do"
		
// Input datasets:
	local FileIn01	"voxA-cr021A_v01a-DecisionsLONG"		// h1-h6, Ggov1-Ggov6, ...

	
	
// Note: These are robustness check results described in footnote 15.

		
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (0) Get the data.
*
use	"`dirIN'/`FileIn01'", clear			




* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (1) "xttobit" RE-Tobit applied to the LO govt provision separately
*			from the HI govt provision.
*
*		Notes:	(1) "xttobit" RE-Tobit results are the random effects parallel 
*						to the FE-Tobit of Alan et al. (2011).
*
*				(2) Because the Tobit is non-linear, these two separate 
*						LO, HI results are not identical to the "LO and HI
*						estimated in one model" results obtained in (2)
*						just below.
*

*** (1.1) Balanced-budget crowd-out at Low
*
preserve

keep if  (Budget==2 | Budget==5)

	quietly xtreg	h Ggov, fe	vce(robust)		// fe . . . just to get the degrees of freedom for the t=test below.
				scalar dfr	= e(df_r)

	xtreg	h Ggov, re	vce(robust)												// re
	xttobit h Ggov, ll(Lower) ul(Upper)	vce(bootstrap, reps(200) seed(303625))	// re-Tobit
		margins , expression(predict(pr(0,.)) * predict(e(0,.))) dydx(*) post	// Low = -.94 (s.e. = .09)
		margins ,  coeflegend 

		lincom	_b[Ggov] - (-1)				// Test for complete crowd-out.
				local sign_coeff = sign(r(estimate))

			test	_b[Ggov] = -1
				scalar F = r(chi2)/r(df)
				
				di "One-sided test - Ho: Ggov >= 1, p-value = " r(p)/2			// Ho: crowd-out is -1 or a bigger negative number.
		
				di "Ho: Ggov >= 1  p-value = " ttail(dfr,`sign_coeff'*sqrt(F))	// using "t" p = .263
				di "Ho: Ggov >= 1  p-value = " 1 - normal(`sign_coeff'*sqrt(F))	// using "z" p = .262
restore



*** (1.2) Balanced-budget crowd-out at High
*
preserve

keep if  (Budget==4 | Budget==6)

	quietly xtreg	h Ggov, fe	vce(robust)		// fe . . . just to get the degrees of freedom for the t=test below.
				scalar dfr	= e(df_r)

	xtreg	h Ggov, re	vce(robust)												// re
	xttobit h Ggov, ll(Lower) ul(Upper)	vce(bootstrap, reps(200) seed(303625))	// re-Tobit
		margins , expression(predict(pr(0,.)) * predict(e(0,.))) dydx(*) post	// High = -.73 (s.e. = .08)
		margins ,  coeflegend 

		lincom	_b[Ggov] - (-1)				// Test for complete crowd-out.
				local sign_coeff = sign(r(estimate))

			test	_b[Ggov] = -1
				scalar F = r(chi2)/r(df)
				
				di "One-sided test - Ho: Ggov >= 1, p-value = " r(p)/2			// Ho: crowd-out is -1 or a bigger negative number.
		
				di "Ho: Ggov >= 1  p-value = " ttail(dfr,`sign_coeff'*sqrt(F))	// using "t" p = .00072707
				di "Ho: Ggov >= 1  p-value = " 1 - normal(`sign_coeff'*sqrt(F))	// using "z" p = .00049657
restore










* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (2) "xttobit" RE-Tobit applied to the LO and HI govt provision at the
*		same time (to estimate the change in crowd-out).
*

preserve

keep if  (Budget==2 | Budget==5 | Budget==4 | Budget==6)	// Balanced-budget crowd-out

	xtreg	h Ggov Ggov_GovtHigh GovtHigh, fe	vce(robust)	// fe . . . just to get the degrees of freedom for the t=test below.
			test Ggov_GovtHigh
				scalar dfr	= r(df_r)
				
	xttobit	h 	c.Ggov 		i.GovtHigh 		i.GovtHigh#c.Ggov, ll(Lower) ul(Upper)	///
									vce(bootstrap, reps(200) seed(303625))			// s.e. = .1145505, .1267951, .1135327 (depends on seed)
									
		margins ,  expression(predict(pr(0,.)) * predict(e(0,.)))	dydx(Ggov)	at(Ggov=(4 28) GovtHigh=(0 1))	post	// Low = -.9527689; High = -.7647187
		margins ,  coeflegend 

			lincom	_b[4._at] - _b[1bn._at]							// Test the difference: .188 (se = .11) . . . . . . Two-sided test.
																	//	p = 0.088 . . . depends on seed.
																	//	Identical: "test	_b[4._at] = _b[1bn._at]"
				return list
					local sign_coeff = sign(r(estimate))

			test	_b[4._at] = _b[1bn._at]							// "test" uses "z" p = 2*.04664772 =.087(5686) = approx= .0876
				return list
					scalar F = r(chi2)/r(df)						// F = 2.918511704347592 / 1
				
				di "One-sided test - Ho: i.GovtHigh#c.Ggov <= 0, p-value = " r(p)/2	// Ho: crowd-out is the same at GovtHigh (as it is at LowGovt),
																					//		 or larger-negative number. p = .0437843
		
				di "Ho: i.GovtHigh#c.Ggov <= 0  p-value = " ttail(dfr,`sign_coeff'*sqrt(F))		// using "t" p = .04563079
				di "Ho: i.GovtHigh#c.Ggov <= 0  p-value = " 1 - normal(`sign_coeff'*sqrt(F))	// using "z" p = .0437843  (approx= .044 . . . use
																								//	"z" because we are assuming normality to begin with here.)
restore



	
	
	log close	log1

	set more on
