//	file: 
		local Namedo "voxA-mod190A-v01a-IndvidEstimatesChecker"
		
//	task:	 Checks to dta-files containing individual-level estimates
//				agains each other.

//	author:
		local aut "\ mow \ 2017-06-20"   

//	Include file to assign directories:
		voxAdir Models
		include "voxA-mod-dirInclude-001A.do"
		
// Input datasets:
	local MatrixIn01	"icd001-2014-02-19"		// *.mmat = 85 x 4 = alpha, beta, sigma. 
	
	
	local FileIn01	"voxA-mod101A-v01a-IndvidEstimates-2017-06-20--16-38"	// alpha_i, beta_i, i = 1, ..., 85.


	
	
	
	
	
	
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* 0. Get the estimates
*
tempfile Old


*** 0.1 First set of estimates.
*
mata:
	mata clear
	mata matuse "`dirMAT'/`MatrixIn01'"
end


* Convert the matrix of estimates into Stata variables.
*
getmata (newid alphaOLD betaOLD sigmaOLD) = icd001


label variable alphaOLD "alpha from Older file"
label variable  betaOLD "beta from Older file"
label variable sigmaOLD "rmse from Older file"	// root mean-squared error

	save `Old'


	
*** 0.2 Second set of estimates, and merge.
*
use "`dirEST'/`FileIn01'" , clear


merge 1:1	newid	using `Old', assert(match) nogenerate




* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* 1. Check.
*
assert (  round(alpha, .001) == round(alphaOLD, .001)  )	if !missing(alphaOLD)

assert (  round( beta, .001) == round( betaOLD, .001)  )	if !missing( betaOLD)



*	list newid alpha alphaOLD	if (round(alpha, .001) ~= round(alphaOLD, .001))






	log close	log1
