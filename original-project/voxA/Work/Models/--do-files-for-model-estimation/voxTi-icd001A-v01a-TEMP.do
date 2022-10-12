//	file: 
		local Namedo "voxTi-icd001A-v01a-TEMP"
		
//	task:	 
//			

//	author:
		local aut "\ mow \ 2017-06-16"   

//	Include file to assign directories:
		voxTdir Models
		include "voxT-mod-dirInclude-001A.do"
		
// Input datasets:
	local FileIn01	"voxT-cr021A_v01a-DecisionsLONG"		// h(b), Ggov(b), ... b = 1, 2, 3, 4, 5, 6.

// Output matrix
	global icdName = "icd001"		// Name the matrix used to store the estimation results.

	
	
// Notes:	(1) 
//					
//
//			(2) 
//
//				
//
//

#delimit ;
*
For each participant:

(1) Get estimates for the impure altruism solution using the positive sign in front of the square-root term (mlcobbPLUScornLN.ado).

	If this works produces sensible estimates of alpha and beta, then check Step (2)

	
(2) Get estimates fof the impure altruism solution using the negative sign in front of the square-root term (mlcobbNEGcornLN.ado).

	If this does not work, go with Step (1) estimates.
	
	If this also produces sensible estimates, choose the Step (1) versus Step (2) estimates based on the largest log-likelihood value.
	
	
(3) If NEITHER Step (1) NOR Step (2) produces satisfactory/sensible results, then:
	
	(a) Estimate the Pure Altruism model (mlcobbPureAltruismcornLN.ado).
	
	(b) Estimate the Pure Warm Glow model (mlcobbPureWGcornLN.ado).
	
	(c) Choose the Step (3a) versus Step (3b) estimates based on the largest log-likelihood value.
;
#delimit cr



		
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (0) Initialization.
*

*** (0.1) Get the data
*
use	"`dirIN'/`FileIn01'", clear			

keep id newid Budget h Ggov z y 

sum y if inrange(Budget, 1, 4)
sum y if inrange(Budget, 5, 6)



*** (0.2) Set up the 85 x 4 matrix into which the estimates will be stored.
*


* Create a timestamped name for the matrix, so that the estimates are
*	not inadvertently overwritten (by a later re-running of this do-file).
*
local date : di %td_CCYY-NN-DD  date("$S_DATE", "DMY")

local hour = substr("$S_TIME", 1, 2)

local minute = substr("$S_TIME", 4, 2)

local timestamp = trim("`date'--`hour'-`minute'")

global mmm = "`dirMAT'/icd-`timestamp'"

* The matrix
mata:

	$icdName = J(85, 4, 99)

	fh = fopen("$mmm", "rw")
	
	fputmatrix(fh, $icdName)
	
	fclose(fh)

end


matrix define B = (9, 8, 7)

global participant = 1

sss

---------------

*** Code to list the matrix.
*

mata:
	
	fh = fopen("$mmm", "r")					//
	
	bbb = fgetmatrix(fh)
	
	bbb[1, .]

	fclose(fh)
end		







-----------------------

	bbb[i, .] = (i, Brow)							// Store the estimates

	fputmatrix(fh, bbb)
mata:
	fh = fopen("$mmm", "rw")
	
	bbb = fgetmatrix(fh)
	
	bbb
	
	fclose(fh)

end





sss

mata matsave 	"$mmm"		$icdName






* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (1) Estimation.
*

*** >>>>>>>>>>>>>>>> Beginning of the loop on i = 1, 2, ..., 85 participants.
***
***
forvalues i = 1(1)1 {

*forvalues i = 10(1)10 {    , difficult iterate(3)

	local participant = `i'

	
di as txt "At the BEGINNING of the estimation for Participant = `participant'" 
di as txt "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" 

	
	
*** *** *** (1) Several participants require estimation options when using their data
*					to estimate the Cobb-Douglas impure altruism model.
*
local EstOptions = " "		// Most participants do not require ML max options

			
	
	
*** *** *** (2) Get starting values
*			
	regress h Ggov z if newid==`participant' , noconstant	// alpha and beta starters

		local  Sstart = e(rmse)
		local  LNSstart = log(`Sstart')
			
		matrix B = e(b)
			matrix x = B[1, "Ggov"]
				local b = 1 + x[1,1]

			matrix x = B[1, "z"]
				local a = x[1,1]
		
			di as result "a = " `a' "  and   b = " `b'


local Astart = `a'			
if `a' < 0 { 
		local Astart = .10
	}
if `a' > 1 {
		local Astart = .45	
	}
	
di "Astart conatins |`Astart'|"


local Bstart = `b'			
if `b' < 0 { 
		local Bstart = .10
	}
if `b' > 1 {
		local Bstart = .45	
	}
	
di as result "Starting values are      Astart = " `Astart' "  and   Bstart = " `Bstart'





*** *** *** (Estimation 1) Impure altruism solution.
*
*version 10
	ml model lf mlcobbPLUScornLN (alpha:) (beta:) (LNsigma:) if newid== `participant'
		ml init alpha:_cons=`Astart'   beta:_cons=`Bstart'    LNsigma:_cons=`LNSstart'
	ml max

	ereturn list

		matrix B = e(b)
			matrix temp = B[1, "alpha:_cons"]
			scalar alphaIMP = temp[1,1]
			di alphaIMP
			
			matrix temp = B[1, "beta:_cons"]
			scalar betaIMP = temp[1,1]
			di betaIMP

			matrix temp = B[1, "LNsigma:_cons"]
			scalar sigmaIMP = exp(temp[1,1])

		scalar loglikeIMP = e(ll)

		scalar convergedIMP = e(converged)

		
*** *** *** (Estimation 2) Pure altruism solution.
*
	ml model lf mlcobbPureAltruismcornLN (alpha:) (LNsigma:) if newid== `participant'
		ml init alpha:_cons=`Astart'       LNsigma:_cons=`LNSstart'
	ml max
	
			matrix B = e(b)
			matrix temp = B[1, "alpha:_cons"]
			scalar alphaALTR = temp[1,1]
			di alphaALTR

			matrix temp = B[1, "LNsigma:_cons"]
			scalar sigmaALTR = exp(temp[1,1])

		scalar loglikeALTR = e(ll)
		
		scalar convergedALTR = e(converged)
	
	
	
*** *** *** (Estimation 3) Pure warm glow solution.
*
	ml model lf mlcobbPureWGcornLN (beta:) (LNsigma:) if newid== `participant'
		ml init beta:_cons=`Bstart'       LNsigma:_cons=`LNSstart'
	ml max
	
		
			matrix B = e(b)
			matrix temp = B[1, "beta:_cons"]
			scalar betaWARM = temp[1,1]
			di betaWARM
			
			matrix temp = B[1, "LNsigma:_cons"]
			scalar sigmaWARM = exp(temp[1,1])

		scalar loglikeWARM = e(ll)

		scalar convergedWARM = e(converged)		
		
		
		
		
*** *** *** (Checking 1) Check whether alphaIMP, betaIMP are sensible.
*
*							Note: "Sensible" means
*									(a) alpha and beta are in "normal" ranges
*											and
*									(b) Estimation converged.
*
	local a		 = alphaIMP
	local b		 =  betaIMP
	
	local aPLUSb = alphaIMP + betaIMP

	
* . . . Check that alpha is in [0, 1].
*
	if (inrange(`a', 0, 1)) {
			local alphaCHECK = "ok"
		}
		else if (`a' < 0) {
			local alphaCHECK = "neg"
		}
		else if (`a' > 1) {
			local alphaCHECK = "biggerONE"
		}

		
* . . . Check that beta is in [0, 1].
*		
	if (inrange(`b', 0, 1)) {
			local betaCHECK = "ok"
		}
		else if (`b' < 0) {
			local betaCHECK = "neg"
		}
		else if (`b' > 1) {
			local betaCHECK = "biggerONE"
		}

	
* . . . Check that   alpha + beta   is in [0, 1].
*		
	if (inrange(`aPLUSb', 0, 1)) {
			local alphaPLUSbetaCHECK = "ok"
		}
		else if (`aPLUSb' < 0) {
			local alphaPLUSbetaCHECK = "neg"
		}
		else if (`aPLUSb' > 1) {
			local alphaPLUSbetaCHECK = "biggerONE"
		}
		

		
* . . . Combine all three into one check: "sensible Impure altruism"
*
*			"Sensible" means:
*
*				(a) The "alpha" and "beta" are in the normal range, and
*
*				(b) ML estimation of the model converged.
*
	if	("`alphaCHECK'" == "ok") & ("`betaCHECK'" == "ok") &	///
		("`alphaPLUSbetaCHECK'" == "ok") {
		
			local rangeIMP = "ok"
		}
		else {
			local rangeIMP = "NOTok"		
		}


	if	("`rangeIMP'" == "ok") & (convergedIMP == 1) {
		
			local sensibleIMP = "yes"
		}
		else {
			local sensibleIMP = "no"		
		}		
		

		
		
		
* . . .	(Checking 2) Now, find the largest log-likelihood value . . .
*
*			(a) . . . among the two single motive models, pure altruism, pure warm glow.
*
*			(b) . . . among hhose two and also including impure alrtuism.
*
	if (loglikeALTR > loglikeWARM) {	/// Biggest loglikelihhod among the single-motive
										//	models.
			local	LogBigSingleL = "ALTR"
			scalar 	LogBigSingle = loglikeALTR
			
		}
		else {
			local	LogBigSingleL = "WARM"		
			scalar 	LogBigSingle = loglikeWARM
		}

		
	if (loglikeIMP > LogBigSingle) {	/// . . . include dual-motive impure altruism
										//	in the comparison.
			local	LogBigDualL = "IMP"
			scalar	LogBigDual = loglikeIMP
			
		}
		else {
			local	LogBigDualL = "`LogBigSingleL'"
			scalar	LogBigDual =   LogBigSingle
		}


			
			

* . . .	(Final step) Select the model that best fits this participant.
*
	if	("`sensibleIMP'" == "yes") &	/// Impure altruism produced sensible estimates. . .
		("`LogBigDualL'" == "IMP") {	//	. . . and the biggest loglikelihood.
	
			foreach x in alpha beta sigma {
				scalar `x' = `x'IMP
			}
			
			scalar ModelSelect = 1
		
		}
		else {							/// Either	(a) Impure altruism NOT sensible, or
										//			(b) NOT the biggest loglikelihood.
		
			if ("`LogBigSingleL'" == "ALTR") {			/// Pure altruism
			
					foreach x in alpha beta sigma {
						scalar `x' = `x'ALTR
					}
					
					scalar ModelSelect = 2
				
			}
			else {										/// Pure warm glow
			
					foreach x in alpha beta sigma {
						scalar `x' = `x'WARM
					}
					
					scalar ModelSelect = 3

			}
			
		}
	


* . . .	and store the estimates in a matrix.
*
	matrix define B = (alpha, beta, sigma)
			
	sss

	
		
		
di _newline(1)
	di as result "Model selection = " ModelSelect
	di "alpha = " alpha
	di "beta  = " beta	
di _newline(1)
	
di as result "At the end of the estimation for Participant = `participant'" 
di as result "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" 
di as result "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" 
di as result "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" 
di _newline(7)		


}
*** <<<<<<<<<<<<<<<<<<<<<< End of the loop on i = 1, 2, ..., 85 participants.
***
***





STOP HERE 2017-06-19 ..............................................




	mata: storeresults
		2 + 3
	end


	mata:
		mata set matastrict on
		
		// storeresults		mow 2017-06-19
		void storeresults(real rowvector var1)
		{	real rowvector x
		
			
			
			
			
			
		mata matuse "`dirMAT'/icd-`timestamp'"				// Get the saved matrix of estimates

		i		= strtoreal(st_local("participant"))		// Participant number
		Brow	= st_matrix("B")							// Estimates: alpha, beta, sigma
	
		$icdName[i, .] = (i, Brow)							// Store the estimates
	
		mata matsave "`dirMAT'/icd-`timestamp'"	$icdName, replace	// Save the matrix of estimates

	end	




		
di in red "Got here .....  1"
di in red "local sensibleIMP conatins |`sensibleIMP'|"		
		
di in red "local LogBigSingleL conatins |`LogBigSingleL'|"		
di in red "local LogBigDualL   conatins |`LogBigDualL'|"		
		
		
			

di in red "Got here .....  2"


			
		
		
	di "alphaCHECK contains |`alphaCHECK'|"
	di "betaCHECK contains |`betaCHECK'|"
	di "alphaPLUSbetaCHECK contains |`alphaPLUSbetaCHECK'|"

	di _newline(1)
	di "sensibleIMP contains |`sensibleIMP'|"

	di _newline(1)
	di "LogBigSingle contains |`LogBigSingle'|"
	di "Scalar LogBigSingle = " LogBigSingle 

	di _newline(1)
	di "LogBigDual   contains |`LogBigDual'|"
	di "Scalar LogBigDual = " LogBigDual 		




----------------------------------------------------------
#delimit ;
*
For each participant:

(1) Get estimates for the impure altruism solution using the positive sign in front of the square-root term (mlcobbPLUScornLN.ado).

	If this works produces sensible estimates of alpha and beta, then check Step (2)

	
(2) Get estimates fof the impure altruism solution using the negative sign in front of the square-root term (mlcobbNEGcornLN.ado).

	If this does not work, go with Step (1) estimates.
	
	If this also produces sensible estimates, choose the Step (1) versus Step (2) estimates based on the largest log-likelihood value.
	
	
(3) If NEITHER Step (1) NOR Step (2) produces satisfactory/sensible results, then:
	
	(a) Estimate the Pure Altruism model (mlcobbPureAltruismcornLN.ado).
	
	(b) Estimate the Pure Warm Glow model (mlcobbPureWGcornLN.ado).
	
	(c) Choose the Step (3a) versus Step (3b) estimates based on the largest log-likelihood value.
;
#delimit cr





















---------------------------------------------------
---------------------------------------------------
---------------------------------------------------
---------------------------------------------------
*** Code to list the matrix.
*
local dirMAT		"/home/r.turquier/media/disk/_Users_remi_Documents_GitHub_robust-micrometrics/original-project/voxT/Work/Models/Matrices" 

mata:
	mata clear
	mata matuse "`dirMAT'/icd-2017-06-19--16-15"	// Get the saved matrix of estimates
	
	$icdName									// List the matrix
end		

 




*	mata matuse "`dirMAT'/icd-2017-06-17--12-59"	// Get the saved matrix of estimates


