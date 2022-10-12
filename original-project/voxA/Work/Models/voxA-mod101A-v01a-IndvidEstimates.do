//	file: 
		local Namedo "voxA-mod101A-v01a-IndvidEstimates"
		
//	task:	 Estimate the Impure altruism, Pure altruism and Pure warm glow 
//				Cobb-Douglas specifications for each individual participant.
//
//			Loop through all the participants.		

//	author:
		local aut "\ mow \ 2017-06-20"   

//	Include file to assign directories:
		voxAdir Models
		include "voxA-mod-dirInclude-001A.do"
		
// Input datasets:
	local FileIn01	"voxA-cr021A_v01a-DecisionsLONG"		// h(b), Ggov(b), ... b = 1, 2, 3, 4, 5, 6.

	local FileIn02	"voxA-cr011A_v01a-BoundaryDecisions"	// The corners, n = 85.

	
// Output dataset:
	local FileOut01	" "		// If you are running this do-file to UPDATE the estimates,
							//	put the filename you want here. Otherwise, leave blank.
	
	*local FileOut01	"voxA-mod101A-v01a-IndvidEstimates-2017-06-20--09-45"
	

// 	Version
	version 14.2			// Version of Stata that produced these estimates.


#delimit ;
*
For each participant:

(1) Get estimates for the impure altruism solution (mlcobbPLUScornLN.ado).

	Check to see if this  (a) produces estimates of alpha and beta in the normal range, that is:
								(i)		alpha in (0, 1)
								(ii)	 beta in (0, 1)
								(iii)	alpha + beta in (0, 1), and then
								
						  (b) check for convergence. 

	
(2) Estimate the Pure Altruism model (mlcobbPureAltruismcornLN.ado).

	Check for alpha in (0, 1) and for convergence.
	
	
(3) Estimate the Pure Warm Glow model (mlcobbPureWGcornLN.ado).

	Check for beta in (0, 1) and for convergence.

	
(4) Choose the estimate (for each participant) based the three checks, and 
		based on the biggest log-likelihood value among the three models.
;
#delimit cr



		
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (0) Initialization.
*

tempfile Decisions Corners


*** (0.1) Get the decision data (Long).
*
use	"`dirIN'\\`FileIn01'", clear			

	keep id newid Budget h Ggov z y 

	sum y if inrange(Budget, 1, 4)
	sum y if inrange(Budget, 5, 6)

	save `Decisions'
	

*** (0.2) Get the corner indicators (n = 85).
*
use	"`dirIN'\\`FileIn02'", clear			

	keep id AlwaysCorn

	tab AlwaysCorn, miss
	
		save `Corners'
	
		
*** (0.3) Merge
*
use `Decisions', clear

merge m:1	id	using	`Corners', assert(match) nogenerate	



*** (0.4) Get the number of distinct indivuduals.
*
	count
		local N = r(N)/6


		
*** (0.5) Set up the dataset into which the estimates will be stored.
*

* Create a timestamped name for the dataset, so that the estimates are
*	not inadvertently overwritten (by a later careless re-running of this do-file).
*
local date : di %td_CCYY-NN-DD  date("$S_DATE", "DMY")

local hour = substr("$S_TIME", 1, 2)

local minute = substr("$S_TIME", 4, 2)

local timestamp = trim("`date'--`hour'-`minute'")

local dirESTNamedoTS = "`dirESTNamedo'-`timestamp'"


if ("`FileOut01'"   ~=   " ") {				// However, if a pre-existing file is
											//	being updated then use that old filename.
	local dirESTNamedoTS = "`dirEST'\\`FileOut01'"
}



* Initialize the dataset that will contain the estimates,
*
*	UNLESS this do-file is being run to UPDATE a pre-existing dta-file.
*
if ("`FileOut01'"   ==   " ") {		// This is a new file (not an update)
	preserve
		keep newid
	
		duplicates drop		// Keep just one record per individual

		foreach x in alpha beta sigma ModelSelect Converged {
			generate `x' = .
		}
	
		label variable alpha "Altruism parameter"
		label variable  beta "Warm glow parameter"
		label variable sigma "RMSE"
	
		label variable ModelSelect "Impure, Altr, WarmGlow"
		label variable Converged   "Convergence report"
		
		note alpha : Estimate of the altruism parameter for this participant \ /*
					*/ ".c" means that this participant was always at a corner /*
					*/ (no estimates are possible) \ /*
					*/ ".s" means that sensible estimates for this participant could not /*
					*/ be obtained (usually means there are too many corners, albeit fewer than six) \ /*
					*/ System missing means the estimates have not yet been done for this participant \ `tag'
				
		note  beta : Estimate of the warm glow parameter for this participant \ /*
					*/ ".c" means that this participant was always at a corner /*
					*/ (no estimates are possible) \ /*
					*/ ".s" means that sensible estimates for this participant could not /*
					*/ be obtained (usually means there are too many corners, albeit fewer than six) \ /*
					*/ System missing means the estimates have not yet been done for this participant \ `tag'
				
		note sigma : Root-mean square error from the estimation for this participant \ /*
					*/ ".c" means that this participant was always at a corner /*
					*/ (no estimates are possible) \ /*
					*/ ".s" means that sensible estimates for this participant could not /*
					*/ be obtained (usually means there are too many corners, albeit fewer than six) \ /*
					*/ System missing means the estimates have not yet been done for this participant \ `tag'
						
		note ModelSelect : Model that generated the estimates for this participant \ /*
					*/ Impure altruism, Altruism (pure), Warm glow (pure) \ /*
					*/ ".c" means that this participant was always at a corner /*
					*/ (no estimates are possible) \ /*
					*/ ".s" means that sensible estimates for this participant could not /*
					*/ be obtained (usually means there are too many corners, albeit fewer than six) \ /*
					*/ System missing means the estimates have not yet been done for this participant \ `tag'
				
		note Converged : 1 means the model selected for this participant did in fact converge \ /*
					*/ ".c" means that this participant was always at a corner /*
					*/ (no estimates are possible) \ /*
					*/ ".s" means that sensible estimates for this participant could not /*
					*/ be obtained (usually means there are too many corners, albeit fewer than six) \ /*
					*/ System missing means the estimates have not yet been done for this participant \ `tag'

		label define ModelSelectLBL 1 "1impure" 2 "2altr" 3 "3warm"
		label value  ModelSelect ModelSelectLBL
	
		label define ConvergedLBL 1 "1yes" 0 "0no"
		label value  Converged ConvergedLBL

		count
			local N = r(n)

		local vvv = c(version)
		label data "VOX expr: Individual-level estimates \ 2017-06-20"
		notes: `Namedo'.dta \ Sample size N = `N' participants \ /*
					*/ Version (during ML estimation) = `vvv' \ `tag'

		save 	"`dirESTNamedoTS'" , replace
	restore
}	
	
	



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (1) Estimation.
*

*** >>>>>>>>>>>>>>>> Beginning of the loop on i = 1, 2, ..., 85 participants.
***
***
forvalues i = 1(1)85 {


local participant = `i'




	
di as txt "At the BEGINNING of the estimation for Participant = `participant'" 
di as txt "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" 



*** *** *** (1) Determine if this participant was always at a corner (i.e., $0,
*					or $40/$46).  If so, then skip estimation for this participant.
*
sum AlwaysCorn if (newid==`participant')

	local Corner6 = r(mean)

if ("`Corner6'" == "0") {		// NOT "6 corner decisions". Proceed with estimation.
								// vvvvvvvvvvvvvvvvvvv Not "6 corners" vvvvvvvvvvvvv





	
*** *** *** (2) Several participants require estimation options when using their data
*					to estimate the Cobb-Douglas impure altruism model.
*
*					These options typically do NOT help estimation problems, because 
*						impure altruism model estimaton is converging toward nonnsese 
*						values, or the loglikelihood function is not concave, or
*						the numerical derivatives cannot be calculated.
*
local EstOptions = " "		// Most participants do not require ML max options

if ("`participant'" == "10") {
	local EstOptions = "difficult log trace iterate(14)"
}

if ("`participant'" == "11") {
	local EstOptions = "difficult iterate(18)"
}

if ("`participant'" == "12") {
	local EstOptions = "difficult iterate(16)"
}	
	
if ("`participant'" == "16") {
	local EstOptions = "iterate(147)"
}	
	
if ("`participant'" == "17") {
	local EstOptions = "iterate(22)"
}	
	
if ("`participant'" == "22") {
	local EstOptions = "difficult iterate(19)"
}		
	
if ("`participant'" == "24") {
	local EstOptions = "difficult log trace iterate(137)"
}		
	
if ("`participant'" == "30") {
	local EstOptions = "difficult log trace iterate(17)"
}	
	
if ("`participant'" == "43") {
	local EstOptions = "difficult log trace iterate(13)"
}	
	
if ("`participant'" == "52") {
	local EstOptions = "difficult log trace iterate(14)"
}	
	
if ("`participant'" == "53") {
	local EstOptions = "difficult log trace iterate(18)"
}	
	
if ("`participant'" == "59") {
	local EstOptions = "difficult log trace iterate(19)"
}	
	
if ("`participant'" == "60") {
	local EstOptions = "difficult"
}	
	
if ("`participant'" == "71") {
	local EstOptions = "difficult log trace iterate(17)"
}	

if ("`participant'" == "78") {
	local EstOptions = "difficult log trace iterate(27)"
}	

if ("`participant'" == "79") {
	local EstOptions = "difficult"
}	

if ("`participant'" == "81") {
	local EstOptions = "difficult"
}	

if ("`participant'" == "83") {
	local EstOptions = "difficult log trace iterate(15)"
}	


	
*** *** *** (3) Get starting values
*			
	regress h Ggov z if (newid==`participant') , noconstant	// alpha and beta starters

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
	
di "Astart contains |`Astart'|"


local Bstart = `b'			
if `b' < 0 { 
		local Bstart = .10
	}
if `b' > 1 {
		local Bstart = .45	
	}


if (`LNSstart' ==.) {
	local LNSstart = 1
}





*** *** *** (Estimation 1) Impure altruism solution.
*
	ml model lf mlcobbPLUScornLN (alpha:) (beta:) (LNsigma:) if (newid== `participant')
		ml init alpha:_cons=`Astart'   beta:_cons=`Bstart'    LNsigma:_cons=`LNSstart'
	ml max, `EstOptions'

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
	ml model lf mlcobbPureAltruismcornLN (alpha:) (LNsigma:) if (newid== `participant')
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
	ml model lf mlcobbPureWGcornLN (beta:) (LNsigma:) if (newid== `participant')
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
*				(b) ML estimation of the impure altruism model converged.
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
		

		
		
		
* . . . Check that the Pure altruism "alpha" is sensible, meaning:
*
*				(a) The "alpha" in (0, 1), and
*
*				(b) ML estimation of the Pure altruism model converged.

	local a		 	 = alphaALTR
	local alphaCHECK = " "				// Precautionary re-set.
	
	
	if (inrange(`a', 0, 1)) {
			local alphaCHECK = "ok"
		}
		else if (`a' < 0) {
			local alphaCHECK = "neg"
		}
		else if (`a' > 1) {
			local alphaCHECK = "biggerONE"
		}
		
		
	if	("`alphaCHECK'" == "ok") & (convergedALTR == 1) {
		
			local sensibleALTR = "yes"
		}
		else {
			local sensibleALTR = "no"		
		}		

		
		
		
		
		
* . . . Check that the Pure warm glow "beta" is sensible, meaning:
*
*				(a) The "beta" in (0, 1), and
*
*				(b) ML estimation of the Pure warm glow model converged.

	local a		 	 = betaWARM
	local betaCHECK = " "				// Precautionary re-set.
	
	
	if (inrange(`a', 0, 1)) {
			local betaCHECK = "ok"
		}
		else if (`a' < 0) {
			local betaCHECK = "neg"
		}
		else if (`a' > 1) {
			local betaCHECK = "biggerONE"
		}
		
		
	if	("`betaCHECK'" == "ok") & (convergedWARM == 1) {
		
			local sensibleWARM = "yes"
		}
		else {
			local sensibleWARM = "no"		
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
				scalar `x'Scalar = `x'IMP
			}
			
			scalar ModelSelectScalar = 1
			scalar ConvergedScalar   = convergedIMP
		
		}
		else if ("`sensibleALTR'"  == "yes") &		/// Pure altruism produced sensible estimates. . .
				("`LogBigSingleL'" == "ALTR") {		//	. . . and the biggest loglikelihood.
		
			
					foreach x in alpha sigma {
						scalar `x'Scalar = `x'ALTR
						}
					
					scalar betaScalar = 0
					
					scalar ModelSelectScalar = 2
					
					scalar ConvergedScalar   = convergedALTR				
			}

		
		else if ("`sensibleWARM'"  == "yes") &		/// Pure warm glow produced sensible estimates. . .
				("`LogBigSingleL'" == "WARM") {		//	. . . and the biggest loglikelihood.
		
			
					foreach x in beta sigma {
						scalar `x'Scalar = `x'WARM
						}
					
					scalar alphaScalar = 0

					scalar ModelSelectScalar = 3
					
					scalar ConvergedScalar   = convergedWARM
			}

		else {										// None of the three models produced
													//	sensible estimates.
			
					foreach x in alpha beta sigma ModelSelect Converged {
						scalar `x'Scalar = .s
					}
		}
			
	


* . . .	and store the estimates in the dta-file.
*
	preserve
	
		use "`dirESTNamedoTS'" , clear

		foreach x in alpha beta sigma ModelSelect Converged {
			replace `x' = `x'Scalar	if (newid == `participant')
		}
		
		tab ModelSelect, miss
		
		mkmat newid alpha beta ModelSelect Converged, matrix(Atemp) // Print the results.
		
		matrix list Atemp, format(%7.2f)		
		
		save "`dirESTNamedoTS'" , replace
		
	restore
	

	
* . . . Also print results to screen.
*		
di _newline(1)
	di as result "Model selection = " ModelSelectScalar
	di as result "Converged (?)   = " ConvergedScalar
	di "alpha = " round(alphaScalar, .01)
	di "beta  = " round( betaScalar, .01)
di _newline(1)







								// End of estimation
}								// ^^^^^^^^^^^^^^^^^^^ Not "6 corners" ^^^^^^^^^^^^^

else {							// This section handles the participants always at a corner.
								// vvvvvvvvvvvvvvvvvvv Yes "6 corners" vvvvvvvvvvvvvvvvvvvv


								
* . . .	Store ".c" in the dta-file.
*
	preserve
	
		use "`dirESTNamedoTS'" , clear

		foreach x in alpha beta sigma ModelSelect Converged {
			replace `x' = .c	if (newid == `participant')
		}

		mkmat newid alpha beta ModelSelect Converged, matrix(Atemp) // Print the results.
		
		matrix list Atemp, format(%7.2f)	
		
		save "`dirESTNamedoTS'" , replace
		
	restore
	
	
* . . . Also print results to screen.
*		
di _newline(1)

di in red "Participant = `participant'  was ALWAYS AT A CORNER." 
di in red "    Estimates are not possible." 
di _newline(1)


								// ^^^^^^^^^^^^^^^^^^^ Yes "6 corners" ^^^^^^^^^^^^^
								// Finished with the Not "6 corners" / Yes "6 corners"
}								//	conditional statement.












	
di as result "At the end of the estimation for Participant = `participant'" 
di as result "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" 
di as result "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" 
di as result "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" 
di _newline(7)		


}
*** <<<<<<<<<<<<<<<<<<<<<< End of the loop on i = 1, 2, ..., 85 participants.
***
***






	log close	log1
