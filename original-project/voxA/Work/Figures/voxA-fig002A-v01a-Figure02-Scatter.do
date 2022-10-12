//	file: 
		local Namedo "voxA-fig002A-v01a-Figure02-Scatter"
	
//	task: Summary statistics and scatter plot for the n = 85 individual-level 
//			estimates from the Cobb-Douglas impure altruism model.

//	author:
		local aut "\ mow \ 2017-06-17"   

//	Include file to assign directories:
		voxAdir Figures
		include "voxA-fig-dirInclude-001A.do"
		
// Input datasets:
	local MatrixIn01	"icd001-2014-02-19"		// *.mmat = 85 x 4 = alpha, beta, 



		
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* 0. Get the estimates
*
mata:
	mata clear
	mata matuse "`dirMAT'/`MatrixIn01'"
end



* Convert the matrix of estimates into Stata variables.
*
getmata (newid alpha beta sigma) = icd001


label variable alpha "alpha"
label variable  beta "beta"
label variable sigma "rmse"		// root mean-squared error




* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* 1. Generate new variables to describe the structural Motives.
*
gen byte PureAltruism	= (alpha >  0) & (beta == 0) & !missing(alpha) & !missing(beta)
gen byte PureWG			= (alpha == 0) & (beta >  0) & !missing(alpha) & !missing(beta)
gen byte Impure			= (alpha >  0) & (beta >  0) & !missing(alpha) & !missing(beta)

replace PureAltruism	=.	if missing(alpha) | missing(beta)
replace PureWG			=.	if missing(alpha) | missing(beta)
replace Impure			=.	if missing(alpha) | missing(beta)

list newid alpha beta if PureAltruism==1		// Check
list newid alpha beta if PureWG==1
list newid alpha beta if Impure==1
list newid PureAltruism PureWG Impure	if missing(alpha) | missing(beta)




gen byte MotivesStr001 = 1	if PureAltruism==1
replace  MotivesStr001 = 2	if Impure==1
replace  MotivesStr001 = 3	if PureWG==1

label define MotivesStrLBL 1 "1altruism"	2 "2impure"		3 "3warm"
label values  MotivesStr001		MotivesStrLBL

tab PureAltruism	MotivesStr001, miss		// Check
tab PureWG			MotivesStr001, miss
tab Impure			MotivesStr001, miss



tab MotivesStr001, miss		// 								Pure altruism .44.	Impure .39.		Pure WG .09		Indeterminate (corners) .08
tab MotivesStr001			// Conditional on NOT missing:	Pure altruism .47.	Impure .42.		Pure WG .10		Indeterminate (corners) .08


label variable	MotivesStr001 "Motives 001 structure"
notes			MotivesStr001 : Motives based on alpha, beta estmates from the "icd001" vector \ 001 version of the individual Cobb-Douglas /*
								*/ estimates \ `tag'


sum PureAltruism Impure PureWG

clonevar alphaZero = PureWG
clonevar  betaZero = PureAltruism

lab var alphaZero "Indicates alpha - zero"
lab var  betaZero "Indicates  beta - zero"

sum alphaZero betaZero





* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* 2. Scatter -- Figure 2 in the AER paper.
*
list alpha beta if beta >= .60 & beta  <= .70

count if !missing(alpha) & !missing(beta)
	local Nobs = r(N)

gen x45 = _n/(`Nobs')			// 45 degree line
gen y45 = _n/(`Nobs')
replace x45 = . if x45 > 1
replace y45 = . if y45 > 1

*** The July 2015 version
*
graph twoway	scatter alpha	beta, mcolor("202 194 126")											///
			||	line	y45		x45	, connect(l) lpattern(shortdash)	lcolor(gs8) lwidth(medium)	///
	graphregion(color(white))																		///
	ylabel(0(.2)1, angle(horizontal))	xlabel(0(.2)1)																	///
	ytitle("Altruism ({it:{&alpha}{sub:i}})") xtitle("Warm glow ({it:{&beta}{sub:i}})")				///
	yline(0, lpattern(solid) lcolor("192 192 192") lwidth(thin) )									///
	xline(0, lpattern(solid) lcolor("192 192 192") lwidth(thin) )									///
	text(.83 .85 "{it:{&alpha}{sub:i}} = {it:{&beta}{sub:i}}"	, placement(e) )					///
	legend(off)																						///
	xsize(4.5) ysize(4.5) 

	
	
	
	
	

	log close	log1
	set more on


