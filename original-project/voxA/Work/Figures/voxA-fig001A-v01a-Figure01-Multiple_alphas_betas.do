//	file: 
		local Namedo "voxA-fig001A-v01a-Figure01-Multiple_alphas_betas"
//
//	task:	Figure 1 for the AER paper.
//
//				For specific (alpha, beta) pairs, solves the Cobb-Douglas
//					Fix q2 as a function of G_{-i}.
//

//	author:
		local aut "\ mow \ 2014-06-20"
	
// Input datasets: None

//	Include file to assign directories:
		voxAdir Figures
		include "voxA-fig-dirInclude-001A.do"

	
	

* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Part 1: The set-up
*
clear all

*** The x-axis
*
local N = 101

set obs `N'

generate Ggov = (_n - 1)		

label var Ggov "Giving by others"



*** Income (or government provision) held constant
*
local y = 40			// Income


*** Cobb-Douglas altruism (alphas), warm glow (betas) parameters
*
local beta1  = .10
local alpha1 = .40

local beta2  = .01
local alpha2 = .27			// .2685698

local beta3  = .35
local alpha3 = .35





* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Part 2: q2 solutions from the Cobb-Douglas:
*
forvalues i = 1/3 {

	gen	S`i' = ( (1 - `alpha`i'')^2 * Ggov^2	 + 2*( (`beta`i'' - `alpha`i'') + (`alpha`i'' + `beta`i'')*`alpha`i'' ) * Ggov * `y'  +  (`alpha`i'' + `beta`i'')^2 * `y'^2		)	

	gen	Sroot`i' = sqrt(S`i')	
	
	gen q2pos`i' = .5 * ( (1 - `beta`i'')			+      ( ((1 - `alpha`i'') - (1 + `alpha`i'')*`beta`i'' ) * Ggov   +  ( (`beta`i'' - `alpha`i'') - (`alpha`i'' + `beta`i'')*`beta`i'' ) * `y') / Sroot`i'	)

	gen KBpos`i'	= abs(-1 + q2pos`i')
	
	label var Sroot`i'	"Square-root term"

	label var q2pos`i'	"Difference between the two income effects"
	
	
	* Optimal giving -- g*
		generate double  B`i' =  (1 - `beta`i'')*Ggov + (`alpha`i'' + `beta`i'')*(Ggov + `y')	
			
		generate double Hplus`i' = .5*(B`i' + sqrt(B`i'^2 - 4*`alpha`i'' * Ggov * (Ggov + `y')) )
		generate double gplus`i' =  Hplus`i' - Ggov

		label variable gplus`i' "Optimal giving"	
}

label var KBpos1	"Altruism strong relative to warm glow"
label var KBpos2	"Altruism super strong relative to warm glow"
label var KBpos3	"Altruism, warm glow similar strength"

list Ggov KBpos1 KBpos2 KBpos3 gplus1 gplus2 gplus3 if Ggov < 10	// Confirm gplus is different

list Ggov KBpos1 KBpos2 KBpos3 gplus1 gplus2 gplus3 if (10 <= Ggov) & (Ggov <= 20)


	
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Part 3: Plot
*

*** Figure 1 Graph of q2, balanced-budget crowd-out
*
*
*		Note: This is Figure 2 in the NBER/September 2014 version, but
*
*				re-christened Figure 1 in the July 2015 version
*


* Alternative version with different axis-label on the vertical axis:
graph twoway	line KBpos1		 Ggov, connect(l) lpattern(solid)			lcolor("220   0   0")	lwidth(thick)		///
			||	line KBpos2		 Ggov, connect(l) lpattern(dash)			lcolor(" 49 105  85")	lwidth(medthick)	///
			||	line KBpos3		 Ggov, connect(l) lpattern(longdash_dot)	lcolor(" 49 105  85")	lwidth(medthick)	///
			graphregion(color(white))	///
	ytitle("Crowd-out, | d{it:g{sub:i}{sup:*}}/d{it:G{sub:-i}} {sub:|d{it:G{sub:-i}} = d{it:w{sub:i}}} | " " ")		///	
	ylabel(0(.2)1, angle(horizontal)) xlabel(0 10 20 40 60 80 100)	///
	xtitle("Giving-by-others ({it:G{sub:-i}})")		legend(position(7))	///
	yline(0,  lstyle(grid) lcolor(gs8) )	///
	xline(10, lstyle(grid) lcolor(gs8) )	///
	note("Notes: {it:U} = {it:{&alpha}} log({it:G}) + {it:{&beta}} log({it:g{sub:i}}) + (1 - {it:{&alpha}} - {it:{&beta}}) log({it:x{sub:i}}). Income is held constant at   {it:w{sub:i}} = \$`y'.") ///
	legend(off)	///
	text(.23 80 "{it:{&alpha}} = `alpha1'0,	{it:{&beta}} = `beta1'0", placement(e) size(small))	///
	text(.07 80 "{it:{&alpha}} = `alpha2',	{it:{&beta}} = `beta2'"	, placement(e) size(small))	///
	text(.57 80 "{it:{&alpha}} = `alpha3',	{it:{&beta}} = `beta3'"	, placement(e) size(small))	



	
	
	log close	log1
	set more on
