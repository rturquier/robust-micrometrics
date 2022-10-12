//	file: 
		local Namedo "voxA-fig019A-v01a-LargeEconomy-q1PLUSq2"
	
//	task:  For the representative Cobb-Douglas parameter estimates,
//			(alpha = .594, beta = .021), find the level of G-i at 
//			which the economy become "large" in the Andreoni (1988) and 
//			Ribar & Wilhelm (2002) sense that motives-at-the-margin are
//			essentially warm glow (only). 
 
//			

//	author:
		local aut "\ mow \ 2017-06-21"   

//	Include file to assign directories:
		voxAdir Figures
		include "voxA-fig-dirInclude-001A.do"
	

//	Note: "Large" is defined to mean q1 + q2 = .99

		
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* The set-up
*

*** The x-axis
*
local N = 20000

set obs `N'

generate Ggov = (_n - 1)/100	

label var Ggov "Giving by others"



*** Income (or government provision) held constant
*
local y = 40			// Income


*** Cobb-Douglas altruism (alpha), warm glow (beta) parameters
*
local alpha = .25		//	.594
local beta  = .25		//	.021

local alpha = .594
local beta  = .021



	
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Optimal giving
*
generate double  B =  (1 - `beta')*Ggov + (`alpha' + `beta')*(Ggov + `y')	
			
generate double Hplus = .5*(B + sqrt(B^2 - 4*`alpha' * Ggov * (Ggov + `y')) )
generate double gplus =  Hplus - Ggov

label variable gplus "+"




* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* q1 and q2 solutions from the Cobb-Douglas:
*
*	One solution (the one with a "+" in front of the square-root term).
*
gen	S = ( (1 - `alpha')^2 * Ggov^2	 + 2*( (`beta' - `alpha') + (`alpha' + `beta')*`alpha' ) * Ggov * `y'  +  (`alpha' + `beta')^2 * `y'^2		)	

gen	Sroot = sqrt(S)	
	
gen q1plusq2POS	= .5 * ( (1 + `alpha')		+      ((1 - `alpha')^2 * Ggov   +   ((`beta' - `alpha') + (`alpha' + `beta')*`alpha') * `y') / Sroot		)

gen q1POS		= .5 * ( (`alpha' + `beta')		+      ( ((`beta' - `alpha') + (`alpha' + `beta')*`alpha' ) * Ggov   +  (`alpha' + `beta')^2 * `y') / Sroot	)

gen q2POS		= .5 * ( (1 - `beta')			+      ( ((1 - `alpha') - (1 + `alpha')*`beta' ) * Ggov   +  ( (`beta' - `alpha') - (`alpha' + `beta')*`beta' ) * `y') / Sroot	)

gen dgdGgovPOS	=   -1 + q1plusq2POS

gen dgdyPOS		=   q1POS

gen Kpos	= abs(-1 + q1plusq2POS)
gen KBpos	= abs(-1 + q2POS)
	
label var Sroot			"Square-root term"
label var dgdGgovPOS	"Deriv wrspt G_{-i} POS"
label var dgdyPOS		"Deriv wrspt y POS"
label var Kpos			"Unfunded crowd-out"
label var KBpos			"Balanced-budget crowd-out"

label var q1POS "q1 positive"
label var q2POS "Difference between the two income effects"

label var q1plusq2POS "Income effect with respect to social income"


gen			Amargin = q1POS / (q1POS + q2POS)
label var	Amargin "Altruism index - marginal"





* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Listing
*
list Ggov gplus q1plusq2POS	in 150



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Convergence at G_{-i} --> Infinity
*

di as result "alpha = `alpha',   beta = `beta'   q1 converges to " round((`beta')/(1 - `alpha'), .001)
di as result "                             q2 converges to " round((1 - (`alpha' + `beta'))/(1 - `alpha'), .001)


sum q*


list Ggov	q1plusq2POS in 15354/15359, sep(3)	// q1 + q2 approx= .99 at Ggov = $153.55
local LrgEcon = 153.55
	
	
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Plots
*

*** Graph of q1, q2, and (q1 + q2)
*
graph twoway	line q1POS Ggov, connect(l) lpattern(solid)	lcolor("202 194 126")	lwidth(thick)				///
			||	line q2POS Ggov, connect(l) lpattern(solid)	lcolor(" 50  50 255")	lwidth(thick)				///
			||	line q1plusq2POS Ggov, connect(l) lpattern(solid)		lcolor("220   0   0")	lwidth(thick)	///
	title("Figure fn19. Cobb-Douglas {it:q{sub:1}} and {it:q{sub:2}} as functions of giving by others.", size(medium))			///
	ytitle("{it:q{sub:1}}, {it:q{sub:2}}")				///	
	xtitle("Giving by others ({it:G{sub:-i}})")			///
	ylabel(0(.2)1)	xlabel(50 100 `LrgEcon' 200)				///
														///
	legend(off)											///
	note("Notes: The Cobb-Douglas parameters are   {it:{&alpha}} = `alpha' and {it:{&beta}} = `beta'."	/// 
		 "            {it:U} = (1 - {it:{&alpha}} - {it:{&beta}}) log({it:c{sub:i}}) + {it:{&alpha}} log({it:G}) + {it:{&beta}} log({it:g{sub:i}}). Income is held constant at   {it:y{sub:i}} = \$`y'.")	///
	text(.68 90 "{it:q{sub:2}}"										, placement(e))		///
	text(.27 90 "{it:q{sub:1}}"										, placement(e))		///
	text(.93 90 "{it:q{sub:1}} + {it:q{sub:2}}"						, placement(e))		///
	xline(`LrgEcon', lpattern(dot) lcolor(gs10) lwidth(medium) )



	
	
	log close	log1
	set more on
