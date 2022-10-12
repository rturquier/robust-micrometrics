//	file:
		local Namedo "voxA-ds101A_v01a-IndividAlphasBetas"
		
//	task:	Summary statistics fr the individual-level alphas and betas.

//	author:
		local aut "\ mow \ 2017-06-20"   

//	Include file to assign directories:
		voxAdir DescStats
		include "voxA-ds-dirInclude-001A.do"
		
// Input datasets:
	local FileIn01	 "voxA-mod101A-v01a-IndvidEstimates-2017-06-20--16-38"
	



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (0) Get the data.
*
use	"`dirEST'\\`FileIn01'", clear

	describe

	
	
	
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (1) Preference-types of the participants:
*
*		Impure altruism, Pure altruism, Pure warm glow.
*
tab ModelSelect, miss

count if (alpha < beta) & (ModelSelect == 1)	//  n = 4

	list newid alpha beta if (alpha < beta) & (ModelSelect == 1)


* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (2) Relative weight of Warm glow (relative to generosity):
*
*				 beta/(alpha + beta)
*
generate g = alpha + beta

generate Bsum = beta/g

sum Bsum, detail	





	log close	log1
	set more on
