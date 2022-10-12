//	file: 
		local Namedo "voxA-mod005A-v01a-CobbDougRepresentative"
		
//	task:	Estimates the Cobb-Douglas alpha and beta, based on the 
//			"positive square root" quadratic solution to the FOC.

//	author:
		local aut "\ mow \ 2016-11-14"   

//	Include file to assign directories:
		voxAdir Models
		include "voxA-mod-dirInclude-001A.do"
		
// Input datasets:
	local FileIn01	"voxA-cr001A_v01a-DecisionsWIDE"		// h1-h6, Ggov1-Ggov6, ...

	
// Notes:	(1) The Random-Effects-Tobit that allows a subject to:
//
//				1. make some choices at the lower corner (e.g., h = $0)
//
//				2. make some choices at the upper corner (e.g., h = $40 and or h = $46)
//
//				3. BUT DOES NOT HANDLE a subject who makes some choices 
//					at the lower corner and other choices at the upper corner.
//
//					[There was no participant who hit $0 for some decisions,
//						then hit the upper corner ($40 or $46) for another 
//						decision.]
//
//			(2) The variables containing the allocations to the child---h1, h2,
//					h3, h4, h5, h6---MUST MUST MUST be ordered so that h5 and h6
//					are the decisions taken with the income = $46.
//


		
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (0.1) Get the data.
*
use	"`dirIN'/`FileIn01'", clear			

keep id h1 h2 h3 h4 h5 h6 Ggov1 Ggov2 Ggov3 Ggov4 Ggov5 Ggov6 ///
		z1 z2 z3 z4 z5 z6 y1 y2 y3 y4 y5 y6 

sum y1 y2 y3 y4 y5 y6, sep(4)






* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (1) Call "tobitCDRElowup.ado"
*
gen sroot = 0


* Note: mdraws will ignore these seeds because I am using Halton sequences and NOT using
*		hrandom and shuffle options.
set seed 56135
gen long seedvar = int((uniform() + 5 - _n)*100000000) in 1/5
local seed1 = seedvar[1]	// 439578138
local seed2 = seedvar[2]	// 326710563
local seed3 = seedvar[3]
local seed4 = seedvar[4]
local seed5 = seedvar[5]

di %15.0g `=`seed1''
di %15.0g `=`seed2''
di %15.0g `=`seed3''
di %15.0g `=`seed4''
di %15.0g `=`seed5''


set more off



*** Impure altruism model.
*
*	Note: When all six decisions are used, h5 and h6 (the decisions made with y = $46)
*			must come last in the list:
*
*			tobitCDRElowup (h1 = Ggov1 z1) (h2 = Ggov2 z2) (h3 = Ggov3 z3) ///
*				(h4 = Ggov4 z4) (h5 = Ggov5 z5) (h6 = Ggov6 z6) , dr(300) an burn(30)
*

cap program drop tobitCDRElowup6			// Six-decision version (standard)

version 8.2

tobitCDRElowup6 (h1 = Ggov1 z1) (h2 = Ggov2 z2) (h3 = Ggov3 z3) ///	dr(300) - takes 26 iterations and a long time to converge
	(h4 = Ggov4 z4) (h5 = Ggov5 z5) (h6 = Ggov6 z6), adoonly 	///
	dr(300) an burn(30)											//	dr( 20) - takes 50 iterations and much less time to converge



	
	log close	log1

	set more on
