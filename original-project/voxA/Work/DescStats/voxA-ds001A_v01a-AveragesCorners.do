//	file:
		local Namedo "voxA-ds001A_v01a-AveragesCorners"
		
//	task:	Summary statistics.

//	author:
		local aut "\ mow \ 2017-06-20"   

//	Include file to assign directories:
		voxAdir DescStats
		include "voxA-ds-dirInclude-001A.do"
		
// Input datasets:
	local FileIn01	 "voxA-cr021A_v01a-DecisionsLONG"
	
	local FileIn02	 "voxA-cr011A_v01a-BoundaryDecisions"
	

	
// Notes:	(1) The work of the present do-file was previously done in
//					"vox-cr001A-v01a-Individuals.do"
//
//					Variables from that previous file that are NOT re-created in 
//						the present file:
//
//						lab var h_avgInc40   "Avg. across the 4 decisions with Income = $40"



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (0) Get the data.
*
tempfile Decisions510	Corners


*** The 85 participants x 6 decisions.
*
use	"`dirIN'/`FileIn01'", clear

	describe
	
	save "`Decisions510'"



*** The corner decisions (n = 85).
*
use	"`dirIN'/`FileIn02'", clear

	describe
	
	save "`Corners'"

	


	
	
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (1) Average giving across the decisions.
*
use "`Decisions510'", clear

sum h



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (2) Corners
*

*** At the 85 x 6 level
*
count if (h == Lower) | (h == Upper)
	scalar Ncorner = r(N)
	
count
	scalar NNN = r(N)

di "Percent of `NNN' decisions at a corner = " 100*round(Ncorner/NNN , .001)




*** At the 85 level (per participant)
*
use "`Corners'", clear

tab AlwaysZero, miss

tab AlwaysMax, miss







	log close	log1
	set more on
