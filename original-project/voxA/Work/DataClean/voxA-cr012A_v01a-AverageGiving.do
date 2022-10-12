//	file:
		local Namedo "voxA-cr012A_v01a-AverageGiving"
		
//	task:	Calculate the average giving per participant.

//	author:
		local aut "\ mow \ 2016-10-03"   

//	Include file to assign directories:
		voxAdir DataClean
		include "voxA-cr-dirInclude-001A.do"
		
// Input datasets:
	local FileIn	 "voxA-cr001A_v01a-DecisionsWIDE"		// h04, ...

	
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
use	"`dirIN'\\`FileIn'", clear

	describe


	
	
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (1) Average giving across the decisions.
*
*

gen		hAvg = (h04 + h10 + h28 + h34 + h04_y46 + h28_y46)/6
lab var hAvg   "Avg giving, six decisions"
note	hAvg : Average giving across the six decisions \ Old name "h_avg" \ `tag'



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Save
*
keep id newid hAvg

order id newid hAvg
	
	
	count	
		local N = r(N)
	
	label data "VOX expr: Average per participant \ 2016-10-03"
	notes: `Namedo'.dta \ Sample size N = `N' participants \ `tag'

	
	save 	"`dirOUTNamedo'" , replace

	
	unab 		 list: _all
	describe 	`list'
	notes 		`list'
	labelbook	


	log close	log1
	set more on
