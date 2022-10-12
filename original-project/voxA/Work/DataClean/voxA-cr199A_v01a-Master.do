//	file: 
		local Namedo "voxA-cr199A_v01a-Master"
		
//
//	task:	Runs multiple data cleaning.
//
//	author:
		local aut "\ mow \ 2017-06-20"
		
//	Include file to assign directories:
		voxAdir DataClean
		include "voxA-cr-dirInclude-001A.do"
		
// Input datasets: 
 

* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Open the master log-file
*
capture log close _all

log using 		"`dirLog'/`Namedo'.log", replace name(master)


* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Do-files
*

do voxA-cr000A_v01a-Summary_File
do voxA-cr001A_v01a-DecisionsWIDE
do voxA-cr011A_v01a-BoundaryDecisions
do voxA-cr012A_v01a-AverageGiving
do voxA-cr021A_v01a-DecisionsLONG

   

	log close	master
	set more on
