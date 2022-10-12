//	file:	"voxA-cr-dirInclude-001A.do"
		
//	task:	Tells do-files in this directory where to find/write *.dta files and *.log files

//	Drive:
	local Drive "C"				// Change drive on which do-files are to be run, if necessary (e.g., "Q").
	
	local Root "`Drive':\Users\mowilhel\Documents"
	
// Specify directories:

*** (1) Log
*
clear all
local dirLog		"`Root'\\voxA\Work\DataClean\Log Files"
	
	if "`Namedo'" ~= "voxA-cr199A_v01a-Master" {
		capture log close log1
		log using 			"`dirLog'\\`Namedo'.log", replace name(log1)
	}
	else {
		capture log close master
		log using 			"`dirLog'\\`Namedo'.log", replace name(master)	
	}


			
*** (2) Data directories
*
local dirSourceIN	"`Root'\\voxA\Work\Datasets\Source"

local dirIN			"`Root'\\voxA\Work\Datasets\Derived" 
local dirOUT		"`Root'\\voxA\Work\Datasets\Derived" 




// Routine opening tasks
	version 13.0
		
	adopath + "`Root'\\ado"

	set maxvar 9999

	set more off

	local dirOUTNamedo "`dirOUT'\\`Namedo'"
	
	local tag "`Namedo'.do  `aut'"
