//	file:	"voxA-mod-dirInclude-001A.do"
		
//	task:	Tells do-files in this directory where to find/write *.dta files and *.log files

//	Drive:
	local Drive "C"				// Change drive on which do-files are to be run, if necessary (e.g., "Q").
	
	local Root "`Drive':\Users\mowilhel\Documents\"

	
// Specify directories:

*** (1) Log
*
clear all
local dirLog		"`Root'\\voxA\Work\Models\Log Files"
	
	if "`Namedo'" ~= "voxA-mod199A_v01a-Master" {
		capture log close log1
		log using 			"`dirLog'\\`Namedo'.log", replace name(log1)
	}
	else {
		capture log close master
		log using 			"`dirLog'\\`Namedo'.log", replace name(master)	
	}


			
*** (2) Data directories
*
local dirIN			"`Root'\\voxA\Work\Datasets\Derived" 
local dirOUT		"`Root'\\voxA\Work\Datasets\Derived" 

local dirBOOT		"`Root'\\voxA\Work\Models\Bootstrap Files" 

local dirMAT		"`Root'\\voxA\Work\Models\Matrices" 

local dirEST		"`Root'\\voxA\Work\Models\Estimates" 



*** (3) ado directory for the Two-side wrappers, marginal effect calculators and
*			do-files that must be pre-run before calling "twoside.ado"
*
local dirTWOSIDE	"`Root'\\voxA\Work\Models\Ado Files\TwoSideWrappers"

local dirCOBB		"`Root'\\voxA\Work\Models\Ado Files"



// Routine opening tasks
	version 13.0
	
	adopath + "`dirTWOSIDE'"

*	adopath + "`Drive':\ado\personal\Honore"

	adopath + "`dirCOBB'"

	adopath + "`Root'\\ado"

	
	
	
	set maxvar 9999

	set more off

	local dirOUTNamedo "`dirOUT'\\`Namedo'"

	local dirESTNamedo "`dirEST'\\`Namedo'"

	local tag "`Namedo'.do  `aut'"
