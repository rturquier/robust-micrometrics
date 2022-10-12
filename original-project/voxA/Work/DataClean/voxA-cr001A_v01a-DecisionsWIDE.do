//	file:
		local Namedo "voxA-cr001A_v01a-DecisionsWIDE"
		
//	task:	Create a full set of Decision data, and the parallel Clone variables
//				with suffixes 1-6 that are easier to "xtset".
//			All in WIDE format.

//	author:
		local aut "\ mow \ 2016-09-22"   

//	Include file to assign directories:
		voxAdir DataClean
		include "voxA-cr-dirInclude-001A.do"
		
// Input datasets:
	local FileIn	 "voxA-cr000A_v01a-Summary_File"

		
// Notes:
//
	
	
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (0) Get the data.
*
use	"`dirIN'\\`FileIn'", clear

	describe

keep id newid h04 h10 h28 h34 h04_y46 h28_y46



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (1) Create the source variables that are otherwise
*		implicit.
*

*** Government provision.
*
generate Ggov04     =  4
generate Ggov10     = 10
generate Ggov28     = 28
generate Ggov34     = 34
generate Ggov04_y46 =  4
generate Ggov28_y46 = 28

*** Own income.
*
generate y04     = 40
generate y10     = 40
generate y28     = 40
generate y34     = 40
generate y04_y46 = 46
generate y28_y46 = 46

*** Social income.
*
generate z04     = y04      +  Ggov04
generate z10     = y10      +  Ggov10
generate z28     = y28      +  Ggov28
generate z34     = y34      +  Ggov34
generate z04_y46 = y04_y46  +  Ggov04_y46
generate z28_y46 = y28_y46  +  Ggov28_y46

*** Decision: Public good level.
*
generate G04     =  4 + h04
generate G10     = 10 + h10
generate G28     = 28 + h28
generate G34     = 34 + h34
generate G04_y46 =  4 + h04_y46
generate G28_y46 = 28 + h28_y46


*** Label and note all the source variables
*
foreach bbb in 04 10 28 34 04_y46 28_y46 {
	label variable Ggov`bbb'	"Government provision"
	label variable z`bbb'		"Social income (y + Ggov)"
	label variable y`bbb'		"Own income"
	label variable G`bbb'		"Public good level (Ggov + h)"
	
	note Ggov`bbb'	: Source \ `tag'
	note z`bbb'		: Source \ `tag'
	note y`bbb'		: Source \ `tag'
	note G`bbb'		: Source \ `tag'
}



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (2) Some parallel clones that we can "xtset"
*
local i = 0

foreach bbb in 04 10 28 34 04_y46 28_y46 {
	local ++i
	
	clonevar Ggov`i'	=	Ggov`bbb'
	clonevar    y`i'	=	   y`bbb'
	clonevar    z`i'	=	   z`bbb'
	clonevar    h`i'	=	   h`bbb'
	clonevar    G`i'	=	   G`bbb'
}


*** Notes indicating which Budget, 1 - 6.
*
foreach x in Ggov y z h G {
	notes replace `x'1 in 1: Budget 1 \ Ggov = 4, y = 40 \ Clone of `x'04 \ `tag'
}

foreach x in Ggov y z h G {
	notes replace `x'2 in 1 : Budget 2 \ Ggov = 10, y = 40 \ Clone of `x'10 \ `tag'
}

foreach x in Ggov y z h G {
	notes replace `x'3 in 1 : Budget 3 \ Ggov = 28, y = 40 \ Clone of `x'28 \ `tag'
}

foreach x in Ggov y z h G {
	notes replace `x'4 in 1 : Budget 4 \ Ggov = 34, y = 40 \ Clone of `x'34 \ `tag'
}

foreach x in Ggov y z h G {
	notes replace `x'5 in 1 : Budget 5 \ Ggov = 4, y = 46 \ Clone of `x'04_y46 \ `tag'
}

foreach x in Ggov y z h G {
	notes replace `x'6 in 1 : Budget 6 \ Ggov = 28, y = 46 \ Clone of `x'28_y46 \ `tag'
}






* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (3) Fix up the order: Source variables first, parallel clones second.
*
order h1 h2 h3 h4 h5 h6, after(G28_y46)
order Ggov1 Ggov2 Ggov3 Ggov4 Ggov5 Ggov6, after(h6)
order y1 y2 y3 y4 y5 y6, after(Ggov6)
order z1 z2 z3 z4 z5 z6, after(y6)
order G1 G2 G3 G4 G5 G6, after(z6)





* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Save
*
	count	
		local N = r(N)
	
	label data "VOX expr: Decision data \ WIDE format \ 2016-09-22"
	notes: `Namedo'.dta \ Sample size N = `N' participants \ `tag'

	
	save 	"`dirOUTNamedo'" , replace

	
	unab 		 list: _all
	describe 	`list'
	notes 		`list'
	labelbook	


	log close	log1
	set more on
