//	file:
		local Namedo "voxA-cr021A_v01a-DecisionsLONG"
		
//	task:	Re-shape the data into LONG (participant x 6-decisions).
//			Create dummies for high govt provision, income.
//			Create interaction variables: Ggov_GovtHigh, y_GovtHigh.
//			Generate "Lower" and "Upper" censoring levels for each Budget.

//	author:
		local aut "\ mow \ 2016-10-04"   

//	Include file to assign directories:
		voxAdir DataClean
		include "voxA-cr-dirInclude-001A.do"
		
// Input datasets:
	local FileIn01	"voxA-cr001A_v01a-DecisionsWIDE"		// h04, ...
	
	local FileIn02	"voxA-cr011A_v01a-BoundaryDecisions"	// For checking purposes (only).
	
	
// Notes:	(1) The work of the present do-file was previously done in
//					"- vox-VariableCreation001A-v01a-Include.do"
//



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (0) Get the data.
*
tempfile Decisions	Boundary

use	"`dirIN'\\`FileIn01'", clear

	describe

	keep id newid	h1-h6 	Ggov1-Ggov6 	y1-y6 	z1-z6 	G1-G6	

	save `Decisions'
	
	
use	"`dirIN'\\`FileIn02'", clear

	describe

	save `Boundary'

	
	
	
	
	
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (1) Reshape for panel models
*
*
use `Decisions', clear

reshape long h Ggov y z G, i(id) j(Budget)

xtset id Budget

label variable Budget "Budget set for this decision"

	note Budget: Table 1 in the VOX (2016) paper \ `tab'

label define BudgetLBL 1 "1-4_40"		2 "2-10_40"		3 "3-28_40"		4 "4-34_40"		5 "5-4_46"		6 "6-28_46"
label values Budget BudgetLBL

* Re-label
	label variable Ggov		"Government provision"
	label variable z		"Social income (y + Ggov)"
	label variable y		"Own income"
	label variable h		"Decision"
	label variable G		"Public good level (Ggov + h)"



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (2) Dummy variables for "high" provision, income.
*
*
gen byte GovtHigh = 0
replace  GovtHigh = 1  if (Budget==3 | Budget==4 | Budget==6)

label variable  GovtHigh "Govt provision is $28 or $34"
note 			GovtHigh : Indicates the three decisions where G-i is high \ Budgets 3, 4, 6 \ `tag'

label define GovtHighLBL 1 "1high"		0 "0low"
label values GovtHigh GovtHighLBL


gen byte IncHigh = 0
replace  IncHigh = 1  if (Budget==5 | Budget==6)

label variable  IncHigh "Income is $40 or $46"
note 			IncHigh : Indicates the two decisions where y_i is high \ Budgets 5, 6 \ `tag'

label values IncHigh GovtHighLBL




* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (3) Interaction variables
*
*
gen y_GovtHigh = y * GovtHigh

label variable  y_GovtHigh "Income interacted w/ Govt = $28 or $34"
note			y_GovtHigh : y * GovtHigh \ `tag'


gen Ggov_GovtHigh = Ggov * GovtHigh

label variable  Ggov_GovtHigh "Govt provision interacted w/ Govt = $28 or $34"
note			Ggov_GovtHigh : Ggov * GovtHigh \ `tag'





* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (4) Lower and upper censoring levels.
*
*
*		Note: These variables are need for some estimation programs,
*				such as random-effects Tobit.
*
generate Lower = 0		// Lower censoring point is always $0
generate Upper = y		// Upper censoring point = own income ($40 or $46)

lab var Lower "Censoring level, lower"
lab var Upper "Censoring level, upper"

note	Lower : Lower level at which this observation was censored ("cornered") \ /*
					*/ Typically $0 \ `tag'
					
note	Upper : Upper level at which this observation was censored ("cornered") \ /*
					*/ Typically $40 or $46 \ `tag'
					
					
tab Budget Upper, miss


tempfile Long
save	`Long'





* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (5) Checking
*
merge m:1	id		using `Boundary', assert(match) nogen

list id Budget GovtHigh IncHigh Ggov y		Ggov_GovtHigh 	y_GovtHigh in 1/24, sep(6) label

list id Budget Ggov y h NcensorZero NcensorMax if NcensorZero==1, sep(6) label

list id Budget Ggov y h NcensorZero NcensorMax if NcensorZero==3, sep(6) label	

list id Budget Ggov y h NcensorZero NcensorMax if NcensorZero==6, sep(6) label	

list id Budget Ggov y h NcensorZero NcensorMax if NcensorMax==1, sep(6) label

list id Budget Ggov y h NcensorZero NcensorMax if NcensorMax==4, sep(6) label	

list id Budget Ggov y h NcensorZero NcensorMax if NcensorMax==6, sep(6) label	




* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Save
*
use `Long', clear

order newid, after(id)
	
	
	count	
		local N = r(N)
	
	label data "VOX expr: Decision data \ LONG format \ 2016-10-04"
	notes: `Namedo'.dta \ Sample size N = `N' participant-decision observations \ `tag'

	
	save 	"`dirOUTNamedo'" , replace

	
	unab 		 list: _all
	describe 	`list'
	notes 		`list'
	labelbook	


	log close	log1
	set more on
