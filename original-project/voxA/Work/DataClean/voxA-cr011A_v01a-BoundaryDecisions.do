//	file:
		local Namedo "voxA-cr011A_v01a-BoundaryDecisions"
		
//	task:	Creates the dummy variable indicators of decisions that hit
//				a corner at $0 or a corner at $40/$46

//	author:
		local aut "\ mow \ 2016-10-03"   

//	Include file to assign directories:
		voxAdir DataClean
		include "voxA-cr-dirInclude-001A.do"
		
// Input datasets:
	local FileIn	 "voxA-cr001A_v01a-DecisionsWIDE"		// h04, ...

		
// Notes:	(1) Some of the work of the present do-file was previously done in
//					"vox-cr001A-v01a-Individuals.do"
//
//					Variables from that previous file that are NOT re-created in 
//						the present file:
//
//						label var NoTrunc69 "Not truncated Budgets 2 and 4"
//
//						label var	NoTruncBBCOLO "No double corners at budgets 2,5"
//						label var	NoTruncBBCOHI "No double corners at budgets 4,6"
//						label var	  NoTruncBBCO "No double corners at budgets 2,5 and 4,6"
//
//						label var	NoTruncUnfndLO "No double corners at budgets 1,2"
//						label var	NoTruncUnfndHI "No double corners at budgets 3,4"
//						label var	NoTruncUnfnd "No double corners at budgets 1,2 and 3,4"
//
//						String variables to summarize the corner decisions (CornerLO, CornerHI)
//
//
//				(2) And some of the work was previously done in
//						"- vox-VariableCreation001A-v01a-Include.do"
//
//								(Ncensorlow, Ncensorup)
//
//					Variables from that previous file that are NOT re-created in 
//						the present file:
//
//						label variable SomeOver40 "This person made some decisions > $40 (strictly >)"
//



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (0) Get the data.
*
use	"`dirIN'\\`FileIn'", clear

	describe


	
	
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (1) Boundary decisions.
*
*

*** NOT truncated at any of the six budgets
*
gen byte NoTrunc = 0
replace	 NoTrunc = 1	if ((0 < h04) & (h04 < 40)) & ((0 < h34) & (h34 < 40))	///
						 & ((0 < h28) & (h28 < 40)) & ((0 < h34) & (h34 < 40))	///
				& ((0 < h04_y46) & (h04_y46 < 46)) & ((0 < h28_y46) & (h28_y46 < 46))
				
label var	NoTrunc "Decision never truncated"
notes 		NoTrunc : Decision never hits $0 or $max Budgets 1-6 \ Old name "NoTrunc66" \ `tag'

label define NoTruncLBL	0 "0YesTrunc"		1 "1NvrTrunc"
label values NoTrunc	NoTruncLBL	

tab NoTrunc, miss



*** At a corner for every decision
*
gen byte	AlwaysZero = (h04 == 0 & h10==0 & h28==0 & h34==0 & h04_y46==0 & h28_y46==0)
label var	AlwaysZero "Every decision = $0"
note		AlwaysZero : Decision always at $0 \ `tag'

gen byte	AlwaysMax = (h04 == 40 & h10==40 & h28==40 & h34==40 & h04_y46==46 & h28_y46==46)
label var	AlwaysMax "Every decision = Maximum ($40, $46)"
note		AlwaysMax : Decision always at $Maximum \ `tag'

tab AlwaysZero, miss
tab AlwaysMax, miss


gen byte	AlwaysCorn = 0
replace		AlwaysCorn = -1 if AlwaysZero == 1
replace		AlwaysCorn = 46 if AlwaysMax  == 1

label var 	AlwaysCorn "Every decision at a corner"
note		AlwaysCorn : Decision always at $0 or $Maximum \ `tag'

	label define AlwaysCornLBL	-1 "-1alwyZero"		46 "46alwyMax"
	label values AlwaysCorn	 AlwaysCornLBL

tab AlwaysCorn, miss



*** Participant hit a Zero corner on one decision, and a Max corner at another decision.
*
gen byte ZeroMaxCorn = 	(h04 ==  0 | h10 ==  0 | h28 ==  0 | h34 ==  0 | h04_y46 ==  0 | h28_y46 ==  0) &	///
						(h04 == 40 | h10 == 40 | h28 == 40 | h34 == 40 | h04_y46 == 46 | h28_y46 == 46)

label var	ZeroMaxCorn "Hit Zero AND Max Corner"
notes		ZeroMaxCorn : Participants who chose $0 for at least one decision AND /*
						*/ chose $Max for at least one other decision \ /*
						*/ The structural Cobb-Douglas program cannot handle a participant /*
						*/ who hits both the lower and upper corners \ `tag'

tab ZeroMaxCorn, miss






* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (2) Indicators of censoring previously created in
*
*			"- vox-VariableCreation001A-v01a-Include"
*

*** At $0
*
gen byte NcensorZero = (h1 == 0)

forvalues i = 2/6 {
	replace NcensorZero = NcensorZero + 1 if h`i' == 0
}

label var	NcensorZero "Num decisions h = 0"

note		NcensorZero : Number of decisions at $0 \ Old name "Ncensorlow" \ `tag'

tab NcensorZero, miss


*** At $40 or $46
*
gen byte NcensorMax = (h1 == 40)

forvalues i = 2/4 {
	replace NcensorMax = NcensorMax + 1 if h`i' == 40
}
forvalues i = 5/6 {
	replace NcensorMax = NcensorMax + 1 if h`i' == 46
}

label variable NcensorMax "Number of decisions h = 40/46"

note	NcensorMax : Number of decisions a $Max \ Old name "Ncensorup" \ `tag'

tab NcensorMax, miss






*** Force a stop is any participant hit both a Zero and a Max corner.
*
tab NcensorZero NcensorMax, miss

tab ZeroMaxCorn, miss

	assert !(NcensorZero >= 1 & NcensorMax >= 1)	// No one has some decisions at $0 and other decisions at $40/$46

	assert !(ZeroMaxCorn == 1)						// Double-check
	




* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Save
*
keep  id newid ZeroMaxCorn NoTrunc AlwaysZero AlwaysMax AlwaysCorn NcensorZero NcensorMax

order id newid ZeroMaxCorn NoTrunc AlwaysZero AlwaysMax AlwaysCorn NcensorZero NcensorMax


	
	count	
		local N = r(N)
	
	label data "VOX expr: Boundary-Corner dummies \ 2016-10-03"
	notes: `Namedo'.dta \ Sample size N = `N' participants \ `tag'

	
	save 	"`dirOUTNamedo'" , replace

	
	unab 		 list: _all
	describe 	`list'
	notes 		`list'
	labelbook	


	log close	log1
	set more on

