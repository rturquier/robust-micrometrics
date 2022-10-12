//	file:
		local Namedo "voxA-cr000A_v01a-Summary_File"
		
//	task:	Read in the source, create a unique participant id, and save as a Stata dataset.

//	author:
		local aut "\ mow \ 2017-06-20"   

//	Include file to assign directories:
		voxAdir DataClean
		include "voxA-cr-dirInclude-001A.do"
		
// Input datasets:
	local FileIn01	"April_data_2007.txt"			// Copy of the original VOX Decision data
	local FileIn02	"April_2007_survey_data.txt"	// Copy of the original VOX Survey response data
	
	local FileIn03	"April_data_order-question-2015.txt"	// Copy of the original VOX "order-of-the-decisions" data
	
	local NewIDin	"newid"							// A simpler id variable.
	
	local FileIn05	"gnpsXXXXXXXXXXXXXXXXXXXXXXX"	// GINPS data


// Notes:
//
	
	
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (0) Get the data.
*
tempfile 	Decisions 		Survey 		Order	smallApril

*** Decisions
*
insheet using "`dirSourceIN'/`FileIn01'", clear names

	save `Decisions'
	
	
*** Survey responses
*
insheet using "`dirSourceIN'/`FileIn02'", clear names
	
	save `Survey'
	

*** Order in which the participant faced the decisions.
*
insheet using "`dirSourceIN'/`FileIn03'", clear names
	
	save `Order'


*** Dataset that contains "newid"
*
*	Note: I have lost track of the provenance of "newid", but 
*		I know that "small_April_2007.dta" has it.
*
*		"newid" is just a sequential 1-to-85 id number that is easier 
*			to work with than is "id".  But ther e is a one-to-one mapping 
*			between "id" and "newid".
*
use "`dirSourceIN'/`NewIDin'", clear

keep id newid

note newid : I got this variable from "small_April_2007.dta" \ /*
		*/ Other than saying that, I have lost track of the /*
		*/ provenance of "newid" \ I made this note on 2017-06-20 \ `tag'
		
label variable newid "Simple sequential numbering"

	save `smallApril'
	


* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (1) Decision data.
*
use `Decisions', clear

	describe

keep session id h4 h10 h28 h34 h4_46 h28_46

rename	h4		h04
rename	h4_46	h04_y46

rename	h28_46	h28_y46


foreach bbb in 04 10 28 34 04_y46 28_y46 {

	label variable	h`bbb' "Decision, amt given"
	note			h`bbb' : SOURCE variable / Amount given by participant /*
						*/ faced with indicated budget (Ggov, y) / `tag'
}


*** Create a unique identifier for the participant.
*
rename id	idOriginal			// Keep a copy of the source-code "id"

label variable	idOriginal	"id created by the experimenters"
note 			idOriginal : SOURCE variable \ Clone of "id" \ `tag'


generate long id = 200700000 + session*1000 + idOriginal

label variable id "unique participant id"
note id : id = 200700000 + session*1000 + idOriginal \ `tag'

label variable session "Session person participated in"
note session : SOURCE variable \ Numerical version of "session" \ `tag'

order id session, first



save `Decisions', replace







* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (2) Survey response data.
*
use `Survey', clear

	describe

	
foreach x	of	 varlist *  {

	note `x' : SOURCE variable / `tag'
}





*** Create a unique identifier for the session, from "date" and "time".
*
generate session = .
replace  session = 1 if date=="4/11/2007" & time=="2:30"
replace  session = 2 if date=="4/11/2007" & time=="4:30"
replace  session = 3 if date=="4/12/2007" & time=="2:30"
replace  session = 4 if date=="4/12/2007" & time=="4:00"
replace  session = 5 if date=="4/13/2007" & time=="2:30"
replace  session = 6 if date=="4/13/2007" & time=="4:00"

assert !missing(session)

label variable session "Session person participated in"
note session : Numerical version of "session" \ Constructed from "date" and "time" \ `tag'




*** Create a unique identifier for the participant.
*
rename id	idOriginal			// Keep a copy of the source-code "id"

label variable	idOriginal	"id created by the experimenters"
note replace	idOriginal in 1 : SOURCE variable \ Clone of "id" \ `tag'


generate long id = 200700000 + session*1000 + idOriginal

label variable id "unique participant id"
note id : id = 200700000 + session*1000 + idOriginal \ `tag'




*** Rename the source variables so that the names are easier to read.
*
forvalues i = 1/21 {
	rename q`i' q1_`i'
}


forvalues i = 1/8 {
	local n = `i' + 25
	rename v`n' q2_`i'
}


forvalues i = 1/12 {
	local n = `i' + 33
	rename v`n' q3_`i'
}








order id session, first

save `Survey', replace





* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (3) Order in which the participant faced the decisions.
*
use `Order', clear

	describe
	

*** Create a unique identifier for the session, from the string version (sessionOriginal)
*
rename session	sessionOriginal 	// Keep a copy of the string-variable "session".

gen byte 	session = .
	replace session = 1 if sessionOriginal == "wedn 2:30"	// The string-to-number mapping is from
	replace session = 2 if sessionOriginal == "wedn 4:30"	// 	"April_data.xls"
	replace session = 6 if sessionOriginal == "fri 4:00"
	replace session = 3 if sessionOriginal == "thur 2:30"
	replace session = 4 if sessionOriginal == "thur 4:00"
	replace session = 5 if sessionOriginal == "fri 2:30"	


assert !missing(session)

label variable session "Session person participated in"
note session : Numerical version of "session" \ Constructed from the string version (sessionOriginal) /*
				*/ that came in "April_data_order-question-2015.txt" \ `tag'


	
*** Create a unique identifier for the participant.
*
rename id	idOriginal			// Keep a copy of the source-code "id"

label variable	idOriginal	"id created by the experimenters"
note replace	idOriginal in 1 : SOURCE variable \ Clone of "id" \ `tag'


generate long id = 200700000 + session*1000 + idOriginal

label variable id "unique participant id"
note id : id = 200700000 + session*1000 + idOriginal \ `tag'


order id session, first

keep id order				// This is all we need, and dropping everything else
							//	will prevent overwriting the survey responses in "Survey"
							//	with the missing data in "Order" (see "pc-vox03A_v02a-PC_Motives_DecisionOrder.do")

save `Order', replace



	



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (4) Merge
*
use `Decisions', clear

merge 1:1	id	using `Survey', assert(match) nogenerate

merge 1:1	id	using `Order', assert(match) nogenerate

order		order, after(h28_y46)



merge 1:1	id	using `smallApril', assert(match) nogenerate	// Just to get "newid"

order		newid, after(id)

	describe newid



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Save
*
	count	
		local N = r(N)
	
	label data "VOX expr: SOURCE data (April 2007) \ 2017-06-20"
	notes: `Namedo'.dta \ Sample size N = `N' participants \ `tag'

	
	save 	"`dirOUTNamedo'" , replace

	
	unab 		 list: _all
	describe 	`list'
	notes 		`list'
	labelbook	


	log close	log1
	set more on
