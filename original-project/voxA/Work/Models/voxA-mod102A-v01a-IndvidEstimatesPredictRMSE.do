//	file: 
		local Namedo "voxA-mod102A-v01a-IndvidEstimatesPredictRMSE"
		
//	task:	 Predicts giving at each of the six budget decisions, 
//				using each individual's alpha, beta (heterogeneous).
//				Cobb-Douglas specifications for each individual participant.
//

//	author:
		local aut "\ mow \ 2017-06-20"   

//	Include file to assign directories:
		voxAdir Models
		include "voxA-mod-dirInclude-001A.do"
		
// Input datasets:
	local FileIn01	"voxA-cr001A_v01a-DecisionsWIDE"		// h(b), Ggov(b), ... b = 1, 2, 3, 4, 5, 6.

	local FileIn02	"voxA-mod101A-v01a-IndvidEstimates-2017-06-20--16-38"	// alpha, beta, i = 1,...85.



		
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (0) Initialization.
*

tempfile DecisionsWIDE 		alpha85


*** (0.1) Get the decision data (Long).
*
use	"`dirIN'/`FileIn01'", clear			

	save `DecisionsWIDE'
	

*** (0.2) Get the estimates (n = 85).
*
use	"`dirEST'/`FileIn02'", clear
	
		save `alpha85'
	
		
*** (0.3) Merge
*
use `DecisionsWIDE', clear

merge 1:1	newid	using	`alpha85', assert(match) nogenerate	







* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Set-up the result matrices, and some vectors
*
tempname Predict Actual

local Nrows = 6
local Ncols = 7
matrix Predict = J(`Nrows',`Ncols',.)

matrix colnames Predict = Ggov y Actual Predict DiffAvgPrdt Mse Rmse
matrix rownames Predict = Bgt1 Bgt2 Bgt3 Bgt4 Bgt5 Bgt6


*** Budgets 1 - 6: Initialize the matrix
*
	forvalues i = 1/6 {
		matrix Predict[`i', 1]	= `i'
		matrix Predict[`i', 2]	= 40
}
matrix Predict[1, 1]	=  4
matrix Predict[2, 1]	= 10
matrix Predict[3, 1]	= 28
matrix Predict[4, 1]	= 34

matrix Predict[5, 1]	=  4
matrix Predict[5, 2]	= 46

matrix Predict[6, 1]	= 28
matrix Predict[6, 2]	= 46

matrix list Predict




* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Actual averages at each budget
*
local i = 0

foreach x in 04 10 28 34 {
	local ++i
		sum h`x'	if !missing(alpha)
			matrix Predict[`i', 3]	= r(mean)		
}

sum	h04_y46			if !missing(alpha)	// Budget 5
	matrix Predict[5, 3]	= r(mean)
	
sum	h28_y46			if !missing(alpha)	// Budget 6
	matrix Predict[6, 3]	= r(mean)
	
	
matrix list Predict,	format(%7.2f)








* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Predictions
*
local i = 0

foreach x in 04 10 28 34 {
	local ++i
	
		predictCDhetero alpha beta Ggov`x' y`x', gen(h`x'hat)

		label var	h`x'hat "Predicted giving at Budget `i'"
		note		h`x'hat : Predictions based on alpha, beta estmates from the "icd001" vector \ 001 version of the individual Cobb-Douglas /*
										*/ estimates \ Missing means that the alpha, beta estimates themselves were missing \ `tag'

		gen			diff`x'sq = (h`x' - h`x'hat)^2
		label var	diff`x'sq "Sqr. diff actual - predict"
		
	sum h`x'hat		if !missing(alpha)
		matrix Predict[`i', 4]	= r(mean)
		
		
	sum diff`x'sq	if !missing(alpha)			// mse for this budget
		matrix Predict[`i', 6]	= r(mean)
}

			
		
			
			
*** Budget 5: Ggov =  $4, y = $46 
*			
predictCDhetero alpha beta Ggov04_y46 y04_y46, gen(h04_y46hat)

label var	h04_y46hat "Predicted giving at Budget 5"
note		h04_y46hat : Predictions based on alpha, beta estmates from the "icd001" vector \ 001 version of the individual Cobb-Douglas /*
								*/ estimates \ Missing means that the alpha, beta estimates themselves were missing \ `tag'

gen			diff04_y46sq = (h04_y46 - h04_y46hat)^2

label var	diff04_y46sq "Sqr. diff actual - predict"

sum	h04_y46hat			if !missing(alpha)	// Budget 5
	matrix Predict[5, 4]	= r(mean)

sum	diff04_y46sq		if !missing(alpha)	// mse for this budget
	matrix Predict[5, 6]	= r(mean)
	
	

*** Budget 6: Ggov = $28, y = $46 
*			
predictCDhetero alpha beta Ggov28_y46 y28_y46, gen(h28_y46hat)

label var	h28_y46hat "Predicted giving at Budget 6"
note		h28_y46hat : Predictions based on alpha, beta estmates from the "icd001" vector \ 001 version of the individual Cobb-Douglas /*
								*/ estimates \ Missing means that the alpha, beta estimates themselves were missing \ `tag'

gen			diff28_y46sq = (h28_y46 - h28_y46hat)^2

label var	diff28_y46sq "Sqr. diff actual - predict"

sum	h28_y46hat			if !missing(alpha)	// Budget 6
	matrix Predict[6, 4]	= r(mean)

sum	diff28_y46sq		if !missing(alpha)	// mse for this budget
	matrix Predict[6, 6]	= r(mean)
	

	
	
* Check	
*gsort - diff04sq
*list h04 h04hat diff04sq



	
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Diff(ActualAvg - PredictedAvg) and RMSE for each Budget
*
matrix temp = 	Predict[1...,3] - Predict[1...,4]	// For each Budget: ActualAvg - PredictedAvg
	matrix list temp

matrix Predict[1,5] = 	temp


forvalues i = 1/`= rowsof(Predict)' { 				// For each Budget: rmse
	matrix Predict[`i', 7] = sqrt(Predict[`i', 6])
} 

matrix list Predict,	format(%7.2f)



*** RMSE taking, for each Budget, the ActualAvg - PredictedAvg as the "error", and then
*
*		taking the average of the SIX squared-errors: (ActualAvg - PredictedAvg)^2
*
*		Note: This is the squared-error considering the Average of the i = 1,...,78 predictions for each Budget as a single prediction,
*				comparing that prediction to the Actual Average of the i = 1,...,78 "h"s (the amounts actually given), squaring that difference,
*				and then taking the average of the SIX squared-differences. And then taking the square-root of that, 
scalar diffSq = 0
forvalues i = 1/6 {
	scalar temp		= Predict[`i', 5]

	scalar diffSq	= diffSq + (temp)^2
}

	scalar mseREP = diffSq/6	// The average of this variable is SUM, b = 1 to 6 of [   (AVG (i = 1,..., 78) h(Budget b) - AVG (i = 1,..., 78) h(Budget b)_HAT)^2   ]

	di "Mean-Square-Error averaged across the SIX Budgets (each Budget taken as one representative individual) = " %7.2f mseREP

	scalar rmseREP = sqrt(mseREP)
	di "RMSE (each Budget taken as one representative individual) = " %7.2f rmseREP
	 




* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* RMSE for the n * 6 individual predictions
*
gen diffsqAllSix = diff04sq + diff10sq + diff28sq + diff34sq + diff04_y46sq + diff28_y46sq

label variable	diffsqAllSix "Sqr. diff for all 6 Budgets"
note			diffsqAllSix : For each participant this is the (Actual h_i - Predicted h_i)-Squared then summed over the six Budgets \ /*
								*/ Predictions based on alpha, beta estmates from the "icd001" vector \ 001 version of the individual Cobb-Douglas /*
								*/ estimates \ Missing means that the alpha, beta estimates themselves were missing \ /*
								*/ Divide the average of this variable by 6 to get the Root-Mean-Square-Error for all N * 6 decisions \ `tag'

sum diffsqAllSix
	scalar mse = r(mean)/6		// The average of this variable is SUM, i = 1 to 78 of [ (h04_i - h04_iHAT)^2 + ... (h28_y46i - h28_y46iHAT)^2]
								//	divided by 78.  Have to divide this by 6 so that the sum is, in the end, divided by 78 * 6

	di "Mean-Square-Error for all 78 * 6 decisions = " %7.2f mse 

	scalar rmse = sqrt(mse)
	di "RMSE for all 78 * 6 decisions = " %7.2f rmse
	

	
	
	
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* Predicted crowd-out
*

*** (1) Just the four "unfunded" budgets.
*
matrix CrowdOut = Predict[1..4,1..7]
matrix colnames CrowdOut = Ggov y Actual Predict COutActual COutPredict COutDiff

matrix CrowdOut[1,5] = 999
matrix CrowdOut[1,6] = 999
matrix CrowdOut[1,7] = 999

forvalues i = 2(1)4 {
	local iLag = `i' - 1
	matrix CrowdOut[`i',5] = CrowdOut[`i',3] - CrowdOut[`iLag',3]
	matrix CrowdOut[`i',6] = CrowdOut[`i',4] - CrowdOut[`iLag',4]
	matrix CrowdOut[`i',7] = CrowdOut[`i',5] - CrowdOut[`i',6]
}

matrix list CrowdOut, format(%7.2f)


* RMSE in the Crowd-out
scalar rmseCO = (	(CrowdOut[2,7])^2	+ (CrowdOut[3,7])^2	+ (CrowdOut[4,7])^2	)/3

di in red "RMSE of crowd-out predictions (unfunded crowd-out only) = " rmseCO




*** (2) Now add the two "balanced-budget" budgets
*
matrix CrowdOut = Predict[1...,1..7]
matrix colnames CrowdOut = Ggov y Actual Predict COutActual COutPredict COutDiff

matrix CrowdOut[1,5] = 999
matrix CrowdOut[1,6] = 999
matrix CrowdOut[1,7] = 999

* The four "unfunded" budgets (just repeating the calculation done about 20 lines earlier)
forvalues i = 2(1)4 {
	local iLag = `i' - 1
	matrix CrowdOut[`i',5] = CrowdOut[`i',3] - CrowdOut[`iLag',3]
	matrix CrowdOut[`i',6] = CrowdOut[`i',4] - CrowdOut[`iLag',4]
	matrix CrowdOut[`i',7] = CrowdOut[`i',5] - CrowdOut[`i',6]
}

* Now, the two "balanced" budgets
matrix CrowdOut[5,5] = CrowdOut[2,3] - CrowdOut[5,3]
matrix CrowdOut[6,5] = CrowdOut[4,3] - CrowdOut[6,3]

matrix CrowdOut[5,6] = CrowdOut[2,4] - CrowdOut[5,4]
matrix CrowdOut[6,6] = CrowdOut[4,4] - CrowdOut[6,4]


* ... and the Actual-Predicted Crowd-outs
matrix CrowdOut[5,7] = CrowdOut[5,5] - CrowdOut[5,6]
matrix CrowdOut[6,7] = CrowdOut[6,5] - CrowdOut[6,6]

matrix list CrowdOut, format(%7.2f)


* RMSE in the Crowd-out
scalar rmseCO = (	(CrowdOut[2,7])^2	+ (CrowdOut[3,7])^2	+ (CrowdOut[4,7])^2	+ (CrowdOut[5,7])^2		+ (CrowdOut[6,7])^2		)/5

di as result "RMSE of crowd-out predictions (this time include the balanced-budgets too) = $" round(rmseCO, .01)




	log close	log1
	set more on


