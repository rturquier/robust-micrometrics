//	file: 
		local Namedo "voxA-mod003A-v01a-IncomeEffects-TwoSide"
		
//	task:	Estimate the change in income effects that we discuss in 
//				VOX (2017), footnote 16.

//	author:
		local aut "\ mow \ 2016-10-08"   

//	Include file to assign directories:
		voxAdir Models
		include "voxA-mod-dirInclude-001A.do"
		
// Input datasets:
	local FileIn01	"voxA-cr021A_v01a-DecisionsLONG"		// h1-h6, Ggov1-Ggov6, ...

	

		
* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (0.1) Get the data.
*
use	"`dirIN'\\`FileIn01'", clear			




* *** *** *** *** *** *** *** *** *** *** *** *** *** ****** *** *** *	
*
* (0.2) Run the preliminary programs that must be run prior to estimation
*
mata							// Evoke mata so that the libraries necessary to run pantob or two_side can be found.
	mata mlib index
end


quietly do "`dirTWOSIDE'\\two_side_mata.do"		// Programs that must be run before calling "two_side"

quietly do "`dirTWOSIDE'\\two_side_prg.do"




* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (1) Income effect LO, Linear-Fixed effects
*
*		Note: Low govt provision:  $4 -> $10 (income $46 -> $40)
*
xtreg h y if (Budget==1 | Budget==5), fe
		
	matrix b = e(b)
	
	matrix list b
	
	matrix IncLow = b["y1", "y"]
	
	scalar IncLow = IncLow[1,1]
	
	di "Income effect at Low Govt-provision ($4) = " round(IncLow, .01)



* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (2) Income effect HI, Linear-Fixed effects
*
*		Note: High govt provision:  $28 -> $34 (income $46 -> $40)
*
xtreg h y if (Budget==3 | Budget==6), fe
		
	matrix b = e(b)
	
	matrix list b
	
	matrix IIncHigh = b["y1", "y"]

	scalar IIncHigh = IIncHigh[1,1]
	
	di "Income effect at High Govt-provision ($28) = " round(IIncHigh, .01)


			
		
			


* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (3) Income effects LO & HI in one model: Change in the income effect
*
*			Two-sided Tobit FE
*

*** (3.1) Keep only the variables needed for estimation, and
*
*		Decide which analysis to run: Income effects ---
*
*					LOW and HIGH Provision.	
*
preserve

keep id Budget h Ggov Ggov_GovtHigh GovtHigh	y y_GovtHigh	Lower Upper

keep if (Budget==1 | Budget==5 | Budget==3 | Budget==6)



*** (3.2) Point estimates
*
two_sideWrap3		 h y y_GovtHigh GovtHigh 	Lower Upper id

	scalar ESTcoLow		= r(coLow)		// Save these for later hypothesis testing.
	scalar ESTcoHigh	= r(coHigh)
	scalar ESTcoChange = r(coChange)



*** (3.3) Standard errors - bootstrapped
*
local dirBOOTNamedo "`dirBOOT'\\`Namedo'-IncomeEffectLO_and_HI-v01a"

xtset , clear

bootstrap coLow = r(coLow)	coHigh = r(coHigh)		coChange = r(coChange), rep(1000) seed(10597219)	/// R = 1000
				 cluster(id) idcluster(newid)	saving("`dirBOOTNamedo'", replace) nowarn : 		///
				 two_sideWrap3		 h y y_GovtHigh GovtHigh 	Lower Upper newid

				 
				 
				 
*** (3.4) Tests
*
use "`dirBOOTNamedo'", clear

sum coLow coHigh coChange	// Income effect at Low, at High, and the Change (ignore the "co" in the name.

	
* Ho: income effect_High == income effect_Low      <==>     kb_High - kb_Low  == 0
	sum coChange, meanonly
	
	gen double SamplDistHo = coChange + ( (0)   -   (r(mean)) )		// Sampling distribution under Ho: kb_High - kb_Low == 0
	sum SamplDistHo, detail
	
	gen byte	Indicator = (SamplDistHo > ESTcoChange) | (SamplDistHo < - ESTcoChange)
	sum			Indicator											// p = .58
	drop		Indicator	SamplDistHo



restore





	
	log close	log1

	set more on
