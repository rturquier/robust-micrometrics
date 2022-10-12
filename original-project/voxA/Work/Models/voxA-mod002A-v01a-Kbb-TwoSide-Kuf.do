//	file: 
		local Namedo "voxA-mod002A-v01a-Kbb-TwoSide-Kuf"
		
//	task:	The results in Table 2 of VOX (2017).
//				Change in balanced-budget crowd-out.
//				Change in unfunded crowd-out.
//				Test of pure warm glow.

//	author:
		local aut "\ mow \ 2016-10-04"   

//	Include file to assign directories:
		voxAdir Models
		include "voxA-mod-dirInclude-001A.do"
		
// Input datasets:
	local FileIn01	"voxA-cr021A_v01a-DecisionsLONG"		// h1-h6, Ggov1-Ggov6, ...

	
// Notes:	(1) VOX (2017), Table 2:
//
//				(col. 1) Kbb LO, Linear-Fixed effects
//				(col. 2) Kbb HI, 		"
//				(col. 3) Kbb LO, Two-sided Tobit FE (account for corners; Alan et al., 2011)
//				(col. 4) Kbb HI,										"
//				(col. 5) Kbb LO & HI in one model: Change in Kbb,	"
//				(col. 6) Kuf LO & HI in one model: Change in Kbb,	"
//
//			(2) Tests pure warm glow: see (6.5) below.
//


		
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
* (1) Kbb LO, Linear-Fixed effects
*
*		Note: Low govt provision:  $4 -> $10 (income $46 -> $40)
*
xtreg h Ggov if (Budget==2 | Budget==5), fe	vce(robust)		// Robust s.e. corresponds to assuming different variances at Low and High.

	test Ggov = -1
	
	di "One-sided test - Ho: Ggov >= 1, p-value = " r(p)/2	// Ho: crowd-out is -1 or a bigger negative number.





* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (2) Kbb HI, Linear-Fixed effects
*
*		Note: High govt provision:  $28 -> $34 (income $46 -> $40)
*
xtreg	h Ggov if (Budget==4 | Budget==6), fe	vce(robust)

	test Ggov = -1
	
	di "One-sided test - Ho: Ggov >= 1, p-value = " r(p)/2	// Ho: crowd-out is -1 or a bigger negative number.





* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (3) Kbb LO:	Two-sided Tobit FE 
*
*				[account for corners:
*						Alan, Honoré, Hu and Leth-Petersen (2014)]
*
preserve


*** (3.1) Keep only the variables needed for estimation, and
*
*		Decide which analysis to run: Balanced budget crowd-out - LOW Provision.	
*
keep id Budget h Ggov Ggov_GovtHigh GovtHigh	y y_GovtHigh	Lower Upper

keep if (Budget==2 | Budget==5)



*** (3.2) Point estimates
*
two_sideWrap1		 h Ggov			Lower Upper id Budget

	scalar ESTme	= r(me)		// Save for later hypothesis testing.



*** (3.3) Standard errors - bootstrapped
*
local dirBOOTNamedo "`dirBOOT'\\`Namedo'-KbbLO-v01a"

xtset , clear

bootstrap me = r(me), rep(1000) seed(81681)	///	R = 1000
				 cluster(id) idcluster(newid)	saving("`dirBOOTNamedo'", replace) nowarn : 		///
				 two_sideWrap1		 h Ggov  	Lower Upper newid Budget

				 
*** (3.4) Tests
*
use "`dirBOOTNamedo'", clear

sum me

* Ho: kb_Low <= -1
	sum me, meanonly
	
	gen double SamplDistHo = me + ( (-1)   -   (r(mean)) )		// Sampling distribution under Ho: kb_Low == -1
	sum SamplDistHo
	
	gen byte	Indicator = (SamplDistHo > ESTme)
	sum			Indicator							// p = .034
	drop		Indicator	SamplDistHo
	
	
	
restore






* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (4) Kbb HI:	Two-sided Tobit FE 
*
*
preserve


*** (4.1) Keep only the variables needed for estimation, and
*
*		Decide which analysis to run: Balanced budget crowd-out - HIGH Provision.	
*

keep id Budget h Ggov Ggov_GovtHigh GovtHigh	y y_GovtHigh	Lower Upper

keep if (Budget==4 | Budget==6)



*** (4.2) Point estimates
*
two_sideWrap1		 h Ggov			Lower Upper id Budget

	scalar ESTme	= r(me)		// Save for later hypothesis testing.



*** (4.3) Standard errors - bootstrapped
*
local dirBOOTNamedo "`dirBOOT'\\`Namedo'-KbbHI-v01a"

xtset , clear

bootstrap me = r(me), rep(1000) seed(15388495)	/// R = 1000
				 cluster(id) idcluster(newid)	saving("`dirBOOTNamedo'", replace) nowarn : 		///
				 two_sideWrap1		 h Ggov  	Lower Upper newid Budget

				 
*** (4.4) Tests
*
use "`dirBOOTNamedo'", clear

sum me

* Ho: kb_High <= -1
	sum me, meanonly
	
	gen double SamplDistHo = me + ( (-1)   -   (r(mean)) )		// Sampling distribution under Ho: kb_Low == -1
	sum SamplDistHo
	
	gen byte	Indicator = (SamplDistHo > ESTme)
	sum			Indicator							// p = .390
	drop		Indicator	SamplDistHo
	
	
	
restore






* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (5) Kbb LO & HI in one model: Change in Kbb				 
*
*			Two-sided Tobit FE
*
preserve


*** (5.1) Keep only the variables needed for estimation, and
*
*		Decide which analysis to run: Balanced budget crowd-out ---
*
*					LOW and HIGH Provision.	
*

keep id Budget h Ggov Ggov_GovtHigh GovtHigh	y y_GovtHigh	Lower Upper

keep if (Budget==2 | Budget==5 | Budget==4 | Budget==6)



*** (5.2) Point estimates
*
two_sideWrap3		 h Ggov Ggov_GovtHigh GovtHigh 	Lower Upper id

	scalar ESTcoLow	= r(coLow)		// Save these for later hypothesis testing.
	scalar ESTcoHigh	= r(coHigh)
	scalar ESTcoChange = r(coChange)



*** (5.3) Standard errors - bootstrapped
*
local dirBOOTNamedo "`dirBOOT'\\`Namedo'-KbbLO_and_HI-v01a"

xtset , clear

bootstrap coLow = r(coLow)	coHigh = r(coHigh)		coChange = r(coChange), rep(1000) seed(901699)	/// R = 1000
				 cluster(id) idcluster(newid)	saving("`dirBOOTNamedo'", replace) nowarn : 		///
				 two_sideWrap3		 h Ggov Ggov_GovtHigh GovtHigh 	Lower Upper newid

*** (5.4) Tests
*
use "`dirBOOTNamedo'", clear

sum coLow coHigh coChange

* (5.4.1) Ho: kb_Low <= -1
	sum coLow, meanonly
	
	gen double SamplDistHo = coLow + ( (-1)   -   (r(mean)) )		// Sampling distribution under Ho: kb_Low == 1
	sum SamplDistHo
	
	gen byte	Indicator = (SamplDistHo > ESTcoLow)
	sum			Indicator											// p = .477
	drop		Indicator	SamplDistHo
	
	
* (5.4.2) Ho: kb_High <= -1
	sum coHigh, meanonly
	
	gen double SamplDistHo = coHigh + ( (-1)   -   (r(mean)) )		// Sampling distribution under Ho: kb_High == 1
	sum SamplDistHo
	
	gen byte	Indicator = (SamplDistHo > ESTcoHigh)
	sum			Indicator											// p = .013
	drop		Indicator	SamplDistHo
	
	
* (5.4.3) Ho: kb_High <= kb_Low      <==>     kb_High - kb_Low  <= 0
	sum coChange, meanonly
	
	gen double SamplDistHo = coChange + ( (0)   -   (r(mean)) )		// Sampling distribution under Ho: kb_High - kb_Low == 0
	sum SamplDistHo, detail
	
	gen byte	Indicator = (SamplDistHo > ESTcoChange)
	sum			Indicator											// p = .07
	drop		Indicator	SamplDistHo

	
	
restore	







* *** *** *** *** *** *** *** *** *** *** *** *** *** *** *	
*
* (6) Kuf LO & HI in one model: Change in Kuf
*
*			Two-sided Tobit FE
*
preserve


*** (6.1) Keep only the variables needed for estimation, and
*
*		Decide which analysis to run: UNFUNDED crowd-out ---
*
*					LOW and HIGH Provision.	

keep id Budget h Ggov Ggov_GovtHigh GovtHigh	y y_GovtHigh	Lower Upper

keep if (Budget==1 | Budget==2 | Budget==3 | Budget==4)



*** (6.2) Point estimates
*
two_sideWrap3		 h Ggov Ggov_GovtHigh GovtHigh 	Lower Upper id

	scalar ESTcoLow		= r(coLow)		// Save these for later hypothesis testing.
	scalar ESTcoHigh	= r(coHigh)
	scalar ESTcoChange	= r(coChange)

	

*** (6.3) Standard errors - bootstrapped
*
local dirBOOTNamedo "`dirBOOT'\\`Namedo'-KufLO_and_HI-v01a"

xtset , clear

bootstrap coLow = r(coLow)	coHigh = r(coHigh)		coChange = r(coChange), rep(1000) seed(3306)	/// R = 1000
				 cluster(id) idcluster(newid)	saving("`dirBOOTNamedo'", replace) nowarn : 		///
				 two_sideWrap3		 h Ggov Ggov_GovtHigh GovtHigh 	Lower Upper newid

				 
*** (6.4) Tests
*
use "`dirBOOTNamedo'", clear

sum coLow coHigh coChange

* (6.4.1) Ho: k_Low <= -1
	sum coLow, meanonly
	
	gen double SamplDistHo = coLow + ( (-1)   -   (r(mean)) )		// Sampling distribution under Ho: k_Low == 1
	sum SamplDistHo
	
	gen byte	Indicator = (SamplDistHo > ESTcoLow)
	sum			Indicator											// p < .001
	drop		Indicator	SamplDistHo
	
	
* (6.4.2) Ho: k_High <= -1
	sum coHigh, meanonly
	
	gen double SamplDistHo = coHigh + ( (-1)   -   (r(mean)) )		// Sampling distribution under Ho: k_High == 1
	sum SamplDistHo
	
	gen byte	Indicator = (SamplDistHo > ESTcoHigh)
	sum			Indicator											// p < .001
	drop		Indicator	SamplDistHo
	
	
* (6.4.3) Ho: k_High <= k_Low      <==>     k_High - k_Low  <= 0
	sum coChange, meanonly
	
	gen double SamplDistHo = coChange + ( (0)   -   (r(mean)) )		// Sampling distribution under Ho: k_High - k_Low == 0
	sum SamplDistHo, detail
	
	gen byte	Indicator = (SamplDistHo > ESTcoChange)
	sum			Indicator											// p = .013
	drop		Indicator	SamplDistHo

	
	
*** (6.5) Testing pure warm glow
*

* (6.5.1) Ho: k_Low == 0
	sum coLow, meanonly
	
	gen double SamplDistHo = coLow + ( (0)   -   (r(mean)) )		// Sampling distribution under Ho: k_Low == 0
	sum SamplDistHo
	
	gen byte	Indicator = (SamplDistHo < ESTcoLow)				// Under Ho, the sampling distribution goes from -.307 to + .222 (centered on zero crowd-out)
	sum			Indicator											// Under this sampling distribution, how many repetitions are at the observed -.64 or more negative?
	drop		Indicator	SamplDistHo								// p < .001
	
	

	
* (6.5.2) Ho: k_High == 0
	sum coHigh, meanonly
	
	gen double SamplDistHo = coHigh + ( (0)   -   (r(mean)) )		// Sampling distribution under Ho: k_High == 0
	sum SamplDistHo
	
	gen byte	Indicator = (SamplDistHo < ESTcoHigh)				// Under this sampling distribution (centered on zero crowd-out), 
	sum			Indicator											// 	how likely to see an estimate of -.42 or more negative?   p < .001
	drop		Indicator	SamplDistHo

restore

	
	log close	log1

	set more on
