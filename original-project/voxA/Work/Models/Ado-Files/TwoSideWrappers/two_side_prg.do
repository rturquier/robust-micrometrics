program two_side,eclass
version 10.1

	if replay() {
		display "Replay not implemented"
	}
	else {
		syntax varlist(min=3) [if] [in] [, absloss details bootstrap vc]
		local nvars: word count of `varlist'
		local depvar: word 1 of `varlist'
		local regs: list varlist - depvar
		local nvars: word count of `regs'


		marksample touse

		local variable_censoring=0
		if ("`vc'"=="vc") local variable_censoring=1
		local estmet=2
		if ("`absloss'"=="absloss") local estmet=1
		local printdetails=1
		if !("`details'"=="details") local printdetails=0
		local boots=0
		if ("`bootstrap'"=="bootstrap") local boots=1


		display " "
		display as text _dup(78) "-"
		display " "
		display "TwoSide: Estimation of Panel Data Regression Models with Two-Sided Censoring"
		display "         Based on Alan, Honoré, Hu and Leth-Petersen and (2011)"
		display "         Version 0.2"
		display " "
		display as text _dup(78) "-"
		display " "

		tempname b V nn ch2 
		mata: two_sidemata("`varlist'","`touse'",`estmet',`printdetails', `boots', `variable_censoring')
		mat b = r(beta)
		mat V= r(V)
		matrix `nn'= r(N)
		matrix `ch2'= r(chi2)
		

 		local nvars=colsof(b)
		local  ddd=""
		local i
		forvalues i=1/`nvars' {
		local expvar: word `i' of `regs'
		local ddd="`ddd' `expvar'" 
		}
		matrix colnames b= `ddd'  
		matrix colnames V= `ddd'   
		matrix rownames V= `ddd' 

		local i =`nvars'+1
		local ddd: word `i' of `regs'		
		local N = `nn'[1,1]

		ereturn clear	
		ereturn post b V, depname(`depvar') obs(`N') esample(`touse')
		ereturn local cmdline `"`0'"'
		ereturn local cmd "bop"

		display " "
		display "Number of Obs  :" %16.0f `nn'[1,1] %22s " " %-10s "Unit ID:" %14s "`ddd'"
		display "Number of Units:" %16.0f `nn'[2,1] %22s " " %-13s "ChiSq:"%11.2f `ch2'[1,1]
		display "Number of Units Used:" %11.0f `nn'[6,1] %22s " " %-18s "Prob > ChiSq:"%6.4f `ch2'[2,1]
		display "Avg No. of Obs per Unit:  " %6.2f `nn'[3,1] %22s " " %-20s "Frac. Cen. Above:"%4.2f `nn'[7,1]
 		display "Min No. of Obs for a Unit:" %6.0f `nn'[4,1] %22s " " %-20s "Frac. Cen. Below:"%4.2f `nn'[8,1]
		display "Max No. of Obs for a Unit:" %6.0f `nn'[5,1]
		display " "
		ereturn display 
		if (`estmet'==1 & `boots'==0) {
		display " " 
		display "WARNING: Std. Err. calculated using a crude automatic bandwidth selection."
		display "         The bootstrap is preferable for final version (but it is slow)."
		}	
		if !(`boots'==0) {
		display " " 
		display "NOTE: Standard errors calculated using bootstrap."
		}	
		display " "
		display as text _dup(78) "="
		display " "
		display " "


	}
end


