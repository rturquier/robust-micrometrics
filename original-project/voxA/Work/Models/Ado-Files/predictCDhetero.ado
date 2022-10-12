*! predictCDhetero	v01a	mow		2013-11-07
program predictCDhetero, rclass
	version 10.1
	
	syntax varlist(min=4 max=4 numeric) [if] [in] [, INCEFF Generate(string)]
	
//	tasks:	For vectors of "alpha" and "beta" (Cobb-Douglas parameters), giving by others, and own income, 
//				this ado-file predicts CD impure altruism:
//					(a) h - level of giving
//					(b) q1 - own income effect
//					(c) q2 - income effect in response to giving by others
//
//	Inputs:	(1) "alpha" - altruism parameter (for individual "i")
//			(2) "beta" - warm glow parameter
//			(3) "Ggov" is the Giving by others.
//			(4) "w" is own income.
//
//
//	Options:	inceff   - request prediction of income effects (q1, q2).
//				generate - select variable names for the predictions (defaults = "ghat", "q1hat", "q2hat").
//
//
//	Notes:	See "predictCDlevel.ado", a similar ado-file that accepts "alpha" and "beta" as scalars and predicts level of giving (h)
//				for several different budgets (Ggov and w values) at a time -- designed to be called for ONE individual at a time.	

marksample touse
	
	quietly count if `touse'
		if `r(N)' == 0 {
			error 2000
		}

		
*** Check if user has defined names for the new variables.
*
tokenize `"`generate'"'
local wc3 : word count `generate'

if "`wc3'" == "1" {
		confirm new var `1'
		local  X `"`1'"'
	}
	else if "`wc3'" == "3" {
		confirm new var `1'
		confirm new var `2'
		confirm new var `3'
		local  X  `"`1'"'
		local  Q1 `"`2'"'
		local  Q2 `"`3'"'
	}
	else if "`wc3'" != "0" {
		di as error "Specify either ONE (g_hat) or THREE (g_hat, q1_hat, q2_hat) variable names"
		error 102
}



* . . . check to see if user did not ask for income effects, but then supplied variable names for the income effects.
if ("`inceff'" != "inceff") & ("`wc3'" == "3") {
		di as error "You did not ask for predictions of income effects, but then supplied names for the predictions of income effects"
		error 197
}
	

	


*** Tokenize the input variable list
*		
tempvar alpha beta Ggov w

tokenize `varlist'

quiet {
	gen double `alpha' = `1'			// Assign meaningful variable names to the (contents of) the positional macros.
	gen double  `beta' = `2'
	gen double  `Ggov' = `3'
	gen double     `w' = `4'
}



*** Predict level of giving - "+" solution
*
quiet {
	tempvar z B Hplus gplus

	gen double `z'		= `Ggov' + `w'		// 	Social income

	gen double `B'		= (1 - `beta')*`Ggov' + (`alpha' + `beta')*`z'	
	gen double `Hplus'	= .5*(`B' + sqrt(`B'^2 - 4*`alpha' * `Ggov' * `z') )
	gen double `gplus'	= `Hplus' - `Ggov'
}
	
	
	
*** Predict income effects (q1, q2) - "+" solution
*
*
quiet {
	if "`inceff'" == "inceff" {
		tempvar Sroot q1POS q2POS

		gen	double `Sroot' = sqrt( (1 - `alpha')^2 * `Ggov'^2	 + 2*( (`beta' - `alpha') + (`alpha' + `beta')*`alpha' ) * `Ggov' * `w'  +  (`alpha' + `beta')^2 * `w'^2		)	
	
		gen double `q1POS' = .5 * ( (`alpha' + `beta')		+      ( ((`beta' - `alpha') + (`alpha' + `beta')*`alpha' ) * `Ggov'   +  (`alpha' + `beta')^2 * `w') / `Sroot'	)

		gen double `q2POS' = .5 * ( (1 - `beta')			+      ( ((1 - `alpha') - (1 + `alpha')*`beta' ) * `Ggov'   +  ( (`beta' - `alpha') - (`alpha' + `beta')*`beta' ) * `w') / `Sroot'	)
	}
}	
	
	


*** Return data vectors
*
quiet {
	if "`wc3'" == "1" {
			gen double `X' = `gplus'
				if "`inceff'" == "inceff" {				// To get here, user supplied a name for predicted giving, and asked for income effects,
					confirm new var q1hat				//	 but did not supply names for the pedicted income effects.
					confirm new var q2hat
					
					gen double q1hat = `q1POS'
					gen  q2hat = `q2POS'
				}
		}
		else if "`wc3'" == "3" {						// To get here, user supplied names for predicted giving and the income effects.
			gen double `X'  = `gplus'
			gen double `Q1' = `q1POS'
			gen double `Q2' = `q2POS'
			}
		else {											// To get here, user did NOT supply names for new variables. . .
			confirm new var ghat
			gen double ghat  = `gplus'
				if "`inceff'" == "inceff" {				// . . . and (if arrived here) also asked for income effects.
					confirm new var q1hat
					confirm new var q2hat
					
					gen double q1hat = `q1POS'
					gen double q2hat = `q2POS'
				}			
		}
}														// End of the quite.



*** Return in r()
*
quiet{
	sum `X'
		return scalar N = r(N)				// Number of valid (not missing) predictions of giving
		
	return local varname `varlist'		// Variables that were in the calling statement
}



end
