*! margCrowdOutTwoSide2hi v01a  MOW 2013-08-31
program margCrowdOutTwoSide2hi, rclass
	version 12.1
	
*	Calculates the own income marginal effect (q1; Budgets 3, 6), unfunded crowd-out-marginal-effect (-1 + q1 + q2; Budgets 3, 4),
*		the warm glow piece (q2), and "q1 + q2"     in a K=2 variable model estimated at high provision (y, high govt. provision).
*																						 ^^^^
*																						 ^^^^
*																						 ^^^^

	args b1 b2 b0 Sigma
	
	tempname y Ggov Lower Upper xbeta

	
	
	
	
	*** Estimate E[h | y = $40, Ggov = $34, Upper = $40] 			Budget 4
	*
	tempname Eh4
	
	scalar `y'			= 40
	scalar `Ggov'		= 34
	scalar `Lower'		=  0
	scalar `Upper'		= 40

	scalar `xbeta'	= `b1' * `y'   +   `b2' * `Ggov'    +   `b0'


	EyTwoSide `Lower' `Upper' `xbeta' `Sigma'

	scalar `Eh4' = r(Ey)


	
	*** Estimate E[h | y = $40, Ggov = $28, Upper = $40] 			Budget 3
	*
	tempname Eh3
	
	scalar `y'			= 40
	scalar `Ggov'		= 28
	scalar `Lower'		=  0
	scalar `Upper'		= 40

	scalar `xbeta'	= `b1' * `y'   +   `b2' * `Ggov'   +   `b0'


	EyTwoSide `Lower' `Upper' `xbeta' `Sigma'

	scalar `Eh3' = r(Ey)
	
	

	*** Estimate E[h | y = $46, Ggov = $28, Upper = $46] 			Budget 6
	*
	tempname Eh6
	
	scalar `y'			= 46
	scalar `Ggov'		= 28
	scalar `Lower'		=  0
	scalar `Upper'		= 46

	scalar `xbeta'	= `b1' * `y'   +   `b2' * `Ggov'   +   `b0'


	EyTwoSide `Lower' `Upper' `xbeta' `Sigma'

	scalar `Eh6' = r(Ey)

	
	
	
	*** Marginal effects
	*
	tempname incHigh coHigh q1plusq2High q2High
	
	scalar `incHigh'		= (`Eh6' - `Eh3') / 6			// Income effect (q1) at HIGH provision
	scalar `coHigh'			= (`Eh4' - `Eh3') / 6			// Unfunded crowd-out (-1 + q1 + q2) at HIGH provision
	scalar `q1plusq2High'	=  `coHigh' + 1					// q1 + q2
	scalar `q2High'			=  `q1plusq2High' - `incHigh'	// q2


	
return scalar incHigh		= `incHigh'
return scalar coHigh		= `coHigh'
return scalar q1plusq2High	= `q1plusq2High'
return scalar q2High		= `q2High'

end
