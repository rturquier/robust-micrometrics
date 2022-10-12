*! margCrowdOutTwoSide3 v01a  MOW 2013-07-27
program margCrowdOutTwoSide3, rclass
	version 12.1
	
*	Calculates the crowd-out-marginal-effect in a K=3 variable model (Low govt. provision, High govt. provision)

	args b1 b2 b3 b0 Sigma
	
	tempname Ggov GovtHigh Lower Upper xbeta

	*** Estimate E[y | Ggov = $4, GovtHigh = 0, Upper = $46] 			Budget 5
	*
	tempname Ey5
	
	scalar `Ggov'		=  4
	scalar `GovtHigh'	=  0
	scalar `Lower'		=  0
	scalar `Upper'		= 46

	scalar `xbeta'	= `b1' * `Ggov'   +   `b2' * `Ggov' * `GovtHigh'   +   `b3' * `GovtHigh'   +   `b0'


	EyTwoSide `Lower' `Upper' `xbeta' `Sigma'

	scalar `Ey5' = r(Ey)

	
	*** Estimate E[y | Ggov = $10, GovtHigh = 0, Upper = $40] 			Budget 2
	*
	tempname Ey2
	
	scalar `Ggov'		= 10
	scalar `GovtHigh'	=  0
	scalar `Lower'		=  0
	scalar `Upper'		= 40

	scalar `xbeta'	= `b1' * `Ggov'   +   `b2' * `Ggov' * `GovtHigh'   +   `b3' * `GovtHigh'   +   `b0'


	EyTwoSide `Lower' `Upper' `xbeta' `Sigma'

	scalar `Ey2' = r(Ey)

	
	
	*** Estimate E[y | Ggov = $28, GovtHigh = 1, Upper = $46] 			Budget 6
	*
	tempname Ey6
	
	scalar `Ggov'		= 28
	scalar `GovtHigh'	=  1
	scalar `Lower'		=  0
	scalar `Upper'		= 46

	scalar `xbeta'	= `b1' * `Ggov'   +   `b2' * `Ggov' * `GovtHigh'   +   `b3' * `GovtHigh'   +   `b0'


	EyTwoSide `Lower' `Upper' `xbeta' `Sigma'

	scalar `Ey6' = r(Ey)
	
	
	
	*** Estimate E[y | Ggov = $34, GovtHigh = 1, Upper = $40] 			Budget 4
	*
	tempname Ey4
	
	scalar `Ggov'		= 34
	scalar `GovtHigh'	=  1
	scalar `Lower'		=  0
	scalar `Upper'		= 40

	scalar `xbeta'	= `b1' * `Ggov'   +   `b2' * `Ggov' * `GovtHigh'   +   `b3' * `GovtHigh'   +   `b0'


	EyTwoSide `Lower' `Upper' `xbeta' `Sigma'

	scalar `Ey4' = r(Ey)

	
	
	*** Marginal effects
	*
	tempname coLow coHigh coChange
	
	scalar `coLow'		= (`Ey2' - `Ey5') / 6		// Crowd-out at LOW  provision
	scalar `coHigh'		= (`Ey4' - `Ey6') / 6		// Crowd-out at HIGH provision
	scalar `coChange'	= `coHigh' - `coLow'		// Change, High to Low


	
return scalar coLow		= `coLow'
return scalar coHigh	= `coHigh'
return scalar coChange	= `coChange'

end
