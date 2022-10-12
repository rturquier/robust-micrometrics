*! margCrowdOutTwoSide1 v01a  MOW 2013-07-30
program margCrowdOutTwoSide1, rclass
	version 12.1
	
*	Calculates the crowd-out-marginal-effect in a K=1 variable model (Low govt. provision or High govt. provision)

	args b1 b0 Sigma G1 Lo1 Up1 G2 Lo2 Up2
	
	tempname xbeta

	di _newline(2)
	di as result "Marginal effect:	x = " `G2' "	Lower = " `Lo2' "		Upper = " `Up2'	"	to . . ."
	di as result "Marginal effect:	x = " `G1' "	Lower = " `Lo1'	"		Upper = " `Up1'
	di _newline(2)	
	
	foreach x in b1 b0 Sigma {
		di as result "`x' = " ``x''
	}
		
	
	
	*** Estimate E[y | Ggov = `G1',  Lower = `Lo1',  Upper = `Up1']
	*
	tempname Ey1
	
	scalar `xbeta'	= `b1' * `G1'     +     `b0'


	EyTwoSide `Lo1' `Up1' `xbeta' `Sigma'

	scalar `Ey1' = r(Ey)


	
	
	*** Estimate E[y | Ggov = `G2',  Lower = `Lo2',  Upper = `Up2']
	*
	tempname Ey2
	
	scalar `xbeta'	= `b1' * `G2'     +     `b0'


	EyTwoSide `Lo2' `Up2' `xbeta' `Sigma'

	scalar `Ey2' = r(Ey)

	
	
	*** Marginal effects
	*
	tempname me
	
	scalar `me'		= (`Ey2' - `Ey1') / 6		// marginal effect

	
return scalar me		= `me'

end
