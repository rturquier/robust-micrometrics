*! EyTwoSide v01a  MOW 2013-07-27
program EyTwoSide, rclass
	version 12.1
	
*	Calculates E(y | x) for a two-limit Tobit model (assuming errors ~ Normal).

	args Lower Upper xbeta Sigma
	
	tempname Ey

	scalar	`Ey' = 		`xbeta' * (   normal((`Upper' - `xbeta')/`Sigma')   -  normal((`Lower' - `xbeta')/`Sigma')   )		///
																															///
					+	`Sigma' * (   normalden((`Lower' - `xbeta')/`Sigma')   -  normalden((`Upper' - `xbeta')/`Sigma')  )	///
																															///
					+	`Upper' * (	  1   -   normal((`Upper' - `xbeta')/`Sigma')   )										///
																															///
					+	`Lower' * (	  normal((`Lower' - `xbeta')/`Sigma')   )
					
	return scalar Ey = `Ey'
end
