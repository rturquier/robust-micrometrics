cap program drop tobitCDRElowup_ll
program define tobitCDRElowup_ll

   version 8.2
   args lnf alpha beta lnsigma1 lnsigma2

/***
	di _newline(2)
	di in red "alpha     = " `alpha'
	di in red "beta      = " `beta'
	di in red "lnsigma1  = " `lnsigma1'
	di in red "lnsigma2  = " `lnsigma2'
**************************************************/







	scalar alpha = `alpha'
	scalar beta  = `beta'


	  	/* -ml- is trying to check a beta < 0 */
	if beta  < 0   {					// If -ml- tries to set beta < 0
		scalar beta  = 0				// must keep beta >= 0 and
									// keep alpha <= min(Ggov/z)
		forvalues i = 1/$X_MVT_NOEQ {
			quietly summ Ggov`i' , meanonly
				scalar Ggov_temp = r(mean)
			quietly summ z`i' , meanonly
				scalar z_temp    = r(mean)
			scalar r`i' = Ggov_temp/z_temp
		}

		scalar alphamin = r1			// Find largest value alpha can take on
		forvalues i = 2/$X_MVT_NOEQ {
			if r`i' < alphamin  scalar alphamin = r`i'
		}

		di in red "beta set to ZERO"

		if alpha > alphamin  { 		// Re-set alpha if necessary
			scalar alpha = alphamin
			di in red "alpha re-set to = " alpha
		}
	}


	  	/* -ml- is trying to check an alpha < 0 */
	if alpha < 0   {
		scalar alpha = 0
		di in red "alpha set to ZERO"
	}



	  	/* -ml- is trying to check an (alpha, beta) pair where  alpha+beta > 1 */
	if alpha+beta > 1   {
		scalar D = alpha+beta
		scalar alpha = (alpha/D)
		scalar beta  =  (beta/D)
		di in red "alpha+beta > 1.  Reset   alpha = " alpha
		di in red "                          beta = " beta		
	}


local deb "*"



quietly {	


	  	/* Generate error terms and "xb_i" terms */
	  forvalues i = 1/$X_MVT_NOEQ {
		 tempvar e`i'  xb`i'
	  	 tempvar B`i'  sroot`i' Hplus`i' gplus`i' res`i'
	  	 generate double `B`i'' = (1 - beta)*Ggov`i' + (alpha + beta)*z`i'	
		 generate double `sroot`i'' = `B`i''^2 - 4*alpha * Ggov`i' * z`i'

*			 noi di "i = |`i'|"
*			 noi di "alpha = " alpha
*			 noi di "beta  = " beta
			 assert `sroot`i'' >= 0

*			 noi di "Now calculating Hplus"
	      generate double `Hplus`i'' = .5*(`B`i'' + sqrt(`sroot`i'') )
			 assert `Hplus`i'' ~= .

		 generate double `gplus`i'' =  `Hplus`i'' - Ggov`i'

		 generate double  `e`i'' = h`i' - `gplus`i''		// residual
		 generate double `xb`i'' =        `gplus`i''		// non-linear function
	  }



	/*	Make covariance matrix for unsorted data - use Cholesky decomposition */
	 tempname Sig
	 matrix `Sig' = J($X_MVT_NOEQ,$X_MVT_NOEQ,0)

		summ `lnsigma1', meanonly
		scalar sigma1temp = exp(r(mean))

		summ `lnsigma2', meanonly
		scalar sigma2temp = exp(r(mean))

		scalar rhotemp = sigma1temp^2/(sigma1temp^2 + sigma2temp^2)


	forvalues i = 1/$X_MVT_NOEQ {
		matrix `Sig'[`i',`i'] = sigma1temp^2 + sigma2temp^2
	}


	forvalues i = 1/$X_MVT_NOEQ {
		 forvalues j = `=`i'+1'/$X_MVT_NOEQ {

			matrix `Sig'[`i',`j'] = rhotemp * sqrt(`Sig'[`i',`i']) * sqrt(`Sig'[`j',`j'])

			matrix `Sig'[`j',`i'] = `Sig'[`i',`j']			
 		 }
	 }		 







	capture mat $X_MVT_C = cholesky(`Sig')
	if _rc != 0 {
		noi di in red "Warning: cannot do Cholesky factorization of covariance matrix &&&&&&&"
		noi matrix list `Sig'

	}




* ****************************************************
*
*	Handle the lower-corners at zero.
*
*		Must also loop on the upper-corner index:
*
*			X_MVT_index_UP = indeksUP0
*
*		the upper-corner censoring schemes that hit
*		NO upper-corners.
*
*		This program CANNOT handle observations that
*		hit a lower corner for some decisions and also
*		hit a upper corner for other decisions.
*
*		Therefore this lower-corner code immediately below
*		must be made to execute ONLY for observations that
*		DO NOT hit upper-corners.  This is ensured by 
*		looping on X_MVT_index_UP = indeksUP0.
*
*		In short, this code takes the program written by
*		Mikkel Barslund and adds the additional looping on
*		X_MVT_index_UP = indeksUP0   to prevent execution
*		on the observations that hit upper-corners.
*		The observations that hit upper-corners will be handled 
*		later (far below).
*

	 
	tempname iSig_a is tempm1 tempvar ip gtmp
	gen double `ip' = 0

	* Observations with ZERO equations censored

	`deb' di _newline(2)
	`deb' di in red "Calculating the ZERO equations censored cases"


	foreach indeks   of global indeks0   {
	foreach indeksUP of global indeksUP0 {	

		local res = ""
		forvalues j = 1/$X_MVT_NOEQ {
			local res `res' `e`j''
		}

		 cap matrix `iSig_a' = syminv(`Sig')
		 if _rc != 0 {
			noi di in red "Warning: cannot invert Sigma matrix"
		 }

		 replace `ip' = 0
		 forvalues i = 1/$X_MVT_NOEQ {
		 	tempvar g`i'
		 	matrix `is' = `iSig_a'[`i',1...]
			matrix colnames `is' = `res'
			matrix score `g`i'' = `is'            if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
			replace `ip' = `ip' + `e`i'' * `g`i'' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
		 }

	  tempvar logfe_1
	  gen double `logfe_1' = -$X_MVT_NOEQ/2 * log(2*_pi) - 1/2*log(det(`Sig')) - .5*`ip' ///
	  						   if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'

	  replace `lnf' = `logfe_1' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'



	  * End zero censoring
	 }  
	 }

			  
			
	/* Observations with ONE equation censored */

	`deb' noi di in red "Calculating the ONE equation censored cases"

	tempname Sig_11 Sig_tmp Sig_1 Sig_21 Sig_22 Sig_22_1 My_2_1

	foreach indeks   of global indeks1   {
	foreach indeksUP of global indeksUP0 {	


		* first the consumed goods
		* find covariance matrix from Sig

		* first reorganize columns
		* start with artifial column
		matrix `Sig_tmp' = `Sig'[1...,1]
		forvalues i = 1/`=$X_MVT_NOEQ-1' {
			matrix `Sig_tmp' = `Sig_tmp', `Sig'[1...,`=word("${ordstr`indeks'}",`i')']
		}
		* remove artificial column
		matrix `Sig_tmp' = `Sig_tmp'[1...,2...]
		* add non-use goods
		matrix `Sig_tmp' = `Sig_tmp', `Sig'[1...,`=word("${unordstr`indeks'}",1)']
	
		* then reorganize rows
		* artifial row to start with
		matrix `Sig_1' = `Sig_tmp'[1,1...]
		forvalues i = 1/`=$X_MVT_NOEQ-1' {
			matrix `Sig_1' = `Sig_1' \ `Sig_tmp'[`=word("${ordstr`indeks'}",`i')',1...]
		}
		* remove artificial row
		matrix `Sig_1' = `Sig_1'[2...,1...]
		* add non-use goods
		matrix `Sig_1' = `Sig_1' \ `Sig_tmp'[`=word("${unordstr`indeks'}",1)',1...]

		matrix `Sig_11' = `Sig_1'[1..`=$X_MVT_NOEQ-1', 1..`=$X_MVT_NOEQ-1']
		matrix `Sig_21' = `Sig_1'[$X_MVT_NOEQ, 1..`=$X_MVT_NOEQ-1']
		matrix `Sig_22' = `Sig_1'[$X_MVT_NOEQ, $X_MVT_NOEQ]
		matrix `Sig_22_1' = `Sig_22' - `Sig_21' * syminv(`Sig_11') * `Sig_21''
	
	
		local res = ""
		local or = ""

		forvalues i = 1/`=$X_MVT_NOEQ-1' {
			local res `res' `e`=word("${ordstr`indeks'}",`i')''
		}

		local or ${ordstr`indeks'}

			
		matrix `iSig_a' = syminv(`Sig_11')
		forvalues i = 1/`=$X_MVT_NOEQ-1' {
		 	tempvar g`i'
		 	matrix `is' = `iSig_a'[`i',1...]
			matrix colnames `is' = `res'
			matrix score `g`i'' = `is'                          if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
			replace `ip' = `ip' + `e`=word("`or'",`i')''*`g`i'' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
		}

		tempvar logfe_1 loghe_1
			  						  
		gen double `logfe_1' = -($X_MVT_NOEQ-1)/2 * log(2*_pi) - 1/2*log(det(`Sig_11')) - .5*`ip' ///
								if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'

		* Then censored equations
		tempvar my_1
		gen double `my_1' = 0 if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
		forvalues i = 1/`=$X_MVT_NOEQ-1' {
			replace `my_1' = `my_1' + `g`i'' * `Sig_21'[1,`i'] if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
		}

		gen double `loghe_1' = log( norm((-`xb`=word("${unordstr`indeks'}",1)'' - 		///
						`my_1')/sqrt(`Sig_22_1'[1,1])) ) if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'

		replace `lnf' = `logfe_1' + `loghe_1' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'


	  * End one censoring
	}  
	}



	/* Observations with TWO equations censored */

	`deb' noi di in red "Calculating the TWO equations censored cases"
	
	foreach indeks   of global indeks2   {
	foreach indeksUP of global indeksUP0 {	

		* first non-censored euations
		* find covariance matrix from Sig
		* reorganize columns

		* start with artifial column
		matrix `Sig_tmp' = `Sig'[1...,1]
		forvalues i = 1/`=$X_MVT_NOEQ-2' {
			matrix `Sig_tmp' = `Sig_tmp',`Sig'[1...,`=word("${ordstr`indeks'}",`i')']
		}
		* remove artificial column
		matrix `Sig_tmp' = `Sig_tmp'[1...,2...]
		* add non-use goods
		forvalues i = 1/2 {
			matrix `Sig_tmp' = `Sig_tmp',`Sig'[1...,`=word("${unordstr`indeks'}",`i')']
		}
		
		* then reorganize rows
		* start with artifial row
		matrix `Sig_1' = `Sig_tmp'[1,1...]
		forvalues i = 1/`=$X_MVT_NOEQ-2' {
			matrix `Sig_1' = `Sig_1' \ `Sig_tmp'[`=word("${ordstr`indeks'}",`i')',1...]
		}
		* remove artificial row
		matrix `Sig_1' = `Sig_1'[2...,1...]
		* add non-use goods
		forvalues i = 1/2 {
			matrix `Sig_1'=`Sig_1' \ `Sig_tmp'[`=word("${unordstr`indeks'}",`i')',1...]
		}

		matrix `Sig_11'   = `Sig_1'[1..`=$X_MVT_NOEQ-2',1..`=$X_MVT_NOEQ-2']
		matrix `Sig_21'   = `Sig_1'[`=$X_MVT_NOEQ-1'...,1..`=$X_MVT_NOEQ-2']
		matrix `Sig_22'   = `Sig_1'[`=$X_MVT_NOEQ-1'...,`=$X_MVT_NOEQ-1'...]
		matrix `Sig_22_1' = `Sig_22'-`Sig_21'*syminv(`Sig_11')*`Sig_21''
		
		
		local res = ""
		local or = ""

		forvalues i=1/`=$X_MVT_NOEQ-2' {
			local res `res' `e`=word("${ordstr`indeks'}",`i')''
		}

		local or ${ordstr`indeks'}

		matrix `iSig_a' = syminv(`Sig_11')
		replace `ip' = 0
		forvalues i = 1/`=$X_MVT_NOEQ-2' {
		 	tempvar g`i'
		 	matrix `is' = `iSig_a'[`i',1...]
			matrix colnames `is' = `res'
			matrix score `g`i'' = `is' 						   if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
			replace `ip' = `ip' + `e`=word("`or'",`i')'' * `g`i'' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
		}

		tempvar logfe_1 loghe_1
				  						  
		gen double `logfe_1' = -($X_MVT_NOEQ-2)/2*log(2*_pi)-1/2*log(det(`Sig_11'))-.5*`ip'  ///
								if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'

		* Then censored equations
		tempname Sig_21_11
			  
		matrix `Sig_21_11' = `Sig_21' * syminv(`Sig_11')
		tempvar my1 my2
		gen double `my1' = 0 if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
		gen double `my2' = 0 if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'

		forvalues i=1/`=$X_MVT_NOEQ-2' {
				replace `my1' = `my1' + `Sig_21_11'[1,`i'] * `e`=word("${ordstr`indeks'}",`i')'' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
				replace `my2' = `my2' + `Sig_21_11'[2,`i'] * `e`=word("${ordstr`indeks'}",`i')'' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
		}

		gen double `loghe_1' = log( binorm( (-`xb`=word("${unordstr`indeks'}",1)''-`my1') / sqrt(`Sig_22_1'[1,1]),	///
						(-`xb`=word("${unordstr`indeks'}",2)''-`my2') / sqrt(`Sig_22_1'[2,2]), `Sig_22_1'[1,2]/		///
						sqrt(`Sig_22_1'[1,1])/sqrt(`Sig_22_1'[2,2]) ) ) if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
	
		replace `lnf' = `logfe_1' + `loghe_1' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'

	* End Two censored equations
	}
	}



	* gen variables for use with mvnp procedure with three or more censored equations
	forvalues z = 1/$X_MVT_NOEQ {
		tempvar u`z'
		gen double `u`z'' = .
	}

				
	* if # eq > 3	
	if $X_MVT_maxcen > 3 | ($X_MVT_maxcen == 3 & $X_MVT_NOEQ > 3) {

	`deb' noi di in red "Calculating the THREE equations censored cases"

		* Observations with THREE or MORE censored equations
		forvalues y = 3/`=$X_MVT_NOEQ-1' {
			foreach indeks   of global indeks`y' {
			foreach indeksUP of global indeksUP0 {					

				* Same procedure as with 2 non-censored equations
				matrix `Sig_tmp' = `Sig'[1...,1]
				forvalues i = 1/`=$X_MVT_NOEQ-`y'' {
					matrix `Sig_tmp' = `Sig_tmp',`Sig'[1...,`=word("${ordstr`indeks'}",`i')']
				}
				* remove artificial column
				matrix `Sig_tmp' = `Sig_tmp'[1...,2...]
				* add non-use goods
				forvalues i = 1/`y' {
					matrix `Sig_tmp' = `Sig_tmp',`Sig'[1...,`=word("${unordstr`indeks'}",`i')']
				}
				
				* Rows
				matrix `Sig_1' = `Sig_tmp'[1,1...]
				forvalues i = 1/`=$X_MVT_NOEQ-`y'' {
					matrix `Sig_1' = `Sig_1' \ `Sig_tmp'[`=word("${ordstr`indeks'}",`i')',1...]
				}
				* remove artificial row
				matrix `Sig_1' = `Sig_1'[2...,1...]
				* add non-use goods
				forvalues i = 1/`y' {
					matrix `Sig_1' = `Sig_1' \ `Sig_tmp'[`=word("${unordstr`indeks'}",`i')',1...]
				}
	
				matrix `Sig_11'   = `Sig_1'[1..`=$X_MVT_NOEQ-`y'',1..`=$X_MVT_NOEQ-`y'']
				matrix `Sig_21'   = `Sig_1'[`=$X_MVT_NOEQ-`y'+1'...,1..`=$X_MVT_NOEQ-`y'']
				matrix `Sig_22'   = `Sig_1'[`=$X_MVT_NOEQ-`y'+1'...,`=$X_MVT_NOEQ-`y'+1'...]
				matrix `Sig_22_1' = `Sig_22'-`Sig_21'*syminv(`Sig_11')*`Sig_21''
					
				local res = ""
				local or = ""
	
				forvalues i = 1/`=$X_MVT_NOEQ-`y'' {
					local res `res' `e`=word("${ordstr`indeks'}",`i')''
				}
	
				local or ${ordstr`indeks'}

				matrix `iSig_a' = syminv(`Sig_11')
				replace `ip' = 0
				forvalues i = 1/`=$X_MVT_NOEQ-`y'' {
				 	tempvar g`i'
				 	matrix `is' = `iSig_a'[`i',1...]
					matrix colnames `is' = `res'
					matrix score `g`i'' = `is'                          if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
					replace `ip' = `ip' + `e`=word("`or'",`i')''*`g`i'' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
				}
		
				tempvar logfe_1 loghe_1
				gen double `logfe_1' = -($X_MVT_NOEQ-`y')/2*log(2*_pi)-1/2*log(det(`Sig_11'))-.5*`ip' ///
				 					if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
						  
				* Then censored equations
				tempname Sig_21_11
				matrix `Sig_21_11' = `Sig_21'*syminv(`Sig_11')
				forvalues z = 1/`y' {
					 tempvar my`z'
					  gen double `my`z'' = 0 if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
				}
				  
				forvalues i = 1/`=$X_MVT_NOEQ-`y'' {
					forvalues z = 1/`y' {
				  		replace `my`z'' = `my`z'' + `Sig_21_11'[`z',`i'] * `e`=word("${ordstr`indeks'}",`i')'' ///
						if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
					}
				}

				* Make cholesky factorization of variance-covariance matrix
				tempname chol_3
				tempvar  prod

				forvalues z = 1/`y' {
					replace `u`z'' = (-`xb`=word("${unordstr`indeks'}",`z')''-`my`z'')   ///
					  				if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
				}
				  

				* make string for argument in mvpn procedure
				local mvpn_str=""
				forvalues z = 1/`y' {
					local mvpn_str = "`mvpn_str'" + "`" + "u`z'" + "' "
				}
				* trim before passing to mvnp
				local mvpn_str = trim("`mvpn_str'")

				matrix `chol_3' = cholesky(`Sig_22_1')
			  	egen `prod' = mvnp(`mvpn_str') if  X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP' ///
							, prefix("$X_MVT_prefix") draws($X_MVT_D) ///
							 chol(`chol_3') $X_MVT_adoonly
	
			  	gen double `loghe_1' = log(`prod')    if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
	
				replace `lnf' = `logfe_1' + `loghe_1' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
			}
			}
		}
	* End if # eq > 3
	}
				
	/* Observations with ALL equations censored */
	if $X_MVT_maxcen == $X_MVT_NOEQ {

		`deb' noi di in red "Calculating the ALL equations censored cases"

		foreach indeks   of global indeks${X_MVT_NOEQ} {
		foreach indeksUP of global indeksUP0           {	
	 		* Make cholesky factorization of variance-covariance matrix
			tempname chol_3
			tempvar prod
			
			matrix `chol_3' = cholesky(`Sig')
			forvalues y = 1/$X_MVT_NOEQ {
		    	replace `u`y'' = (-`xb`=word("${unordstr`indeks'}",`y')'') if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
			}					  
	
			* make string for argument in mvpn procedure
			local mvpn_str = ""
			forvalues z = 1/$X_MVT_NOEQ {
				local mvpn_str = "`mvpn_str'" + "`" + "u`z'" + "' "
			}
	
			* trim before passing to mvnp
			local mvpn_str = trim("`mvpn_str'")
	
		    egen `prod' = mvnp(`mvpn_str') if  X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP' ///
							, prefix("$X_MVT_prefix")   ///
							draws($X_MVT_D) chol(`chol_3') $X_MVT_adoonly
	
		    tempvar loghe_1_15
			gen double `loghe_1_15' = log(`prod') if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
	
			replace `lnf' = `loghe_1_15'          if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
	 	}
		}
    * end all eq censored
	}








*di as error "Inside tobitCDRElowup_ll....made it to Upper Corner"




* ****************************************************
*
*	Handle the upper-corners at $40 and $46.
*
*		Must also loop on the lower-corner index:
*
*			X_MVT_index = indeks0
*
*		the lower-corner censoring schemes that hit
*		NO lower-corners.
*
*		This program CANNOT handle observations that
*		hit a lower corner for some decisions and also
*		hit a upper corner for other decisions.
*
*		Therefore this upper-corner code immediately below
*		must be made to execute ONLY for observations that
*		DO NOT hit lower-corners.  This is ensured by 
*		looping on X_MVT_index = indeks0.
*
*		In short, this code re-writes the program written by
*		Mikkel Barslund to handle upper-corners, but also
*		must adds additional looping on
*		X_MVT_index_UP = indeks0   to prevent execution
*		on the observations that hit lower-corners.
*		The observations that hit lower-corners were handled 
*		above.
*

	 
	tempname iSig_a is tempm1 tempvar ip gtmp
	gen double `ip' = 0


	* Observations with ZERO equations censored

	*	Observations that hit NO upper-corners were already
	*	processed above in the code that handles the 
	*	observations that hit NO LOWER-corners.
 
			
	/* Observations with ONE equation upper-corner */

	`deb' noi di in red "Calculating the ONE equation upper-corner cases"

	tempname Sig_11 Sig_tmp Sig_1 Sig_21 Sig_22 Sig_22_1 My_2_1

	foreach indeksUP of global indeksUP1 {
	foreach indeks   of global indeks0   {

		* first the consumed goods
		* find covariance matrix from Sig

		* first reorganize columns
		* start with artifial column
		matrix `Sig_tmp' = `Sig'[1...,1]
		forvalues i = 1/`=$X_MVT_NOEQ-1' {
			matrix `Sig_tmp' = `Sig_tmp', `Sig'[1...,`=word("${ordstrUP`indeksUP'}",`i')']
		}
		* remove artificial column
		matrix `Sig_tmp' = `Sig_tmp'[1...,2...]
		* add non-use goods
		matrix `Sig_tmp' = `Sig_tmp', `Sig'[1...,`=word("${unordstrUP`indeksUP'}",1)']
	
		* then reorganize rows
		* artifial row to start with
		matrix `Sig_1' = `Sig_tmp'[1,1...]
		forvalues i = 1/`=$X_MVT_NOEQ-1' {
			matrix `Sig_1' = `Sig_1' \ `Sig_tmp'[`=word("${ordstrUP`indeksUP'}",`i')',1...]
		}
		* remove artificial row
		matrix `Sig_1' = `Sig_1'[2...,1...]
		* add non-use goods
		matrix `Sig_1' = `Sig_1' \ `Sig_tmp'[`=word("${unordstrUP`indeksUP'}",1)',1...]

		matrix `Sig_11' = `Sig_1'[1..`=$X_MVT_NOEQ-1', 1..`=$X_MVT_NOEQ-1']
		matrix `Sig_21' = `Sig_1'[$X_MVT_NOEQ, 1..`=$X_MVT_NOEQ-1']
		matrix `Sig_22' = `Sig_1'[$X_MVT_NOEQ, $X_MVT_NOEQ]
		matrix `Sig_22_1' = `Sig_22' - `Sig_21' * syminv(`Sig_11') * `Sig_21''
	
		local res = ""
		local or = ""

		forvalues i = 1/`=$X_MVT_NOEQ-1' {
			local res `res' `e`=word("${ordstrUP`indeksUP'}",`i')''
		}

		local or ${ordstrUP`indeksUP'}

			
		matrix `iSig_a' = syminv(`Sig_11')
		forvalues i = 1/`=$X_MVT_NOEQ-1' {
		 	tempvar g`i'
		 	matrix `is' = `iSig_a'[`i',1...]
			matrix colnames `is' = `res'
			matrix score `g`i'' = `is'                          if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
			replace `ip' = `ip' + `e`=word("`or'",`i')''*`g`i'' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
		}

		tempvar logfe_1 loghe_1
			  						  
		gen double `logfe_1' = -($X_MVT_NOEQ-1)/2 * log(2*_pi) - 1/2*log(det(`Sig_11')) - .5*`ip' ///
								if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'


		* Then censored equations
		tempvar my_1
		gen double `my_1' = 0 if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
		forvalues i = 1/`=$X_MVT_NOEQ-1' {
			replace `my_1' = `my_1' + `g`i'' * `Sig_21'[1,`i'] if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
		}

		* There are two ways to have one upper-corner:
		*	(1) the upper-corner is at $40 or
		*	(2) the upper-corner is at $46.
		gen double `loghe_1' = log( norm(-(40 -`xb`=word("${unordstrUP`indeksUP'}",1)'' - 	///
						`my_1')/sqrt(`Sig_22_1'[1,1])) ) 							///
						if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP' &		///
						   X_MVT_up1 == 1
		replace    `loghe_1' = log( norm(-(46 -`xb`=word("${unordstrUP`indeksUP'}",1)'' - 	///
						`my_1')/sqrt(`Sig_22_1'[1,1])) ) 							///
						if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP' &		///
						   X_MVT_up2 == 1

		replace `lnf' = `logfe_1' + `loghe_1' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'


	  * End one upper-corner
	}  
	}



	/* Observations with TWO equations upper-corner */

	`deb' noi di in red "Calculating the TWO equations upper-corner cases"
	
	foreach indeksUP of global indeksUP2 {
	foreach indeks   of global indeks0   {	

		* first non-censored euations
		* find covariance matrix from Sig
		* reorganize columns

		* start with artifial column
		matrix `Sig_tmp' = `Sig'[1...,1]
		forvalues i = 1/`=$X_MVT_NOEQ-2' {
			matrix `Sig_tmp' = `Sig_tmp',`Sig'[1...,`=word("${ordstrUP`indeksUP'}",`i')']
		}
		* remove artificial column
		matrix `Sig_tmp' = `Sig_tmp'[1...,2...]
		* add non-use goods
		forvalues i = 1/2 {
			matrix `Sig_tmp' = `Sig_tmp',`Sig'[1...,`=word("${unordstrUP`indeksUP'}",`i')']
		}
		
		* then reorganize rows
		* start with artifial row
		matrix `Sig_1' = `Sig_tmp'[1,1...]
		forvalues i = 1/`=$X_MVT_NOEQ-2' {
			matrix `Sig_1' = `Sig_1' \ `Sig_tmp'[`=word("${ordstrUP`indeksUP'}",`i')',1...]
		}
		* remove artificial row
		matrix `Sig_1' = `Sig_1'[2...,1...]
		* add non-use goods
		forvalues i = 1/2 {
			matrix `Sig_1'=`Sig_1' \ `Sig_tmp'[`=word("${unordstrUP`indeksUP'}",`i')',1...]
		}

		matrix `Sig_11'   = `Sig_1'[1..`=$X_MVT_NOEQ-2',1..`=$X_MVT_NOEQ-2']
		matrix `Sig_21'   = `Sig_1'[`=$X_MVT_NOEQ-1'...,1..`=$X_MVT_NOEQ-2']
		matrix `Sig_22'   = `Sig_1'[`=$X_MVT_NOEQ-1'...,`=$X_MVT_NOEQ-1'...]
		matrix `Sig_22_1' = `Sig_22'-`Sig_21'*syminv(`Sig_11')*`Sig_21''
		
		
		local res = ""
		local or = ""

		forvalues i=1/`=$X_MVT_NOEQ-2' {
			local res `res' `e`=word("${ordstrUP`indeksUP'}",`i')''
		}

		local or ${ordstrUP`indeksUP'}

		matrix `iSig_a' = syminv(`Sig_11')
		replace `ip' = 0
		forvalues i = 1/`=$X_MVT_NOEQ-2' {
		 	tempvar g`i'
		 	matrix `is' = `iSig_a'[`i',1...]
			matrix colnames `is' = `res'
			matrix score `g`i'' = `is' 					    if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
			replace `ip' = `ip' + `e`=word("`or'",`i')'' * `g`i'' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
		}

		tempvar logfe_1 loghe_1
				  						  
		gen double `logfe_1' = -($X_MVT_NOEQ-2)/2*log(2*_pi)-1/2*log(det(`Sig_11'))-.5*`ip'  ///
								if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'

		* Then censored equations
		tempname Sig_21_11
			  
		matrix `Sig_21_11' = `Sig_21' * syminv(`Sig_11')
		tempvar my1 my2
		gen double `my1' = 0 if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
		gen double `my2' = 0 if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'

		forvalues i=1/`=$X_MVT_NOEQ-2' {
				replace `my1' = `my1' + `Sig_21_11'[1,`i'] * `e`=word("${ordstrUP`indeksUP'}",`i')'' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
				replace `my2' = `my2' + `Sig_21_11'[2,`i'] * `e`=word("${ordstrUP`indeksUP'}",`i')'' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
		}


		* There are three ways to have TWO upper-corners:
		*	(1) two upper-corners at $40, or
		*	(2) one upper-corner  at $40 and one upper-corner at $46, or
		*	(3) two upper-corners at $46.
		gen double `loghe_1' = log( binorm( -(40 -`xb`=word("${unordstrUP`indeksUP'}",1)''-`my1') / sqrt(`Sig_22_1'[1,1]),	///
						-(40 -`xb`=word("${unordstrUP`indeksUP'}",2)''-`my2') / sqrt(`Sig_22_1'[2,2]), `Sig_22_1'[1,2]/		///
						sqrt(`Sig_22_1'[1,1])/sqrt(`Sig_22_1'[2,2]) ) ) 			///
						if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP' &	///
						   X_MVT_up1 == 2

		replace    `loghe_1' = log( binorm( -(40 -`xb`=word("${unordstrUP`indeksUP'}",1)''-`my1') / sqrt(`Sig_22_1'[1,1]),	///
						-(46 -`xb`=word("${unordstrUP`indeksUP'}",2)''-`my2') / sqrt(`Sig_22_1'[2,2]), `Sig_22_1'[1,2]/		///
						sqrt(`Sig_22_1'[1,1])/sqrt(`Sig_22_1'[2,2]) ) ) 			///
						if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP' &	///
						   X_MVT_up1 == 1 & X_MVT_up2 == 1

		replace    `loghe_1' = log( binorm( -(46 -`xb`=word("${unordstrUP`indeksUP'}",1)''-`my1') / sqrt(`Sig_22_1'[1,1]),	///
						-(46 -`xb`=word("${unordstrUP`indeksUP'}",2)''-`my2') / sqrt(`Sig_22_1'[2,2]), `Sig_22_1'[1,2]/		///
						sqrt(`Sig_22_1'[1,1])/sqrt(`Sig_22_1'[2,2]) ) ) 			///
						if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP' &	///
						   X_MVT_up2 == 2

		replace `lnf' = `logfe_1' + `loghe_1' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'

	* End Two upper-corner equations
	}
	}






	* gen variables for use with mvnp procedure with three or more upper-corner equations
	forvalues z = 1/$X_MVT_NOEQ {
		tempvar u`z'
		gen double `u`z'' = .
	}

				
	* if # eq > 3	
	if $X_MVT_maxcen_UP > 3 | ($X_MVT_maxcen_UP == 3 & $X_MVT_NOEQ > 3) {


	`deb' noi di in red "Calculating the THREE equations upper-corner cases"

		* Observations with THREE or MORE upper-corner equations
		forvalues y = 3/`=$X_MVT_NOEQ-1' {
			foreach indeksUP of global indeksUP`y' {
			foreach indeks   of global indeks0     {					

				* Same procedure as with 2 non-censored equations
				matrix `Sig_tmp' = `Sig'[1...,1]
				forvalues i = 1/`=$X_MVT_NOEQ-`y'' {
					matrix `Sig_tmp' = `Sig_tmp',`Sig'[1...,`=word("${ordstrUP`indeksUP'}",`i')']
				}
				* remove artificial column
				matrix `Sig_tmp' = `Sig_tmp'[1...,2...]
				* add non-use goods
				forvalues i = 1/`y' {
					matrix `Sig_tmp' = `Sig_tmp',`Sig'[1...,`=word("${unordstrUP`indeksUP'}",`i')']
				}

				* Rows
				matrix `Sig_1' = `Sig_tmp'[1,1...]
				forvalues i = 1/`=$X_MVT_NOEQ-`y'' {
					matrix `Sig_1' = `Sig_1' \ `Sig_tmp'[`=word("${ordstrUP`indeksUP'}",`i')',1...]
				}
				* remove artificial row
				matrix `Sig_1' = `Sig_1'[2...,1...]

				* add non-use goods
				forvalues i = 1/`y' {
					matrix `Sig_1' = `Sig_1' \ `Sig_tmp'[`=word("${unordstrUP`indeksUP'}",`i')',1...]
				}

				matrix `Sig_11'   = `Sig_1'[1..`=$X_MVT_NOEQ-`y'',1..`=$X_MVT_NOEQ-`y'']
				matrix `Sig_21'   = `Sig_1'[`=$X_MVT_NOEQ-`y'+1'...,1..`=$X_MVT_NOEQ-`y'']
				matrix `Sig_22'   = `Sig_1'[`=$X_MVT_NOEQ-`y'+1'...,`=$X_MVT_NOEQ-`y'+1'...]
				matrix `Sig_22_1' = `Sig_22'-`Sig_21'*syminv(`Sig_11')*`Sig_21''

					
				local res = ""
				local or = ""
	
				forvalues i = 1/`=$X_MVT_NOEQ-`y'' {
					local res `res' `e`=word("${ordstrUP`indeksUP'}",`i')''
				}
	
				local or ${ordstrUP`indeksUP'}

				matrix `iSig_a' = syminv(`Sig_11')
				replace `ip' = 0
				forvalues i = 1/`=$X_MVT_NOEQ-`y'' {
				 	tempvar g`i'
				 	matrix `is' = `iSig_a'[`i',1...]
					matrix colnames `is' = `res'
					matrix score `g`i'' = `is'                          if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
					replace `ip' = `ip' + `e`=word("`or'",`i')''*`g`i'' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
				}

				tempvar logfe_1 loghe_1
				gen double `logfe_1' = -($X_MVT_NOEQ-`y')/2*log(2*_pi)-1/2*log(det(`Sig_11'))-.5*`ip' ///
				 					if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
						  
				* Then censored equations
				tempname Sig_21_11
				matrix `Sig_21_11' = `Sig_21'*syminv(`Sig_11')
				forvalues z = 1/`y' {
					 tempvar my`z'
					  gen double `my`z'' = 0 if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
				}
				  
				forvalues i = 1/`=$X_MVT_NOEQ-`y'' {
					forvalues z = 1/`y' {
				  		replace `my`z'' = `my`z'' + `Sig_21_11'[`z',`i'] * `e`=word("${ordstrUP`indeksUP'}",`i')'' ///
						if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
					}
				}

				* Make cholesky factorization of variance-covariance matrix
				tempname chol_3
				tempvar  prod

				* Handle all the upper-corners, except the last 2.
				* Only the last 2 upper-corner eqns (decisions 5-6) could possibly
				*	hit a corner at $46.
				* Therefore, the previous corners must all be at $40
				*
				* Remember, in this part of the code there are 3+ upper-corner eqns, so
				* there will be at least one upper-corner eqn at $40.
				forvalues z = 1/`=`y'-2' {
					replace `u`z'' = -(40 -`xb`=word("${unordstrUP`indeksUP'}",`z')''-`my`z'')   		///
					  				if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
				}


				* Handle the last 2 upper-corners.
				*
				* Remember, in this part of the code there are 3+ upper-corner eqns, so
				* the last 2 eqns are necessarily at upper-corners.  We just need to figure 
				* out the the last 2 upper-corners are at $40 or $46.
				*
				* There are only three cases.  The last 2 upper-corners can be:
				*	(1) two upper-corners at $40 (X_MVT_up2 = 0), or
				*	(2) one upper-corner  at $40 and one upper-corner at $46 (X_MVT_up2 = 1), or
				*	(3) two upper-corners at $46 (X_MVT_up2 = 2).

				* (1)
				forvalues z = `=`y'-1'/`y' {
					replace `u`z'' = -(40 -`xb`=word("${unordstrUP`indeksUP'}",`z')''-`my`z'')   		///
					  				if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP' &	///
										X_MVT_up2 == 0
				} 

				* (2)
				forvalues z = `=`y'-1'/`=`y'-1' {
					replace `u`z'' = -(40 -`xb`=word("${unordstrUP`indeksUP'}",`z')''-`my`z'')   		///
					  				if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP' &	///
										X_MVT_up2 == 1
				}
				forvalues z = `y'/`y' {
					replace `u`z'' = -(46 -`xb`=word("${unordstrUP`indeksUP'}",`z')''-`my`z'')   		///
					  				if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP' &	///
										X_MVT_up2 == 1
				}


				* (3)
				forvalues z = `=`y'-1'/`y' {
					replace `u`z'' = -(46 -`xb`=word("${unordstrUP`indeksUP'}",`z')''-`my`z'')   		///
					  				if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP' &	///
										X_MVT_up2 == 2
				}


				* make string for argument in mvpn procedure
				local mvpn_str=""
				forvalues z = 1/`y' {
					local mvpn_str = "`mvpn_str'" + "`" + "u`z'" + "' "
				}
				* trim before passing to mvnp
				local mvpn_str = trim("`mvpn_str'")


				matrix `chol_3' = cholesky(`Sig_22_1')
			  	egen `prod' = mvnp(`mvpn_str') if  X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP' ///
							, prefix("$X_MVT_prefix") draws($X_MVT_D) ///
							 chol(`chol_3') $X_MVT_adoonly
	
			  	gen double `loghe_1' = log(`prod')    if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
	
				replace `lnf' = `logfe_1' + `loghe_1' if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
			}
			}
		}
	* End if # eq > 3
	}
	



			
	/* Observations with ALL equations upper-corner */
	if $X_MVT_maxcen_UP == $X_MVT_NOEQ {

		`deb' noi di in red "Calculating the ALL equations upper-corner cases"

		foreach indeksUP of global indeksUP${X_MVT_NOEQ} {
		foreach indeks   of global indeks0               {	
	 		* Make cholesky factorization of variance-covariance matrix
			tempname chol_3
			tempvar prod
			
			matrix `chol_3' = cholesky(`Sig')

			* Handle all upper-corners, except the last 2.
			*	These upper-corners hit $40.
			forvalues y = 1/`=$X_MVT_NOEQ-2' {
		    	replace `u`y'' = -(40 -`xb`=word("${unordstrUP`indeksUP'}",`y')'') if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
			}					  
	
			* Handle the last 2 all upper-corners.
			*	These upper-corners hit $46.
			forvalues y = `=$X_MVT_NOEQ-1'/$X_MVT_NOEQ {
		    	replace `u`y'' = -(46 -`xb`=word("${unordstrUP`indeksUP'}",`y')'') if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
			}


			* make string for argument in mvpn procedure
			local mvpn_str = ""
			forvalues z = 1/$X_MVT_NOEQ {
				local mvpn_str = "`mvpn_str'" + "`" + "u`z'" + "' "
			}
	
			* trim before passing to mvnp
			local mvpn_str = trim("`mvpn_str'")
	
		    egen `prod' = mvnp(`mvpn_str') if  X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP' ///
							, prefix("$X_MVT_prefix")   ///
							draws($X_MVT_D) chol(`chol_3') $X_MVT_adoonly
	
		    tempvar loghe_1_15
			gen double `loghe_1_15' = log(`prod') if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
	
			replace `lnf' = `loghe_1_15'              if X_MVT_index == `indeks' & X_MVT_index_UP == `indeksUP'
	 	}
		}
    * end all eq upper-corner
	}





	* End quietly and evaluator			
	 }




end
