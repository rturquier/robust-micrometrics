*! Random Effects Tobit - for a Cobb-Douglass Impure Altruism model
*!
*!
*!    Constrains  sigmas in each equation to be equal, and
*!    Constrains  rhos   between pairs of equations to be equal.
*!
*!
*!	  Handles corners at
*!
*!			1. $0
*!		or 
*!			2. $40/$46
*!		but
*!			3. NOT a subject that sometimes chooses $0 and
*!					other times chooses $40/$46.
*!
*!	This program borrowed heavily (well, hacked):
*! 		mvtobit.ado starting point was:
*! 		version 1.0  August 2007 by Mikkel Barslund
*! 		Multivariate tobit by method of MSL.  
cap program drop tobitCDRElowup6
program define tobitCDRElowup6, eclass byable(onecall) sortpreserve
    version 8.2
    if replay() {
        if "`e(cmd)'" != "tobitCDRElowup6" {
            di as error "results for tobitCDRElowup6 not found"
            exit 301
        }
        if _by() { 
            error 190 
        } 
        Display `0'
        exit `rc'
    }
    if _by() {
        by `_byvars'`_byrc0': Estimate `0'
    }
    else    Estimate `0'
end

program define Estimate, eclass byable(recall)

        /* First parse the equation stuff into constituent 
           components and thereby get the number of equations.
           Completely hacked from -mvprobit-.
        */

    loc i = 1

    loc paren "("

    while "`paren'" == "(" {    
        gettoken e`i' 0:0, parse(" ,[") match(paren)
        loc left "`0'"
        loc junk: subinstr loc e`i' ":" ":", count(loc number)
        if "`number'" == "1" {
            gettoken dep`i'n  e`i':  e`i', parse(":")
            gettoken junk  e`i':  e`i', parse(":")
        }
        loc  e`i' : subinstr loc  e`i' "=" " "
        gettoken dep`i' 0:  e`i', parse(" ,[") 
        unab dep`i': `dep`i'' 

            /* collecting together -- for e.g. marking sample */
        loc deps "`deps' `dep`i''"

        confirm variable `dep`i''
        if "`dep`i'n'" == "" {
            loc dep`i'n "`dep`i''"
        }
        syntax [varlist(default=none)] [, noCONstant]
        loc ind`i' `varlist'        
        loc ninds`i' : word count `ind`i''
        if "`constant'" == "" {
            loc ninds`i' = `ninds`i'' + 1
        }
        loc inds "`inds' `ind`i''"
        loc nc`i' `constant'
        loc  0 "`left'"
            /* reset `paren' to empty when thru eqns */
        gettoken check : left, parse(" ,[") match(paren)
        loc i = `i' + 1
    }


            /* Clear (permanent) variables for future use. */
    cap drop  X_MVT_index
    cap drop  X_MVT_index_UP
    cap drop  X_MVT_low1
    cap drop  X_MVT_up1
    cap drop  X_MVT_up2


            /* Ensure globals used later are already clear */
            /* Using horrible prefix reduces chance of overwriting */
    foreach g in X_MVT_i X_MVT_D X_MVT_atrho X_MVT_slatrho X_MVT_std ///
                 X_MVT_eqs X_MVT_slstd X_MVT_prefix S_MVT_rho { 
            glo `g' 
    } 

            /* number of equations */
    glo X_MVT_NOEQ = `i' - 1

    if $X_MVT_NOEQ < 2 { 
        di as error "More than 1 equation required (Use -Tobit- for single equation)." 
        exit 198
    }


	/* Create globals with the names of the dep. vars, ind. variables, and betas */
    forvalues i = 1/$X_MVT_NOEQ  {
	      global X_TRE_dep`i'  `dep`i''
		global X_TRE_ind`i'  `ind`i''
    }

    forvalues i = 2/$X_MVT_NOEQ  {
		local ilast `i'-1
	      if "`ninds`i''" != "`ninds`ilast''" {
			di as error "Tobit RE must have same number of xs in each equation"
			exit 301
		}
    }

    foreach g in X_TRE_b X_TRE_slb X_TRE_cmb X_TRE_beq { 	// Clear globals
            glo `g' 
    } 

    local Ktemp = `ninds1' - 1		// ninds1 = K + 1 (includes the constant term)
    forvalues k = 1/`Ktemp'  {
		glo X_TRE_b   "$X_TRE_b b`k' "
		glo X_TRE_slb "$X_TRE_slb /b`k' "
		glo X_TRE_beq "$X_TRE_beq  b`k':_cons "
	}
		glo X_TRE_beq "$X_TRE_beq  b0:_cons lnsigma1:_cons"	// Add const. term and ln sigma1



            /* remaining options in command line */
    loc 0 "`left'"
    syntax [if] [in] [pw fw iw aw] [, DRaws(integer 5)  Robust Cluster(varname)      ///
           Level(integer $S_level) Beta0  Seed(integer 123456789)                    ///
           ATRho0(string) AN PREfix(string) BUrn(integer 0) RANdom HRANdom SHuffle   ///
           ADOONly PRIMes(string) INIT(string) noLOG MLOpts(string) * ]

            /* Set various options */
    loc draws "`draws'"
    glo X_MVT_D = `draws'

    glo X_MVT_prefix "X_MVT" 
    if "`prefix'" != ""  glo X_MVT_prefix "`prefix'"

    glo X_MVT_adoonly = "`adoonly'"

    
    if "`primes'" != "" loc primes "primes(`primes')"

    set seed `seed'

    loc option0 `options'
    marksample touse
    markout `touse' `deps' `inds'   

    loc wtype `weight'
    loc wtexp `"`exp'"'
    if "`weight'" != "" { 
        loc wgt `"[`weight'`exp']"'  
    }
    if "`weight'" == "pweight" | "`cluster'" != "" {
            loc robust "robust"
    }

    if "`cluster'" ! = "" { 
        loc clopt "cluster(`cluster')" 
    }
    mlopts stdopts, `option0'

    if "`level'" != "" {
        loc level "level(`level')"
    }
        

    if "`log'" == "" {
                loc log "noisily"
        }
        else   {
         loc log "quietly"
    }

    loc log2 = cond("`beta0'" == "", "quietly", "noisily")


            /*  Checking of depvars etc. */
    quietly {
        count if `touse' 
        if r(N) == 0 { 
            di as error "no valid observations"
            error 2000
        }
        loc N = r(N)
        foreach var of loc deps {
            capture assert (`var' >= 0) if `touse'
                if _rc==9 {
                    di as error "depvar `var' should be greater or equal to zero (>= 0)"
                    exit 450
                }
            count if `var' == 0 & `touse'
            loc d0 = r(N)
            if `d0' == 0 {
                di as error "`var' is never zero"
***                exit 2000
            }
            else if `d0' == `N' {
                di as error "`var' is always zero"
                exit 2000 
            }
        }


										// Check the upper corner.
        forvalues d = 1/4 {
            capture assert (`dep`d'' <= 40) if `touse'
                if _rc==9 {
                    di as error "depvar `dep`d'' should be less than or equal to 40 (<= 40)"
                    exit 450
                }
		}
        forvalues d = 5/6 {
            capture assert (`dep`d'' <= 46) if `touse'
                if _rc==9 {
                    di as error "depvar `dep`d'' should be less than or equal to 46 (<= 46)"
                    exit 450
                }
		}


    }							// End of quietly



    tempname C
    glo X_MVT_C "`C'"  /* used for matrix name in evaluation program */

        /*  DO NOT Get starting values from marginal univariate tobits 
            DO NOT check collinearities among RHS vbles 
            DO     create macros containing bits of syntax for parsing to -ml- evaluator
            DO     Fit univariate Tobits regressed on JUST A CONSTANT TERM to get
            comparison likelihood value                                                 */
        

    forval i = 1/$X_MVT_NOEQ {

        tempvar y0  
        _rmcoll `ind`i'' `wgt' if `touse'
        loc ind`i' "`r(varlist)'"
        gen double `y0' = `dep`i''
        qui replace `y0' = . if `dep`i''==0
        `log2' intreg `y0' `dep`i''             `wgt' if `touse'

        if `i' == 1 {
            loc ll0 = e(ll)
        }

        if `i' > 1  {
            loc ll0 = e(ll) + `ll0'   /* logL for comparison model */
        }

        glo X_MVT_eqs "$X_MVT_eqs (`dep`i'n': `dep`i'' = `ind`i'') "
        glo X_MVT_i "$X_MVT_i xb`i' "

        glo X_MVT_std "$X_MVT_std lnsigma`i' "
        glo X_MVT_slstd "$X_MVT_slstd /lnsigma`i' "
        
        forval j = `=`i'+1'/$X_MVT_NOEQ {
            glo X_MVT_atrho "$X_MVT_atrho atrho`i'`j'"
            glo X_MVT_slatrho "$X_MVT_slatrho /atrho`i'`j'"
            glo S_MVT_rho "$S_MVT_rho rho`j'`i' ="
        }
        drop `y0'
    }




*** *** *** ***   BOOKKEEPING  starts  here   *** *** *** ***

            /* New stuff - Lower corner.  Mikkel Barslund's original code for y <= 0 and y > 0. */
    qui {
    

            /* Find max censored equations */
        loc consstr = ""
        forval i=1/$X_MVT_NOEQ {
            tempvar i`i'
            gen `i`i'' = (`dep`i''>0) if `touse'
            loc consstr = "`consstr' " + "`i`i'' "
        }
        tempvar cons noncons
        egen `cons' = rsum(`consstr') if `touse'
    
        gen `noncons' = $X_MVT_NOEQ-`cons' if `touse'
        gen X_MVT_low1 = `noncons' if `touse'  // perm. var. for the log-like

        sum `cons'  if `touse'
        glo X_MVT_maxcen = $X_MVT_NOEQ - r(min)
    
        /* General principle:
           A permanent variable (X_MVT_index) holds a unique number for each 
				lower-corner censoring scheme.
           For each censoring scheme - i.e. each unique value of X_MVT_index there are 3 global string
           variables: ordstr,  unordstr  and resstr.
           * ordstr(X_MVT_index) - contains ordered (numbers for ) non-censored equations for indeks value = indeks
           * resstr(X_MVT_index) - contains ordered residuals for non-censored equations   - :: -
           * unordstr(X_MVT_index) - contains ordered (numbers for) censored equations   - :: -
    
           Another permanent variable (X_MVT_low1; generated above) holds the number of
              lower-corners; created for the log-likelihood program.

		   See below for censoring schemes for the upper-corners.

		   This program CANNOT handle an observation that has both a lower-corner,
		   and an upper-corner.  An observation with both lower- and upper- corners
		   will cause an error-interrupt below. 									*/



                /* Make strings with lower-corner censoring pattern */
        tempvar indexstr 
        gen `indexstr'=""  if `touse'
        forval j=1/$X_MVT_NOEQ {
            replace `indexstr' = `indexstr' + " 1" if `dep`j''==0 & `touse'
            replace `indexstr' = `indexstr' + " 0" if `dep`j''>0  & `touse'
        }
        replace `indexstr' = trim(`indexstr') if `touse'
    
                /* Assign unique value to each censoring pattern */
        egen X_MVT_index = group(`indexstr') if `touse'
        sum X_MVT_index if `touse'
        loc i_min = r(min)
        loc i_max = r(max)
    
                /* Create: ordstr,  unordstr  and resstr */
        forval indeks = `i_min'/`i_max' {
            tempvar tempstr tempstr2
            glo ordstr`indeks' = ""
            glo unordstr`indeks' = ""
            glo resstr`indeks' = ""
    
            gen `tempstr2' = `indexstr' if X_MVT_index == `indeks' &  `touse'
            egen `tempstr' = mode(`tempstr2')
            loc seq = `tempstr'
            drop `tempstr' `tempstr2'
    
            forval j=1/$X_MVT_NOEQ {
                if word("`seq'",`j') == "1" {
                    glo unordstr`indeks' = "${unordstr`indeks'} " + "`j'"
                }
                else {
                     glo ordstr`indeks' = "${ordstr`indeks'} " + "`j'"
                     glo resstr`indeks' = "${resstr`indeks'} " + "e`j'"
                }
            }
        }
    }				// end of qui(etly)


            /* Create meta indices - indices constaining information on which
                index belongs to which censoring scheme  */

    forval y=0/$X_MVT_NOEQ {
        glo indeks`y' = ""
        forval indeks = `i_min'/`i_max' {
            qui count if X_MVT_index == `indeks' & `cons' == $X_MVT_NOEQ-`y' & `touse'
            if r(N)~=0 {
                glo indeks`y' = "${indeks`y'} " + "`indeks'"
            }
        }
    }




			/* Upper corner.  My code for y >= 40 or  y < 40, and	*/
			/*							y >= 46 or  y < 46.		*/
    qui {
    
            /* Find max censored equations */
        loc consstrUP  = ""
        loc consstrUP1 = ""
        loc consstrUP2 = ""

        forval i=1/4 {
            tempvar d`i'
            gen `d`i'' = (`dep`i''<40) if `touse'
            loc consstrUP  = "`consstrUP' "  + "`d`i'' "
            loc consstrUP1 = "`consstrUP1' " + "`d`i'' "
			sum `dep`i''
             if r(max)>40 {
                 di as error "An observation(s) in eqns. 1-4 has max > $40.  Not allowed"
                 exit 459
             }
        }

        forval i=5/$X_MVT_NOEQ {
            tempvar d`i'
            gen `d`i'' = (`dep`i''<46) if `touse'
            loc consstrUP  = "`consstrUP' "  + "`d`i'' "
            loc consstrUP2 = "`consstrUP2' " + "`d`i'' "
			sum `dep`i''
             if r(max)>46 {
                 di as error "An observation(s) in eqns. 5-6 has max > $46.  Not allowed"
                 exit 459
             }
        }

        tempvar consUP nonconsUP consUP1 consUP2
        egen `consUP'  = rsum(`consstrUP')  if `touse'
        egen `consUP1' = rsum(`consstrUP1') if `touse'
        egen `consUP2' = rsum(`consstrUP2') if `touse'
    
        gen `nonconsUP' = $X_MVT_NOEQ-`consUP'  if `touse'
        gen X_MVT_up1   = 4          -`consUP1' if `touse'
        gen X_MVT_up2   = 2          -`consUP2' if `touse'

        sum `consUP'  if `touse'
        glo X_MVT_maxcen_UP = $X_MVT_NOEQ - r(min)
    
        /* General principle:
           A permanent variable (X_MVT_index_UP) holds a unique number for each 
				upper-corner censoring scheme.
           For each censoring scheme - i.e. each unique value of X_MVT_index_UP there are 3 global string
           variables: ordstr,  unordstr  and resstr.
           * ordstrUP(X_MVT_index_UP) - contains ordered (numbers for ) non-censored equations for indeks value = indeks
           * resstrUP(X_MVT_index_UP) - contains ordered residuals for non-censored equations   - :: -
           * unordstrUP(X_MVT_index_UP) - contains ordered (numbers for) censored equations   - :: -
    
           Other permanent variables (X_MVT_up1 and X_MVT_up2; generated above) hold the number of
              upper-corners; created for the log-likelihood program.
		   The are two types of upper-corners: 
				(1) at $40 that can happen at decisions 1-4, and
				(2) at $46 that can happen at decisions 5-6.

		   See above for censoring schemes for the lower-corners.

		   This program CANNOT handle an observation that has both a lower-corner,
		   and an upper-corner.  An observation with both lower- and upper- corners
		   will cause an error-interrupt right now. 									*/

         capture assert (X_MVT_low1==0 | (X_MVT_up1==0 & X_MVT_up2==0)) if `touse'
             if _rc==9 {
                 di as error "An observation(s) has BOTH lower- and upper- corners.  Not allowed"
                 exit 459
             }

             if $X_MVT_NOEQ~=6 {
                 di as error "Number of equations must be 6."
                 di as error "Also, equations 1-4: upper-corner must be $40."
                 di as error "Also, equations 5-6: upper-corner must be $46."
                 exit 459
             }
    
                /* Make strings with upper-corner censoring pattern */
        tempvar indexstrUP
        gen `indexstrUP'=""  if `touse'
        forval j=1/4 {
            replace `indexstrUP' = `indexstrUP' + " 1" if `dep`j''==40 & `touse'
            replace `indexstrUP' = `indexstrUP' + " 0" if `dep`j''<40  & `touse'
        }
        forval j=5/$X_MVT_NOEQ {
            replace `indexstrUP' = `indexstrUP' + " 1" if `dep`j''==46 & `touse'
            replace `indexstrUP' = `indexstrUP' + " 0" if `dep`j''<46  & `touse'
        }

        replace `indexstrUP' = trim(`indexstrUP') if `touse'
    
                /* Assign unique value to each censoring pattern */
        egen X_MVT_index_UP = group(`indexstrUP') if `touse'		// Here is the new permanent variable.
        sum X_MVT_index_UP if `touse'
        loc d_min = r(min)
        loc d_max = r(max)
    
                /* Create: ordstr,  unordstr  and resstr */
        forval indeks = `d_min'/`d_max' {
            tempvar tempstrUP tempstr2UP
            glo ordstrUP`indeks' = ""
            glo unordstrUP`indeks' = ""
            glo resstrUP`indeks' = ""
    
            gen `tempstr2UP' = `indexstrUP' if X_MVT_index_UP == `indeks' &  `touse'
            egen `tempstrUP' = mode(`tempstr2UP')
            loc seq = `tempstrUP'
            drop `tempstrUP' `tempstr2UP'
    
            forval j=1/$X_MVT_NOEQ {
                if word("`seq'",`j') == "1" {
                    glo unordstrUP`indeks' = "${unordstrUP`indeks'} " + "`j'"
                }
                else {
                     glo ordstrUP`indeks' = "${ordstrUP`indeks'} " + "`j'"
                     glo resstrUP`indeks' = "${resstrUP`indeks'} " + "e`j'"
                }
            }
        }
    }				// end of qui(etly)


            /* Create meta indices - indices constaining information on which
                index belongs to which censoring scheme  */

    forval y=0/$X_MVT_NOEQ {
        glo indeksUP`y' = ""
        forval indeks = `d_min'/`d_max' {
            qui count if X_MVT_index_UP == `indeks' & `consUP' == $X_MVT_NOEQ-`y' & `touse'
            if r(N)~=0 {
                glo indeksUP`y' = "${indeksUP`y'} " + "`indeks'"
            }
        }
    }


*list id `dep1' `dep2' `dep3' `dep4' `dep5' `dep6' `indexstrUP' 


*** *** *** ***   End   of   BOOKKEEPING   *** *** *** ***




            /* Random drawing using -mdraws- */

    if $X_MVT_NOEQ==2 {
        di in green "Maximum # of equations is $X_MVT_NOEQ == 2. Simulations are not needed." _n
        di in green "Continues with conventional ML" _n
    }
    else {
        di in green "Maximum # of censored equations is $X_MVT_maxcen" _n
        di in green "Draws (Halton/pseudo random) are being made:" _n
	   di _newline(3)
	   di in green "Local |an|      contains |`an'|"    
	   di in green "Local |burn|    contains |`burn'|"    
	   di in green "Local |random|  contains |`random'|"    
	   di in green "Local |hrandom| contains |`hrandom'|"    
	   di in green "Local |shuffle| contains |`shuffle'|"    
	   di in green "Local |seed|    contains |`seed'|"    


        mdraws if `touse', draws($X_MVT_D) neq($X_MVT_maxcen) prefix($X_MVT_prefix) `an' burn(`burn') ///
            `random' `primes' `hrandom' `shuffle' replace seed(`seed')

        * Double number of draws if antithetic is on
        if "`an'" != ""  glo X_MVT_D = 2*$X_MVT_D
    }

    if $X_MVT_NOEQ==2 {
        loc title "Random Effects Tobit from mvtobit hack, NOEQ == 2"
    }
    else {
        loc title "Random Effects Tobit using MSL (# draws = $X_MVT_D)"
    }






            /* Estimate model - I changed "$_MVT_maxcen" to "$X_MVT_NOEQ"  */
    if $X_MVT_NOEQ==2 {
        ml model lf tobitCDRElowup_twoll /alpha /beta /lnsigma1 /lnsigma2 `wgt' if `touse', 	///
			maximize  difficult ///
               collinear wald(0) 									///
			init(alpha:_cons=.50      	   beta:_cons=.10	 		///
			  lnsigma1:_cons= 2.3025851  lnsigma2:_cons=2.8518912)		///
			title(`title') `robust'                 				///
               search(off) `clopt' `level' `mlopts' `stdopts'


    }
    else {
        ml model lf tobitCDRElowup_ll /alpha /beta /lnsigma1 /lnsigma2 `wgt' if `touse', 	///
			maximize  difficult ///
               collinear wald(0) 									///
			init(alpha:_cons=.50      	   beta:_cons=.10	 		///
			  lnsigma1:_cons= 2.3025851  lnsigma2:_cons=2.8518912)		///
			title(`title') `robust'                 				///
               search(off) `clopt' `level' `mlopts' `stdopts'
    }
            
/*****************************
For alpha = 0, beta = 0, rho = 0, sigma2 = 13.62937 (the value of sigma2 determined by the total 
	sigma1^2 + sigma2^2  	from the run we report):
			init(alpha:_cons=.0      	   beta:_cons=.0	 		///
			  lnsigma1:_cons= -16.118096  lnsigma2:_cons=2.612227)	///
			iterate(0)	///

For alpha = 0, beta = 0, rho = .902872, sigma2 = 4.247643 (the values of sigma1 and sigma2 determined 
	estimates from the run we report. . .hence rho = .902872):
			init(alpha:_cons=.0      	   beta:_cons=.0	 		///
			  lnsigma1:_cons= 2.56114  lnsigma2:_cons=1.446364)	///
			iterate(0)	///


For the run we report:
			init(alpha:_cons=.50      	   beta:_cons=.10	 		///
			  lnsigma1:_cons= 2.3025851  lnsigma2:_cons=2.8518912)	///


			init(alpha:_cons=.50      	   beta:_cons=.10	 		///

difficult


For the actual data (n = 80):
			init(alpha:_cons=.5622838       beta:_cons=.0286859   	 	///
			  lnsigma1:_cons=1.113418  lnsigma2:_cons=.24338796)		///

For simulated data:
			init(alpha:_cons=.50      	   beta:_cons=.10	 		///
			  lnsigma1:_cons= 2.3025851  lnsigma2:_cons=2.8518912)	///



wald(-$X_MVT_NOEQ): The Wald test for 6 equations tries to do a test of a "[]" variable
. . . do not understand why.


****************************/


            /* drop variables used for random draws */
    cap drop ${X_MVT_prefix}*

            /* Prepare output - again mostly taken from -mvprobit */
    eret scalar neqs = $X_MVT_NOEQ
    eret scalar draws = $X_MVT_D
    eret scalar seed = `seed'
    eret scalar ll0 = `ll0'
    eret scalar chi2_c = abs(-2*(e(ll0)-e(ll)))
    eret scalar nrho = (e(neqs)-1)*e(neqs)/2
    eret loc an no
    if "`an'" != "" {
        eret loc an yes
    }

    eret loc cmd "tobitRE"

    forval i = 1/$X_MVT_NOEQ {
        eret loc rhs`i' "`ind`i''"
        eret loc nrhs`i' "`ninds`i''"

        loc t = [lnsigma1]_b[_cons]
        loc tse = [lnsigma1]_se[_cons]
        eret scalar sigma1 =  exp(`t')
        eret scalar sesigma1 = `tse'*exp(`t')

        loc t = [lnsigma2]_b[_cons]
        loc tse = [lnsigma2]_se[_cons]
        eret scalar sigma2 =  exp(`t')
        eret scalar sesigma2 = `tse'*exp(`t')


    }


    Display, `level' 

    /* Taken entirely from mvprobit:
        don't report LR test of "every rho = 0" for models 
       with constraints imposed: (i) constraints on coeffs can't be 
       imposed on the initial -probit- models. (ii) test has no sense 
       if constraints are placed on rhos 
       ==> report test only if no constraints (in which case "Cns" exists)
    */

    tempname c
    capture mat `c' = get(Cns)
    if (_rc != 0 & $X_MVT_NOEQ > 1 ) {
                    /* LR test of models without & without rhos */
        di as txt "Likelihood ratio test of $S_MVT_rho 0:  " 
        di as txt "             chi2(" as res "`e(nrho)'" as txt ") = " /*
            */ as res %8.0g e(chi2_c) _c
        di as txt "   Prob > chi2 = " as res %6.4f /*
            */ chiprob(`e(nrho)', e(chi2_c))
    }

                    /* clear up globals no longer needed */
    foreach g in X_MVT_i X_MVT_C S_MLE_I S_MLE_tvar S_MLE_rho S_MLE_atrho X_MVT_slatrho  { 
            glo `g' 
    } 
    capture macro drop S_MLE_z* 
end


program define Display
    syntax [,Level(int $S_level)]
    ml display, level(`level') neq($X_MVT_NOEQ) plus

			DispA lnsigma1 /lnsigma1 `level'            
			DispA lnsigma2 /lnsigma2 `level'
    
			DispB lnsigma1 sigma1 `level'
			DispB lnsigma2 sigma2 `level'

			DispR lnsigma1 lnsigma2 rho `level'
end


program define DispA
    loc level = `3'
    _diparm `1', level(`level') label("`2'") prob
    di in smcl as txt "{hline 13}{c BT}{hline 64}"
end

program define DispB
    loc level = `3'
    _diparm `1', level(`level') exp label("`2'") prob
    di in smcl as txt "{hline 13}{c +}{hline 64}"
end

program define DispR
    loc level = `4'
    _diparm lnsigma1 lnsigma2, level(95) 						///
		func(exp(@1*2)/(exp(@1*2)+exp(@2*2)))				///
		der(2*exp((@1+@2)*2)/((exp(@1*2)+exp(@2*2)))^2 -2*exp((@1+@2)*2)/((exp(@1*2)+exp(@2*2)))^2)	///
		ci(probit) label("rho") prob
    di in smcl as txt "{hline 13}{c +}{hline 64}"
end
