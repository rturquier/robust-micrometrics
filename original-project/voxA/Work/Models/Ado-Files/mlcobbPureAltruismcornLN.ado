        program mlcobbPureAltruismcornLN
                version 10
                args lnf alpha LNsigma

        // Some temporary variables
		    tempvar g res sigma
			
			quietly generate double `sigma' = exp(`LNsigma')
			
		    quietly generate double  `g' = (`alpha')*z - Ggov

		    quietly generate double  `res' = h - `g'


		    quietly replace `lnf' = -0.5*ln(2*_pi)-ln(`sigma')-0.5*`res'^2/`sigma'^2

		    quietly replace `lnf' =     ln(normal(-`g'/`sigma'))        if h==0
		    quietly replace `lnf' = ln( 1 - normal((40 -`g')/`sigma') ) if h==40 & y==40
		    quietly replace `lnf' = ln( 1 - normal((46 -`g')/`sigma') ) if h==46 & y==46

            end
