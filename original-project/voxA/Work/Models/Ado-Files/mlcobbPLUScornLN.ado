        program mlcobbPLUScornLN
                version 10
                args lnf alpha beta LNsigma

        // Some temporary variables
		    tempvar B Hplus gplus res sigma
			
			quietly generate double `sigma' = exp(`LNsigma')
			
		    quietly generate double `B' = (1 - `beta')*Ggov + (`alpha' + `beta')*z	
                quietly generate double `Hplus' = .5*(`B' + sqrt(`B'^2 - 4*`alpha' * Ggov * z) )
		    quietly generate double  `gplus' =  `Hplus' - Ggov

		    quietly generate double  `res' = h - `gplus'


		    quietly replace `lnf' = -0.5*ln(2*_pi)-ln(`sigma')-0.5*`res'^2/`sigma'^2

		    quietly replace `lnf' =     ln(normal(-`gplus'/`sigma'))        if h==0
		    quietly replace `lnf' = ln( 1 - normal((40 -`gplus')/`sigma') ) if h==40 & y==40
		    quietly replace `lnf' = ln( 1 - normal((46 -`gplus')/`sigma') ) if h==46 & y==46

            end
