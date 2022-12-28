
/*Ado File for updated command*/	
run "C:\Users\tinl\Documents\GitHub\dca.stata\dca.ado"

/*Format for graphs to make it a bit clearer*/
	global format=`"legend(size(small) cols(1) textwidth(100)) scheme(s1color) ylabel(, format("%9.2f")) xlabel(, format("%9.2f"))"'
/*Using *Updated* dataset*/
	use "O:\Outcomes\Andrew\Methodology Work\Vickers dca function updates\Example Dataset (Updated).dta", clear

**********************
*** UNIVARIATE DCA ***
**********************
logit cancer famhistory
printmodel, text

dca cancer famhistory, ${format}
dca cancer famhistory, xstop(0.35) xlabel(0(0.05)0.35) ${format} //limiting x-axis

************************
*** MULTIVARIATE DCA ***
************************

/* Evaluation of New Model */
logit cancer marker age famhistory
// predict cancerpredmarker 	//variable already exists in dataset
dca cancer cancerpredmarker famhistory, xstop(0.35) xlabel(0(0.05)0.35) ${format}

/* Evaluation of Published models */
gen logodds_Brown =  0.75*(famhistory)+ 0.26*(age) - 17.5
gen phat_Brown = invlogit(logodds_Brown)
label var phat_Brown "Brown Model"
dca cancer phat_Brown, xstop(0.35) ${format} xlabel(0(0.05)0.35)

/*Joint or Conditional Tests */
gen high_risk = risk_group == "high"
label var high_risk "Treat Only High Risk Group"
gen joint = risk_group == "high" | cancerpredmarker > .15
label var joint "Treat via Joint Approach"
gen conditional = risk_group == "high" | risk_group == "intermediate" & cancerpredmarker > .15
label var conditional "Treat via Conditional Approach"

dca cancer high_risk joint conditional, xstop(0.35) xlabel(0(0.05)0.35) ${format}

/* Incorporating Harms into Model Assessment */
local harm_marker = 0.0333
gen intermediate_risk = risk_group =="intermediate"
sum intermediate_risk
local harm_conditional = r(mean)*`harm_marker'
dca cancer high_risk joint conditional, harm(0 `harm_marker' `harm_conditional') xstop(0.35) xlabel(0(0.05)0.35) ${format}

/*Saving out Net Benefit Values*/
dca cancer marker, prob(no) xstart(0.05) xstop(0.35) xby(0.05) saving("DCA marker (dca).dta")

/*Interventions Avoided*/
dca cancer marker, prob(no) inter xstart(0.05) xstop(0.35) xby(0.05) xlabel(0.05(0.05)0.35) ${format} ylabel(,format("%9.0f")) saving("intervention", replace)

*****************************
*** Survival Outcomes DCA ***
*****************************
/*Do File for "stdca" command*/	
	run "O:\Outcomes\Andrew\Methodology Work\Vickers stdca function updates\Stata\stdca.do"

/*Format for graphs to make it a bit clearer*/
	global format=`"legend(size(small) cols(1) textwidth(100)) scheme(s1color) ylabel(, format("%9.2f")) xlabel(, format("%9.2f"))"'

/*Using *Updated* (Clean) dataset (so created variables for our Binary Outcome section are erased*/ 
	use "O:\Outcomes\Andrew\Methodology Work\Vickers dca function updates\Example Dataset (Updated).dta", clear

/*Basic Set-up */
	stset ttcancer, f(cancer)

	stcox age famhistory marker, basesurv(surv_func)	//saved baseline survival function
	predict xb, xb										//linear prediction

	sum surv_func if _t <= 1.5							//taking all survival times less than 1.5
	local base = r(min)									//lowest survival probability is the survival probability closes to 1.5 year

	g pr_failure18 = 1 - `base'^exp(xb)					//probability of failure
	label var pr_failure18 "Probility of Failure at 18 months"

/*Multivariable DCA*/
	stdca pr_failure18, timepoint(1.5) xstop(.5) smooth ${format} saving("survival mult", replace)

/*DCA with Competing Risks*/
	*'Creating' status = cancercr	
	gen status = 0
	replace status = 1 if cancer == 1
	replace status = 2 if cancer==0 & dead == 1
	
stset ttcancer, f(status=1)

stdca pr_failure18, timepoint(1.5) compet1(2) xstop(.5) smooth ${format} saving("survival multi compete", replace)

/*Showing both (Kaplan Meier & Competing Risk) together on the same graph*/

	* start with the standard Kaplan Meier model, saving the results to a temporary file
	stset ttcancer, f(cancer)
	tempfile km
	stdca pr_failure18, timepoint(1.5) xstop(.5) nograph saving(`km')  
	
	* now do the competing risk model, again saving the results to a temporary file
	stset ttcancer, f(status=1)
	tempfile cr
	stdca pr_failure18, timepoint(1.5) compet1(2) xstop(.5) nograph saving(`cr') 
	
	* using the temp files, sort by threshold (for merging later) & rename all the models so that we can distinguish them *
	use `km', clear
	sort threshold
	rename pr_failure18 kmmodel
	label var kmmodel "Kaplan Meier: Pr(Failure) at 1.5 years"
	rename all kmall
	label var kmall "Kaplan-Meier: Treat All"
	tempfile kmsort
	save `kmsort'
	
	use `cr', clear
	sort threshold
	rename pr_failure18 crmodel
	label var crmodel "Competing Risk: Pr(Failure) at 1.5 years"
	rename all crall
	label var crall "Competing Risk: Treat All"
	
	* Merging File*
	merge 1:1 threshold using `kmsort'
	
	* Create Graph *
	twoway (line kmall crall threshold if kmall>-0.05 & crall > -0.05, sort) || (line kmmodel crmodel none threshold, sort ${format} ytitle("Net Benefit"))

************************
*** Case Control DCA ***
************************
/*Do File for "dca" command*/	
	run "O:\Outcomes\Andrew\Methodology Work\Vickers dca function updates\Stata Code\dca.do"
/*Format for graphs to make it a bit clearer*/
	global format=`"legend(size(small) cols(1) textwidth(100)) scheme(s1color) ylabel(, format("%9.2f")) xlabel(, format("%9.2f"))"'
/*Using *Updated* dataset*/
	use "O:\Outcomes\Andrew\Methodology Work\Vickers dca function updates\Example Dataset (Updated).dta", clear

	drop if casecontrol == 0		//only use the patients that were in the casecontrol study
	logit cancer famhistory age
	predict xb, xb					//log odds of diaseas

	local true = 0.05
	sum cancer
	local design = r(mean)
	local Bayes=log((`true'/(1-`true'))/(`design'/(1-`design')))
	replace xb=xb+`Bayes'
	
	g phat=invlogit(xb)
	dca phat phat, xstop(0.35) xlabel(0(0.05)0.35) ${format} 

*********************************
*** DCA Correction for Overfit***
*********************************
/*Do File for "dca" command*/	
	run "O:\Outcomes\Andrew\Methodology Work\Vickers dca function updates\Stata Code\dca.do"
/*Format for graphs to make it a bit clearer*/
	global format=`"legend(size(small) cols(1) textwidth(100)) scheme(s1color) ylabel(, format("%9.2f")) xlabel(, format("%9.2f"))"'
/*Using *Updated* dataset*/
	use "O:\Outcomes\Andrew\Methodology Work\Vickers dca function updates\Example Dataset (Updated).dta", clear

	/*Overall command to do loop 200 times*/
	forvalues i=1(1)200 {
		/* Local macros to store event of interest, predictors, and names of model. */
		local event="cancer"
		local predictors1 = "age famhistory"
		local predictors2 = "age famhistory marker"
		local prediction1 = "base"
		local prediction2 = "full"

		/* Variable which will store probabilities from prediction model*/
		quietly g `prediction1'=.
		quietly g `prediction2'=.
		
		/*Assign Patients to one of 10 groups*/
		quietly g u = uniform()
		sort `event' u		//sorted such that we will have an equal number of events in each group
		g group = mod(_n, 10) + 1

			forvalues num = 1/2 { 	//this loop is just so we can do the base and full model together
				forvalues j=1(1)10 {
					quietly logit `event' `predictors`num'' if group!=`j'
					quietly predict ptemp if group==`j'
					quietly replace `prediction`num'' = ptemp if group==`j'
					drop ptemp
				}
			}
		
	/* Work out decision curve for data set created by Cross Validation */
	tempfile dca`i'
	quietly dca `event' `prediction1' `prediction2', nograph xstop(.5) saving("`dca`i''")

	drop u group `prediction1' `prediction2'
}

/* Using the tempfiles saved out, combine all 200 iterations */
	use "`dca1'", clear
	forvalues i=2(1)200 {
		append using "`dca`i''"
		}

	/* Corrected Net Benefits: Mean of the 200 Iterations */
	collapse all none base full base_i full_i, by(threshold)
	
	/* Labeling the variables so the legend looks better for the graph */
	label var all "(Mean) Net Benefit: Treat All"
	label var none "(Mean) Net Benefit: Treat None"
	label var base "(Mean) Net Benefit: Base Model"
	label var full "(Mean) Net Benefit: Full Model"
	label var base_i "(Mean) Intervention: Base Model"
	label var full_i "{Mean) Intervention: Full Model"
	
	/* Saving the Corrected Net Benefits */
	save "Cross Validation DCA output.dta", replace

	/*Creating the Graph*/
	twoway (line all threshold if all>-0.05, sort) || (line none base full threshold, sort ${format} ytitle("Net Benefit"))



