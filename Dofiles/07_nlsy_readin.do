import delimited "$data/NLSY97/base/iigr.csv", clear
ren *, upper
do "$data/NLSY97/base/iigr-value-labels"
ren *, lower

//renaming
ren sex sex
ren bdate_y birthyr
drop bdate_m
ren pubid uniqid
ren race race
ren cv_hh_income parent_income
ren cv_hh_net parent_wealth
drop cv_highest_degree_ever_edt_1998
ren cv_sample sample
drop asvab_prev

//coding things to missing
ds asvab*, not
foreach var in `r(varlist)'{
	replace `var' = . if `var' == -1
	replace `var' = . if `var' == -2
	replace `var' = . if `var' == -3
	replace `var' = . if `var' == -4
	replace `var' = . if `var' == -5	
}


//renaming
forval y = 1997/2017{
	cap drop cv_hgc_ever_`y'
	cap drop cv_highest_degree_ever_`y'
	
	
	cap ren cv_enroll*`y' enroll_`y'
	
	if _rc{
		cap drop cv_enrollstat_`y'
		cap ren cv_enroll*`y' enroll_`y'
	}
	
	cap ren cv_hgc*`y' hgc_`y'
	cap ren cv_highest*`y' hideg_`y'
	cap ren yinc_1700_`y' income_`y'	
}

//hours/weeks
forval y = 1997/2017{
    local stub = substr("`y'", -2, 2)
	cap ren cvc_wkswk_yr_all_`stub'_xrnd weeks_worked_`y'
	cap ren cvc_hours_wk_yr_all_`stub'_xrnd hours_`y'
}
//drop weeks_worked*
drop cvc*

//student loan stuff
ren yast25_5012_comb_xrnd educ_loan_personal
ren yast25_5016_comb_xrnd educ_loan_govt

replace educ_loan_p = 0  if educ_loan_p<0
replace educ_loan_g = 0  if educ_loan_g<0
gen educ_loans_25 = educ_loan_g + educ_loan_p
drop yast*

drop parent_income
ren parent_wealth parent_wealth_1997
//order uniqid sample sex race birthyr ability parent_wealth hgc* enroll* hideg* income* hours* educ_loan*
save "$temp/nlsy_noweights", replace

*******add on cumulative weights*******
import delimited "$data/NLSY97/weights/weights.csv", clear
ren *, upper
do "$data/NLSY97/weights/weights-value-labels"
ren *, lower

ren pubid uniqid
drop key* cv* *panel*

//renaming
forval y = 1997/2017{
	cap ren sampling_weight_cc_`y' weight_`y'
}
merge 1:1 uniqid using "$temp/nlsy_noweights", keep(match) nogen

//ASVAB stuff
local asvabs `"gs ar wk pc no cs ai si mk mc ei ao"'
foreach asvab in `asvabs'{
	replace asvab_`asvab'_ability_est_pos_1999 = 0 if asvab_`asvab'_ability_est_pos_1999<0
	replace asvab_`asvab'_ability_est_neg_1999 = 0 if asvab_`asvab'_ability_est_neg_1999<0
	gen asvab_`asvab' = asvab_`asvab'_ability_est_pos_1999 - asvab_`asvab'_ability_est_neg_1999
	drop asvab_`asvab'_ability_est_pos_1999 asvab_`asvab'_ability_est_neg_1999
}
drop if asvab_ar == 0 & asvab_pc == 0
gen age_at_asvab = 1997 - birthyr
tab age_at_asvab


//convert to SD measure within age groups
qui{
forval i = 13/17{
	foreach a in `asvabs'{
		su asvab_`a' [fw = weight_1997] if age_at_asvab == `i'
		replace asvab_`a' = (asvab_`a' - `r(mean)') / (`r(sd)') if age_at_asvab == `i'
	}
}
}

pca asvab*
predict pca1, score
xtile ability = pca1 [fw = weight_1997], nq(2) //binary ability
drop asvab* pca* age_at_asvab
compress
order uniqid sample sex race birthyr ability parent_wealth hgc* enroll* hideg* income* hours* educ_loan* weight*
save "$temp/nlsy_weights", replace

********add on 1997 family income*******
import delimited "$data/NLSY97/parent_income/parent_income.csv", clear
ren *, upper
do "$data/NLSY97/parent_income/parent_income-value-labels"
ren *, lower

ren pubid uniqid
drop key* cv_sample
ren cv_income parent_income_1997

//missing codings
ds uniqid*, not
foreach var in `r(varlist)'{
	replace `var' = . if `var' == -1
	replace `var' = . if `var' == -2
	replace `var' = . if `var' == -3
	replace `var' = . if `var' == -4
	replace `var' = . if `var' == -5	
}

merge 1:1 uniqid using "$temp/nlsy_weights", keep(match) nogen
compress
order uniqid sample sex race birthyr ability parent_wealth parent_income educ_loan* hgc* enroll* hideg* income* hours* weight*

//final cleaning
replace parent_income = 0 if parent_income<0

//clean up and consolidate student loan variables
ds educ*
foreach var in `r(varlist)'{
	replace `var' = 0 if `var' == .
}
replace educ_loans_25 = educ_loan_g + educ_loan_p if educ_loans_25 == 0 & (educ_loan_g!=0 | educ_loan_p!=0)
drop educ_loan_g educ_loan_p
save "$temp/nlsy_base_sample", replace



*****************add on 1997 family housing wealth**************
import delimited "$data/NLSY97/housing_wealth/housing_wealth.csv", clear
ren *, upper
do "$data/NLSY97/housing_wealth/housing_wealth-value-labels"
ren *, lower
drop key* cv_sample*
ren pubid uniqid
ren p5_101 ownershp
ren p5_112 house_value
drop p5*
replace house_value = 0 if ownershp == 2 | house_value == -4
replace house_value = . if house_value == -2 | house_value == -1
ren house house_value_parent
ren ownershp ownership_parent
replace ownership_parent = . if ownership_parent<-0
save "$temp/nlsy_parent_house", replace
//p112: value of house
//p113: large range of value of house
//115/116: owed/range owed on mortgage


*****************add on actual wage data**************
import delimited "$data/NLSY97/hourly_pay/hourly_pay.csv", clear
ren *, upper
do "$data/NLSY97/hourly_pay/hourly_pay-value-labels"
ren *, lower
drop key* cv_sample*
ren pubid uniqid

//missing codings
ds uniqid*, not
foreach var in `r(varlist)'{
	replace `var' = . if `var' == -1
	replace `var' = . if `var' == -2
	replace `var' = . if `var' == -3
	replace `var' = . if `var' == -4
	replace `var' = . if `var' == -5	
}

//hours/weeks
forval y = 1997/2017{
	cap ren cv*`y' wage_`y'
	cap replace wage_`y' = wage_`y' / 100
	
	cap qui su wage_`y', d //windsorizing happens here
	cap replace wage_`y' = `r(p99)' if wage_`y' > `r(p99)' & wage_`y'!=. //windsorize at p99
}

save "$temp/nlsy_wages", replace

//all together
merge 1:1 uniqid using "$temp/nlsy_parent_house", keep(match) nogen
merge 1:1 uniqid using "$temp/nlsy_base_sample", keep(match) nogen
drop parent_wealth ownership_parent
ren house_value_parent parent_house_value
drop income* hours*

order uniqid sample sex race birthyr ability parent_house_value parent_income educ_loan* hgc* enroll* hideg* wage* weight*
drop parent* enroll* educ_loan*
save "$temp/nlsy_base_sample", replace //re-write


*****occupation/fields
import delimited "$data/NLSY97/occ_field/tso.csv", clear
ren *, upper
do "$data/NLSY97/occ_field/tso-value-labels"
ren *, lower
drop key* cv_sample* cv_hgc* cv_highest*
ren pubid uniqid

//missing codings
ds uniqid*, not
foreach var in `r(varlist)'{
	replace `var' = . if `var' == -1
	replace `var' = . if `var' == -2
	replace `var' = . if `var' == -3
	replace `var' = . if `var' == -4
	replace `var' = . if `var' == -5	
}

forval y = 1997/2017{
	cap ren ysch*`y' degfield_`y'
	cap ren yemp*`y' occ_`y'
}

ren cvc_trn_cert_xrnd has_license
drop cvc*
merge 1:1 uniqid using "$temp/nlsy_base_sample", keep(match) nogen //merge on other variables

//construct education field and teaching occupation dummies
forval i = 1997/2017{
	cap gen educ_`i' = (degfield_`i' == 12)
	
	//compressed occupation codes
	gen comp_occ_`i' = .
	cap replace comp_occ_`i' = 2 if occ_`i'>=2200 & occ_`i'<=2340 //teachers
	if _rc{
		drop comp_occ_`i'
		continue
	}
	cap replace comp_occ_`i' = 1 if occ_`i'!=. & comp_occ_`i'!=2 //assigned a non-teaching occupation
	cap replace comp_occ_`i' = 0 if weeks_worked_`i' < 10 //if didn't work enough, code as not-workign	
	cap drop occ_`i'
	ren comp_occ_`i' occ_`i'
}

egen degree_educ_ever = rowmax(educ_*)
drop degfield* educ*
order uniqid sample sex race birthyr ability degree_educ_ever hgc* hideg* has_license weeks* wage* occ* weight* 

//reshape to long
reshape long hgc_ hideg_ weeks_worked_ wage_  occ_ weight_, i(uniqid) j(year)

******deflate stuff***********
replace year = year-1 //wages reported are for preceding year
merge m:1 year using "$data/GDP/gdp_pce_deflator", keep(match) nogen
replace deflator = deflator/100 //normalize
replace wage = wage/deflator //deflate labor income
drop deflator
replace year = year+1
gen age = year - birthyr
compress
save "$temp/nlsy_final_long", replace


//end of dofile



//