//import and format psid data
clear all
do "$data/PSID/J287545"
do "$data/PSID/J287545_formats"

gen uniqid = ER30001*1000 + ER30002 //generation of uniqid identifier variable
duplicates report uniqid //no duplicates!
drop ER30002 ER30000 
ren ER32000 sex
drop V3529

//drop variables we don't care about
//loop over variables and drop what we can
ds uniqid, not
foreach var in `r(varlist)'{
	
	//fetch label
	local lab: variable label `var'
	
	//kill sequence/release/interview numbers
	if strpos("`lab'", "SEQUENCE") | strpos("`lab'", "RELEASE") | strpos("`lab'", "INTERVIEW")  | strpos("`lab'",  "WEALTH") | strpos("`lab'",  "HOME") | strpos("`lab'",  "EMPL"){
		drop `var'
		continue
	}
	
	//get year for individual files
	local num = word("`lab'", -1)
	
	//skip variables with labels that i don't like
	if length("`num'")!=2 | "`var'" == "V4373" {
		continue
	}	
	
	if `num'<20{ //2000s variable
		local year = 2000 + `num'
	}
	if `num'>20{ //1900s variable
		local year = 1900 + `num'
	}
	
	****rename according to label contents
	//age
	if strpos("`lab'", "AGE OF"){
		ren `var' age_`year'
	}
	
	//relation to head
	if strpos("`lab'", "RELATION"){
		ren `var' relate_`year'
	} 
	
	//occupation/industry OF wife/head
	if strpos("`lab'", "OCCUPATION OF HEAD"){
		ren `var' occ_head_`year'
	} 
	
	if strpos("`lab'", "INDUSTRY   OF HEAD"){
		ren `var' ind_head_`year'
	} 
	
	//
	if strpos("`lab'", "OCCUPATION OF WIFE"){
		ren `var' occ_sp_`year'
	} 
	
	if strpos("`lab'", "INDUSTRY   OF WIFE"){
		ren `var' ind_sp_`year'
	}
}

ds sex occ* ind* relate* age* uniqid, not
local vars "`r(varlist)'"
local counter 0
forval y = 1981/2001{
	if `y' == 1998 | `y' == 2000 {
		continue
	}
	local start = 4 * `counter'
	local start1 = `start' + 1
	local start2 = `start' + 2
	local start3 = `start' + 3
	local start4 = `start' + 4
	local counter `++counter'
	
	//rename
	local var = word("`vars'", `start1')
	ren `var' occ_head_`y'
	
	local var = word("`vars'", `start2')
	ren `var' ind_head_`y'
	
	local var = word("`vars'", `start3')
	ren `var' occ_sp_`y'
	
	local var = word("`vars'", `start4')
	ren `var' ind_sp_`y'
}

//rename 2003 onward
ds sex occ* ind* relate* age* uniqid, not
local vars "`r(varlist)'"
local counter 0
local years `"2003 2005 2007 2009 2011 2013 2015 2017"'

foreach y in `years'{
	local start = 4 * `counter'
	local start1 = `start' + 1
	local start2 = `start' + 2
	local start3 = `start' + 3
	local start4 = `start' + 4
	local counter `++counter'
	
	//rename
	local var = word("`vars'", `start1')
	ren `var' occ_sp_`y'
	
	local var = word("`vars'", `start2')
	ren `var' ind_sp_`y'
	
	local var = word("`vars'", `start3')
	ren `var' occ_head_`y'
	
	local var = word("`vars'", `start4')
	ren `var' ind_head_`y'
}

//reshape 
reshape long age_ occ_sp_ ind_sp_ ind_head_ occ_head_ relate_, i(uniqid) j(year)
replace relate = relate*10 if relate<10
keep if relate == 10 | relate == 20
gen occ = . 
replace occ = occ_head if relate == 10
replace occ = occ_sp if relate == 20

gen ind = . 
replace ind = ind_head if relate == 10
replace ind = ind_sp if relate == 20
ren age age 
ren relate relate
keep uniqid year age relate occ ind sex
save "$temp/psid_occ_long", replace

//code teachers
use "$temp/psid_occ_long", clear
drop if year == 2017 //weird occupation codes
gen unemp = (occ == 0) //simple hack for now; can validate later


//coding
gen teacher = 0
replace teacher = 1 if occ>=142 & occ<=145 & year<=2001 //1970 codes
replace teacher = 1 if oc>=230 & occ<=234 & year>=2003 //2000 codes

//simplified tabulation of occupation
tab teacher
gen occ_simple = .
replace occ_simple = 0 if unemp
replace occ_simple = 1 if !teacher & !unemp
replace occ_simple = 2 if teacher
tab sex if teacher //sanity check

//run analyses
keep if age>=20 & age<=60

****tabulations of switching based on 10-year age bins

qui{
//whole sample
preserve
gen temp = (occ_simple == 1)
bys uniqid: egen ever_other = max(temp)
drop temp
gen temp = (occ_simple == 2)
bys uniqid: egen ever_teach = max(temp)
gen switcher = (ever_other & ever_teach) //individual observed both teaching and non-teaching
collapse (mean) switcher, by(uniqid)
noi tab switcher
restore

//20-30
preserve
keep if age>=20 & age<=30
gen temp = (occ_simple == 1)
bys uniqid: egen ever_other = max(temp)
drop temp
gen temp = (occ_simple == 2)
bys uniqid: egen ever_teach = max(temp)
gen switcher = (ever_other & ever_teach) //individual observed both teaching and non-teaching
collapse (mean) switcher, by(uniqid)
noi tab switcher
restore

//30-40
preserve
keep if age>=30 & age<=40
gen temp = (occ_simple == 1)
bys uniqid: egen ever_other = max(temp)
drop temp
gen temp = (occ_simple == 2)
bys uniqid: egen ever_teach = max(temp)
gen switcher = (ever_other & ever_teach) //individual observed both teaching and non-teaching
collapse (mean) switcher, by(uniqid)
noi tab switcher
restore

//40-50
preserve
keep if age>=40 & age<=50
gen temp = (occ_simple == 1)
bys uniqid: egen ever_other = max(temp)
drop temp
gen temp = (occ_simple == 2)
bys uniqid: egen ever_teach = max(temp)
gen switcher = (ever_other & ever_teach) //individual observed both teaching and non-teaching
collapse (mean) switcher, by(uniqid)
noi tab switcher
restore

//50-60
preserve
keep if age>=50 & age<=60
gen temp = (occ_simple == 1)
bys uniqid: egen ever_other = max(temp)
drop temp
gen temp = (occ_simple == 2)
bys uniqid: egen ever_teach = max(temp)
gen switcher = (ever_other & ever_teach) //individual observed both teaching and non-teaching
collapse (mean) switcher, by(uniqid)
noi tab switcher
restore

}





