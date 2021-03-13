//read in ACS data; get ratio of teacher salaries vs. non-teachers by state and year

/*
use "$data/ACS/usa_00030", clear
gen teacher = (occ>=2300 & occ<=2340)*100 //teacher
gen teacher_d = 0
replace teacher_d = 1 if occ == 2310 //elementary/middle school
replace teacher_d = 2 if occ == 2320 //secondary
replace teacher_d = 3 if teacher & teacher_d == 0   //other teacher
gen teacher_elem = (teacher_d == 1)*100
gen teacher_sec = (teacher_d == 2)*100
gen teacher_oth = (teacher_d == 3)*100

//get state fips
merge 1:1 serial pernum year using "$data/ACS/acs_0518_bpl", keep(match) nogen
save "$temp/acs_0518_statefip", replace
*/


//clean
use "$temp/acs_0518_statefip", clear
keep if age>=23 & age<=55 //cohort restriction
keep if educ>=10
gen birthyr = year-age

//deflate earnings
merge m:1 year using "$data/GDP/gdp_pce_deflator", keep(match) nogen
replace inctot = inctot/deflator * 100
replace inctot = 0 if inctot<0
su inctot, d
replace inctot = `r(p99)' if inctot>`r(p99)'
gen count = 1

//labor supply
gen lab_supply = .
replace lab_supply = 0 if empstat>1 //not working
replace lab_supply = 2 if empstat == 1 & wkswork2>=4 & uhrswork>30 //full time: at least 40 weeks and 30 hours
replace lab_supply = 1 if empstat == 1 & lab_supply == . //part time: working but don't fullfill above requiremnets
tab lab_supply
gen nowork = (lab_supply == 0)*100
gen partime = (lab_supply == 1)*100
gen fulltime = (lab_supply == 2)*100
//keep if lab_supply==2 //keep full-time workers
keep if empstat == 1

//collapse by year, state, 23-
preserve
collapse (sum) count (mean) inctot [w = perwt], by(teacher statefip year)
drop count
save "$temp/wages_all_0518", replace
reshape wide inctot, i(year statefip) j(teacher)
gen ratio = inctot0 / inctot1
keep year statefip ratio
ren ratio ratio_all
save "$temp/ratios_all", replace
restore

//23-35
keep if age>=23 & age<=35
collapse (sum) count (mean) inctot [w = perwt], by(teacher statefip year)
drop count
save "$temp/wages_2335_0518", replace
reshape wide inctot, i(year statefip) j(teacher)
gen ratio = inctot0 / inctot1
keep year statefip ratio
ren ratio ratio_2335
save "$temp/ratios_2335", replace


/*
//clean
use "$temp/acs_degfield", clear
keep if age>=23 & age<=55 //cohort restriction
gen birthyr = year-age
merge 1:1 serial pernum year using "$data/ACS/usa_00036", keep(match) nogen
merge 1:1 serial pernum year using "$data/ACS/usa_00035", keep(match) nogen //get bpl

***code up some variables***

//female dummy
replace sex = sex - 1
ren sex female 

//black/hispanic
gen black_hispan = 0
replace black_hispan = 1 if race == 2 | hispan > 0  

//education degree and detailed version
gen educ_degree = (degfield == 23 | degfield2 == 23)*100
gen educ_degree_d = 0
replace educ_degree_d = 1 if degfieldd == 2300 | (degfield2d == 2300) //general
replace educ_degree_d = 2 if degfieldd == 2304 | (degfield2d == 2304) //elementary
replace educ_degree_d = 2 if degfield2d == 2304 & educ_degree_d == 0
replace educ_degree_d = 3 if educ_degree_d == 0 & educ_degree //all other
gen educ_general = (educ_degree_d == 1)*100
gen educ_elem = (educ_degree_d == 2)*100
gen educ_other = (educ_degree_d == 3)*100

//occupation
gen teacher = (occ>=2300 & occ<=2340)*100 //teacher
gen teacher_d = 0
replace teacher_d = 1 if occ == 2310 //elementary/middle school
replace teacher_d = 2 if occ == 2320 //secondary
replace teacher_d = 3 if teacher & teacher_d == 0   //other teacher
gen teacher_elem = (teacher_d == 1)*100
gen teacher_sec = (teacher_d == 2)*100
gen teacher_oth = (teacher_d == 3)*100

//education: has a higher degre
gen hideg = 0
replace hideg = 100 if educ == 11
gen masters = (educd == 114)*100 //specifically a masters

//earnings
merge m:1 year using "$data/GDP/gdp_pce_deflator", keep(match) nogen
replace inctot = inctot/deflator * 100
replace inctot = 0 if inctot<0
su inctot, d
replace inctot = `r(p99)' if inctot>`r(p99)'
su inctot, d
replace inctot = inctot/1000

//labor supply
gen lab_supply = .
replace lab_supply = 0 if empstat>1 //not working
replace lab_supply = 2 if empstat == 1 & wkswork2>=4 & uhrswork>30 //full time: at least 40 weeks and 30 hours
replace lab_supply = 1 if empstat == 1 & lab_supply == . //part time: working but don't fullfill above requiremnets
tab lab_supply
gen nowork = (lab_supply == 0)*100
gen partime = (lab_supply == 1)*100
gen fulltime = (lab_supply == 2)*100

save "$temp/acs_degfield_cleaned", replace
*/

//do the thing
use "$temp/acs_degfield_cleaned", clear
ren statefip state_live
ren bpl statefip
gen year_16 = year - (age-16) //year when respondent was 16
gen year_18 = year - (age-18) //year when respondent was 18
gen year_20 = year - (age-20) //year when respondent was 20
gen year_22 = year - (age-22) //year when respondent was 22
ren year year_current

local years `"16 18 20 22"'
foreach year in `years'{
		ren year_`year' year 
		
		//year variables
		merge m:1 year statefip using "$temp/ratios_all", keep (1 3) nogen
		ren ratio_all ratio_all_`year'
		
		merge m:1 year statefip using "$temp/ratios_2335", keep (1 3) nogen
		ren ratio_2335 ratio_2335_`year'

		//rename
		ren year year_`year'
}

//start playing with regressions
lab var female "Female"
lab var black_hispan "Black or Hispanic"
ds ratio*
foreach var in `r(varlist)'{
	lab var `var' "Ratio of NT Earnings to T"
}

local years `"16 18 20 22"'
local counter = 0
local replace `"replace"'
local controls `"female black_hispan i.statefip i.birthyr i.year_current"'

//main loop
foreach year in `years'{
	local counter `++counter'
	if `counter'>1{
		local replace `""'
	}
	
	//23-35 salary ratio
	ren ratio_2335_`year' ratio
	reg educ_degree ratio `controls' [w=perwt], cl(statefip)
	outreg2 using "$output/acs_DiD", `replace' tex lab keep(ratio female black_hispan)
	ren ratio ratio_2335_`year'
	
	//full life-cycle salary ratio
	ren ratio_all_`year' ratio
	reg educ_degree ratio `controls' [w=perwt], cl(statefip)
	outreg2 using "$output/acs_DiD", tex lab keep(ratio female black_hispan)
	ren ratio ratio_all_`year'	
}




local years `"16 18 20 22"'
local counter = 0
local replace `"replace"'
local controls `"female black_hispan i.statefip i.birthyr i.year_current"'

//main loop
foreach year in `years'{
	local counter `++counter'
	if `counter'>1{
		local replace `""'
	}
	
	//23-35 salary ratio
	ren ratio_2335_`year' ratio
	reg teacher ratio `controls' [w=perwt], cl(statefip)
	outreg2 using "$output/acs_DiD_teach", `replace' tex lab keep(ratio female black_hispan)
	ren ratio ratio_2335_`year'
	
	//full life-cycle salary ratio
	ren ratio_all_`year' ratio
	reg teacher ratio `controls' [w=perwt], cl(statefip)
	outreg2 using "$output/acs_DiD_teach", tex lab keep(ratio female black_hispan)
	ren ratio ratio_all_`year'	
}

//note: ratio is non-teacher over teacher. So higher level means that teaching is worse in comparison

/*
So I'm starting to go back through this teacher literature. I will write something up soon, and likely ask Jane to help me compile articles. 

Some of the papers use district level data. In the interest of moving on this, let me just write this. Really only the conceivable thing we can do with the ACS and Census data is something like this:

a) at state (or msa or whatever), construct ratio of avg. teacher salaries to non-teacher, college educated salaries

b) regress fraction of college cohort choosing education major on this ratio, include time and state fixed effects

The tricky part is i) timing -- which year for salaries to use, how forward looking are college cohort? and ii) which moments of the life-cycle earnings to use--could use starting earnings (at ages 23-25, say), could also include later earnings, or some function of it.

We can also look at fraction teaching (using occupation) rather than education major as well. And we could also divide by elementary vs. secondary, or public vs. private (using sector).

How does this sound as a start? I don't take this seriously as a "causal effect"--this is descriptive. We can also consider (?) instruments for salary--like school age population as demand for teachers shifter.
*/






//end of dofile