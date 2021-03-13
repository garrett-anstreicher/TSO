//do stuff for 05-18
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
keep if lab_supply==2 //keep full-time workers
//keep if empstat == 1

//collapse by year, state, 23-
preserve
collapse (sum) count (mean) inctot [w = perwt], by(teacher statefip year)
drop count
save "$temp/wages_all_0518", replace
restore

//23-35
keep if age>=23 & age<=35
collapse (sum) count (mean) inctot [w = perwt], by(teacher statefip year)
drop count
save "$temp/wages_2335_0518", replace

//1990/2000
use "$data/ACS/census_occ", clear
keep if age>=23 & age<=55 //cohort restriction
keep if educ>=10
gen teacher = 0
replace teacher = 100 if year == 1990 & (occ>=155 & occ<=159)
replace teacher = 100 if year == 2000 & (occ>=230 & occ<=234)


merge m:1 year using "$data/GDP/gdp_pce_deflator", keep(match) nogen
replace inctot = inctot/deflator * 100
replace inctot = 0 if inctot<0
su inctot, d
replace inctot = `r(p99)' if inctot>`r(p99)'

//labor supply
gen lab_supply = .
replace lab_supply = 0 if empstat>1 //not working
replace lab_supply = 2 if empstat == 1 & wkswork1>=40 & uhrswork>30 //full time: at least 40 weeks and 30 hours
replace lab_supply = 1 if empstat == 1 & lab_supply == . //part time: working but don't fullfill above requiremnets
tab lab_supply
gen nowork = (lab_supply == 0)*100
gen partime = (lab_supply == 1)*100
gen fulltime = (lab_supply == 2)*100
keep if lab_supply==2 //keep full-time workers

//collapsing
preserve
collapse (mean) inctot [w = perwt], by(teacher statefip year)
save "$temp/wages_all_9000", replace
restore

//23-35
keep if age>=23 & age<=35
collapse (mean) inctot [w = perwt], by(teacher statefip year)
save "$temp/wages_2335_9000", replace



//end of odfile