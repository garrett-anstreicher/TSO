//read in big data and shrink
/*
use "$data/ACS/usa_00030", clear
keep if age>=20 & age<=60 //age restriction
drop if degfield == 0
keep if educ>=10
compress
save "$temp/ACS_smaller", replace
*/


use "$temp/ACS_smaller", clear
drop if year == 2018 //occupation codes swithced he4e

//RHS variables
gen teacher = (occ >=2200 & occ<=2340)
gen educ_degree = (degfield == 23) //education degree dummy
tab teacher if educ_degree
tab educ_degree if teacher
gen ma = (educd==114) //ma dummy
drop if educd>114 //nix PHDs
replace sex = sex - 1 //sex dummy
gen exp = age-22 //years of potential experience
drop if exp<0

gen black_hispan = 0
replace black_hispan = 1 if (race == 2 | hispan>0)

//work requirements
gen nwork = 0
replace nwork = 1 if empstat > 1
replace nwork = 1 if uhrswork<30
replace nwork = 1 if wkswork2<4
replace nwork = 1 if inctot<0

//deflate and windsorize earnings
merge m:1 year using "$data/GDP/gdp_pce_deflator", keep(match) nogen
replace deflator = deflator/100 
replace inctot = inctot/deflator //deflate earnings
su inctot, d
replace inctot = `r(p99)'*1.5 if inctot>`r(p99)'
replace inctot = inctot / 2080 //convert to hourly wages as though worked full-time
replace inctot = log(inctot)

//get estimates
reg inctot ma c.exp##c.exp##c.exp i.year [aw=perwt] if teacher, robust
predict resid, resid 
su resid if teacher
drop resid

exit

reg inctot sex black_hispan educ_degree ma c.exp##c.exp##c.exp i.year [aw=perwt] if !teacher, robust
predict resid, resid 
su resid if !teacher

//