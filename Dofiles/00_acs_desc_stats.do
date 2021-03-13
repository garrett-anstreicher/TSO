//get statefips
/*
global tempdir "C:\Users\Garrett\Documents\Grad_School\Papers\FRHC\Data\ACS"
use "$tempdir/usa_00025", clear
keep serial pernum year statefip
save "$temp/acs_statefip", replace


use "$data/ACS/usa_00030", clear
merge 1:1 serial pernum year using "$temp/acs_statefip", keep(match) nogen

//sample restrictions
keep if year>=2009 //years with degfield
keep if age>=20 & age<=60

//keep only if individual has a degree
drop if degfield==0 | degfield==.
compress
save "$temp/acs_degfield", replace
*/

use "$temp/acs_degfield", clear
keep if age>=23 & age<=55 //cohort restriction
gen birthyr = year-age
merge 1:1 serial pernum year using "$data/ACS/usa_00036", keep(match) nogen

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

*****begin producing descriptive stats
//summary stats by demographics
eststo wmen: estpost su educ_degree educ_gen educ_elem educ_other teacher teacher_elem ///
teacher_sec teacher_oth hideg masters inctot nowork partime fulltime if !female & !black_hispan [w=perwt]

eststo wwomen: estpost su educ_degree educ_gen educ_elem educ_other teacher teacher_elem ///
teacher_sec teacher_oth hideg masters inctot nowork partime fulltime if female & !black_hispan [w=perwt]

eststo nwmen: estpost su educ_degree educ_gen educ_elem educ_other teacher teacher_elem ///
teacher_sec teacher_oth hideg masters inctot nowork partime fulltime if !female & black_hispan [w=perwt]

eststo nwwomen: estpost su educ_degree educ_gen educ_elem educ_other teacher teacher_elem ///
teacher_sec teacher_oth hideg masters inctot nowork partime fulltime if female & black_hispan [w=perwt]

esttab wmen wwomen nwmen nwwomen using "$output/summary_by_demog.tex", replace cell(mean(fmt(2)) sd(par fmt(2)))

//fraction of teachers with educ degrees, masters degrees
eststo teacher: estpost su educ_degree educ_gen educ_elem educ_other hideg masters if teacher [w=perwt]
eststo teacher_elem: estpost su educ_degree educ_gen educ_elem educ_other hideg masters if teacher_elem [w=perwt]
eststo teacher_sec: estpost su educ_degree educ_gen educ_elem educ_other hideg masters if teacher_sec [w=perwt]
eststo teacher_oth: estpost su educ_degree educ_gen educ_elem educ_other hideg masters if teacher_oth [w=perwt]
esttab teacher teacher_elem teacher_sec teacher_oth using "$output/teacher_educ.tex", replace cell(mean(fmt(2)) sd(par fmt(2)))

//earnings of teachers by age, masters, degree, by occupation
preserve
keep if fulltime==100
collapse (p10) inc_p10 = inctot (p25) inc_p25 = inctot (p50) inc_p50 = inctot ///
(p75) inc_p75 = inctot (p90) inc_p90 = inctot [w=perwt], by(age teacher)

twoway line inc_p10 age if teacher, graphregion(color(white)) bgcolor(white) ///
|| line inc_p25 age if teacher ///
|| line inc_p50 age if teacher ///
|| line inc_p75 age if teacher ///
|| line inc_p90 age if teacher, ///
legend(lab(1 "p10") lab(2 "p25") lab(3 "p50") lab(4 "p75") lab(5 "p90")) xtitle("Age") ytitle("Earnings ($1,000 2012)") ylabel(0(25)200, angle(0) grid)
graph export "$output/earnings_teacher.png", replace

twoway line inc_p10 age if !teacher, graphregion(color(white)) bgcolor(white) ///
|| line inc_p25 age if !teacher ///
|| line inc_p50 age if !teacher ///
|| line inc_p75 age if !teacher ///
|| line inc_p90 age if !teacher, ///
legend(lab(1 "p10") lab(2 "p25") lab(3 "p50") lab(4 "p75") lab(5 "p90")) xtitle("Age") ytitle("Earnings ($1,000 2012)") ylabel(0(25)200, angle(0) grid)
graph export "$output/earnings_nonteacher.png", replace
restore


preserve
keep if fulltime==100
collapse ///
(p10) inc_p10 = inctot ///
(p20) inc_p20 = inctot ///
(p30) inc_p30 = inctot ///
(p40) inc_p40 = inctot ///
(p50) inc_p50 = inctot ///
(p60) inc_p60 = inctot ///
(p70) inc_p70 = inctot ///
(p80) inc_p80 = inctot ///
(p90) inc_p90 = inctot [w=perwt], by(teacher hideg educ_degree)

reshape long inc_p, i(teacher hideg educ_degree) j(decile)

line inc_p decile if teacher & !hideg & educ_degree, graphregion(color(white)) bgcolor(white) ///
|| line inc_p decile if teacher & hideg & educ_degree ///
|| line inc_p decile if teacher & !hideg & !educ_degree ///
|| line inc_p decile if teacher & hideg & !educ_degree, ///
legend(lab(1 "Coll, Edu") lab(2 ">Coll, Edu") lab(3 "Coll, Non-Edu") lab(4 ">Coll, Non-Edu")) xtitle("Decile") ytitle("Earnings ($1,000 2012)") ylabel(0(25)200, angle(0) grid)
graph export "$output/earnings_by_deg_teacher.png", replace

line inc_p decile if !teacher & !hideg & educ_degree, graphregion(color(white)) bgcolor(white) ///
|| line inc_p decile if !teacher & hideg & educ_degree ///
|| line inc_p decile if !teacher & !hideg & !educ_degree ///
|| line inc_p decile if !teacher & hideg & !educ_degree, ///
legend(lab(1 "Coll, Edu") lab(2 ">Coll, Edu") lab(3 "Coll, Non-Edu") lab(4 ">Coll, Non-Edu")) xtitle("Decile") ytitle("Earnings ($1,000 2012)") ylabel(0(25)200, angle(0) grid)
graph export "$output/earnings_by_deg_nonteacher.png", replace
restore

//fraction of educ majors by year of birth
preserve
collapse (mean) educ_degree [w = perwt], by(birthyr female black_hispan)
drop if birthyr<1980 | birthyr>1995
twoway connected educ_degree birthyr if !female & !black_hispan, graphregion(color(white)) bgcolor(white) ///
|| connected educ_degree birthyr if female & !black_hispan ///
|| connected educ_degree birthyr if !female & black_hispan ///
|| connected educ_degree birthyr if female & black_hispan, ///
legend(lab(1 "White Men") lab(2 "White Women") lab(3 "Nonwhite Men") lab(4 "Nonwhite Women")) xtitle("Birth Year") ytitle("% Education Majors") 
graph export "$output/teachers_by_birthyr.png", replace
restore

//employment probabilities by age, educ major, sex, race
preserve
//collapse work categories
gen work_n = (nowork==100)*100
gen work_teach = (nowork==0 & teacher==100) * 100
gen work_nonteach = (nowork==0 & teacher!=100) * 100
collapse (mean) work_n work_teach work_nonteach [w = perwt], by(age educ_degree female black_hispan)
sort female black_hispan educ_degree age 
save "$temp/employment_probabilities", replace
restore

//age profiles of occupation choice
gen work_n = (nowork==100)*100
gen work_teach = (nowork==0 & teacher==100) * 100
gen work_nonteach = (nowork==0 & teacher!=100) * 100
collapse (mean) work_n work_teach work_nonteach [w = perwt], by(age female)
twoway line work_n age if female|| line work_teach age if female || line work_nonteach age if female, ///
bgcolor(white) graphregion(color(white)) ylabel(0(20)100, angle(0) grid) legend(lab(1 "Not Working") lab(2 "Teacher") lab(3 "Non-Teacher"))
graph export "$output/occupation_age_profile_female.png", replace

twoway line work_n age if !female|| line work_teach age if !female || line work_nonteach age if !female, ///
bgcolor(white) graphregion(color(white)) ylabel(0(20)100, angle(0) grid) legend(lab(1 "Not Working") lab(2 "Teacher") lab(3 "Non-Teacher"))
graph export "$output/occupation_age_profile_male.png", replace




//end of dofile