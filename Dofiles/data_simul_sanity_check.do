clear all
//read in and rename
qui{
cap log close
log using "$temp/sanity_check", replace
import delimited "$dir/Model/test_data.csv", clear
ren v1 id
ren v2 period
ren v3 cq
ren v4 gender
ren v5 race
ren v6 ability
ren v7 major
ren v8 ma
ren v9 license
ren v10 e1
ren v11 e2
ren v12 va_type
ren v13 job_last
ren v14 job
ren v15 wage
ren v16 va

//labels
lab def gen_lab 0 "Male" 1 "Female"
lab val gender gen_lab

lab def race_lab 0 "Non Black/Hispanic" 1 "Black/Hispanic"
lab val race race_lab

lab def ability_lab 0 "Low" 1 "High"
lab val ability ability_lab

lab def major_lab 0 "Non-Education" 1 "Education"
lab val major major_lab

lab def ma_lab 0 "No Educ MA" 1 "Educ MA"
lab val ma ma_lab

lab def license_lab 0 "No License" 1 "License"
lab val license license_lab

lab def va_type_lab -2 "Low" 0 "Middle" 2 "High"
replace va_type = va_type*10
gen int va_int = round(va_type)
lab val va_int va_type_lab

lab def job_lab 0 "Not Working" 1 "Non-Teaching" 2 "Teaching"
lab val job job_lab

//currently agenty types are sampled entirely randomly, so we have equal proportions of gender, race, etc.
noi di "################"
noi di "TABULATIONS OF TYPES (CURRENTLY, SIMULATION DRAWS EVENLY)"
noi tab gender //0 = man. 1 = woman.
noi di ""
noi tab race //0 = not black/hispanic. 1 = black/hispanic
noi di ""
noi tab ability //0=low. 1 = high
noi di ""
noi di ""
noi di "################"
noi di ""

//major
noi di "TABULATION OF MAJOR CHOICES"
noi tab major //0 = not teaching. 1 = teaching
noi di ""

noi di "TABULATION OF MAJOR CHOICES BY ABILITY"
noi di "LOW"
noi tab major if ability == 0
noi di ""

noi di "HIGH"
noi tab major if ability == 1 //higher ability types more frequently get non-teaching majors.
noi di ""

noi di "TABULATION OF MAJOR CHOICES BY GENDER"
noi di "MEN"
noi tab major if gender == 0
noi di ""

noi di "WOMEN"
noi tab major if gender == 1 //higher ability types more frequently get non-teaching majors.

noi di ""
noi di "################"
noi di ""

//MA
noi di "TABULATIONS OF EDUC MA/LICENSES"
noi di "MA"
noi tab ma
noi di ""
noi di "MA, NO TEACHING UNDERGRAD"
noi tab ma if major == 0
noi di ""
noi di "MA, TEACHING UNDERGRAD"
noi tab ma if major == 1 //people who get teaching MAs are overwhelmingly people who majored in education
noi di ""

//license
noi di "LICENSE"
noi tab license
noi di ""
noi di "LICENSE, NO TEACHING UNDERGRAD"
noi tab license if major == 0
noi di ""
noi di "LICENSE, TEACHING UNDERGRAD"
noi tab license if major == 1 //people who get teaching licenses are overwhelmingly people who majored in education

noi di ""
noi di "################"
noi di ""

//VA
noi di "TABULATION OF VA TYPES"
noi tab va_int
noi di ""
noi di "AMONG TEACHERS"
noi tab va_int if job == 2 //people working in teaching on average are better teachers.
noi di ""
noi di "AMONG NON-TEACHERS"
noi tab va_int if job!=2 

noi di ""
noi di "################"
noi di ""

//wages
noi di "AVERAGE WAGES"
noi su wage //average earnings, normalized to an hourly wage measure. reasonable

noi di ""
noi di "################"
noi di ""

//jobs held

noi di "TABULATIONS OF JOBS HELD"
noi tab job //0 = not working. 1 = non-teaching. 2 = teaching.
noi di ""

noi di "AMONG PEOPLE WITH NON-EDUCATION MAJORS"
noi tab job if major == 0
noi di ""

noi di "AMONG EDUATION MAJORS"
noi tab job if major == 1 //people generally work in jobs corresponding to their major
noi di ""

noi di "LOW ABILITY"
noi tab job if ability == 0
noi di ""

noi di "HIGH ABILITY"
noi tab job if ability == 1 //higher ability people work in non-teaching jobs more and are more likely to be employed in general.
noi di ""

noi di "MEN"
noi tab job if gender == 0
noi di ""

noi di "WOMEN"
noi tab job if gender == 1 //1 = woman. Women work in teaching more and are employed less.
noi di ""

noi di "NOT BLACK/HISPANIC"
noi tab job if race == 0
noi di ""

noi di "BLACK/HISPANIC"
noi tab job if race == 1 //1 = black or hispanic. Minorities work in teaching more and are employed less.
noi di ""

noi di ""
noi di "################"
noi di ""

//job switching
gen switch = (job!=job_last)
noi di "JOB SWITCHING RATES"
noi tab switch //switching rate of around 20%


log close
}
//end of dofile