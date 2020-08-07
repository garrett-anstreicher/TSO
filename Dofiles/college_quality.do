clear all
do "$data/CQ/STATA_RV_7192020-802"
save "$temp/cq_1", replace

do "$data/CQ/STATA_RV_7192020-567"
merge 1:1 unitid using "$temp/cq_1", keep(match) nogen


//keep schools that offer a bachelors
keep if level5 == 1
keep if carnegie<50 //nix specialty schools
replace satmt75 = satmt75/10
replace saltotl = saltotl/1000


replace dvic01 = 100 - dvic01 //flip to rejection rate
replace stufacr = 1/stufacr*1000


//pca: ratio, satv75, satm75, actcm75 (sat verb, math, and act composite)
pca stufacr saltotl dvic01 satmt75
//predict pca1, score

//keep if pca!=.
//xtile qual_pct = pca[fw=totenrl], nq(100)
//keep unitid instnm qual_pct
//save "$temp/qual_pct", replace

/*
first eigenvector: 
-stufacr: 0.3633
-saltotl: 0.5533
-dvic01: 0.4760
-satmt75: 0.5790
*/

local vars `"stufacr salprof dvic01 satmt75"'
foreach var in `vars'{
    gen `var'_miss = (`var' == .)
}
egen miss_total = rowtotal(*_miss)
drop if miss_total>2 //at least two quality proxies must be present

//fill in
gen qual = stufacr * 0.3633 if stufacr!=.
replace qual = qual + saltotl*0.5533 if saltotl!=.
replace qual = qual + dvic01*0.4760 if dvic01!=.
replace qual = qual + satmt75*0.5790 if satmt75!=.

//other adjustments
replace qual = qual/4 if miss_total == 0
replace qual = qual/3 if miss_total == 1
replace qual = qual/2 if miss_total == 2

//percentiles
xtile qual_pct = qual [fw=totenrl], nq(100) //pctiles
keep unitid instnm qual_pct
save "$temp/college_qualities", replace


//