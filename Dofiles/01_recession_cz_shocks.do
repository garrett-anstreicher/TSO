//create plot of CZ recession shocks

//forval y = 07/09


//read in county unemployment data
forval y = 7/9{
	import excel using "$data/BLS/laucnty0`y'.xlsx", clear 
	gen fips = B + C
	keep fips J
	ren J unemp_rate
	drop if _n<7
	destring *, replace
	ren unemp_rate unemp_rate_`y'
	duplicates tag fips, gen(dup)
	drop if dup
	drop dup
	save "$temp/county_unemp_0`y'", replace
}




//read in county-cz crosswalk
import excel "$data/crosswalks/cz00_eqv_v1.xls", clear 
ren A fips
ren C cz 
ren G pop
keep fips cz pop
drop if _n == 1
destring *, replace

//merge
forval y = 7/9{
	merge 1:1 fips using "$temp/county_unemp_0`y'", keep(match) nogen
}


tostring fips, replace
replace fips = "0" + fips if length(fips) == 4
gen statefip = substr(fips, 1, 2)
destring fips statefip, replace
collapse (mean) unemp* [fw = pop], by(statefip) //collapse by cz, weighting by population
gen delta = unemp_rate_9 - unemp_rate_7
keep statefip delta
save "$temp/gr_state_shocks", replace
