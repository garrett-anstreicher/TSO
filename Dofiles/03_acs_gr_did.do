use "$temp/acs_degfield_cleaned", clear
ren statefip state_live
ren bpl statefip
merge m:1 statefip using "$temp/gr_state_shocks", keep(match) nogen //merge GR shock to state of birth

//basic spec: educ_degree = demo + stFE + cohortFE + yearFE + after + after*shock
gen after = (birthyr>=1988) //age 20 in 2008; perhaps time to change degree. Can play with this.
gen after_shock = after*delta

//do the regression
local controls `"after i.birthyr i.year female black_hispan after"'
areg educ_degree after_shock `controls' [w=perwt], a(statefip) cl(statefip)






//exit