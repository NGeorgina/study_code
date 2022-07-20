
set more off

*participants with three vaccinations 
local num "three"
 use vac_safety\AIRD_cohort\covidvacWideFlareWide_airdEdited, clear

tab num_vac
gen num_vac2=num_vac
replace num_vac=3 if num_vac==4  
tab num_vac

keep if num_vac==3

 *start of observation period 
gen cutp1=max(crd, study_start)

*end of observation period 
gen study_end=td(31/12/2021)
format study_end %td

gen cutp2=min(study_end, tod, dod,lcd,date_vac4)  
reshape long date_flare, i(patid) j(epi)
drop if date_flare==.

drop epi

//create exposure category cut points 
*exposure 1

generate cutp3 = date_vac1-7  //pre-exposure period 
generate cutp4 = date_vac1 
generate cutp5 = date_vac1 + 21

*exposure 2
generate cutp6 = date_vac2-7
generate cutp7 = date_vac2 
generate cutp8 = date_vac2 + 21

*exposure 3
generate cutp9 = date_vac3-7
generate cutp10 = date_vac3 
generate cutp11 = date_vac3 + 21


*ensure exposure groups will not overlap
foreach i of numlist 3/10{
local j = `i'+1
replace cutp`i' = cutp`j' if cutp`i' > cutp`j' 

}


foreach var of varlist cutp*{
replace `var' = cutp1 if `var' < cutp1
replace `var' = cutp2 if `var' > cutp2
}

compress

 
sort patid date_flare 
reshape long cutp, i(patid date_flare) j(type )
sort patid date_flare cutp type

*number of adverse events within each interval
by patid: generate int nevents = 1 if date_flare > cutp[_n-1]+0.5 & date_flare <= cutp[_n]+0.5

order cutp, after(type)
format cutp %td


collapse (sum) nevents num_vac num_flare date_vac1 date_vac2 date_vac3 firstvac_brand1 pfizer1 pfizer2 pfizer3 astrazeneca1 astrazeneca2 astrazeneca3 moderna1 moderna2 moderna3, by(patid date_flare cutp type)

*intervals
by patid: generate int interval = cutp[_n] - cutp[_n-1] if type!=1
tab num_flare
order interval, after(cutp)


*exposure groups
 
generate exgr = type-3 if type > 2
count if exgr>=.
local nmiss = r(N)
local nchange = 1
while `nchange'>0{
by patid: replace exgr = exgr[_n+1] if exgr>=.
count if exgr>=.
local nchange = `nmiss'-r(N)
local nmiss = r(N)
}
replace exgr = 0 if exgr==.

save vac_safety\AIRD_cohort\observnPeriods_only`num'Airdvac, replace

*recode the exposure groups:three doses combined 
recode exgr (0=0) (1=1) (2=2) (3=0) (4=1) (5=2) (6=0) (7=1) (8=2)


drop if interval ==0 | interval==.

generate loginterval = log(interval)


//season 
gen winter_start=td(01/12/2020)
format winter_start %td

gen winter_end=td(28/02/2021)
format winter_end %td

gen spring_start=td(01/03/2021)
format spring_start %td

gen spring_end=td(31/05/2021)
format spring_end %td

gen summer_start=td(01/06/2021)
format summer_start %td

gen summer_end=td(31/08/2021)
format summer_end %td

gen autumn_start=td(01/09/2021)
format autumn_start %td

gen autumn_end=td(30/11/2021)
format autumn_end %td

gen winter_start2=td(01/12/2021)
format winter_start2 %td

gen winter_end2=td(31/12/2021)
format winter_end2 %td


bysort patid: generate int season = 1 if winter_start <= cutp[_n] & cutp[_n] <= winter_end
bysort patid: replace season = 1 if winter_start2 < cutp[_n] & cutp[_n] <= winter_end2

order season, after(cutp*)
sort cutp*

bysort patid: replace season = 2 if spring_start <= cutp[_n] & cutp[_n] <= spring_end

bysort patid: replace season = 3 if summer_start <= cutp[_n] & cutp[_n] <= summer_end

bysort patid: replace season = 4 if autumn_start <= cutp[_n] & cutp[_n] <= autumn_end
sort id type
save vac_safety\AIRD_cohort\vacSCCS_airds`num'21new, replace  


/*append SCCS data for single, two and three vaccinations together to create a 
SCCS analysis dataset for all C-19 vaccine doses per individual   
*/
use vac_safety\AIRD_cohort\vacSCCS_airdssingle21, clear  
     
gen only_vac=1
 append using vac_safety\AIRD_cohort\onlytwoVacSCCS_airds21 
 replace only_vac =2 if only_vac==.
append using vac_safety\AIRD_cohort\vacSCCS_airdsthree21new 
replace only_vac=3 if only_vac==.
tab only_vac, miss 

save vac_safety\AIRD_cohort\RepeatVacSCCS_airdsNew21new, replace  
 
//create numeric patid 

keep patid
sort patid
by patid: gen epi=_n
keep if epi==1
count
drop epi
gen id=_n
save vac_safety\AIRD_cohort\RepeatVacSCCS_airdspatidsNew21new, replace

merge 1:m patid using vac_safety\AIRD_cohort\RepeatVacSCCS_airdsNew21new
drop _m

save vac_safety\AIRD_cohort\RepeatVacSCCS_airdsEdited21new, replace 

*fit model
xi: xtpoisson nevents i.exgr, fe i(id) offset(loginterval) irr
xi: xtpoisson nevents i.exgr i.season, fe i(id) offset(loginterval) irr

*stratified analysis according to disease type

use vac_safety\AIRD_cohort\airdstudypon_covidvacFlareDemographics, clear 

tab aird
label define gg 1"RA" 2"SpA" 3"CTD/vasculitis" 4"GCA/PMR"
label values aird gg
tab aird
keep patid aird
merge 1:m patid using vac_safety\AIRD_cohort\RepeatVacSCCS_airdsEdited21new
drop _m

forvalues i=1/4{
xi: xtpoisson nevents i.exgr if aird==`i' , fe i(id) offset(loginterval) irr
xi: xtpoisson nevents i.exgr i.season if aird==`i' , fe i(id) offset(loginterval) irr
}


//Stratified analysis according to 1st, 2nd or 3rd vaccine dose

local num "three"
use vac_safety\AIRD_cohort\observnPeriods_only`num'Airdvac, clear 

*recode exposure groups  
recode exgr (0=0) (1=1) (2=2) (3=0) (4=3) (5=4) (6=0) (7=5) (8=6)


drop if interval ==0 | interval==.

generate loginterval = log(interval)


//season 
gen winter_start=td(01/12/2020)
format winter_start %td

gen winter_end=td(28/02/2021)
format winter_end %td

gen spring_start=td(01/03/2021)
format spring_start %td

gen spring_end=td(31/05/2021)
format spring_end %td

gen summer_start=td(01/06/2021)
format summer_start %td

gen summer_end=td(31/08/2021)
format summer_end %td

gen autumn_start=td(01/09/2021)
format autumn_start %td

gen autumn_end=td(30/11/2021)
format autumn_end %td

gen winter_start2=td(01/12/2021)
format winter_start2 %td

gen winter_end2=td(31/12/2021)
format winter_end2 %td


bysort patid: generate int season = 1 if winter_start <= cutp[_n] & cutp[_n] <= winter_end
bysort patid: replace season = 1 if winter_start2 <= cutp[_n] & cutp[_n] <= winter_end2

order season, after(cutp*)
sort cutp*

bysort patid: replace season = 2 if spring_start <= cutp[_n] & cutp[_n] <= spring_end

bysort patid: replace season = 3 if summer_start <= cutp[_n] & cutp[_n] <= summer_end

bysort patid: replace season = 4 if autumn_start <= cutp[_n] & cutp[_n] <= autumn_end
sort id type
save vac_safety\AIRD_cohort\vacSCCS_airds`num'21doseEffect, replace  


//append SCCS data for single, two and three vaccinations 
use vac_safety\AIRD_cohort\vacSCCS_airdssingle21, clear  
     
gen only_vac=1
 append using vac_safety\AIRD_cohort\onlytwoVacSCCS_airds21doseEffect 
 replace only_vac =2 if only_vac==.
append using vac_safety\AIRD_cohort\vacSCCS_airdsthree21doseEffect 
replace only_vac=3 if only_vac==.
tab only_vac, miss 

save vac_safety\AIRD_cohort\RepeatVacSCCS_airdsNew21doseEffect, replace  
 
//create numeric patid 
drop id
keep patid
sort patid
by patid: gen epi=_n
keep if epi==1
count
drop epi
gen id=_n
save vac_safety\AIRD_cohort\RepeatVacSCCS_airdspatidsNew21doseEffect, replace

merge 1:m patid using vac_safety\AIRD_cohort\RepeatVacSCCS_airdsNew21doseEffect
drop _m

save vac_safety\AIRD_cohort\RepeatVacSCCS_airdsEdited21doseEffect, replace 

order patid id only_vac exgr

*fit model
xi: xtpoisson nevents i.exgr, fe i(id) offset(loginterval) irr
xi: xtpoisson nevents i.exgr i.season, fe i(id) offset(loginterval) irr



