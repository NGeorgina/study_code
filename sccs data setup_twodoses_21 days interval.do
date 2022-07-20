
set more off 

*participants with two vaccinations 
local num "two"
 use vac_safety\AIRD_cohort\covidvacWideFlareWide_airdEdited, clear

tab num_vac

keep if num_vac==2

 *start of observation period 
gen cutp1=max(crd, study_start)

*end of observation period 
gen study_end=td(31/12/2021)
format study_end %td

gen cutp2=min(study_end, tod, dod,lcd) 

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


*ensure exposure groups will not overlap
foreach i of numlist 3/7{
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
collapse (sum) nevents date_firstvac1 firstvac_brand1 date_vac1 pfizer1 astrazeneca1 moderna1 date_vac2 pfizer2 astrazeneca2 moderna2 num_vac num_flare, by(patid date_flare cutp type)

*intervals
by patid: generate int interval = cutp[_n] - cutp[_n-1] if type!=1


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
*recode the exposure groups: two doses combined 
recode exgr (0=0) (1=1) (2=2) (3=0) (4=1) (5=2) 

drop if interval ==0 | interval==.

generate loginterval = log(interval)
save vac_safety\AIRD_cohort\observnPeriods_only`num'VacsAirds21, replace 
 
//string patid is not compartible with xtpossion syntax below 

keep patid
by patid: gen epi=_n
keep if epi==1
gen id=_n
drop epi
save vac_safety\AIRD_cohort\only`num'VacsAirds_airdsSCCSpatids21, replace 

merge 1:m patid using vac_safety\AIRD_cohort\observnPeriods_only`num'VacsAirds21
drop _m 
sort id type

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

save vac_safety\AIRD_cohort\only`num'VacSCCS_airds21, replace 


//For stratified analysis according to 1st and 2nd vaccine dose

local num "two"
use vac_safety\AIRD_cohort\observnPeriods_only`num'Airdvac, clear 

*recode exposure groups  
recode exgr (0=0) (1=1) (2=2) (3=0) (4=3) (5=4) 
save vac_safety\AIRD_cohort\onlytwoVacSCCS_airds21doseEffect, replace 





