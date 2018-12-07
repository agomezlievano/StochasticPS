// EComm_0002_fromR2Stata.do

*use "C:\Users\agomez\Dropbox\Harvard\LittleProjects\StochasticPS\outputdata\EComm_0002_fromR2Stata.dta", clear
import delimited "C:\Users\agomez\Dropbox\Harvard\LittleProjects\StochasticPS\outputdata\EComm_0002_fromR2Stata.csv", clear case(preserve)

local period 10

cap drop rY`period'yr
cap drop rX`period'yr
bysort exporter (year): gen rY`period'yr = (1+(GDPpcconst2010USD[_n + `period'] - GDPpcconst2010USD[_n])/GDPpcconst2010USD[_n])^(1.0/`period') - 1 if _n <= _N - `period'
bysort exporter (year): gen rX`period'yr = (1+(TotalExports[_n + `period'] - TotalExports[_n])/TotalExports[_n])^(1.0/`period')-1 if _n <= _N - `period'

* Remove outliers
qui sum rY`period'yr, d
drop if rY`period'yr < `r(p1)' | rY`period'yr > `r(p99)'

gen logpop = log(population)
gen loggdppc = log(GDPpcconst2010USD)

export delimited "C:\Users\agomez\Dropbox\Harvard\LittleProjects\StochasticPS\outputdata\EComm_0002_2_fromStata2R.csv", replace quote


*cor eci Distance2Origin_cp Distance2Origin_ccp

* levels
reg loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp

reg loggdppc logpop c.eci##c.Distance2Origin_cp c.eci##c.Distance2Origin_ccp SimDirection2MostComplex_cp SimDirection2MostComplex_ccp Distance2MostComplex_cp Distance2MostComplex_ccp
reg loggdppc logpop eci Distance2Origin_cp Distance2Origin_ccp SimDirection2MostComplex_cp SimDirection2MostComplex_ccp Distance2MostComplex_cp Distance2MostComplex_ccp

* growth
reg rY10yr loggdppc logpop eci leftDistance2Origin_cp rightDistance2Origin_cp leftSimDirection2MostComplex_cp rightSimDirection2MostComplex_cc leftDistance2MostComplex_cp




reg rY10yr loggdppc logpop eci Distance2Origin_ccp   SimDirection2MostComplex_cp        ///
              c.Distance2Origin_ccp#c.loggdppc   ///
             c.logpop#c.loggdppc   c.SimDirection2MostComplex_cp#c.logpop   c.SimDirection2MostComplex_cp#c.loggdppc   ///
             c.Distance2Origin_ccp#c.logpop   c.eci#c.SimDirection2MostComplex_cp   ///
             c.eci#c.Distance2MostComplex_cp   c.Distance2Origin_ccp#c.Distance2MostComplex_cp


			 



* ==============================================================================
* Ricardo pointed out that regressions across many years give t-statistics
* that are inflated because observations are not independent. Hence, do the 
* regressions on a single year (maybe compare the stability of the coefficients
* by repeating the regression separately across all years).
* ------------------------------------------------------------------------------
encode region, generate(region_encode)

* levels
eststo clear
eststo: reg loggdppc logpop                                                                                                           if year==2005
eststo: reg loggdppc        eci                                                                                                       if year==2005
eststo: reg loggdppc logpop eci                                                                                                       if year==2005
eststo: reg loggdppc logpop     rightDistance2Origin_ccp                                                                              if year==2005
eststo: reg loggdppc logpop                              rightSimDirection2MostComplex_cp                                             if year==2005
eststo: reg loggdppc logpop                                                               leftDistance2MostComplex_cp                 if year==2005
eststo: reg loggdppc logpop eci rightDistance2Origin_ccp                                                                              if year==2005 
eststo: reg loggdppc logpop eci                          rightSimDirection2MostComplex_cp                                             if year==2005 
eststo: reg loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp                 if year==2005
eststo: reg loggdppc logpop     rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp                 if year==2005
eststo: reg loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode if year==2005
eststo: reg loggdppc logpop     rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode if year==2005
eststo: reg loggdppc logpop                              rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode if year==2005
esttab , ar2 nogaps compress

* levels
eststo clear
eststo: reg loggdppc logpop                                                                                                           , vce(cluster exporter)
eststo: reg loggdppc        eci                                                                                                       , vce(cluster exporter)
eststo: reg loggdppc logpop eci                                                                                                       , vce(cluster exporter)
eststo: reg loggdppc logpop     rightDistance2Origin_ccp                                                                              , vce(cluster exporter)
eststo: reg loggdppc logpop                              rightSimDirection2MostComplex_cp                                             , vce(cluster exporter)
eststo: reg loggdppc logpop                                                               leftDistance2MostComplex_cp                 , vce(cluster exporter)
eststo: reg loggdppc logpop eci rightDistance2Origin_ccp                                                                              , vce(cluster exporter) 
eststo: reg loggdppc logpop eci                          rightSimDirection2MostComplex_cp                                             , vce(cluster exporter) 
eststo: reg loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp                 , vce(cluster exporter)
eststo: reg loggdppc logpop     rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp                 , vce(cluster exporter)
eststo: reg loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode , vce(cluster exporter)
eststo: reg loggdppc logpop     rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode , vce(cluster exporter)
eststo: reg loggdppc logpop                              rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode , vce(cluster exporter)
esttab , ar2 nogaps compress

* growth
eststo clear
eststo: reg rY10yr loggdppc logpop                                                                                                           if year==2005
eststo: reg rY10yr loggdppc        eci                                                                                                       if year==2005
eststo: reg rY10yr loggdppc logpop eci                                                                                                       if year==2005
eststo: reg rY10yr loggdppc logpop     rightDistance2Origin_ccp                                                                              if year==2005
eststo: reg rY10yr loggdppc logpop                              rightSimDirection2MostComplex_cp                                             if year==2005
eststo: reg rY10yr loggdppc logpop                                                               leftDistance2MostComplex_cp                 if year==2005
eststo: reg rY10yr loggdppc logpop eci rightDistance2Origin_ccp                                                                              if year==2005 
eststo: reg rY10yr loggdppc logpop eci                          rightSimDirection2MostComplex_cp                                             if year==2005 
eststo: reg rY10yr loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp                 if year==2005
eststo: reg rY10yr loggdppc logpop     rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp                 if year==2005
eststo: reg rY10yr loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode if year==2005
eststo: reg rY10yr loggdppc logpop     rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode if year==2005
eststo: reg rY10yr loggdppc logpop                              rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode if year==2005
esttab , ar2 nogaps compress


* growth
eststo clear
eststo: reg rY10yr loggdppc logpop                                                                                                           , vce(cluster exporter)
eststo: reg rY10yr loggdppc        eci                                                                                                       , vce(cluster exporter)
eststo: reg rY10yr loggdppc logpop eci                                                                                                       , vce(cluster exporter)
eststo: reg rY10yr loggdppc logpop     rightDistance2Origin_ccp                                                                              , vce(cluster exporter)
eststo: reg rY10yr loggdppc logpop                              rightSimDirection2MostComplex_cp                                             , vce(cluster exporter)
eststo: reg rY10yr loggdppc logpop                                                               leftDistance2MostComplex_cp                 , vce(cluster exporter)
eststo: reg rY10yr loggdppc logpop eci rightDistance2Origin_ccp                                                                              , vce(cluster exporter) 
eststo: reg rY10yr loggdppc logpop eci                          rightSimDirection2MostComplex_cp                                             , vce(cluster exporter) 
eststo: reg rY10yr loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp                 , vce(cluster exporter)
eststo: reg rY10yr loggdppc logpop     rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp                 , vce(cluster exporter)
eststo: reg rY10yr loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode , vce(cluster exporter)
eststo: reg rY10yr loggdppc logpop     rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode , vce(cluster exporter)
eststo: reg rY10yr loggdppc logpop                              rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode , vce(cluster exporter)
esttab , ar2 nogaps compress

* ------------------------------------------------------------------------------




			 
eststo clear
eststo: reg loggdppc logpop  
eststo: reg loggdppc        eci 
eststo: reg loggdppc logpop eci 
eststo: reg loggdppc logpop     rightDistance2Origin_ccp 
eststo: reg loggdppc logpop                              rightSimDirection2MostComplex_cp 
eststo: reg loggdppc logpop                                                               leftDistance2MostComplex_cp
eststo: reg loggdppc logpop eci rightDistance2Origin_ccp 
eststo: reg loggdppc logpop eci                          rightSimDirection2MostComplex_cp 
eststo: reg loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp
eststo: reg loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode
esttab , ar2

eststo clear
eststo: reg rY10yr loggdppc                                                                                                   
eststo: reg rY10yr loggdppc logpop                                                                                           
eststo: reg rY10yr loggdppc logpop eci                                                                                      
eststo: reg rY10yr loggdppc logpop     rightDistance2Origin_ccp                                                             
eststo: reg rY10yr loggdppc logpop                              rightSimDirection2MostComplex_cp                             
eststo: reg rY10yr loggdppc logpop                                                               leftDistance2MostComplex_cp 
eststo: reg rY10yr loggdppc logpop eci rightDistance2Origin_ccp                                                              
eststo: reg rY10yr loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp 
eststo: reg rY10yr loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode
esttab , ar2


eststo clear
eststo: reg loggdppc logpop                              leftSimDirection2MostComplex_cp 
eststo: reg loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp
eststo: reg loggdppc logpop eci rightSimDirection2MostComplex_cp rightSimDirection2MostComplex_cc rightDistance2MostComplex_cp rightDistance2MostComplex_ccp
eststo: reg loggdppc logpop eci rightDistance2Origin_cp rightDistance2Origin_cc rightSimDirection2MostComplex_cp rightSimDirection2MostComplex_cc rightDistance2MostComplex_cp rightDistance2MostComplex_cc
esttab , r2 ar2


reg loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp


reg loggdppc logpop eci rightDistance2Origin_ccp rightSimDirection2MostComplex_cp leftDistance2MostComplex_cp i.region_encode
