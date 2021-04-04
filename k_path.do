* Full name: Minh Tri Hoang (Tri Hoang)
* MSc student in Economics, Department of Economics, University of Bonn.

* COURSE: MACROECONOMICS (PhD-level)


*_________________________PROBLEM SET (Deadline: 16/11/2020)____________________*
* Question 1.13

clear all
local A = 20
local alpha = 0.3
local beta = 0.6

local alpha_0 = 1/(1-`beta')*(log(1-`alpha'*`beta') + (`alpha'*`beta')/(1-`alpha'*`beta')*log(`alpha'*`beta') + log(`A')/(1-`alpha'*`beta'))
local alpha_1 = `alpha'/(1-`alpha'*`beta')

local k_set

local k_0 = 2

forval i = 1/20 {
local k_optimal = `beta'*`alpha_1'*`A'*(`k_0'^`alpha')/(1+`beta'*`alpha_1')
local k_set `k_set' `k_optimal'
local k_0 = `k_optimal'
}

matrix input k = (`k_set')
svmat k
