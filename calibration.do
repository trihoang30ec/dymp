* Full name: Minh Tri Hoang (Tri Hoang)
* MSc student in Economics, Department of Economics, University of Bonn

* COURSE: MACROECONOMICS (PhD-level)


*_________________________PROBLEM SET (Deadline: 16/11/2020)____________________*
* Question 1.11 and 1.12

/*
Utility function: U(c) = log(c) = log(f(k) - k')
A = 20
alpha = 0.3
beta = 0.6
K = {2, 4, 6, 8, 10}
v(k) = max{u(f(k)-k') + beta*v'(k')} on k' (0<k'<f(k)). We have to prove that v(k) = v'(k) for any k
stopping criterion: |v_{n+1} - v_{n}| < 10^-6
*/

*___________________________________RUN ALL_____________________________________*
clear all

local obs 100
* Create 100 observations on u from 2 to 10
range u 2 10 100
sort u

* STEP 1: Suppose v'(k) = 0 for any k in K, then k' = 2. We have: v1 = log(20*k^0.3 - 2) + 0.6*0

local A = 20
local alpha = 0.3
local beta = 0.6
local k_set 

forval z = 1/`obs' {
local e = u[`z']
local k_set `k_set' `e'
}

di `k_set'


* Create a sequence of k'
local k_prime
foreach k in `k_set' {
local k_prime `k_prime' `k'
}


* STEP 2: Calculate v1(k) for any k
local v1
local k_index = 0
foreach i in `k_set' {
local ++k_index
local v1_k`k_index' = log(`A'*(`i'^`alpha') - 2) + `beta'*0 // value function
di "v1_k`k_index': " `v1_k`k_index''
local v1 `v1' `v1_k`k_index''
}


* STEP 3: Calculate |v_1 - v_0| (Remark: |.| is the Euclidean norm)

* Create a column vector v1.
matrix input v1 = (`v1')
matrix v1 = v1'
* Calculate norm of v = v_1 - v_0 (remark: v_0 = 0 for any k in K)
svmat v1 // convert a vector to a variable
replace v1 = v1^2
egen norm_v1 = total(v1)
replace norm_v1 = sqrt(norm_v1)
local norm_v = norm_v1[1] // Remark: norm_v = |v1 - v0| = |v1 - 0| because v0(k) = 0 for any k. 
clear // Delete old dataset (IMPORTANT!!!)
* Calculate v2(k), v3(k),..., vn(k)
* Define v_{n+1} = v1(k), and v_{n} = v0

local n = 1
while `norm_v' >= 10^-6 /* stopping criterion */ {
local ++n
local seq_v
local seq_k

local k_index = 0
foreach i in `k_set' {
local ++k_index
local vec_v  // for each k, create a vector v(k) corresponding to k'

local k_index_p = 0
foreach j in `k_set' {
local ++k_index_p
local n_1 = `n'-1
di "v`n_1'_k`k_index_p'"
local v0 = `v`n_1'_k`k_index_p'' // v'(k'). For instance, v0(k') is known, find v1(k)
local v1 = log(`A'*(`i'^`alpha') - `j') + `beta'*`v0' // i = k, j = k', v1 = v(k), v0 = v'(k')
* We will calculate: v1, v2, v3,..., vn(k). For each k, impute k' and find maximum value of v(k) (k' is given)
local v`n'_k`k_index'_`k_index_p' = `v1' // vn_k(k'). 
local vec_v `vec_v' `v`n'_k`k_index'_`k_index_p''
}

di "`vec_v'"
local x : subinstr local vec_v " " ",", all
di "`x'"
local v`n'_k`k_index' = max(`x') // maximum value of v(k) on k'
di "v`n'_k`k_index': " `v`n'_k`k_index''

* Create sequence of value function v(k) for each k
local seq_v `seq_v' `v`n'_k`k_index'' 
di "`seq_v'"

* Trace out k': 
local index = 0
foreach v_k in `vec_v' {
local ++index 
if `v`n'_k`k_index'' == `v_k' {
local k_star_`n'_k`k_index' : word `index' of `k_prime'
}
}
}
local seq_v_k`k_index' `seq_v' // Create sequence of value function v(k) for each k

* Create a column vector
matrix input v`n' = (`seq_v_k`k_index'')
matrix v`n' = v`n''

* Calculate norm of v = v_{n} - v_{n-1}
matrix v = v`n' - v`n_1'
svmat v // convert a vector to a variable
replace v = v^2
egen norm_v = total(v)
replace norm_v = sqrt(norm_v)
local norm_v = norm_v[1]
di "|v_{n} - v_{n-1} = " `norm_v'
clear // Delete old dataset (IMPORTANT!!!)
}


* STEP 4: Result

local z = `n' - 1 
set obs `z'

local k_index = 0
foreach i in `k_set' {
local ++k_index
gen v_k`k_index' = .
gen k_k`k_index' = .
gen k_`k_index' = .
forval t = 2/`n' {
local j = `t' - 1
replace v_k`k_index' = `v`t'_k`k_index'' in `j'
replace k_k`k_index' = `k_star_`t'_k`k_index'' in `j'
replace k_`k_index' = `i'
}
}

order v*

gen n = _n
local t_converged = _N
reshape long v_k k_k k_, i(n) j(i)

la var v_k "V_n(k) (Numerical Solution)"
la var k_ "Initial Capital (K)"
la var n "index n in V_n(k)"
la var k_k "Next period's Capital (K') (Numerical Solution)"
gen k_p_ana = .
gen k_p_gap = .
gen v_ana = .
gen v_gap = .
la var k_p_ana "Next period's Capital (K') (Analytic Solution)"
la var v_ana "V_n(k) (Analytic Solution)"
la var v_gap "Compare v(k)"
la var k_p_gap "Compare g(k)"


* STEP 5: Compare analytic solution and numerical solution

local alpha_0 = 1/(1-`beta')*(log(1-`alpha'*`beta') + (`alpha'*`beta')/(1-`alpha'*`beta')*log(`alpha'*`beta') + log(`A')/(1-`alpha'*`beta'))
local alpha_1 = `alpha'/(1-`alpha'*`beta')

replace k_p_ana = `beta'*`alpha_1'*`A'*(k_^`alpha')/(1+`beta'*`alpha_1')
replace k_p_gap = k_p_ana - k_k if n == `t_converged'
replace v_ana = `alpha_0' + `alpha_1'*log(k_) if n == `t_converged'
replace v_gap = v_ana - v_k
gsort -n i

local a_0 = round(`alpha_0', 0.01)
local a_1 = round(`alpha_1', 0.01)

* Parameters
di "A = " `A'
di "alpha = " `alpha'
di "beta = " `beta'
di "alpha_0 = " `alpha_0'
di "alpha_1 = " `alpha_1' 
di "t_converged = " `t_converged'

* Convergence of Value Function
graph twoway (scatter v_k k_, mcolor(blue) msymbol(x) connect(none) yvarlab("v_n(k) for n > 1") connect(ascending) xtitle("Captial Stock k Today") ytitle("Value Function v(k)")) ///
(scatter v_k k_ if n == `t_converged', yvarlab("Converged v(k)") msymbol(X) connect(none) mcolor(black) connect(ascending)) ///
(function y = `alpha_0' + `alpha_1'*log(x), range(2 10) color(red) plotregion(style(none)) yvarlab("v(k) = `a_0' + `a_1'*log(k) (analytic solution)"))

graph save value_function, replace


* Convergence of Policy Function
graph twoway (scatter k_k k_, mcolor(blue) msymbol(x) connect(none) yvarlab("g(k)") xtitle("Captial Stock k Today") ytitle("Policy Function g(k)")) ///
(scatter k_k k_ if n == `t_converged', yvarlab("Converged g(k)") msymbol(X) connect(none) mcolor(black) connect(ascending)) ///
(function y = `beta'*`alpha_1'*`A'*(x^`alpha')/(1+`beta'*`alpha_1'), range(2 10) color(red) plotregion(style(none)) yvarlab("k' = g(k) (analytic solution)"))

graph save policy_function, replace

* Recall graphs
graph use policy_function // Policy Function
graph use value_function // Value Function 











