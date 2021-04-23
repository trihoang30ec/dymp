# Dynamic Programming

Use the following parametrization of the model: ```A = 20```, ```α = 0.3```, and ```β = 0.6```.
Set up a grid for possible capital values in the economy ```K = {2, 4, 6, 8, 10}```. 
Write a computer code to solve the problem numerically. 

Use: ```max k∈K |vn+1(k) − vn(k)| < 10−6``` as stopping criterion.

1. Set up grids for possible capital values in the economy ``K′`` with 100 (1,000) points between 2 and 10. Write a computer code to solve the problem numerically. (Solution: see ```numso.do```)
2. Compare the value and policy functions from the analytic solution and the solutions with discretized choices. (Solution: see ```numso.do```)
3. Determine the path for ``{kt}`` under the optimal solution starting from ``k0 = 2``. (Solution: see ```kpath.do```)
