# routing-and-scheduling-for-offshore-wind-farms
Using Dantzig-Wolfe Decomposition to solve routing and scheduling problem for offshore wind farms

Optimization of Maintenance Routing and Scheduling Problem for offshore wind farms
---
A solution to the routing&scheduling problem using Dantzig-Wolfe Decomposition. Implementation with python and SCIP.
The algorithm based on the paper: Optimisation of maintenance routing and scheduling for offshore wind farms[^1]. 
[^1]: Irawan, Chandra Ade, Djamila Ouelhadj, Dylan Jones, Magnus Stålhane and Iver Bakken Sperstad. “Optimisation of maintenance routing and scheduling for offshore wind farms.” Eur. J. Oper. Res. 256 (2017): 76-89.

Dependencies
---
- Python 3.7
- gurobi 9.0 (license needed)
- python package: numpy, itertools


GB: an optimization model for routing & scheduling, 3-7 days, multi bases & multi wind farms
---
####  Model:
v = vessels, d = days, f = wind farms, p = type of technician, J = set of turbines
$dei=e$
decision variables:  
- $ y_{v,d,f,i,j} $ binary variable
- $t_{v,d,f,i}$ continous variable
- $q_{v,d,f,p,i}$ integer variable
- $x_{f,i}$: integer variable, the number of delayed days for maintenance task on turbine i in wind farm f

objective function:
$\min\ Cost=c^{q}+c^{t}+c^{l}$
$c^{q}=\sum\limits_{v\in V} \sum\limits_{d\in D} \sum\limits_{p\in P}(q_{v,d,p,0} \cdot \hat{c}_{p})$ Personnel Cost
$c^{t}=\sum\limits_{v\in V} \sum\limits_{d\in D} \sum\limits_{f\in F} \sum\limits_{i\in J^{*}} \sum\limits_{j\in J^{*}} (c_{i,j}\cdot y_{v,d,f,i,j})$ Travel Cost
$c^{l} = \sum\limits_{f\in F}\sum\limits_{j\in J^{-}_{f}}x_{f,j}\cdot \tilde{c}_{f,j}$ Penalty Cost

constrains:
$\sum\limits_{j \in J^{-}_{f}}y_{v,d,f,0,j} \leq 1\ , \forall f\in F,\ v \in V,\ d \in D$, leave base at most once
$\sum\limits_{i \in J^{+}_{f}}y_{v,d,f,i,(2n+1)} \leq 1\ , \forall f\in F,\ v \in V,\ d \in D$, return base at most once
$\sum\limits_{v \in V} \sum\limits_{d \in D} \sum\limits_{i \in J^{*}_{f}}y_{v,d,f,i,j} = 1\ , \forall j \in J^{-}_{f},\ f \in F$, turbine must be visited once
$\sum\limits_{j \in J^{*}_{f}}y_{v,d,f,i,j} = \sum\limits_{j \in J^{*}_{f}}y_{v,d,f,j,i}\ , \forall i\in J^{*}_{f},\ f \in F,\ d\in D,\ v\in V$, flow conservation
$\sum\limits_{j \in J^{*}_{f}}y_{v,d,f,i,j} = \sum\limits_{j \in J^{*}_{f}}y_{v,d,f,(n+i),j}\ , \forall i \in J^{-}_{f},\ f \in F,\ d\in D,\ v\in V$, ensure drop off and pick up are done
$\sum\limits_{v\in V}\sum\limits_{d\in D}y_{v,d,f,i,(n+i)} = 1,\ \forall i \in J^{V}_{f},\ f \in F$
$t_{v,d,f,(2n+1)} \leq \psi_{v,d,f}\ , \forall f \in F,\ d\in D,\ v\in V$
$t_{v,d,f,0}=0\ , \forall f \in F,\ d\in D,\ v\in V$
$t_{v,d,f,(n+i)}-t_{v,d,f,i}-\sum\limits_{j \in J^{-}_{f}}y_{v,d,f,0,j}\cdot (\hat{\tau}_{i}+\tilde{\tau}) \geq 0\ , \forall i\in J^{-}_{f},\ d\in D,\ v\in V$, ensure travelling time compatibility
$y_{v,d,f,i,j}\cdot (t_{v,d,f,i}+\tau_{v,i,j}+\tilde{\tau}-t_{v,d,f,j})\leq 0,\ \forall i,j\in J_{f},\ j\neq i+n,\ f \in F,\ d \in D,\ v\in V$
$0\leq \sum\limits_{p\in P} q_{v,d,f,p,i} \leq \tilde{\rho}_{v}\ ,\forall i\in J^{*}_{f},\ f \in F,\ d\in D,\ v\in V$, personnel capacity
$\sum\limits_{i \in J^{-}_{f}} \sum\limits_{j \in J^{*}_{f}}(y_{v,t,f,i,j}\cdot w_{i}) \leq \hat{w}_{v}\ ,\forall f\in F,\ t\in T,\ v\in V$, weight capacity
$q_{v,d,f,p,i}\leq q_{v,d,f,p,0}\ ,\forall p\in P,\ i\in J^{*}_{f},\ f \in F$
$\sum\limits_{v\in V} \sum\limits_{f\in F} q_{v,d,f,p,0}\leq \hat{\rho}_{p}\ ,\forall p\in P,\ i\in J^{*}_{f},\ f \in F$, personnel availability at base
$y_{v,d,f,i,j}\cdot (q_{v,d,f,p,i}-\rho_{p,j}-q_{v,d,f,p,j})=0\ ,\forall p\in P,\ j\in J^{-}_{f},\ i\in J^{*}_{f},\ f \in F$
$y_{v,d,f,i,j}\cdot (q_{v,d,f,p,i}+\rho_{p,j}-q_{v,d,f,p,j})=0\ ,\forall p\in P,\ j\in J^{+}_{f},\ i\in J^{*}_{f},\ f \in F$
$\sum\limits_{i\in J^{*}_{f}} \sum\limits_{v\in V} \sum\limits_{d\in D}(d\cdot y_{v,d,f,i,j} - x_{f,j}) \leq \delta_{f,j}\ ,\forall j\in J^{-}_{f},\ f \in F$, soft constraint for lateness of maintenance
