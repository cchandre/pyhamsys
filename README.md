# pyHamSys
pyHamSys is a Python package for scientific computing involving Hamiltonian systems

![PyPI](https://img.shields.io/pypi/v/pyhamsys)
![License](https://img.shields.io/badge/license-BSD-lightgray)

Installation within a Python virtual environment: 
```
python3 -m pip install pyhamsys
```
For more information on creating a Python virtual environment, click [here](https://realpython.com/python-virtual-environments-a-primer/).

## Symplectic Integrators
pyHamSys includes a class SymplecticIntegrator containing the following symplectic splitting integrators:

- `Verlet` (order 2, all purpose), also referred to as Strang or Störmer-Verlet splitting 
- From [Forest, Ruth, Physica D 43, 105 (1990)](https://doi.org/10.1016/0167-2789(90)90019-L): 
    - `FR` (order 4, all purpose)
- From [Yoshida, Phys. Lett. A 150, 262 (1990)](https://doi.org/10.1016/0375-9601(90)90092-3):
    - `Yo#`: # should be replaced by an even integer, e.g., `Yo6` for 6th order symplectic integrator (all purpose)
    - `Yos6`: (order 6, all purpose) optimized symplectic integrator (solution A from Table 1)
- From [McLachlan, SIAM J. Sci. Comp. 16, 151 (1995)](https://doi.org/10.1137/0916010):
    - `M2` (order 2, all purpose)
    - `M4` (order 4, all purpose)
- From [Omelyan, Mryglod, Folk, Comput. Phys. Commun. 146, 188 (2002)](https://doi.org/10.1016/S0010-4655(02)00451-4): 
    - `EFRL` (order 4) optimized for *H* = *A* + *B*
    - `PEFRL` and `VEFRL` (order 4) optimized for *H* = *A*(*p*) + *B*(*q*). For `PEFRL`, *chi* should be exp(*h* X<sub>A</sub>)exp(*h* X<sub>B</sub>). For `VEFRL`, *chi* should be exp(*h* X<sub>B</sub>)exp(*h* X<sub>A</sub>).
- From [Blanes, Moan, J. Comput. Appl. Math. 142, 313 (2002)](https://doi.org/10.1016/S0377-0427(01)00492-7):
    - `BM4` (order 4, all purpose) refers to S<sub>6</sub> 
    - `BM6` (order 6, all purpose) refers to S<sub>10</sub>
    - `RKN4b` (order 4) refers to SRKN<sub>6</sub><sup>*b*</sup> optimized for *H* = *A*(*p*) + *B*(*q*). Here *chi* should be exp(*h* X<sub>B</sub>)exp(*h* X<sub>A</sub>).
    - `RKN6b` (order 6) refers to SRKN<sub>11</sub><sup>*b*</sup> optimized for *H* = *A*(*p*) + *B*(*q*). Here *chi* should be exp(*h* X<sub>B</sub>)exp(*h* X<sub>A</sub>).
    - `RKN6a` (order 6) refers to SRKN<sub>14</sub><sup>*a*</sup> optimized for *H* = *A*(*p*) + *B*(*q*). Here *chi* should be exp(*h* X<sub>A</sub>)exp(*h* X<sub>B</sub>).
- From [Blanes, Casas, Farrés, Laskar, Makazaga, Murua, Appl. Numer. Math. 68, 58 (2013)](http://dx.doi.org/10.1016/j.apnum.2013.01.003):
    - `ABA104` (order (10,4)) optimized for *H* = *A* + &epsilon; *B*. Here *chi* should be exp(*h* X<sub>A</sub>)exp(*h* X<sub>B</sub>).
    - `ABA864` (order (8,6,4)) optimized for *H* = *A* + &epsilon; *B*. Here *chi* should be exp(*h* X<sub>A</sub>)exp(*h* X<sub>B</sub>).
    - `ABA1064` (order (10,6,4)) optimized for *H* = *A* + &epsilon; *B*. Here *chi* should be exp(*h* X<sub>A</sub>)exp(*h* X<sub>B</sub>).
    
All purpose integrators are for any splitting of the Hamiltonian *H*=&sum;<sub>*k*</sub> *A*<sub>*k*</sub> in any order of the functions *A*<sub>*k*</sub>. Otherwise, the order of the operators is specified for each integrator. These integrators are used in the functions `solve_ivp_symp` and `solve_ivp_sympext` by specifying the entry `method` (default is `BM4`). 

----
## HamSys class

### Parameters
- `ndof` : number of degrees of freedom of the Hamiltonian system
       	'ndof' should be an integer or half an integer. Half integers denote an explicit time dependence.

### Attributes
- `hamiltonian` : callable  
	A function of (*t*, *y*) which returns the Hamiltonian *H*(*t*,*y*) where *y* is the state vector.
- `y_dot` : callable  
  	A function of (*t*, *y*) which returns {*y*,*H*(*t*,*y*)} where *y* is the state vector and *H* is the Hamiltonian. In canonical coordinates (used, e.g., in `solve_ivp_sympext`) where *y* = (*q*, *p*), this function returns (&part;*H*/&part;*p*, -&part;*H*/&part;*q*).
- `k_dot` : callable  
	A function of (*t*, *y*) which returns {*k*,*H*(*t*,*y*)} = -&part;*H*/&part;*t* where *k* is canonically conjugate to *t* and *H* is the Hamiltonian.

### Functions
- `compute_vector_field` : from a callable function (Hamiltonian in canonical coordinates) written with symbolic variables (SymPy), computes the vector fields, `y_dot` and `k_dot`.

	Determine Hamilton's equations of motion from a given scalar function &ndash;the Hamiltonian&ndash; *H*(*q*, *p*, *t*) where *q* and *p* are respectively positions and momenta.

	#### Parameters
	- `hamiltonian` : callable
		Function *H*(*q*, *p*, *t*) &ndash;the Hamiltonian expressed in symbolic variables&ndash;, expressed using [SymPy](https://www.sympy.org/en/index.html) functions.
	- `output` : bool, optional
		If True, displays the equations of motion. Default is False.
	
	The function `compute_vector_field` determines the HamSys function attributes `y_dot` and `k_dot` to be used in `solve_ivp_sympext`. The derivatives are computed symbolically using SymPy.

- `compute_energy` : callable
  	A function of `sol` &ndash;a solution provided by `solve_ivp_sympext`&ndash; and `maxerror`, a boolean indicating whether the maximum error in total energy is given (if True) or all the values of the total energy (if False). 
	#### Parameters
	- `sol` : OdeSolution  
   		Solution provided by `solve_ivp_sympext`. 
 	- `maxerror` : bool, optional  
    		Default is True.

---
## solve_ivp_symp and solve_ivp_sympext

The functions `solve_ivp_symp` and `solve_ivp_sympext` solve an initial value problem for a Hamiltonian system using an element of the class SymplecticIntegrator, an explicit symplectic splitting scheme (see [1]). These functions numerically integrate a system of ordinary differential equations given an initial value:  
	&nbsp; d*y* / d*t* = {*y*, *H*(*t*, *y*)}  
	&nbsp; *y*(*t*<sub>0</sub>) = *y*<sub>0</sub>  
Here *t* is a 1-D independent variable (time), *y*(*t*) is an N-D vector-valued function (state). A Hamiltonian *H*(*t*, *y*) and a Poisson bracket {. , .} determine the differential equations. The goal is to find *y*(*t*) approximately satisfying the differential equations, given an initial value *y*(*t*<sub>0</sub>) = *y*<sub>0</sub>. 

The function `solve_ivp_symp` solves an initial value problem using an explicit symplectic integration. The Hamiltonian flow is defined by two functions `chi` and `chi_star` of (*h*, *t*, *y*) (see [2]). This function works for any set of coordinates, canonical or non-canonical, provided that the splitting *H*=&sum;<sub>*k*</sub> *A*<sub>*k*</sub> leads to facilitated expressions for the operators exp(*h* X<sub>*k*</sub>) where X<sub>*k*</sub> = {*A*<sub>*k*</sub> , &centerdot;}.

The function `solve_ivp_sympext` solves an initial value problem using an explicit symplectic approximation obtained by an extension in phase space (see [3]). This symplectic approximation works for canonical Poisson brackets, and the state vector should be of the form *y* = (*q*, *p*). 

### Parameters:  

  - `chi` (for `solve_ivp_symp`) : callable  
	Function of (*h*, *t*, *y*) returning exp(*h* X<sub>*n*</sub>)...exp(*h* X<sub>1</sub>) *y* at time *t*. If the selected integrator is not all purpose, refer to the list above for the specific ordering of the operators. The operator X<sub>*k*</sub> is the Liouville operator associated with the function *A*<sub>*k*</sub>, i.e., for Hamiltonian flows X<sub>*k*</sub> = {*A*<sub>*k*</sub> , &centerdot;} where {&centerdot; , &centerdot;} is the Poisson bracket.
	`chi` must return an array of the same shape as `y`.
  - `chi_star` (for `solve_ivp_symp`) : callable   
	Function of (*h*, *t*, *y*) returning exp(*h* X<sub>1</sub>)...exp(*h* X<sub>*n*</sub>) *y* at time *t*.
	`chi_star` must return an array of the same shape as `y`.
  - `hs` (for `solve_ivp_sympext`) : element of class HamSys  
	The attributes `y_dot` of `hs` should be defined. If `check_energy` is True. It the Hamiltonian system has an explicit time dependence (i.e., the parameter `ndof` of `hs`  is half an integer), the attribute `k_dot` of `hs` should be specified. 
  - `t_span` : 2-member sequence  
	Interval of integration (*t*<sub>0</sub>, *t*<sub>f</sub>). The solver starts with *t*=*t*<sub>0</sub> and integrates until it reaches *t*=*t*<sub>f</sub>. Both *t*<sub>0</sub> and *t*<sub>f</sub> must be floats or values interpretable by the float conversion function.	
  - `y0` : array_like, shape (n,)  
	Initial state.
  - `step` : float   
	Step size.
  - `t_eval` : array_like or None, optional  
	Times at which to store the computed solution, must be sorted and equally spaced, and lie within `t_span`. If None (default), use points selected by the solver.
  - `method` : string, optional  
 	Integration methods are listed on [pyhamsys](https://pypi.org/project/pyhamsys/).   
	'BM4' is the default.
  - `omega` (for `solve_ivp_sympext`) : float, optional  
   	Coupling parameter in the extended phase space (see [3]). Default = 10.
  - `command` : function of (*t*, *y*)  
	Function to be run at each step size (e.g., plotting an observable associated with the state vector *y*, or register specific events).
  - `check_energy` (for `solve_ivp_sympext`) : bool, optional  
	If True, the attribute `hamiltonian` of `hs` should be defined. Default is False. 

### Returns:  
&nbsp; Bunch object with the following fields defined:
   - `t` : ndarray, shape (n_points,)  
	Time points.
   - `y` : ndarray, shape (n, n_points)  
	Values of the solution `y` at `t`.
   - `k` (for `solve_ivp_sympext`) : ndarray, shape (n//2, n_points)
     	Values of `k` at `t`. Only for `solve_ivp_sympext` and if `check_energy` is True for a Hamiltonian system with an explicit time dependence (i.e., the parameter `ndof` of `hs`  is half an integer).
   - `err` (for `solve_ivp_sympext`) : float
     	Error in the computation of the total energy. Only for `solve_ivp_sympext` and if `check_energy` is True.
   - `step` : step size used in the computation.

### Remarks:   
  - Use `solve_ivp_symp` is the Hamiltonian can be split and if each partial operator exp(*h* X<sub>*k*</sub>) can be easily expressed/computed. Otherwise use `solve_ivp_sympext` if your coordinates are canonical.  
  - If `t_eval` is a linearly spaced list or array, or if `t_eval` is None (default), the step size is slightly readjusted so that the output times contain the values in `t_eval`, or the final time *t*<sub>f</sub> corresponds to an integer number of step sizes. The step size used in the computation is recorded in the solution as `sol.step`.  

### References:  
  - [1] Hairer, Lubich, Wanner, 2003, *Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations* (Springer)  
  - [2] McLachlan, *Tuning symplectic integrators is easy and worthwhile*, Commun. Comput. Phys. 31, 987 (2022); [arxiv:2104.10269](https://arxiv.org/abs/2104.10269)  
  - [3] Tao, M., *Explicit symplectic approximation of nonseparable Hamiltonians: Algorithm and long time performance*, Phys. Rev. E 94, 043303 (2016)

### Example

```python
import numpy as xp
import sympy as sp
import matplotlib.pyplot as plt
from pyhamsys import HamSys, solve_ivp_sympext
hs = HamSys()
hamiltonian = lambda q, p, t: p**2 / 2 - sp.cos(q)
hs.compute_vector_field(hamiltonian, output=True)
sol = solve_ivp_sympext(hs, (0, 20), xp.asarray([3, 0]), step=1e-1, check_energy=True)
print(f"Error in energy : {sol.err}")
plt.plot(sol.y[0], sol.y[1])
plt.show()
```
---
For more information: <cristel.chandre@cnrs.fr>
