# pyHamSys
pyHamSys is a Python package for scientific computing involving Hamiltonian systems

![PyPI](https://img.shields.io/pypi/v/pyhamsys)
![License](https://img.shields.io/badge/license-BSD-lightgray)

Installation: 
```
pip install pyhamsys
```

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
    
All purpose integrators are for any splitting of the Hamiltonian *H*=&sum;<sub>*k*</sub> *A*<sub>*k*</sub> in any order of the functions *A*<sub>*k*</sub>. Otherwise, the order of the operators is specified for each integrator.

Usage: *integrator* = SymplecticIntegrator(*name*) where *name* is one of the names listed above. 

The functions `solve_ivp_symp` and `solve_ivp_sympext` solve an initial value problem for a Hamiltonian system using an element of the class SymplecticIntegrator, an explicit symplectic splitting scheme (see [1]). These functions numerically integrate a system of ordinary differential equations given an initial value:  
	&nbsp; d*y* / d*t* = {*y*, *H*(*t*, *y*)}  
	&nbsp; *y*(*t*<sub>0</sub>) = *y*<sub>0</sub>  
Here *t* is a 1-D independent variable (time), *y*(*t*) is an N-D vector-valued function (state). A Hamiltonian *H*(*t*, *y*) and a Poisson bracket {. , .} determine the differential equations. The goal is to find *y*(*t*) approximately satisfying the differential equations, given an initial value *y*(*t*<sub>0</sub>) = *y*<sub>0</sub>. 

The function `solve_ivp_symp` solves an initial value problem using an explicit symplectic integration. The Hamiltonian flow is defined by two functions `chi` and `chi_star` of (*h*, *t*, *y*) (see [2]). 

The function `solve_ivp_sympext` solves an initial value problem using an explicit symplectic approximation obtained by an extension in phase space (see [3]). The Hamiltonian flow is defined by one function `fun` of (*t*, *y*) and one coupling parameter `omega`. 

### Parameters:  

  - `chi` (for `solve_ivp_symp`) : callable  
	Function of (*h*, *t*, *y*) returning exp(*h* X<sub>*n*</sub>)...exp(*h* X<sub>1</sub>) *y* at time *t*. If the selected integrator is not all purpose, refer to the list above for the specific ordering of the operators. The operator X<sub>*k*</sub> is the Liouville operator associated with the function *A*<sub>*k*</sub>, i.e., for Hamiltonian flows X<sub>*k*</sub> = {*A*<sub>*k*</sub> , &centerdot;} where {&centerdot; , &centerdot;} is the Poisson bracket.
	`chi` must return an array of the same shape as `y`.
  - `chi_star` (for `solve_ivp_symp`) : callable   
	Function of (*h*, *t*, *y*) returning exp(*h* X<sub>1</sub>)...exp(*h* X<sub>*n*</sub>) *y* at time *t*.
	`chi_star` must return an array of the same shape as `y`.
  - `fun` (for `solve_ivp_sympext`) : callable  
	Right-hand side of the system: the time derivative of the state *y* at time *t*. i.e., {*y*, *H*(*t*, *y*)}. The calling signature is `fun(t, y)`, where `t` is a scalar and `y` is an ndarray with `len(y) = len(y0)`. `fun` must return an array of the same shape as `y`. 
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

### Remarks:   
  - Use `solve_ivp_symp` is the Hamiltonian can be split and if each partial flow exp(*h* X<sub>*k*</sub>) can be easily computed. Otherwise use `solve_ivp_sympext`.  
  - If `t_eval` is a linearly spaced list or array, or if `t_eval` is None (default), the step size is slightly readjusted so that the output times contain the values in `t_eval`, or the final time *t*<sub>f</sub> corresponds to an integer number of step sizes.  

### Returns:  
&nbsp; Bunch object with the following fields defined:
   - `t` : ndarray, shape (n_points,)  
	Time points.
   - `y` : ndarray, shape (n, n_points)  
	Values of the solution at `t`.
   - `step` : step size used in the computation.

### References:  
  - [1] Hairer, Lubich, Wanner, 2003, *Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations* (Springer)  
  - [2] McLachlan, *Tuning symplectic integrators is easy and worthwhile*, Commun. Comput. Phys. 31, 987 (2022); [arxiv:2104.10269](https://arxiv.org/abs/2104.10269)  
  - [3] Tao, M., *Explicit symplectic approximation of nonseparable Hamiltonians: Algorithm and long time performance*, Phys. Rev. E 94, 043303 (2016)

### Example

```python
>>> import numpy as xp
>>> import matplotlib.pyplot as plt
>>> from pyhamsys import solve_ivp_sympext
>>> def fun(t,y):
	x, p = xp.split(y, 2)
	return xp.concatenate((p, -xp.sin(x)), axis=None)
>>> sol = solve_ivp_sympext(fun, (0, 20), xp.asarray([3, 0]), step=1e-1)
>>> plt.plot(sol.y[0], sol.y[1])
>>> plt.show()
```
