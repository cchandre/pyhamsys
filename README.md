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

Usage: *integrator* = SymplecticIntegrator(*name*, *step*) where *name* is one of the names listed above and *step* is the time step of the integrator (float). 

The functions `solve_ivp_symp` and `solve_ivp_sympext` solve an initial value problem for a Hamiltonian system using an explicit symplectic splitting scheme (see [1]). These functions numerically integrate a system of ordinary differential equations given an initial value::
	d*y* / d*t* = {y, H(t, y)}
	y(t0) = y0

		Here t is a 1-D independent variable (time), y(t) is an N-D 
		vector-valued function (state), and a Hamiltonian H(t, y) and a 
		Poisson bracket {. , .} determine the differential equations.
		The goal is to find y(t) approximately satisfying the differential
		equations, given an initial value y(t0)=y0. The Hamiltonian flow
		is defined by two functions chi and chi_star (see [2]).

The function *integrator*.`integrate` integrates the Hamiltonian flow from the initial conditions specified by the initial state vector *y* using *integrator*, one of the selected symplectic splitting integrators. It returns the value of *y* at times defines by the float, list or numpy array *times*.
The function *integrator*.`integrate` integrates the Hamiltonian flow from the initial conditions specified by the initial state vector *y* using *integrator*, one of the selected symplectic splitting integrators. It returns the value of *y* at times defines by the integer, float, list or numpy array t_span.



		

Parameters:  
  - *chi* : function of (*h*, *y*), *y* being the state vector.
    Function returning exp(*h* X<sub>*n*</sub>)...exp(*h* X<sub>1</sub>) *y*. If the selected integrator is not all purpose, refer to the list above for the specific ordering of the operators. The operator X<sub>*k*</sub> is the Liouville operator associated with the function *A*<sub>*k*</sub>, i.e., for Hamiltonian flows X<sub>*k*</sub> = {*A*<sub>*k*</sub> , &centerdot;} where {&centerdot; , &centerdot;} is the Poisson bracket.
  - *chi_star* : function of (*h*, *y*).
    Function returning exp(*h* X<sub>1</sub>)...exp(*h* X<sub>*n*</sub>) *y*.
  - *y* : initial state vector (numpy array)
  - *times* : times at which the values of the state vector are computed
  - *command* : function of (*t*, *y*).
    Function to be run at each time step (e.g., plotting an observable associated with the state vector, or register specific events). 

Remarks:  
    - If the vector field is explicitly time dependent, it should be first autonomized by adding time as an extra variable  
    - If *times* is a linearly spaced list or array, or if *times* is an integer or a float, the time step is slightly readjusted so that the output times contain the values in *times*  
    - If *times* is not linearly spaced, a linear interpolation of the solution is performed; the accuracy of the integrator might be lost   

Returns:  
&nbsp; Bunch object with the following fields defined:
   - t : final integration time if *times* is a float or integer  
         &nbsp; *times* if *times* is a list or a numpy array  
         &nbsp; all computed times if *times* is a list or numpy array with a single element
   - y : state vector at times t
   - time_step : time step used in the computation

References:
  - Hairer, Lubich, Wanner, 2003, *Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations* (Springer)
  - McLachlan, *Tuning symplectic integrators is easy and worthwhile*, Commun. Comput. Phys. 31, 987 (2022); [arxiv:2104.10269](https://arxiv.org/abs/2104.10269)
