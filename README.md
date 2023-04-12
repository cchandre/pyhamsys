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
    - `PEFRL` and `VEFRL` (order 4) optimized for *H* = *A*(*p*) + *B*(*q*). For `PEFRL`, *chi* should be exp(*h*A)exp(*h*B). For `VEFRL`, *chi* should be exp(*h*B)exp(*h*A).
- From [Blanes, Moan, J. Comput. Appl. Math. 142, 313 (2002)](https://doi.org/10.1016/S0377-0427(01)00492-7):
    - `BM4` (order 4, all purpose) refers to S<sub>6</sub> 
    - `BM6` (order 6, all purpose) refers to S<sub>10</sub>
    - `RKN4b` (order 4) refers to SRKN<sub>6</sub><sup>*b*</sup> optimized for *H* = *A*(*p*) + *B*(*q*). Here *chi* should be exp(*h*B)exp(*h*A).
    - `RKN6b` (order 6) refers to SRKN<sub>11</sub><sup>*b*</sup> optimized for *H* = *A*(*p*) + *B*(*q*). Here *chi* should be exp(*h*B)exp(*h*A).
    - `RKN6a` (order 6) refers to SRKN<sub>14</sub><sup>*a*</sup> optimized for *H* = *A*(*p*) + *B*(*q*). Here *chi* should be exp(*h*A)exp(*h*B).
- From [Blanes, Casas, Farrés, Laskar, Makazaga, Murua, Appl. Numer. Math. 68, 58 (2013)](http://dx.doi.org/10.1016/j.apnum.2013.01.003):
    - `ABA104` (order (10,4)) optimized for *H* = *A* + &epsilon; *B*. Here *chi* should be exp(*h*A)exp(*h*B).
    - `ABA864` (order (8,6,4)) optimized for *H* = *A* + &epsilon; *B*. Here *chi* should be exp(*h*A)exp(*h*B).
    - `ABA1064` (order (10,6,4)) optimized for *H* = *A* + &epsilon; *B*. Here *chi* should be exp(*h*A)exp(*h*B).
    
All purpose integrators are for any splitting of the Hamiltonian *H*=&sum;<sub>*k*</sub> *X*<sub>*k*</sub> in any order of the functions *X*<sub>*k*</sub>. Otherwise, the order of the operators is specified for each integrator.

Usage: *integrator* = SymplecticIntegrator(*name*, *step*)
where *name* is one of the names listed above and *step* is the time step of the integrator (float). 

The function *integrator*.`_integrate` integrates the Hamiltonian flow by one step.

The function *integrator*.`integrate` integrates the Hamiltonian flow from the initial conditions specified by the initial state vector *y* using *integrator*, one of the selected symplectic splitting integrators. It returns the value of *y* at times defines by the float, list or array *times*.

Parameters:
  - *chi* : function of (*h*, *y*), y being the state vector.
    Function returning exp(*h* X<sub>*n*</sub>)...exp(*h* X<sub>1</sub>) *y*. If the selected integrator is not all purpose, refer to the list above for the specific ordering of the operators. The operator X<sub>*k*</sub> is the Liouville operator associated with the function *X*<sub>*k*</sub>, i.e., X<sub>*k*</sub> = {*X*<sub>*k*</sub>, &centerdot;}.
  - *chi_star* : function of (*h*, *y*).
    Function returning exp(*h* X<sub>1</sub>)...exp(*h* X<sub>*n*</sub>) *y*.
  - *y* : initial state vector (numpy array)
  - *times* : times at which the values of the state vector are computed
  - *command* : function of (*t*, *y*).
    Function to be run at each time step (e.g., plotting an observable associated with the state vector, or register specific events). 
  - *autonomous* : boolean.
    If autonomous is False, the state vector y should be of the form *y* = [*t*, *x*], where the first coordinate is time. 

Returns:
   - If *times* is a float of integer, the output is a tuple (*t*, *y* or *x*) where *y* is the value of the state vector and *y* = [*t*, *x*] if autonomous is False.
   - If *times* is a list or array, returns the times and values of *y* or *x* at *times*. 
   - If *times* is a list or array with a single element, returns the times and values of *y* or *x* at all computed times. 

References:
  - Hairer, Lubich, Wanner, 2003, *Geometric-Preseving Algorithms for Ordinary Differential Equations* (Springer)
  - McLachlan, 2022, *Tuning symplectic integrators is easy and worthwhile*, [arxiv:2104.10269](https://arxiv.org/abs/2104.10269)
