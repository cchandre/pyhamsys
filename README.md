# pyHamSys
pyHamSys (short for Python Hamiltonian Systems) is a Python library for scientific computing involving Hamiltonian systems. It provides tools to model, analyze, and simulate Hamiltonian systems. In particular, the library offers a streamlined and user-friendly environment for implementing and running symplectic-split integrators. These specialized numerical methods are designed to preserve the geometric structure of Hamiltonian systems, ensuring accurate and stable simulations of their dynamics over long time periods. 

![PyPI](https://img.shields.io/pypi/v/pyhamsys)
![License](https://img.shields.io/badge/license-BSD-lightgray)
![PyPI - Downloads](https://img.shields.io/pypi/dm/pyhamsys.svg?label=PyPI%20downloads)

## Installation 
Installation within a Python virtual environment: 
```
python3 -m pip install pyhamsys
```
For more information on creating a Python virtual environment, click [here](https://realpython.com/python-virtual-environments-a-primer/). For a summary with the main steps, click [here](https://github.com/cchandre/HamLorenz/wiki/Python-Virtual-Environment-Primer).


## Symplectic Integrators

pyHamSys features a dedicated class, `SymplecticIntegrator`, which provides a comprehensive implementation of symplectic-split integrators. These integrators are designed specifically for the numerical integration of Hamiltonian systems, ensuring the preservation of the symplectic structure of phase space—a key property that underpins the stability and accuracy of long-term simulations of such systems.
Symplectic-split integrators decompose the Hamiltonian into subcomponents that are individually solvable, allowing for efficient and accurate integration. This decomposition is particularly effective for complex or high-dimensional systems, as it minimizes numerical drift and conserves critical invariants like energy over extended time intervals.
The `SymplecticIntegrator` class offers a variety of splitting methods, enabling users to select the most appropriate scheme for their specific Hamiltonian system and computational requirements. Each integrator is implemented to handle both autonomous and non-autonomous systems, supporting applications in classical mechanics, molecular dynamics, astrophysics, and quantum mechanics.


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
The `HamSys` class provides a robust framework for defining and integrating Hamiltonian systems. It allows users to specify the number of degrees of freedom, coordinate representations, and key attributes like the Hamiltonian and associated equations of motion.

### Parameters
- `ndof` : integer or half-integer, optional   
       	The number of degrees of freedom in the Hamiltonian system. Half integers denote an explicit time dependence. Default is 1.  
- `btype` : str, optional
  	Information on the Poisson bracket used in the equations of motion. For btype='pq', a canonical Poisson bracket in (p,q) is used, i.e., the dynamical variables (*q*, *p*) are real and canonically conjugate. If btype='psi', the dynamical variables are (&psi;, &psi;<sup>*</sup>) where $\psi=(q + i p)/\sqrt{2}$. Default is 'pq'. For other btype, the function `coupling` should be specified for the element of the class `HamSys`. 

### Parameters and Attributes
- `y_dot` : callable, optional   
  	A function of (*t*, *y*) which returns {*y*,*H*(*t*,*y*)} where *y* is the state vector and *H* is the Hamiltonian. In (real) canonical coordinates (used, e.g., in `solve_ivp_sympext`) where *y* = (*q*, *p*), this function returns (&part;*H*/&part;*p*, -&part;*H*/&part;*q*). In complex coordinate &psi;, this function returns -i &part;*H*/&part;&psi;<sup>*</sup>. For practical implementation, the state vector *y* should be represented as a one-dimensional array with a shape of (*n*,), where *n* denotes the total number of dynamical variables in the system. This ensures compatibility with numerical solvers and facilitates efficient computation of the system's evolution.  
- `k_dot` : callable, optional   
	A function of (*t*, *y*) which returns {*k*,*H*(*t*,*y*)} = -&part;*H*/&part;*t* where *k* is canonically conjugate to *t* and *H* is the Hamiltonian.
- `hamiltonian` : callable, optional   
	A function of (*t*, *y*) which returns the Hamiltonian *H*(*t*,*y*) where *y* is the state vector.
- `coupling` : callable, optional
  	A function of (*h*, *y*, &omega;) which advances *y* from time *t* to *t*+*h* for the coupling Hamiltonian $\omega (y - \bar{y})^2/2$. This function is already computed for the types btype='pq' and 'psi'. For any other type, it should be provided. 

### Functions
- `compute_vector_field` : from a callable function (Hamiltonian in canonical coordinates) written with symbolic variables ([SymPy](https://www.sympy.org/en/index.html)), computes the vector fields, `y_dot` and `k_dot`.

	Determine Hamilton's equations of motion from a given scalar function &ndash;the Hamiltonian&ndash; *H*(*q*, *p*, *t*) where *q* and *p* are respectively positions and momenta. However, it is preferrable to code explicitly and optimize `y_dot` and `k_dot`.

	#### Parameters
	- `hamiltonian` : callable   
		Function *H*(*q*, *p*, *t*) &ndash;the Hamiltonian expressed in symbolic variables&ndash;, expressed using [SymPy](https://www.sympy.org/en/index.html) functions.
	- `output` : bool, optional   
		If True, displays the equations of motion. Default is False.
	
	The function `compute_vector_field` determines the HamSys function attributes `y_dot` and `k_dot` to be used in `solve_ivp_sympext`. The derivatives are computed symbolically using SymPy.

- `compute_energy` : callable   
  	A function of `sol` &ndash;a solution provided by `solve_ivp_sympext`&ndash; and a boolean `maxerror`. If `maxerror`is `True` the function `compute_energy` returns the maximum error in total energy; otherwise, it returns all the values of the total energy. 
	#### Parameters
	- `sol` : OdeSolution  
   		Solution provided by `solve_ivp_sympext`. 
 	- `maxerror` : bool, optional  
    		Default is True.

- `integrate` : callable   
    Integrate the Hamiltonian system using either a **symplectic solver** (see above for a complete list) or a **standard IVP solver** ('RK23', 'RK45', 'DOP853', 'Radau', 'BDF', 'LSODA'). Supports optional *symplectic extension* and *energy conservation checks*.

   #### Parameters
   - **z0** (`array_like`)  
  Initial condition(s) of the system.
   - **t_eval** (`array_like`)  
  Times at which the solution is evaluated and stored.
   - **timestep** (`float`)  
  Fixed integration time step (used in symplectic solvers and to bound steps in IVP solvers).
   - **solver** (`str`, optional, default=`"BM4"`)  
  Solver method. Must be a member of `METHODS` (symplectic solvers), or  `IVP_METHODS` (classical IVP solvers).  
   - **extension** (`bool`, optional, default=`False`)  
  If `True`, use a symplectic extension method in phase space.  
   - **check_energy** (`bool`, optional, default=`False`)  
  If `True`, appends an auxiliary variable to track the Hamiltonian. Requires `hamiltonian` and `k_dot` to be defined.  
   - **omega** (`float`, optional, default=`10`)  
  Restrain parameter for symplectic extension solvers.  
   - **tol** (`float`, optional, default=`1e-8`)  
  Absolute and relative tolerance for IVP solvers.  
   - **display** (`bool`, optional, default=`True`)  
  If `True`, prints runtime information such as CPU time, error in energy, and copy distance (if available).  

    #### Returns
   - **sol** (`object`)  
    Solution object. Its attributes depend on the solver used:
     - `y` : state trajectory  
     - `t` : time points  
     - `step` : integration time step  
     - `cpu_time` : total CPU time used  
     - `err` : (if `check_energy=True`) maximum error in energy  
     - `dist_copy` : (if applicable) distance between trajectory copies  

    #### Notes
    - **Symplectic solvers (`METHODS`)**  
    Require `chi` and `chi_star` to be defined in the class. Preserves geometric properties of Hamiltonian flows.  
    - **IVP solvers (`IVP_METHODS`)**  
    Require `y_dot` (and `k_dot` if `check_energy=True`). Allow adaptive step sizes bounded by `timestep`.  
    - **Energy checking**  
    When `check_energy=True`, an auxiliary variable is added and the error in Hamiltonian is computed relative to its initial value. 

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
	The attribute `y_dot` of `hs` should be defined. If `check_energy` is True, the attribute `hamiltonian`  and if the Hamiltonian system has an explicit time dependence (i.e., the parameter `ndof` of `hs`  is a half-integer), the attribute `k_dot` of `hs` should be specified. 
  - `t_span` : 2-member sequence  
	Interval of integration (*t*<sub>0</sub>, *t*<sub>f</sub>). The solver starts with *t*=*t*<sub>0</sub> and integrates until it reaches *t*=*t*<sub>f</sub>. Both *t*<sub>0</sub> and *t*<sub>f</sub> must be floats or values interpretable by the float conversion function.	
  - `y0` : array_like  
	Initial state. For `solve_ivp_sympext`, the vector `y0` should be with shape (n,).
  - `step` : float   
	Step size.
  - `t_eval` : array_like or None, optional  
	Times at which to store the computed solution, must be sorted, and lie within `t_span`. If None (default), use points selected by the solver.
  - `method` : string, optional  
 	Integration methods are listed on [pyhamsys](https://pypi.org/project/pyhamsys/). Default is 'BM4'.
  - `omega` (for `solve_ivp_sympext`) : float, optional  
   	Coupling parameter in the extended phase space (see [3]). Default is 10.
  - `command` : void function of (*t*, *y*), optional    
	Void function to be run at each step size (e.g., plotting an observable associated with the state vector *y*, modify global or mutable variables, or register specific events).
  - `check_energy` (for `solve_ivp_sympext`) : bool, optional  
	If True, the attribute `hamiltonian` of `hs` should be defined. Default is False. 

### Returns:  
&nbsp; Bunch object with the following fields defined:
   - `t` : ndarray, shape (n_points,)  
	Time points.
   - `y` : ndarray, shape y0.shape + (n_points,)  
	Values of the solution `y` at `t`.
   - `k` (for `solve_ivp_sympext`) : ndarray, shape (n_points,)   
     	Values of `k` at `t`. Only for `solve_ivp_sympext` and if `check_energy` is True for a Hamiltonian system with an explicit time dependence (i.e., the parameter `ndof` of `hs`  is half an integer).
   - `dist_copy` (for `solve_ivp_sympext`) : float   
        Maximum distance between the two copies of the state in the extended phase space.   
   - `err` (for `solve_ivp_sympext`) : float   
     	Error in the computation of the total energy. Only for `solve_ivp_sympext` and if `check_energy` is True.
   - `step` : step size used in the computation.

### Remarks:   
  - Use `solve_ivp_symp` if the Hamiltonian can be split and if each partial operator exp(*h* X<sub>*k*</sub>) can be easily and explicitly expressed/computed. Otherwise use `solve_ivp_sympext` if your coordinates are canonical, i.e., in $(q,p)$ or $(\psi,\psi^*)$ variables.  
  - The step size is slightly readjusted so that the final time *t*<sub>f</sub> corresponds to an integer number of step sizes. The step size used in the computation is recorded in the solution as `sol.step`.
  - For integrating multiple trajectories at the same time, extend phase space and define a state vector y = (y<sub>1</sub>, y<sub>2</sub>,...y<sub>N</sub>) where N is the number of trajectories. The Hamiltonian is given by $H(t,\mathbf{y})=\sum_{i=1}^N h(t, y_i)$.

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
For more examples, see the folder [Examples](https://github.com/cchandre/pyhamsys/tree/main/Examples)

---

This work has been carried out within the framework of the French Federation for Magnetic Fusion Studies (FR-FCM).

---
For more information: <cristel.chandre@cnrs.fr>
