# pyHamSys
pyHamSys is a Python package for scientific computations involving Hamiltonian systems

## Symplectic Integrators
PyHamSys includes a class SymplecticIntegrator containing the following integrators:

Names of integrators include:
- 'Verlet' (order 2)
- 'Forest-Ruth' (order 4) from [Forest, Ruth, Physica D 43, 105 (1990)](https://doi.org/10.1016/0167-2789(90)90019-L)
- 'EFRL', 'PEFRL' or 'VEFRL' (order 4) from [Omelyan, Mryglod, Folk, Comput. Phys. Commun. 146, 188 (2002)](https://doi.org/10.1016/S0010-4655(02)00451-4). Optimized for *H* = *A* + *B*: 
- 'BM4' (order 4) refers to PRK<sub>6</sub>4 from [Blanes, Moan, J. Comput. Appl. Math. 142, 313 (2002)](https://doi.org/10.1016/S0377-0427(01)00492-7)
- 'BM6' (order 6) refers to PRK<sub>10</sub>6 from [Blanes, Moan, J. Comput. Appl. Math. 142, 313 (2002)](https://doi.org/10.1016/S0377-0427(01)00492-7)

Usage: *integrator* = SymplecticIntegrator(*name*, *step*, *order*)
where *name* is one of the names listed above, *step* is the time step of the integrator (float), and *order* is the order of the splitting, so 1 or -1 depending on the order AB or BA of the splitting. 

The function *integrator*.`_integrate` integrates the Hamiltonian flow by one step.

The function *integrator*.`integrate` integrates the Hamiltonian flow from the initial conditions specified by the initial state vector *y* using *integrator*, one of the selected symplectic splitting integrators. It returns the value of *y* at times defines by the float, list or array *times*.

Parameters:
  - chi : function of (*h*, *y*), y being the state vector.
    Function returning exp(*h* X<sub>*n*</sub>)...exp(*h* X<sub>1</sub>) *y*.
  - chi_star : function of (*h*, *y*).
    Function returning exp(*h* X<sub>1</sub>)...exp(*h* X_<sub>*n*</sub>) *y*.
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
  - McLachlan, R.I, 2022, *Tuning symplectic integrators is easy and worthwhile*, [arxiv:2104.10269](https://arxiv.org/abs/2104.10269)


Installation: 
```
pip install pyhamsys
```
