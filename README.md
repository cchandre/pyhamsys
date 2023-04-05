# pyHamSys
pyHamSys is a Python package for scientific computations involving Hamiltonian systems

## Symplectic Integrators
PyHamSys includes a class SymplecticIntegrator containing the following integrators:

Names of integrators include:
- 'Verlet' (order 2)
- 'Forest-Ruth' (order 4) from [Forest, Ruth, Physica D 43, 105 (1990)](https://doi.org/10.1016/0167-2789(90)90019-L)
- 'PEFRL' (order 4) from [Omelyan, Mryglod, Folk, Comput. Phys. Commun. 146, 188 (2002)](https://doi.org/10.1016/S0010-4655(02)00451-4)
- 'BM4' (order 4) and 'BM6' (order 6) refer respectively to BM<sub>6</sub>4 and BM<sub>10</sub>6 from [Blanes, Moan, J. Comput. Appl. Math. 142, 313 (2002)](https://doi.org/10.1016/S0377-0427(01)00492-7)

Usage: integrator = SymplecticIntegrator(*name*, *step*, *order*)
where *name* is one of the names listed above, *step* is the time step of the integrator (float), and *order* is the order of the splitting, so 1 or -1 depending on the order AB or BA of the splitting. 

The function `integrate` integrates the Hamiltonian flow from the initial conditions specified by the initial state vector y using one of the selected symplectic splitting integrators. It returns the value of y at times defines by the float, list or array times.

    Parameters
      chi : function of (h, y), y being the state vector
        function returning exp(h X_n)...exp(h X_1) y.
      chi_star : function of (h, y)
        function returning exp(h X_1)...exp(h X_n) y.
      command : function of (t, y) 
        function to be run at each time step. 
      autonomous : boolean
        if autonomous is False, the state vector y should be of the form 
        y = [t, x]. 

    Returns
      out : array of times, and values of y at times
        If times is a float of integer, the output is a tuple (t, y or x) where
        y is the value of the state vector and y = [t, x] if autonomous is False.
        If times is a list or array, returns the times and values of y or x at 
        times. If times is a list or array with a single element, returns the 
        times and values of y or x at all computed times. 

    References
        McLachlan, R.I, 2022, "Tuning symplectic integrators is easy and 
        worthwhile", *arxiv:2104.10269*. 
