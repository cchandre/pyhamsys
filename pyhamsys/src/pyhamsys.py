#
# BSD 2-Clause License
#
# Copyright (c) 2023, Cristel Chandre
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as xp
from scipy.fft import rfft, irfft, rfftfreq
from typing import Callable, Union, Tuple
from scipy.optimize import OptimizeResult
import sympy as sp
	
def eqns_of_motion(hamiltonian:Callable, ndof:int=1, output:bool=False) -> Callable:
	"""
	Determine Hamilton's equations of motion from a given Hamiltonian 
	H(q, p, t) where q and p are N-D vector (resp., positions and momenta).
	The output is a NumPy function providing the equations of motion ready to
	use in solve_ivp and solve_ivp_sympext. 

	Parameters
	----------
	hamiltonian : callable
		Function H(q, p, t) expressed with SymPy functions.
		`hamiltonian` must return a scalar.
	ndof : int, optional
		Number of degrees of freedom, i.e., number of positions. Default is 1.
	output : bool, optional
		If True, displays the equations of motion. Default is False.

	Returns
	-------
	Function of (t, y) where y = (q, p) returning (dH/dp, -dH/dq). If there is
	an explicit dependence on time, this function returns 
	(dH/dp, -dH/dq, -dH/dt). The input y and the output are ndarrays.  
	"""
	q = sp.symbols('q0:%d'%ndof) if ndof>=2 else sp.Symbol('q')
	p = sp.symbols('p0:%d'%ndof) if ndof>=2 else sp.Symbol('p')
	t = sp.Symbol('t')
	eqn = sp.simplify(sp.derive_by_array(hamiltonian(q, p, t), [q, p]))
	eqn = sp.flatten([eqn[1], -eqn[0]])
	if sp.diff(hamiltonian(q, p, t), t)!=0:
		eqn += [-sp.simplify(sp.diff(hamiltonian(q, p, t), t))]	
	if output:
		print(eqn)
	eqn = sp.lambdify([q, p, t], eqn)
	def fun(t, y):
		y_ = xp.split(y.flatten(), 2)
		return xp.asarray(eqn(y_[0], y_[1], t)).flatten()
	return fun

def antiderivative(vec:xp.ndarray, N:int=2**10) -> xp.ndarray:
	nu = rfftfreq(N, d=1/N)
	div = xp.divide(1, 1j * nu, where=nu!=0)
	div[0] = 0
	return irfft(rfft(vec, axis=-1) * div, axis=-1)

def padwrap(vec:xp.ndarray) -> xp.ndarray:
	return xp.concatenate((vec, vec[..., 0][..., xp.newaxis]), axis=-1)

def get_last_elements(a:xp.ndarray, axis:int=0) -> None:
    shape = list(a.shape)
    shape[axis] = 1
    return xp.take(a, -1, axis=axis).reshape(tuple(shape))

def cart2sph(state:xp.ndarray, dim:int=2) -> xp.ndarray:
    if dim == 2:
        x, y, px, py = state
        r, phi = xp.hypot(x, y), xp.arctan2(y, x)
        p_r = (x * px + y * py) / r
        p_phi = x * py - y * px
        return xp.concatenate(([r], [phi], [p_r], [p_phi]), axis=0)
    elif dim == 3:
        x, y, z, px, py, pz = state
        xy, phi = xp.hypot(x, y), xp.arctan2(y, x)
        r, theta = xp.hypot(xy, z), xp.arctan2(xy, z)
        p_r = (x * px + y * py + z * pz) / r
        p_theta = ((x * px + y * py) * z - pz * xy**2) / xy
        p_phi = x * py - y * px
        return xp.concatenate(([r], [theta], [phi], [p_r], [p_theta], [p_phi]), axis=0)

def sph2cart(state:xp.ndarray, dim:int=2) -> xp.ndarray:
    if dim == 2:
        r, phi, p_r, p_phi = state
        x, y = r * xp.cos(phi), r * xp.sin(phi)
        px = p_r * xp.cos(phi) - p_phi * xp.sin(phi) / r
        py = p_r * xp.sin(phi) + p_phi * xp.cos(phi) / r
        return xp.concatenate(([x], [y], [px], [py]), axis=0)
    elif dim == 3:
        r, theta, phi, p_r, p_theta, p_phi = state
        x, y, z = r * xp.sin(theta) * xp.cos(phi), r * xp.sin(theta) * xp.sin(phi), r * xp.cos(theta)
        px = p_r * xp.sin(theta) * xp.cos(phi) + p_theta * xp.cos(theta) * xp.cos(phi) / r - p_phi * xp.sin(phi) / (r * xp.sin(theta))
        py = p_r * xp.sin(theta) * xp.sin(phi) + p_theta * xp.cos(theta) * xp.sin(phi) / r + p_phi * xp.cos(phi) / (r * xp.sin(theta))
        pz = p_r * xp.cos(theta) - p_theta * xp.sin(theta) / r
        return xp.concatenate(([x], [y], [z], [px], [py], [pz]), axis=0)
    
def rotating(state:xp.ndarray, angle:float, type:str='cartesian', dim:int=2) -> xp.ndarray:
    state_ = state.copy()
    if type == 'spherical':
        state_  = sph2cart(state_)
    if dim == 2:
        x, y, px, py = state_
    elif dim == 3:
        x, y, z, px, py, pz = state_
    xr = x * xp.cos(angle) + y * xp.sin(angle)
    yr = -x * xp.sin(angle) + y * xp.cos(angle)
    pxr = px * xp.cos(angle) + py * xp.sin(angle)
    pyr = -px * xp.sin(angle) + py * xp.cos(angle)
    if dim == 2:
        state_ = xp.concatenate(([xr], [yr], [pxr], [pyr]), axis=0)
    elif dim == 3:
        state_ = xp.concatenate(([xr], [yr], [z], [pxr], [pyr], [pz]), axis=0)
    return cart2sph(state_) if type == 'spherical' else state_

def field_envelope(t:float, te_au:xp.ndarray, envelope:str='sinus') -> float:
    te = xp.cumsum(te_au)
    if envelope == 'sinus':
        return xp.where(t<=0, 0, xp.where(t<=te[0], xp.sin(xp.pi * t / (2 * te[0]))**2, xp.where(t<=te[1], 1, xp.where(t<=te[2], xp.sin(xp.pi * (te[2] - t) / (2 * te_au[2]))**2, 0))))
    elif envelope == 'const':
        return 1
    elif envelope == 'trapez':
        return xp.where(t<=0, 0, xp.where(t<=te[0], t / te[0], xp.where(t<=te[1], 1, xp.where(t<=te[2], (te[2] - t) / te_au[2], 0))))
    else:
	    raise NameError(f'{envelope} envelope not defined')
	
METHODS = ['Verlet', 'FR', 'Yo# with # any integer', 'Yos6', 'M2', 'M4', 'EFRL', 'PEFRL', 'VEFRL', 'BM4', 'BM6', 'RKN4b', 'RKN6b', 'RKN6a', 'ABA104', 'ABA864', 'ABA1064']

class OdeSolution(OptimizeResult):
    pass

class SymplecticIntegrator:
	"""
    Some explicit symplectic splitting integrators in Python

    Attributes
    ----------
    name : str
        Name of the symplectic integrator. 
		Integration methods are listed in https://pypi.org/project/pyhamsys/ 
    """
	def __repr__(self) -> str:
		return f'{self.__class__.__name__}({self.name})'
	
	def __str__(self) -> str:
		return f'{self.name}'

	def __init__(self, name:str) -> None:
		self.name = name
		if (self.name not in METHODS) and (self.name[:2] != 'Yo'):
			raise ValueError(f"The chosen integrator must be one of {METHODS}.")
		if self.name == 'Verlet':
			alpha_s = [0.5]
		elif self.name == 'FR':
			theta = 1 / (2 - 2**(1/3))
			alpha_s = [theta / 2, theta / 2, 0.5 - theta]
		elif self.name == 'Yos6':
			a = [0.784513610477560, 0.235573213359357, -1.17767998417887, 1.31518632068390]
			alpha_s = [a[0]/2,  a[0]/2, a[1]/2, a[1]/2, a[2]/2, a[2]/2, a[3]/2]
		elif self.name[:2] == 'Yo':
			try:
				alpha_s = xp.asarray([0.5])
				for n in range(1, int(self.name[2:]) // 2):
					x1 = 1 / (2 - 2**(1/(2*n+1)))
					x0 = 1 - 2 * x1
					alpha_ = xp.concatenate((alpha_s, xp.flip(alpha_s)))
					alpha_ = xp.concatenate((x1 * alpha_, x0 * alpha_s))
					alpha_s = alpha_.copy()
			except:
				raise NameError(f'{self.name} integrator not defined') 
		elif self.name.endswith('EFRL'):
			if self.name.startswith('V'):
				xi, lam, chi = 0.1644986515575760, -0.02094333910398989, 1.235692651138917
			elif self.name.startswith('P'):
				xi, lam, chi = 0.1786178958448091, -0.2123418310626054, -0.06626458266981849
			else:
				xi, lam, chi = 0.1720865590295143, -0.09156203075515678, -0.1616217622107222
			alpha_s = [xi, 0.5 - lam - xi, lam + xi + chi -0.5, 0.5 - chi - xi]
		elif self.name == 'M2':
			y = (2*xp.sqrt(326)-36)**(1/3)
			z = (y**2 + 6 * y -2) / (12 * y)
			alpha_s = [z, 0.5 - z]
		elif self.name == 'M4':
			alpha_s = [(14-xp.sqrt(19))/108, (146+5*xp.sqrt(19))/540, (-23-20*xp.sqrt(19))/270, (-2+10*xp.sqrt(19))/135, 1/5]
		elif self.name == 'BM4':
			alpha_s = [0.0792036964311957, 0.1303114101821663, 0.2228614958676077, -0.3667132690474257, 0.3246481886897062, 0.1096884778767498]
		elif self.name == 'BM6':
			alpha_s = [0.050262764400392, 0.098553683500650, 0.314960616927694, -0.447346482695478, 0.492426372489876, -0.425118767797691, 0.237063913978122, 0.195602488600053, 0.346358189850727, -0.362762779254345]
		elif self.name == 'RKN4b':
			alpha_s = [0.082984406417405, 0.162314550766866, 0.233995250731502, 0.370877414979578, -0.409933719901926, 0.059762097006575]
		elif self.name == 'RKN6b':
			alpha_s = [0.041464998518262, 0.081764777428009, 0.116363894490058, 0.174189903309500, -0.214196095413653, 0.087146882788236, -0.011892898486655, -0.234438862575420, 0.222927475154732, 0.134281397641196, 0.102388527145735]
		elif self.name == 'RKN6a':
			alpha_s = [0.0378593198406116,  0.053859832783850, 0.048775800318585, 0.135207369686421, -0.161075257952980, 0.104540892120091, 0.209700510951356, -0.204785822176643, 0.074641362659228, 0.069119764509130, 0.037297935860413, 0.291269757886391, -0.300064001014902, 0.103652534528448]
		elif self.name == 'ABA104':
			alpha_s = [0.04706710064597251, 0.07181481672222451, 0.1129421186948636, 0.128108341856638, 0.1545976638231982, -0.4278843305285221, 0.4133542887856252]
		elif self.name == 'ABA864':
			alpha_s = [0.07113342649822312, 0.1119502609739741, 0.129203166982666, 0.1815796929159088, 0.3398320688569059, -0.3663966873688647, 0.03269807114118675]
		elif self.name == 'ABA1064':
			alpha_s = [0.03809449742241219, 0.05776438341466301, 0.08753433270225074, 0.116911820440748, 0.0907158752847932, 0.1263544726941979, 0.3095552309573282, -0.3269306129163933]
		elif self.name != 'Simple':
			raise NameError(f'{self.name} integrator not defined')
		if self.name == 'Simple':
			self.alpha_s = [1]
			self.alpha_o = [1]
		else:
			self.alpha_s = xp.concatenate((alpha_s, xp.flip(alpha_s)))
			self.alpha_o = xp.tile([1, 0], len(alpha_s))
	
def solve_ivp_symp(chi:Callable, chi_star:Callable, t_span:tuple, y0:xp.ndarray, step:float, t_eval:Union[list, xp.ndarray]=None,
				   method:str='BM4', command:Callable=None) -> OdeSolution:
	"""
	Solve an initial value problem for a Hamiltonian system using an explicit 
	symplectic splitting scheme (see [1]).

	This function numerically integrates a system of ordinary differential
	equations given an initial value:

	dy / dt = {y, H(t, y)}
	y(t0) = y0

	Here t is a 1-D independent variable (time), y(t) is an N-D vector-valued 
	function (state), and a Hamiltonian H(t, y) and a Poisson bracket {. , .} 
	determine the differential equations. The goal is to find y(t) 
	approximately satisfying the differential equations, given an initial value 
	y(t0)=y0. The Hamiltonian flow is defined by two functions `chi` and 
	`chi_star` (see [2]). 

	Parameters
	----------
	chi : callable
		Function of (h, t, y) returning exp(h X_n)...exp(h X_1) y at time t.
		`chi` must return an array of the same shape as y.
	chi_star : callable 
		Function of (h, t, y) returning exp(h X_1)...exp(h X_n) y at time t.
		`chi_star` must return an array of the same shape as y.
	t_span : 2-member sequence
		Interval of integration (t0, tf). The solver starts with t=t0 and
		integrates until it reaches t=tf. Both t0 and tf must be floats or 
		values interpretable by the float conversion function.	
	y0 : array_like, shape (n,)
		Initial state.
	step : float
		Step size.
	t_eval : array_like or None, optional
		Times at which to store the computed solution, must be sorted and 
		equally spaced, and lie within `t_span`. If None (default), use points 
		selected by the solver.
	method : string, optional
        Integration methods are listed on https://pypi.org/project/pyhamsys/ 
		'BM4' is the default.
	command : function of (t, y) 
		Function to be run at each step size.   

	Returns
	-------
	Bunch object with the following fields defined:
	t : ndarray, shape (n_points,)  
		Time points.
	y : ndarray, shape (n, n_points)  
		Values of the solution at `t`.
	step : step size used in the computation

	References
	----------
		[1] Hairer, Lubich, Wanner, 2003, Geometric Numerical Integration: 
		Structure-Preserving Algorithms for Ordinary Differential Equations
		(Springer)
		[2] McLachlan, R.I, 2022, "Tuning symplectic integrators is easy and 
		worthwhile", Commun. Comput. Phys. 31, 987 (2022)
	"""
	integrator = SymplecticIntegrator(method)
	t0, tf = map(float, t_span)

	if t_eval is None or xp.isclose(t_eval[0], t0):
		t_vec, y_vec = [t0], y0[..., xp.newaxis] 
	else:
		t_vec, y_vec = [], []
	
	if t0 > tf:
		raise ValueError("Values in `t_span` are not properly sorted.")
	if t_eval is not None:
		t_eval = xp.asarray(t_eval)
		if t_eval.ndim != 1:
			raise ValueError("`t_eval` must be 1-dimensional.")
		if xp.any(t_eval < t0) or xp.any(t_eval > tf):
			raise ValueError("Values in `t_eval` are not within `t_span`.")
		if xp.any(xp.diff(t_eval) <= 0):
			raise ValueError("Values in `t_eval` are not properly sorted.")
		if not xp.allclose(xp.diff(t_eval), xp.diff(t_eval)[0]):
			raise ValueError("Values in `t_eval` are not equally spaced.")
		if not xp.isclose(t_eval[0], t0):
			t_eval = xp.insert(t_eval, 0, t0)

	if t_eval is not None:
		step = ((t_eval[1] - t_eval[0])) / xp.ceil((t_eval[1] - t_eval[0]) / step)
		spacing = int(xp.ceil((t_eval[1] - t_eval[0]) / step))
	else:
		step = xp.abs((tf - t0)) / xp.ceil(xp.abs((tf - t0)) / step)
	alpha_s = integrator.alpha_s * step

	def _integrate(t:float, y:xp.ndarray) -> Tuple[float, xp.ndarray]:
		for h, st in zip(alpha_s, integrator.alpha_o):
			y = chi(h, t + h, y) if st==0 else chi_star(h, t, y)
			t += h
		return t, y

	count = 0
	t, y_ = t0, y0.copy()
	while t < tf:
		t, y_ = _integrate(t, y_)
		if t_eval is not None and t > t_eval.max():
			break
		if t_eval is not None:
			count += 1
			if count % spacing == 0:
				count = 0
				t_vec.append(t)
				y_vec = xp.concatenate((y_vec, y_[..., xp.newaxis]), axis=-1)
		else:
			t_vec.append(t)
			y_vec = xp.concatenate((y_vec, y_[..., xp.newaxis]), axis=-1)
		if command is not None:
			command(t, y_)
	t_vec = xp.asarray(t_vec)
	return OdeSolution(t=t_vec, y=y_vec, step=step)

def solve_ivp_sympext(fun:Callable, t_span:tuple, y0:xp.ndarray, step:float, t_eval:Union[list, xp.ndarray]=None, 
					  method:str='BM4', omega:float=10, command:Callable=None, check_trajs:int=None) -> OdeSolution:
	"""
	Solve an initial value problem for a Hamiltonian system using an explicit 
	symplectic approximation obtained by an extension in phase space (see [1]).

	This function numerically integrates a system of ordinary differential 
	equations given an initial value:

	dy / dt = {y, H(t, y)}
	y(t0) = y0

	Here t is a 1-D independent variable (time), y(t) is an N-D vector-valued 
	function (state), and a Hamiltonian H(t, y) and a canonical Poisson bracket
	{. , .} determine the differential equations. The goal is to find y(t) 
	approximately satisfying the differential equations, given an initial value 
	y(t0)=y0.
	
	The state y(t) could be of the form (q(t), p(t)) or (q(t), p(t), k(t)). 
	Here k is a canonically conjugate variable to time, for checking the 
	conservation of energy (optional). 

	Parameters
	----------
	fun : callable
		Right-hand side of the system: the time derivative of the state y at  
		time t. i.e., {y, H(t, y)}. The calling signature is `fun(t, y)`, where  
		`t` is a scalar and `y` is an ndarray with `len(y) = len(y0)`.  
		y = (q, p) or (q, p, k) where k is conjugate to time (if energy needs  
		to be ckecked). In that case, specify the number of trajectories  
		check_trajs.
		`fun` must return an array of the same shape as `y`.   
	t_span : 2-member sequence  
		Interval of integration (t0, tf). The solver starts with t=t0 and  
		integrates until it reaches t=tf. Both t0 and tf must be floats or   
		values interpretable by the float conversion function.	 
	y0 : array_like, shape (n,)
		Initial state.
	step : float
		Step size.
	t_eval : array_like or None, optional
		Times at which to store the computed solution, must be sorted and   
		equally spaced, and lie within `t_span`. If None (default), use points  
		selected by the solver.
	method : string, optional
        Integration methods are listed on https://pypi.org/project/pyhamsys/  
		'BM4' is the default.
	omega : float, optional
		Coupling parameter in the extended phase space (see [1])
	command : function of (t, y) or None, optional
		Function to be run at each step size. 
	check_trajs: None (default) or int, optional
		Number of trajectories.
		If y = (q, p, k), the number of trajectories needs to be specified. 
		This is used to check the conservation of energy.

	Returns
	-------
	Bunch object with the following fields defined:
	t : ndarray, shape (n_points,)  
		Time points.
	y : ndarray, shape (n, n_points)  
		Values of the solution at `t`.
	step : step size used in the computation

	References
	----------
		[1] Tao, M., 2016, "Explicit symplectic approximation of nonseparable 
		Hamiltonians: Algorithm and long time performance", 
		Phys. Rev. E 94, 043303
	"""
	if len(y0) // check_trajs % 2 == 0:
		check_trajs = None
	if check_trajs is not None: 
		ny, ce0, ce1 = len(y0), len(y0)-check_trajs, 2*len(y0)-check_trajs
		slices = xp.r_[0:ce0, ny:ce1]
	
	def _coupling(h:float) -> xp.ndarray:
		return (xp.array([[1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]])\
			  + xp.cos(2 * omega * h) * xp.array([[1, 0, -1, 0], [0, 1, 0, -1], [-1, 0, 1, 0], [0, -1, 0, 1]])\
			  + xp.sin(2 * omega * h) * xp.array([[0, -1, 0, 1], [1, 0, -1, 0], [0, 1, 0, -1], [-1, 0, 1, 0]])) / 2

	def _chi_ext(h:float, t:float, y:xp.ndarray) -> xp.ndarray:
		y_ = xp.split(y, 2)
		y_[0] += h * fun(t, y_[1])
		y_[1] += h * fun(t, y_[0])
		y_ = xp.concatenate((y_[0], y_[1]), axis=None)
		if check_trajs is None:
			return xp.einsum('ij,j...->i...', _coupling(h), xp.split(y_, 4)).flatten()
		yr = xp.split(xp.einsum('ij,j...->i...', _coupling(h), xp.split(y_[slices], 4)).flatten(), 2)
		return xp.concatenate((yr[0], y_[ce0:ny], yr[1], y_[ce1:]), axis=None) 
		
	def _chi_ext_star(h:float, t:float, y:xp.ndarray) -> xp.ndarray:
		if check_trajs is None:
			yr = xp.einsum('ij,j...->i...', _coupling(h), xp.split(y, 4)).flatten()
		else:
			yr = xp.split(xp.einsum('ij,j...->i...', _coupling(h), xp.split(y[slices], 4)).flatten(), 2)
			y_ = xp.concatenate((yr[0], y[ce0:ny], yr[1], y[ce1:]), axis=None) 
		y_ = xp.split(y_, 2)
		y_[1] += h * fun(t, y_[0])
		y_[0] += h * fun(t, y_[1])
		return xp.concatenate((y_[0], y_[1]), axis=None)
	
	y_ = xp.tile(y0, 2)
	sol = solve_ivp_symp(_chi_ext, _chi_ext_star, t_span, y_, method=method, step=step, t_eval=t_eval, command=command)
	if check_trajs is None:
		y_ = xp.split(sol.y, 4, axis=0)
		sol.y = xp.concatenate(((y_[0] + y_[2]) / 2, (y_[1] + y_[3]) / 2), axis=0)
	else:
		y_ = xp.split(sol.y[slices], 4, axis=0)
		sol.y = xp.concatenate(((y_[0] + y_[2]) / 2, (y_[1] + y_[3]) / 2, (sol.y[ce0:ny] + sol.y[ce1:]) / 2), axis=0)
	return sol
