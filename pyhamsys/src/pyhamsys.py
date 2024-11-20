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
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
from sklearn.metrics import r2_score
from scipy.fft import rfft, irfft, rfftfreq
from typing import Callable, Union, Tuple
from scipy.optimize import OptimizeResult
import sympy as sp
from functools import partial 

class OdeSolution(OptimizeResult):
    pass

class HamSys:
	def __init__(self, ndof:float=1, complex:bool=False) -> None:
		if str(ndof) != str(int(ndof)) + '.5' * bool(str(ndof).count('.5')):
			raise ValueError('Number of degrees of freedom should be an integer or half an integer.')
		self._ndof = int(ndof)
		self._time_dependent = bool(str(ndof).count('.5'))
		self._complex = complex

	def _split(self, y:xp.ndarray, n:int, check_energy:bool=False):
		yr = y[:-1] if check_energy else y
		ys = xp.split(yr, n)
		if not check_energy:
			return ys
		return ys.append(y[-1])
	
	def _create_function(self, t:float, y:xp.ndarray, eqn:Callable) -> xp.ndarray:
		y_ = xp.split(y, 2)
		return xp.asarray(eqn(y_[0], y_[1], t)).flatten()
	
	def rectify_sol(self, sol:OdeSolution, check_energy:bool=False) -> OdeSolution:
		if not check_energy:
			return sol
		if self._time_dependent:
			sol.y, sol.k = sol.y[:-1, :], sol.y[-1, :]
		sol.err = self.compute_energy(sol)
		return sol
	
	def compute_vector_field(self, hamiltonian:Callable, output:bool=False) -> None:
		if self.complex:
			raise ValueError("Computation of vector fields not implemented for complex Hamiltonians")
		q = sp.symbols('q0:%d'%self._ndof) if self._ndof>=2 else sp.Symbol('q')
		p = sp.symbols('p0:%d'%self._ndof) if self._ndof>=2 else sp.Symbol('p')
		t = sp.Symbol('t')
		energy = sp.lambdify([q, p, t], hamiltonian(q, p, t))
		self.hamiltonian = partial(self._create_function, eqn=energy)
		eqn = sp.simplify(sp.derive_by_array(hamiltonian(q, p, t), [q, p]).doit())
		eqn = sp.flatten([eqn[1], -eqn[0]])
		if output:
			print('y_dot : ', eqn)
		eqn = sp.lambdify([q, p, t], eqn)
		self.y_dot = partial(self._create_function, eqn=eqn)
		eqn_t = -sp.simplify(sp.diff(hamiltonian(q, p, t), t))
		if output and eqn_t!=0:
			print('k_dot : ', eqn_t)
		eqn_t = sp.lambdify([q, p, t], eqn_t)
		self.k_dot = partial(self._create_function, eqn=eqn_t)

	def compute_energy(self, sol:OdeSolution, maxerror:bool=True) -> xp.ndarray:
		if not hasattr(self, 'hamiltonian'):
			raise ValueError("In order to check energy, the attribute 'hamiltonian' must be provided.")
		val_h = xp.empty_like(sol.t)
		for _, t in enumerate(sol.t):
			val_h[_] = self.hamiltonian(t, sol.y[:, _])
		if self._time_dependent:
			val_h += sol.k
		return xp.max(xp.abs(val_h - val_h[0])) if maxerror else val_h

def compute_msd(sol:OdeSolution, plot_data:bool=False, output_r2:bool=False):
	x, y = xp.split(sol.y, 2)
	nt = len(sol.t)
	r2 = xp.zeros(nt)
	for _ in range(nt):
		r2[_] = ((x[:, _:] - x[:, :-_ if _ else None])**2 + (y[:, _:] - y[:, :-_ if _ else None])**2).mean()
	t_win, r2_win = sol.t[nt//8:7*nt//8], r2[nt//8:7*nt//8]
	res = linregress(t_win, r2_win)
	diff_data = [res.slope, res.intercept, res.rvalue**2]
	func_fit = lambda t, a, b: (a * t)**b
	popt = curve_fit(func_fit, t_win, r2_win, bounds=((0, 0.25), (xp.inf, 3)))[0]
	r2_fit = func_fit(t_win, *popt)
	interp_data = [*popt, r2_score(r2_win, r2_fit)]
	if plot_data:
		plt.plot(sol.t, r2, ':', color='r', lw=1)
		plt.plot(t_win, r2_win, '-', color='r', lw=2)
		plt.plot(t_win, r2_fit, '-.', color='r', lw=2)
		plt.xlabel('$t$')
		plt.ylabel('$r^2$')
		plt.show()
	if output_r2:
		return sol.t, r2, diff_data, interp_data
	return diff_data, interp_data

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
		if t_eval is not None and t >= t_eval.max() + step / 2:
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

def solve_ivp_sympext(hs:HamSys, t_span:tuple, y0:xp.ndarray, step:float, t_eval:Union[list, xp.ndarray]=None, method:str='BM4', omega:float=10, command:Callable=None, check_energy:bool=False) -> OdeSolution:
	"""
	Solve an initial value problem for a Hamiltonian system using an explicit 
	symplectic approximation obtained by an extension in phase space (see [1]).

	This function numerically integrates a system of ordinary differential 
	equations in canonical coordinates given an initial value:

	dy / dt = {y, H(t, y)}
	y(t0) = y0

	Here t is a 1-D independent variable (time), y(t) is an N-D vector-valued 
	function (state), and a Hamiltonian H(t, y) and a canonical Poisson bracket
	{. , .} determine the differential equations. The goal is to find y(t) 
	approximately satisfying the differential equations, given an initial value 
	y(t0)=y0. The state y(t) is of the form (q(t), p(t)) for a real Hamiltonian 
	system, and (q(t) + i p(t)) / sqrt(2) for complex Hamiltonian systems. 
	
	If the Hamiltonian has an explicit time dependence, there is the 
	possibility to check energy by computing k(t) where k is a canonically 
	conjugate variable to time. Its evolution is given by

	dk / dt = -dH / dt
	k(t0)=0

	Parameters
	----------
	hs : HamSys
		Hamiltonian system containing the Hamiltonian vector field. The 
		attributes `y_dot` (for dy / dt) should be specified. If there is an 
		explicit time dependence and `check_energy` is True, the attribute 
		`k_dot` (for dk / dt) should be specified. 
	t_span : 2-member sequence  
		Interval of integration (t0, tf). The solver starts with t=t0 and  
		integrates until it reaches t=tf. Both t0 and tf must be floats or   
		values interpretable by the float conversion function.	 
	y0 : array_like, shape (2n,)
		Initial state y0. If hs is complex y0 = (q0 + i p0) / sqrt(2) where q0 are the initial positions and p0 the initial momenta. Otherwise y0 = (q0, p0). 
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
	check_energy : bool, optional
		If True, computes the total energy. Default is False.  	

	Returns
	-------
	Bunch object with the following fields defined:
	t : ndarray, shape (n_points,)  
		Time points.
	y : ndarray, shape (2n, n_points) real array or  (n, n_points) complex array
		If hs real, y(t) = (q(t), p(t)) at `t`. If hs complex, y(t) = (q(t) + i p(t)) / sqrt(2)
	k : ndarray, shape (n_points,)
		Values of k(t) at `t` if `check_energy` is True and if the Hamiltonian
		system has an explicit time dependence.   
	err : float
		Error in the computation of the total energy, computed only if
		`check_energy` is True. 
	step : float
		Step size used in the computation.

	References
	----------
		[1] Tao, M., 2016, "Explicit symplectic approximation of nonseparable 
		Hamiltonians: Algorithm and long time performance", 
		Phys. Rev. E 94, 043303
	"""
	check_energy_ = check_energy * hs._time_dependent
	
	def _coupling(h:float) -> xp.ndarray:
		if hs._complex:
			return 	(xp.array([[1, 1], [1, 1]]) + xp.exp(2j * omega * h) * xp.array([[1, -1], [-1, 1]])) / 2
		return (xp.array([[1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]])\
			  + xp.cos(2 * omega * h) * xp.array([[1, 0, -1, 0], [0, 1, 0, -1], [-1, 0, 1, 0], [0, -1, 0, 1]])\
			  + xp.sin(2 * omega * h) * xp.array([[0, -1, 0, 1], [1, 0, -1, 0], [0, 1, 0, -1], [-1, 0, 1, 0]])) / 2

	def _chi_ext(h:float, t:float, y:xp.ndarray) -> xp.ndarray:
		y_ = hs._split(y, 2, check_energy=check_energy_)
		y_[1] += h * hs.y_dot(t, y_[0])
		if check_energy_:
			y_[-1] += h * hs.k_dot(t, y_[0])
		y_[0] += h * hs.y_dot(t, y_[1])
		if check_energy_:
			y_[-1] += h * hs.k_dot(t, y_[1])
		yr = xp.concatenate((y_[0], y_[1]), axis=None)
		yr = xp.einsum('ij,j...->i...', _coupling(h), hs._split(yr, 2 if hs._complex else 4)).flatten()
		if not check_energy_:
			return yr
		return xp.concatenate((yr, y_[-1]), axis=None) 
		
	def _chi_ext_star(h:float, t:float, y:xp.ndarray) -> xp.ndarray:
		y_ = hs._split(y, 2 if hs._complex else 4, check_energy=check_energy_)
		yr = y_ if not check_energy_ else y_[:-1]
		yr = xp.einsum('ij,j...->i...', _coupling(h), yr).flatten()
		if check_energy_:
			yr = xp.concatenate((yr, y_[-1]), axis=None) 
		y_ = hs._split(yr, 2, check_energy=check_energy_)
		y_[0] += h * hs.y_dot(t, y_[1])
		if check_energy_:
			y_[-1] += h * hs.k_dot(t, y_[1])
		y_[1] += h * hs.y_dot(t, y_[0])
		if check_energy_:
			y_[-1] += h * hs.k_dot(t, y_[0])
		return xp.concatenate([_ for _ in y_], axis=None)
	
	if not hasattr(hs, 'y_dot'):
		raise ValueError("The attribute 'y_dot' must be provided.")
	if check_energy_ and not hasattr(hs, 'k_dot'):
		raise ValueError("In order to check energy for a time-dependent system, the attribute 'k_dot' must be provided.")
	y_ = xp.tile(y0, 2)
	if check_energy_:
		y_ = xp.pad(y_, (0, 1), 'constant')
	sol = solve_ivp_symp(_chi_ext, _chi_ext_star, t_span, y_, method=method, step=step, t_eval=t_eval, command=command)
	y_ = hs._split(sol.y, 2 if hs._complex else 4, check_energy=check_energy_)
	if not hs._complex:
		sol.y = xp.concatenate(((y_[0] + y_[2]) / 2, (y_[1] + y_[3]) / 2), axis=0)
	else:
		sol.y = (y_[0] + y_[1]) / 2
	if check_energy_:
		sol.k = y_[-1] / 2
	if check_energy:
		sol.err = hs.compute_energy(sol)
	return sol
