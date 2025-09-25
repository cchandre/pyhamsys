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
from scipy.integrate import solve_ivp
from scipy.integrate._ivp.ivp import METHODS as IVP_METHODS
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
from sklearn.metrics import r2_score
from scipy.fft import rfft, irfft, rfftfreq
from typing import Callable, Union, Tuple
from scipy.optimize import OptimizeResult
from scipy.io import savemat
import sympy as sp
from functools import partial 
import time
from datetime import datetime

METHODS = ['Verlet', 'FR', 'Yo# with # any integer', 'Yos6', 'M2', 'M4', 'EFRL', 'PEFRL', 'VEFRL', 'BM4', 'BM6', 'RKN4b', 'RKN6b', 'RKN6a', 'ABA104', 'ABA864', 'ABA1064']

IVP_METHODS = list(IVP_METHODS.keys())
ALL_METHODS = METHODS + IVP_METHODS

class OdeSolution(OptimizeResult):
    pass

class HamSys:
	def __init__(self, ndof:float=1, btype:str='pq', y_dot:Callable=None,\
			   k_dot:Callable=None, hamiltonian:Callable=None) -> None:
		if str(ndof) != str(int(ndof)) + '.5' * bool(str(ndof).count('.5')):
			raise ValueError('Number of degrees of freedom should be an integer or half an integer.')
		self._ndof = int(ndof)
		self._time_dependent = bool(str(ndof).count('.5'))
		self.btype = btype
		self._ysplit = 4 if btype=='pq' else 2
		self._complex = True if btype=='psi' else False
		if y_dot is not None:
			self.y_dot = y_dot
		if k_dot is not None:
			self.k_dot = k_dot
		if hamiltonian is not None:
			self.hamiltonian = hamiltonian

	def _split(self, y:xp.ndarray, n:int, check_energy:bool=False):
		ys = xp.split(y[:-1] if check_energy else y, n)
		if check_energy:
			ys.append(y[-1])
		return ys
	
	def _create_function(self, t:float, y:xp.ndarray, eqn:Callable) -> xp.ndarray:
		q, p = xp.split(y, 2)
		return xp.asarray(eqn(q, p, t)).flatten()
	
	def rectify_sol(self, sol:OdeSolution, check_energy:bool=False) -> OdeSolution:
		if not check_energy:
			return sol
		if self._time_dependent:
			sol.y, sol.k = sol.y[:-1, :], sol.y[-1, :]
		if not hasattr(self, 'hamiltonian'):
			raise ValueError("In order to check energy, the attribute 'hamiltonian' must be provided.")
		sol.err = self.compute_energy(sol)
		return sol
	
	def compute_vector_field(self, hamiltonian:Callable, output:bool=False, check_energy:bool=False) -> None:
		if self.btype != 'pq':
			raise ValueError("Computation of vector fields not implemented for these variables")
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
		if self._time_dependent and check_energy:
			eqn_t = -sp.simplify(sp.diff(hamiltonian(q, p, t), t))
			if output and eqn_t!=0:
				print('k_dot : ', eqn_t)
			eqn_t = sp.lambdify([q, p, t], eqn_t)
			self.k_dot = partial(self._create_function, eqn=eqn_t)

	def compute_energy(self, sol:OdeSolution, maxerror:bool=True) -> xp.ndarray:
		val_h = xp.empty_like(sol.t)
		for _, t in enumerate(sol.t):
			val_h[_] = self.hamiltonian(t, sol.y[:, _])
		if self._time_dependent:
			val_h += sol.k
		return xp.max(xp.abs(val_h - val_h[0])) if maxerror else val_h
	
	def _y_dot_ext(self, t, z):
		return xp.concatenate((self.y_dot(t, z[:-1]), self.k_dot(t, z[:-1])), axis=None)
	
	def integrate(self, z0, t_eval, timestep, solver="BM4", extension=False, check_energy=False, omega=10, diss=0, tol=1e-8, display=True):
		"""
		Integrate the system using either an IVP solver or a symplectic solver.

		Parameters
		----------
		z0 : array_like
			Initial condition(s).
		t_eval : array_like
			Times at which to store the solution.
		timestep : float
			Integration time step.
		solver : str, optional
			Solver method. Must be in METHODS or IVP_METHODS.
		extension : bool, optional
			Whether to use symplectic extension.
		check_energy : bool, optional
			If True, adds an auxiliary variable to check energy conservation.
		omega : float, optional
			Frequency parameter for symplectic extension solvers.
		tol : float, optional
			Absolute and relative tolerance for IVP solvers.
		display : bool, optional
			Whether to print runtime information.

		Returns
		-------
		sol : object
			Solution object with attributes depending on solver used.
		"""
		
		if solver not in ALL_METHODS:
			raise ValueError(f"Solver '{solver}' not recognized. "
                 f"Valid solvers are {ALL_METHODS}.")
		if solver in IVP_METHODS or (solver in METHODS and extension):
			if not hasattr(self, 'y_dot'):
				raise ValueError("In order to use an IVP solver or an extension in phase space, 'y_dot' must be provided.")
		if solver in METHODS and not extension:
			if not hasattr(self, 'chi') or not hasattr(self, 'chi_star'):
				raise ValueError("In order to use a symplectic integrator, 'chi' and 'chi_star' must be provided.")
		if solver in IVP_METHODS and check_energy:
			if not hasattr(self, 'k_dot'):
				raise ValueError("In order to check energy with an IVP solver, 'k_dot' must be provided.")
			z0 = xp.concatenate([z0, xp.zeros(1, dtype=z0.dtype)])
		start = time.process_time()
		if solver in IVP_METHODS:
			sol = solve_ivp(self._y_dot_ext if check_energy else self.y_dot, (t_eval[0], t_eval[-1]), z0, t_eval=t_eval, method=solver, atol=tol, rtol=tol, max_step=timestep)
			sol = self.rectify_sol(sol, check_energy=check_energy)
			sol.step = timestep
		elif extension:
			sol = solve_ivp_sympext(self, (t_eval[0], t_eval[-1]), z0, step=timestep, t_eval=t_eval, method=solver, check_energy=check_energy, omega=omega, diss=diss)
		else:
			sol = solve_ivp_symp(self.chi, self.chi_star, (t_eval[0], t_eval[-1]), z0, step=timestep, t_eval=t_eval, method=solver)
		sol.cpu_time = time.process_time() - start
		if display:
			print(f'\033[90m        Computation finished in {int(sol.cpu_time)} seconds \033[00m')
			if hasattr(sol, 'err'):
				print(f'\033[90m           with error in energy = {sol.err:.2e} \033[00m')
			if hasattr(sol, 'dist_copy'):
				print(f'\033[90m           with distance in copies = {sol.dist_copy:.2e}\033[00m')
		return sol
	
	def compute_lyapunov(self, tf, z0, reortho_dt, tol=1e-8, solver='RK45', display=True):
		if solver not in IVP_METHODS:
			raise ValueError(f"Solver {solver} is not recognized for Lyapunov exponent computation."
							 f"Available solvers are {IVP_METHODS}.")
		if not hasattr(self, 'y_dot_lyap'):
			raise ValueError("In order to compute the Lyapunov spectrum, 'y_dot_lyap' must be provided.")
		start = time.time()
		n = len(z0) // 2
		lyap_sum = xp.zeros((2, n), dtype=xp.float64)
		t, z = 0, xp.concatenate((z0, xp.ones(n), xp.zeros(n), xp.zeros(n), xp.ones(n)), axis=None)
		for _ in range(int(tf / reortho_dt)):
			sol = solve_ivp(self.y_dot_lyap, (t, t + reortho_dt), z, method=solver, t_eval=[t + reortho_dt], atol=tol, rtol=tol)
			z, Q = sol.y[:2 * n, -1], xp.moveaxis(sol.y[2 * n:, -1].reshape((2, 2, n)), -1, 0)
			for i in range(n):
				q, r = xp.linalg.qr(Q[i])
				lyap_sum[:, i] += xp.log(xp.abs(xp.diag(r)))
				Q[i] = q
			t += reortho_dt
			z = xp.concatenate((z, xp.moveaxis(Q, 0, -1)), axis=None)
		if display:
			print(f'\033[90m        Computation finished in {int(time.time() - start)} seconds \033[00m')
		return xp.sort(lyap_sum / tf)
	
	def save_data(self, data, params, filename=None, author=None, display=True):
		params.update({'data': data})
		params.update({'date': datetime.now().strftime(" %B %d, %Y\n"), 'author': author})
		filename += '_' + datetime.now().strftime("%Y%m%d_%H%M%S") + '.mat'
		savemat(filename, params)
		if display:
			print(f'\033[90m        Results saved in {filename} \033[00m')

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

def antiderivative(vec:xp.ndarray, N:int=2**10, axis=-1) -> xp.ndarray:
    nu = rfftfreq(N, d=1/N) 
    div = xp.zeros_like(nu, dtype=complex)
    div[1:] = 1 / (1j * nu[1:]) 
    shape = [1] * vec.ndim
    shape[axis] = -1
    div = div.reshape(shape)
    return irfft(rfft(vec, axis=axis) * div, axis=axis)

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
	y0 : array_like
		Initial state.
	step : float
		Step size.
	t_eval : array_like or None, optional
		Times at which to store the computed solution, must be sorted and, 
		lie within `t_span`. If None (default), use points selected by the 
		solver.
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
		n_eval = len(t_eval) - 1
	else:
		n_eval = int((tf - t0) / step) + 1
	nstep = (int(xp.ceil((tf - t0) / step)) // n_eval) * n_eval + n_eval
	step = (tf - t0) / nstep
	alpha_s = integrator.alpha_s * step

	times = xp.linspace(t0, tf, nstep + 1)
	if t_eval is None:
		t_eval = times.copy()
	t_vec = xp.empty_like(t_eval)
	y_vec = xp.empty(y0.shape + t_vec.shape, dtype=y0.dtype)
	y_vec[:] = xp.nan

	def _integrate(t:float, y:xp.ndarray) -> Tuple[float, xp.ndarray]:
		for h, st in zip(alpha_s, integrator.alpha_o):
			y = chi(h, t + h, y) if st==0 else chi_star(h, t, y)
			t += h
		return t, y
	
	count, y_ = 0, y0.copy()
	for _, t in enumerate(times):
		if (count <= len(t_eval) - 1) and (xp.abs(times - t_eval[count]).argmin() == _):
			t_vec[count] = t
			y_vec[..., count] = y_
			count += 1
		if command is not None:
			command(t, y_)
		if t != times[-1]:
			y_ = _integrate(t, y_)[1]
	return OdeSolution(t=t_vec, y=y_vec, step=step)

def solve_ivp_sympext(hs:HamSys, t_span:tuple, y0:xp.ndarray, step:float, t_eval:Union[list, xp.ndarray]=None, method:str='BM4', omega:float=10, diss:float=0, command:Callable=None, check_energy:bool=False) -> OdeSolution:
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
		Initial state y0. If hs is complex y0 = (q0 + i p0) / sqrt(2) where q0
		are the initial positions and p0 the initial momenta. 
		Otherwise y0 = (q0, p0). 
	step : float
		Step size.
	t_eval : array_like or None, optional
		Times at which to store the computed solution, must be sorted, and lie 
		within `t_span`. If None (default), use points selected by the solver.
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
		If hs real, y(t) = (q(t), p(t)) at `t`. 
		If hs complex, y(t) = (q(t) + i p(t)) / sqrt(2)
	k : ndarray, shape (n_points,)
		Values of k(t) at `t` if `check_energy` is True and if the Hamiltonian
		system has an explicit time dependence.  
	dist_copy : float
		Maximum distance between the two copies of the state in the extended 
		phase space.
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

	if hs.btype not in ['pq', 'psi'] and hasattr(hs, 'coupling')==False:
		raise ValueError("The attribute 'coupling' should be defined")

	J20, J22 = xp.array([[1, 1], [1, 1]]) / 2, xp.array([[1, -1], [-1, 1]]) / 2
	J40, J42c, J42s = xp.array([[1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]) / 2, \
		xp.array([[1, 0, -1, 0], [0, 1, 0, -1], [-1, 0, 1, 0], [0, -1, 0, 1]]) / 2, \
		xp.array([[0, -1, 0, 1], [1, 0, -1, 0], [0, 1, 0, -1], [-1, 0, 1, 0]]) / 2
	
	def _coupling(h:float, y:xp.ndarray) -> xp.ndarray:
		if hasattr(hs, 'coupling'):
			return hs.coupling(h, y, omega).flatten()
		if hs.btype == 'psi':
			J = J20 + xp.exp(2j * omega * h) * J22
		elif hs.btype == 'pq': 
			J = J40 + xp.cos(2 * omega * h) * J42c + xp.sin(2 * omega * h) * J42s
		return xp.einsum('ij,j...->i...', J, xp.split(y, hs._ysplit)).flatten()
	
	def _dissipate(h:float, y:xp.ndarray) -> xp.ndarray:
		u = xp.exp(-h)
		up, um = (1 + u) / 2, (1 - u) / 2
		z1, z2 = xp.split(y, 2)
		return xp.concatenate((up * z1 + um * z2, um * z1 + up * z2), axis=None)
	
	def _command(t:float, y:xp.ndarray):
		return command(t, hs._split(y, 2, check_energy=check_energy_)[0])

	def _chi_ext(h:float, t:float, y:xp.ndarray) -> xp.ndarray:
		y_ = hs._split(y, 2, check_energy=check_energy_)
		y_[1] += h * hs.y_dot(t, y_[0])
		if check_energy_:
			y_[-1] += h * hs.k_dot(t, y_[0])
		y_[0] += h * hs.y_dot(t, y_[1])
		if check_energy_:
			y_[-1] += h * hs.k_dot(t, y_[1])
		yr = xp.concatenate((y_[0], y_[1]), axis=None)
		yr = _coupling(h, yr)
		if diss != 0:
			yr = _dissipate(h * diss, yr)
		if not check_energy_:
			return yr
		return xp.concatenate((yr, y_[-1]), axis=None) 
		
	def _chi_ext_star(h:float, t:float, y:xp.ndarray) -> xp.ndarray:
		yr = y if not check_energy_ else y[:-1]
		if diss != 0:
			yr = _dissipate(h * diss, yr)
		y_ = xp.split(_coupling(h, yr), 2)
		y_[0] += h * hs.y_dot(t, y_[1])
		if check_energy_:
			y[-1] += h * hs.k_dot(t, y_[1])
		y_[1] += h * hs.y_dot(t, y_[0])
		if check_energy_:
			y[-1] += h * hs.k_dot(t, y_[0])
		if not check_energy_:
			return xp.concatenate((y_[0], y_[1]), axis=None)
		return xp.concatenate((y_[0], y_[1], y[-1]), axis=None)
	
	if not hasattr(hs, 'y_dot'):
		raise ValueError("The attribute 'y_dot' must be provided.")
	if check_energy_ and not hasattr(hs, 'k_dot'):
		raise ValueError("In order to check energy for a time-dependent system, the attribute 'k_dot' must be provided.")
	if check_energy and not hasattr(hs, 'hamiltonian'):
		raise ValueError("In order to check energy, the attribute 'hamiltonian' must be provided.")
	
	y_ = xp.tile(y0, 2).astype(xp.complex128 if hs._complex else xp.float64)
	if check_energy_:
		y_ = xp.append(y_, 0)
	sol = solve_ivp_symp(_chi_ext, _chi_ext_star, t_span, y_, method=method, step=step, t_eval=t_eval, command=_command if command is not None else None)
	y_ = hs._split(sol.y, hs._ysplit, check_energy=check_energy_)
	if hs._ysplit == 2:
		sol.y = (y_[0] + y_[1]) / 2
		sol.dist_copy = xp.amax(xp.abs(y_[0] - y_[1]))
	else:
		sol.y = xp.concatenate(((y_[0] + y_[2]) / 2, (y_[1] + y_[3]) / 2), axis=0)
		sol.dist_copy = max(xp.amax(xp.abs(y_[0] - y_[2])), xp.amax(xp.abs(y_[1] - y_[3])))
	if check_energy_:
		sol.k = y_[-1] / 2
	if check_energy:
		sol.err = hs.compute_energy(sol)
	return sol

