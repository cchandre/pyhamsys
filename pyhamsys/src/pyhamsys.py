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
from scipy.interpolate import interp1d
from typing import Callable, Union
from scipy.optimize import OptimizeResult
from inspect import signature
import warnings

warnings.simplefilter('once', UserWarning)

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
    Some symplectic splitting integrators in Python
    .. version:: 0.2.0

    Attributes
    ----------
    name : str
        the name of the symplectic integrator
    step : float
        the time step for the integrator

    Methods
    ------- 
	integrate : integrate the flow from initial to final time 
    """
	def __repr__(self) -> str:
		return f'{self.__class__.__name__}({self.name}, {self.step})'
	
	def __str__(self) -> str:
		return f'{self.name}'

	def __init__(self, name:str, step:float) -> None:
		self.name = name
		self.step = step
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

	def _integrate(self, chi_:Callable, chi_star_:Callable, t, y):
		for h, st in zip(self.alpha_s_, self.alpha_o):
			t += h
			y = chi_(h, t, y) if st==0 else chi_star_(h, t, y)
		return t, y
	
	def integrate(self, chi:Callable, chi_star:Callable, y:xp.ndarray, t_span:tuple, t_eval:Union[int, float, list, xp.ndarray]=None, command:Callable=None) -> OdeSolution:
		"""
		Integrate the Hamiltonian flow from the initial conditions 
		specified by the initial state vector y using one of the selected 
		symplectic splitting integrators.
		Returns the value of y at times defines by the integer, float, list 
		or numpy array times.

		Parameters
		----------
		chi : function of (h, t, y) or (h, y), y being the state vector
			function returning exp(h X_n)...exp(h X_1) y.
		chi_star : function of (h, t, y) or (h, y)
			function returning exp(h X_1)...exp(h X_n) y.
		y : initial state vector (numpy array)
		t_span : tuple of floats or integers; (initial time, final time)
		t_eval : array_like or None, optional  
			times at which the values of the state vector are computed, must  
			be sorted and lie within t_span. If None (default), use points  
			selected by the solver.
		command : function of (t, y) 
			function to be run at each time step.   

		Returns
		-------
		Bunch object with the following fields defined:
		t : ndarray, shape (n_points,)  
        	Time points.
		y : ndarray, shape (n, n_points)  
        	Values of the solution at `t`.
		time_step : time step used in the computation

		References
		----------
			Hairer, Lubich, Wanner, 2003, Geometric Numerical Integration: 
			Structure-Preserving Algorithms for Ordinary Differential Equations
			(Springer)
			McLachlan, R.I, 2022, "Tuning symplectic integrators is easy and 
			worthwhile", Commun. Comput. Phys. 31, 987 (2022)
		"""
		number_vars = len(signature(chi).parameters) - 1
		if number_vars==2:
			chi_, chi_star_ = chi, chi_star
		elif number_vars==1:
			chi_ = lambda h, t, y: chi(h, y)
			chi_star_ = lambda h, t, y: chi_star(h, y)
		t, y_ = t_span[0], y.copy()
		if t_eval is None or xp.isclose(t_eval[0], t_span[0], rtol=1e-12, atol=1e-12):
			t_vec, y_vec = [t_span[0]], y_[..., xp.newaxis] 
		else:
			t_vec, y_vec = [], []
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
			if not xp.isclose(t_eval[0], t_span[0], rtol=1e-12, atol=1e-12):
				t_eval = xp.insert(t_eval, 0, t_span[0])
		evenly_spaced = True if (t_eval is not None and len(t_eval)>=2 and xp.allclose(xp.diff(t_eval)-xp.diff(t_eval)[0], 0, rtol=1e-12, atol=1e-12)) else False
		timestep = self.step
		if evenly_spaced:
			timestep = ((t_eval[1] - t_eval[0])) / xp.ceil((t_eval[1] - t_eval[0]) / self.step)
			spacing = int(xp.ceil((t_eval[1] - t_eval[0]) / timestep))
		elif t_eval is None:
			timestep = xp.abs((tf - t0)) / xp.ceil(xp.abs((tf - t0)) / self.step)
			spacing = int(xp.ceil(xp.abs((tf - t0)) / timestep))
		if xp.abs(timestep - self.step) >= 1e-12:
			print(f"\033[91m        The time step is redefined: old ({self.step}) -> new ({timestep}) \033[00m")
		self.step = timestep
		self.alpha_s_ = self.alpha_s * self.step
		count = 0
		while t < t_span[-1]:
			t, y_ = self._integrate(chi_, chi_star_, t, y_)
			if evenly_spaced or t_eval is None:
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
		if evenly_spaced or t_eval is None:
			return OdeSolution(t=t_vec, y=y_vec, time_step=self.step)
		else:
			return OdeSolution(t=t_eval, y=interp1d(t_vec, y_vec, assume_sorted=True)(t_eval), time_step=self.step)
		

def solve_ivp_sympext(fun:Callable, t_span:tuple, y0:xp.ndarray, step:float, t_eval:Union[int, float, list, xp.ndarray]=None, 
					  method:str='BM4', omega:float=10, command:Callable=None) -> OdeSolution:
	"""
	Solve an initial value problem for a Hamiltonian system using an 
	explicit symplectic approximation obtained by an extension in 
	phase space (see [1]).

	This function numerically integrates a system of ordinary differential
	equations given an initial value::

	dy / dt = {y, H(t, y)}
	y(t0) = y0

	Here t is a 1-D independent variable (time), y(t) is an
	N-D vector-valued function (state), and a Hamiltonian H(t, y)
	and a Poisson bracket {. , .} determine the differential equations.
	The goal is to find y(t) approximately satisfying the differential
	equations, given an initial value y(t0)=y0.

	Parameters
	----------
	fun : callable
		Right-hand side of the system: the time derivative of the state ``y``
		at time ``t``. i.e., {y, H(t, y)}. The calling signature is 
		``fun(t, y)``, where ``t`` is a scalar and ``y`` is an ndarray with 
		``len(y) = len(y0)``. ``fun`` must return an array of the same shape 
		as ``y``. 
	t_span : 2-member sequence
		Interval of integration (t0, tf). The solver starts with t=t0 and
		integrates until it reaches t=tf. Both t0 and tf must be floats
		or values interpretable by the float conversion function.	
	y0 : array_like, shape (n,)
		Initial state.
	step : float
		Step size.
	t_eval : array_like or None, optional
		Times at which to store the computed solution, must be sorted and lie
		within `t_span`. If None (default), use points selected by the solver.	
	method : string, optional
        Integration methods are listed on https://pypi.org/project/pyhamsys/ 
		'BM4' is the default.
	omega : float, optional
		Coupling parameter in the extended phase space (see [1])
	command : function of (t, y) or None, optional
		function to be run at each time step.   

	Returns
	-------
	Bunch object with the following fields defined:
	t : ndarray, shape (n_points,)  
		Time points.
	y : ndarray, shape (n, n_points)  
		Values of the solution at `t`.
	time_step : time step used in the computation

	References
	----------
		[1] Tao, M., 2016, "Explicit symplectic approximation of nonseparable 
		Hamiltonians: Algorithm and long time performance", Phys. Rev. E 94, 043303
	"""
	integrator = SymplecticIntegrator(method, step)
	rotation_e = lambda h: (xp.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [0, 0, 1, 1]])\
			  + xp.cos(2 * omega * h) * xp.array([[1, -1, 0, 0], [-1, 1, 0, 0], [0, 0, 1, -1], [0, 0, -1, 1]])\
			  + xp.sin(2 * omega * h) * xp.array([[0, 0, -1, 1], [0, 0, 1, -1], [1, -1, 0, 0], [-1, 1, 0, 0]])) / 2

	def chi_ext(h, t, y, fun:Callable):
		y_ = xp.split(y, 2)
		y_[0] += h * fun(t, y_[1])
		y_[1] += h * fun(t, y_[0])
		y_ = xp.concatenate((y_[0], y_[1]), axis=None)
		y_ = xp.einsum('ij,j...->i...', rotation_e(h), xp.split(y_, 4)).flatten()
		return t + h, y_
		
	def chi_ext_star(h, t, y, fun:Callable):
		t += h
		y_ = xp.einsum('ij,j...->i...', rotation_e(h), xp.split(y, 4)).flatten()
		y_ = xp.split(y_, 2)
		y_[1] += h * fun(t, y_[0])
		y_[0] += h * fun(t, y_[1])
		return t, xp.concatenate((y_[0], y_[1]), axis=None)
	
	chi = lambda h, t, y: chi_ext(h, t, y, fun)
	chi_star = lambda h, t, y: chi_ext_star(h, t, y, fun)
	y_ = xp.tile(y0, 2)
	sol = integrator.integrator(step).integrate(chi, chi_star, t_span, y_, t_eval=t_eval, command=command)
	y_ = xp.split(sol.y[1], 4, axis=0)
	sol.y = xp.concatenate(((y_[0] + y_[2]) / 2, (y_[1] + y_[3]) / 2), axis=0)
	return sol
