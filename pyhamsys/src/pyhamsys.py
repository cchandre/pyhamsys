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
from typing import Callable, Tuple, Union

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

class SymplecticIntegrator:
	"""
    Some symplectic splitting integrators in Python
    .. version:: 0.1.0

    Attributes
    ----------
    name : str
        the name of the symplectic integrator
    step : float
        the time step for the integrator

    Methods
    -------
    _integrate : integrate the state vector y by one step.
	integrate : integrate the state vector y from 0 to tf.
        if times is int or float : tf = times, and returns only the final value of the state vector y. 
		if times is a list or numpy array : tf = times[-1]. 
		if len(times)=1 : returns all intermediate steps. 
		if len(times)>=2 : returns only the times 
		specified in times.
    """
	def __repr__(self) -> str:
		return f'{self.__class__.__name__}({self.name}, {self.step})'
	
	def __str__(self) -> str:
		return f'{self.name}'

	def __init__(self, name:str, step:float) -> None:
		self.name = name
		self.step = step
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
				N = int(self.name[2:]) // 2
				alpha_s = xp.asarray([0.5])
				for n in range(1, N):
					x1 = 1 / (2 - 2**(1/(2*n+1)))
					x0 = 1 - 2 * x1
					alpha_ = xp.concatenate((alpha_s, xp.flip(alpha_s)))
					alpha_ = xp.concatenate((x1 * alpha_, x0 * alpha_s))
					alpha_s = alpha_.copy()
			except:
				raise NameError(f'{self.name} integrator not defined') 
		elif self.name.endswith('EFRL'):
			if self.name.startswith('V'):
				xi = 0.1644986515575760
				lam = -0.02094333910398989
				chi = 1.235692651138917
			elif self.name.startswith('P'):
				xi = 0.1786178958448091
				lam = -0.2123418310626054
				chi = -0.06626458266981849
			else:
				xi = 0.1720865590295143
				lam = -0.09156203075515678
				chi = -0.1616217622107222
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
		else:
			raise NameError(f'{self.name} integrator not defined')
		self.alpha_s = xp.concatenate((alpha_s, xp.flip(alpha_s)))
		self.alpha_s *= self.step
		self.alpha_o = xp.tile([1, 0], len(alpha_s))
		self.alpha_s = xp.concatenate((alpha_s, xp.flip(alpha_s)))
		self.alpha_s *= self.step
		self.alpha_o = xp.tile([1, 0], len(alpha_s))
	
	def _integrate(self, chi:Callable, chi_star:Callable, y):
		for h, st in zip(self.alpha_s, self.alpha_o):
			y = chi(h, y) if st==0 else chi_star(h, y)
		return y
	
	def integrate(self, chi:Callable, chi_star:Callable, y:xp.ndarray, times:Union[int, float, list, xp.ndarray], command:Callable=None, autonomous:bool=True) -> Tuple[Union[float, xp.ndarray], xp.ndarray]:
		"""
		Integrate the Hamiltonian flow from the initial conditions specified 
		by the initial state vector y using one of the selected symplectic 
		splitting integrators.
		Returns the value of y at times defines by the float, list or array times.
		.. versionadded:: 0.1.0

		Parameters
		----------
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
		-------
		out : array of times, and values of y at times
			If times is a float of integer, the output is a tuple (t, y or x) where
			y is the value of the state vector and y = [t, x] if autonomous is False.
			If times is a list or array, returns the times and values of y or x at 
			times. If times is a list or array with a single element, returns the 
			times and values of y or x at all computed times. 

		References
		----------
			McLachlan, R.I, 2022, "Tuning symplectic integrators is easy and 
			worthwhile", *arxiv:2104.10269*.
		"""
		t, y_ = 0, y.copy()
		t_vec, y_vec = [0], y_.copy()[..., xp.newaxis] if autonomous else y_[1][..., xp.newaxis]
		times = xp.asarray(times) if isinstance(times, list) else times
		while t < (times if isinstance(times, (int, float)) else times.max()):
			y_ = self._integrate(chi, chi_star, y_)
			t += self.step
			if not isinstance(times, (int, float)):
				t_vec.append(t)
				yt = y_[..., xp.newaxis] if autonomous else y_[1][..., xp.newaxis]
				y_vec = xp.concatenate((y_vec, yt), axis=-1)
			if command is not None:
				command(t, y_)
		t_vec = xp.asarray(t_vec)
		if isinstance(times, (int, float)):
			return t, y_ if autonomous else y_[1]
		elif len(times) == 1:
			return t_vec, y_vec
		else:
			times.sort()
			return xp.asarray(times[times>=0]), interp1d(t_vec, y_vec, assume_sorted=True)(times[times>=0])