import numpy as np
import matplotlib.pyplot as plt
from pyhamsys import HamSys, solve_ivp_sympext

## Parameters
epsilon = 0.027				# parameter of the Hamiltonian system
tstep = 2 * np.pi / 11		# integration timestep 
N = 200  					# number of trajectories
nf = 1500 					# number of points on the Poincaré section per trajectory

## Hamiltonian systems and equations of motion 
def eqn(t, y):
	y_ = np.split(y, 2)
	return np.concatenate((y_[1], epsilon * (np.sin(y_[0] - t) + np.sin(y_[0]))), axis=None)

def hamiltonian(t, y):
	y_ = np.split(y, 2)
	return np.sum(y_[1]**2 / 2 + epsilon * (np.cos(y_[0] - t) + np.cos(y_[0])))

def keqn(t, y):
	y_ = np.split(y, 2)
	return -epsilon * np.sum(np.sin(y_[0] - t))

hs = HamSys(ndof=N + 0.5, hamilton_eqn=eqn, hamilton_keqn=keqn, hamiltonian=hamiltonian)
	
x0 = 2 * np.pi * np.random.random(N)
p0 = np.random.random(N)
y0 = np.concatenate((x0, p0), axis=None)

## Integration
sol = solve_ivp_sympext(hs, (0, 2 * np.pi * nf), y0, step=tstep, t_eval=2 * np.pi * np.arange(nf + 1), check_energy=True)
print(sol.err)

## Plot of the Poincaré section
fig, ax = plt.subplots()
x, p = np.split(sol.y, 2)
ax.plot(x % (2 * np.pi), p, '.b')
ax.set_xlim((0, 2 * np.pi))
plt.show()