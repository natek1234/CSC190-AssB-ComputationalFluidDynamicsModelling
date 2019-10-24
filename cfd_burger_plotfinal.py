import numpy as numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot, cm
from matplotlib.colors import Normalize

def equation_of_motion(u, v, …):
   # generate the next state as a function of the old state
   un = u.copy()
   vn = v.copy()
   u[1:-1, 1:-1] = # see reference
   v[1:-1, 1:-1] = # see reference
   return (u,v)

def boundary(u, v, nozzle_u, nozzle_v, nx, ny, t_step):
   u[0, :] = 0
   u[-1, :] = 0
   u[:, 0] = 0
   u[:, -1] = 0

   v[0, :] = 0
   v[-1, :] = 0
   v[:, 0] = 0
   v[:, -1] = 0

   # special nozzle BC
   u[ny//2-2:ny//2+2, 0] = nozzle_u[t_step]
   v[ny//2-2:ny//2+2, 0] = nozzle_v[t_step]

   return (u, v)

def evolve(u, v, …, steps):
   for i in range(steps):
      (u, v) = equation_of_motion(u, v, …)
      (u, v) = boundary(u, v, …, i)
   return (u, v)

def main():
   nx = 41
   ny = 41
   dx = 2 / float(nx - 1)
   dy = 2 / float(ny - 1)
   u = numpy.ones((ny, nx))
   v = numpy.ones((ny, nx))
   sigma = .001
   nu = 0.01
   dt = sigma * dx * dy / nu
   nt = 2510 # the number of steps we’re simulating
   ###Assign initial conditions
   initial_u = numpy.zeros((nx, ny))
   initial_v = numpy.zeros((nx, ny))
   ### Special BC for nozzle
   # located at (0, 1)
   nozzle_u = numpy.append(10*numpy.ones(1000), numpy.zeros(nt))
   nozzle_v = numpy.append(10*numpy.ones(1000), numpy.zeros(nt))
   (final_u, final_v) = evolve(initial_u, initial_v, …, nt)

main()
