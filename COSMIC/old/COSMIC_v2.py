"""
* COSMIC : Cluster Orbital SysteM Integration Code
* Version 2 - 2024-02-27
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (Intern)
"""
import numpy as np
import amuse.units.units as u
import datetime


from amuse.units import nbody_system
from amuse.ic.kingmodel import new_king_model
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.datamodel.particles import Channels
from amuse.community.ph4.interface import Ph4
from amuse.community.seba.interface import SeBa
from amuse.couple.bridge import Bridge


# * Parameters
N_part = 10  # Number of particle
N_iter = 30  # Number of iteration
dt = 1|u.Myr   # Time increment

R0 = 1|u.parsec       # Radius of the system
W0 = 9                # King parameter
z = 0.02                 # Metallicity

# * Variables
date = datetime.datetime.now()
date = "{:04d}-{:02d}-{:02d}_{:02d}h{:02d}".format(date.year, date.month, date.day, date.hour, date.minute)

directory = "./output/"
filename_cluster = "{}_cluster.csv".format(date)
filename_stars = "{}_stars.csv".format(date)

# * Function
def get_total_momentum(M, X, V): 
   """
   Compute the total momentum of particles with mass M, positions X and velocities V.
   """
   l = 0
   for j in range(len(X)):
      m = M[j].value_in(u.MSun)
      x = X[j].value_in(u.pc)
      v = V[j].value_in(u.pc/u.Myr)
      l += m*np.cross(x,v)
   return np.sqrt(np.sum(l**2))|(u.pc**2 * u.MSun / u.Myr)

# * Initialization
def create_cluster(N:int, R=1|u.pc, W=9):
   """
   Create cluster of N stars of radius R (default: 1 pc) following King model of parameter W (default: 9).
   """
   print("\033[94m"+"Info: Simulation with "+"\033[0m"+"{}".format(N)+"\033[94m"+" stars"+"\033[0m")
   M_stars = new_salpeter_mass_distribution(N)
   converter = nbody_system.nbody_to_si(M_stars.sum(), R)
   stars = new_king_model(N, W, converter)
   stars.mass = M_stars
   stars.zams_mass = M_stars
   return stars, converter

# * Evolution
def stellar_gravity_evolution(stars, converter, t_max=16|u.Myr, dt=1|u.Myr, out_dt=1|u.Myr):
   """
   Evolve a model with both gravity and stellar evolution. dt is the simulation time step, while out_dt is the delay between printing in file.
   """
   print("\033[94m"+"Info: Simulation duration of "+"\033[0m"+"{:.2f} Myr".format(t_max.value_in(u.Myr)))
   init_files()

   gravity = Ph4(converter) # gravity model
   stellar = SeBa()         # stellar model
   bridge = Bridge()        # bridge between models

   # TODO gravity.initialize_code()

   gravity.particles.add_particles(stars)
   stellar.particles.add_particles(stars)

   stellar.set_metallicity(z)
   stellar.commit_particles()
   stellar.update_particles(stars) # TODO pourquoi ?
   
   bridge.add_system(gravity)
   bridge.add_system(stellar)

   bridge.synchronize_model()

   bridge.channels.add_channel(stellar.particles.new_channel_to(gravity.particles, attributes=["mass", "radius"])) 
   # Sync the mass and radius between the two models

   out_channel = Channels()
   out_channel.add_channel(stellar.particles.new_channel_to(stars))
   out_channel.add_channel(gravity.particles.new_channel_to(stars))

   bridge.timestep = dt # TODO timestep dynamique
   print(gravity.get_initial_timestep_median()) # TODO ????
   i = 0
   t = 0 | u.Myr
   E0 = gravity.get_total_energy()
   L0 = get_total_momentum(stars.mass, stars.position, stars.velocity)
   while t <= t_max:
      print("\033[36m"+"t: "+"\033[0m"+"{:02.2f} Myr/{:02.2f} Myr".format(t.value_in(u.Myr),t_max.value_in(u.Myr)))
      bridge.evolve_model(t)
      out_channel.copy()
      E = gravity.get_total_energy()
      L = get_total_momentum(stars.mass, stars.position, stars.velocity)
      write_files(stars, i, t, E, E0, L, L0)
      i += 1
      t += out_dt
   
   print("\033[36m"+"t: "+"\033[0m"+"{:02.2f} Myr/{:02.2f} Myr".format(t_max.value_in(u.Myr),t_max.value_in(u.Myr)))
   print("\033[36m"+"Simulation completed"+"\033[0m")
   return 0

# * Files initialization and write 
def init_files():
   file = open(directory+filename_cluster, "w")
   file.write("#i (time), t [Myr], E [pc2 MSun Myr-2], (E-E0)/E0, L [pc2 MSun Myr-1], (L-L0)/L0\n")
   file.close()

   file = open(directory+filename_stars, "w")
   file.write("#i (time), j (star), m [MSun], x [pc], y [pc], z [pc], vx [km s-1], vy [km s-1], vz [km s-1], L [LSun], R [RSun], T [K]\n")
   file.close()

   print("\033[94m"+"Info: Created two output files with date " + "\033[0m" + date.replace("h", ":").replace("_", " "))
   return 0

def write_files(stars, i, t, E, E0, L, L0):
   file = open(directory+filename_cluster, "a")
   file.write("{:d}, ".format(i))
   file.write("{:+e}, ".format(t.value_in(u.Myr)))
   file.write("{:+e}, ".format(E.value_in(u.pc**2 * u.MSun * u.Myr**-2)))
   file.write("{:+e}, ".format((E-E0)/E0))
   file.write("{:+e}, ".format(L.value_in(u.pc**2 * u.MSun * u.Myr**-1)))
   file.write("{:+e}\n".format((L-L0)/L0))
   file.close()

   file = open(directory+filename_stars, "a")
   for j in range(N_part):
      file.write("{:d}, ".format(i))
      file.write("{:d}, ".format(j))
      file.write("{:+e}, ".format(stars.mass[j].value_in(u.MSun)))
      file.write("{:+e}, ".format(stars.position[j,0].value_in(u.pc)))
      file.write("{:+e}, ".format(stars.position[j,1].value_in(u.pc)))
      file.write("{:+e}, ".format(stars.position[j,2].value_in(u.pc)))
      file.write("{:+e}, ".format(stars.velocity[j,0].value_in(u.km * u.s**-1)))
      file.write("{:+e}, ".format(stars.velocity[j,1].value_in(u.km * u.s**-1)))
      file.write("{:+e}, ".format(stars.velocity[j,2].value_in(u.km * u.s**-1)))
      file.write("{:+e}, ".format(stars.luminosity[j].value_in(u.LSun)))
      file.write("{:+e}, ".format(stars.radius[j].value_in(u.RSun)))
      file.write("{:+e}\n".format(stars.temperature[j].value_in(u.K)))
   file.close()
   return 0

# * Main
cluster, nbody_converter = create_cluster(N_part, R0, W0)
stellar_gravity_evolution(cluster, nbody_converter, N_iter*dt, dt, dt)




