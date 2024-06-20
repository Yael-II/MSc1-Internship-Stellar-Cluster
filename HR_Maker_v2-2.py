"""
* HR Maker: Building HR diagram from BVR.csv files and fit the main sequence
* Version 2.2 - 2024-04-16
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (Intern)
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import LinearNDInterpolator
from amuse.community.seba.interface import SeBa
from amuse.datamodel import Particles
from amuse.units import units as u
from astropy import units as u2
import scipy.constants as const
from astropy.io import fits
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarQuery

# * Parameters

if "YII" in plt.style.available: plt.style.use("YII")

# * Variables
in_directory = "./Cluster_files/"
fig_directory = "../Figs/HR_Maker_figs/"

Cluster_prefix = "Cluster_"
MS_prefix = "MS_"
BVR_suffix = "_BVR.csv"

metallicity = 0
ZSun = 0.0134
bessel_directory = "YBC_tables/ubvrijhklm/regrid/" # Bessel
atlas9_prefix = "Avodonnell94Rv3.1f" # atlas9
atlas9_suffix = "k2odfnew.BC.fits" # atlas9

phoenix_filename = "Avodonnell94Rv3.1BT-Settl_M-0.0_a+0.0.BC.fits" # phoenix
atlas12_filename = "Avodonnell94Rv3.147TucP1.BC.fits" # atlas12

if metallicity >= 0: atlas9_filename = atlas9_prefix + "p{:02d}".format(metallicity) + atlas9_suffix
else: atlas9_filename = atlas9_prefix + "m{:02d}".format(metallicity) + atlas9_suffix

M_bol_sun = 4.74 # Solar bolometric magnitude
sigma = const.Stefan_Boltzmann | u.W * u.m**-2 * u.K**-4 # Stefan-Boltzmann constant

YES = ["Y", "YES", "1", "OUI"]

stellar = SeBa()
t_min = 0 # Myr
t_max = 2e5 # Myr
N_t = int(t_max//10+1)
N_color = 1000
N_stars = 500

print("\033[94m"+"Info: loading extinction map..."+"\033[0m")
ext_query = BayestarQuery(version='bayestar2019')
print("\033[94m"+"Info: done!"+"\033[0m")

coeff_ext_BV = 3.626-2.742
coeff_ext_V = 2.742
R_ext = 3.1

# * Observational data import
Cluster_filelist = np.sort(os.listdir(in_directory))
Cluster_filelist = [f for f in Cluster_filelist if Cluster_prefix in f]
MS_filelist = np.sort(os.listdir(in_directory))
MS_filelist = [f for f in MS_filelist if MS_prefix in f]
Cluster_sources = []
MS_sources = []
sources = []
for i in range(len(Cluster_filelist)):
   k = Cluster_filelist[i][len(Cluster_prefix):].find("_BVR")
   source = Cluster_filelist[i][len(Cluster_prefix):k+len(Cluster_prefix)]
   if source not in Cluster_sources: Cluster_sources.append(source)
for i in range(len(MS_filelist)):
   k = MS_filelist[i][len(MS_prefix):].find("_BVR")
   source = MS_filelist[i][len(MS_prefix):k+len(MS_prefix)]
   if source not in MS_sources: MS_sources.append(source)

sources = np.intersect1d(Cluster_sources, MS_sources)

print("\033[36m"+"Available source:"+"\033[0m")
for i in range(len(sources)):
   print("\033[36m"+"\t{}: ".format(i+1) + "\033[0m" + "{}".format(sources[i]))
k = int(input("\033[32m"+"Select source number: "+"\033[0m"))-1
print("")
source = sources[k]

BVR_Cluster = np.loadtxt(in_directory+Cluster_prefix+source+BVR_suffix, skiprows=1, usecols=(0,1,2,3,4,5,6,7), comments="#", delimiter=",")
BVR_MS = np.loadtxt(in_directory+MS_prefix+source+BVR_suffix, skiprows=1, usecols=(0,1,2,3,4,5,6,7), comments="#", delimiter=",")

RA_Cluster = BVR_Cluster[:,0]
DE_Cluster = BVR_Cluster[:,1]
B_Cluster = BVR_Cluster[:,2]
V_Cluster = BVR_Cluster[:,3]
R_Cluster = BVR_Cluster[:,4]
dB_Cluster = BVR_Cluster[:,5]
dV_Cluster = BVR_Cluster[:,6]
dR_Cluster = BVR_Cluster[:,7]

B_MS = BVR_MS[:,2]
V_MS = BVR_MS[:,3]
R_MS = BVR_MS[:,4]
dB_MS = BVR_MS[:,5]
dV_MS = BVR_MS[:,6]
dR_MS = BVR_MS[:,7]

N_MS = len(B_MS)
N_Cluster = len(B_Cluster)

ra = np.mean(RA_Cluster)
de = np.mean(DE_Cluster)

Coords_Cluster = SkyCoord(ra*u2.deg, de*u2.deg, distance=0*u2.pc, frame='icrs')
print("\033[36m"+"Cluster coordinates: "+"\033[0m" + Coords_Cluster.to_string('hmsdms'))
# * Interpolation data import

hdul = fits.open(bessel_directory+atlas9_filename)
# hdul = fits.open(bessel_directory+phoenix_filename)
# hdul = fits.open(bessel_directory+atlas12_filename)
interp_data = hdul[1].data

# * Main sequence generation
print("\033[94m"+"Info: generating and evolving stars..."+"\033[0m")


def stars_MS(M_min = .5, M_max = 5, z=ZSun, N = 100):
   M_stars = np.linspace(M_min, M_max, N) | u.MSun
   stars = Particles(len(M_stars), mass=M_stars)
   stellar = SeBa()
   stellar.particles.add_particles(stars)
   stellar.set_metallicity(z)
   stellar.commit_particles()
   stellar.update_particles(stars)
   return stars

# * Interpolation functions

interp_logTeff = interp_data.field(0)
interp_logg = interp_data.field(1)
BC_B_bessel = interp_data.field(3)
BC_V_bessel = interp_data.field(4)
BC_R_bessel = interp_data.field(5)

"""
Interpolation using YBC tables (Yang, 2019)
@ https://ui.adsabs.harvard.edu/abs/2019A%26A...632A.105C/abstract
@ https://sec.center/YBC/
"""
interp_BC_B = LinearNDInterpolator(list(zip(interp_logTeff, interp_logg)), BC_B_bessel)
interp_BC_V = LinearNDInterpolator(list(zip(interp_logTeff, interp_logg)), BC_V_bessel) 
interp_BC_R = LinearNDInterpolator(list(zip(interp_logTeff, interp_logg)), BC_R_bessel)

# * Interpolation quantities computation

def extract_BVR_stars(stars):
   L_stars = stars.luminosity.value_in(u.LSun)
   Teff_stars = ((stars.luminosity / (4*np.pi * stars.radius**2 * sigma))**(1/4)).value_in(u.K)
   G_stars = (u.constants.G * stars.mass / stars.radius**2).value_in(u.m * u.s**(-2))
   # type_stars = stars.stellar_type.value_in(u.stellar_type)
   BC_B_stars = interp_BC_B(np.log10(Teff_stars),np.log10(G_stars))
   BC_V_stars = interp_BC_V(np.log10(Teff_stars),np.log10(G_stars))
   BC_R_stars = interp_BC_R(np.log10(Teff_stars),np.log10(G_stars))

   M_bol_stars = M_bol_sun - 2.5*np.log10(L_stars)

   B_abs_stars = M_bol_stars-BC_B_stars
   V_abs_stars = M_bol_stars-BC_V_stars
   R_abs_stars = M_bol_stars-BC_R_stars

   # w1_not = np.where(type_stars != 1)

   # B_abs_stars[w1_not] = np.nan
   # V_abs_stars[w1_not] = np.nan
   # R_abs_stars[w1_not] = np.nan

   return B_abs_stars, V_abs_stars, R_abs_stars


stars = stars_MS(N = N_stars, z=0.2)
N_stars = len(stars)
stellar.particles.add_particles(stars)
stellar.commit_particles()

Time_stars = np.linspace(t_min, t_max, N_t)

def BVR_abs_stars(stellar, stars):
   B_abs_stars = np.empty((N_t,N_stars))
   V_abs_stars = np.empty((N_t,N_stars))
   R_abs_stars = np.empty((N_t,N_stars))
   for i in range(N_t):
      t = Time_stars[i]
      stellar.evolve_model(t|u.Myr)
      stellar.update_particles(stars)
      B_abs_stars[i], V_abs_stars[i], R_abs_stars[i] = extract_BVR_stars(stars)
      
   return B_abs_stars, V_abs_stars, R_abs_stars
 

B_abs_stars, V_abs_stars, R_abs_stars = BVR_abs_stars(stellar, stars)

def V(t, mu):
   i = np.argmin(np.abs(Time_stars-t))
   return V_abs_stars[i]+mu

# * Theoretical values
BV_MS = B_MS-V_MS

BV_Cluster = B_Cluster-V_Cluster

print("\033[94m"+"Info: done!"+"\033[0m")

# * Figure and fitting

def print_results(t_fit, dt_fit, d_fit, dd_fit, ext, ext_BV, ext_V):
   print("\033[36m"+"t = "+"\033[0m" + "{:.2f}±{:.2f}".format(t_fit, dt_fit)+"\033[36m"+" Myr"+"\033[0m")
   print("\033[36m"+"d = "+"\033[0m" + "{:.2f}±{:.2f}".format(d_fit/1000, dd_fit/1000)+"\033[36m"+" kpc"+"\033[0m")
   print("\033[36m"+"ext = "+"\033[0m" + "{:.2f}".format(ext)+"\033[36m"+" mag"+"\033[0m")
   print("\033[36m"+"E(B-V) = "+"\033[0m" + "{:.2f}".format(ext_BV)+"\033[36m"+" mag"+"\033[0m")
   print("\033[36m"+"A_V = "+"\033[0m" + "{:.2f}".format(ext_V)+"\033[36m"+" mag"+"\033[0m")
   return None

def mu(d): return 5*np.log10(d)-5
fig1, ax1 = plt.subplots(1)
fig4, ax4 = plt.subplots(1)

w = np.argwhere(ext_query.distances.to(u2.pc).value < 5000)
ext_tot = ext_query(SkyCoord(ra*u2.deg, de*u2.deg, frame='icrs'), mode='mean')
ax4.plot(ext_query.distances.to(u2.pc).value[w], ext_tot[w], color="C0")

sct = ax4.plot(0,0, color="C4", marker="+", linestyle="", alpha=1, markersize=5)
sctp = ax4.plot(0,0, color="C1", marker="+", linestyle="", alpha=0.65, markersize=3)
sctm = ax4.plot(0,0, color="C1", marker="+", linestyle="", alpha=0.65, markersize=3)


ax4.set_xlabel("Distance [pc]")
ax4.set_ylabel("Reddening")


t_fit = 0
dt_fit = 0
d_fit = 10
dd_fit = 0
ext = 0
ext_BV = 0
ext_V = 0

continue_fit = True
auto_ext = True

ax1.scatter(BV_Cluster, V_Cluster, s=1, marker='.', color="C3")
sct_BV = ax1.scatter(BV_MS, V_MS, s=1, marker='.', color="C0")


X_BV = B_abs_stars[0] - V_abs_stars[0]

Y = V(t_fit, mu(d_fit))
Yp = V(t_fit+dt_fit, mu(d_fit-dd_fit))
Ym = V(t_fit-dt_fit, mu(d_fit+dd_fit))

fit_BV = ax1.plot(X_BV, Y, color="C4", marker=".", linestyle="", markersize=2)
err_BVp = ax1.plot(X_BV, Yp, color="C1", marker=".", linestyle="", alpha=0.65, markersize=1.5)
err_BVm = ax1.plot(X_BV, Ym, color="C1", marker=".", linestyle="", alpha=0.65, markersize=1.5)

ax1.set_xlabel("Color $B-V$")
ax1.set_ylabel("Magnitude $V$")
ax1.invert_yaxis()
ax1.legend([sct_BV, fit_BV[0]], ["Observational data (main sequence)", "$t =$ {:.2f} Myr, $d =$ {:.2f} kpc".format(t_fit, d_fit/1000)], loc="upper right")
fig1.tight_layout()

plt.show(block=False)

print("\033[36m"+"Actions:"+"\033[0m")
print("\033[36m"+"\t- t+/t- [int/float, Myr] (time increase/decrease)"+"\033[0m")
print("\033[36m"+"\t- d+/d- [int/float, pc] (distance increase/decrease)"+"\033[0m")
print("\033[36m"+"\t- vt+/vt- [int/float, Myr] (time uncertainty increase/decrease)"+"\033[0m")
print("\033[36m"+"\t- vd+/vd- [int/float, pc] (distance uncertainty increase/decrease)"+"\033[0m")
print("\033[36m"+"\t- ext+/ext- [int/float] (extinction value)"+"\033[0m")
print("\033[36m"+"\t- ext auto [int/float] (auto extinction value)"+"\033[0m")
print("\033[36m"+"\t- adjust (adjust plots)"+"\033[0m")
print("\033[36m"+"\t- print (show results)"+"\033[0m")
print("\033[36m"+"\t- ok (end)"+"\033[0m")
print("\033[36m"+"\t- cancel (end)"+"\033[0m")


while continue_fit:
   var = input(("\033[32m"+"Enter action: "+"\033[0m"))
   if var == "":
      continue
   elif var == "ok":
      continue_fit = False
   elif var == "cancel":
      continue_fit = False
      plt.close()
      exit()
   elif var == "print":
      print_results(t_fit, dt_fit, d_fit, dd_fit, ext, ext_BV, ext_V)
   elif var == "adjust":
      ax1.set_ylim(np.max(V_Cluster), np.min(V_Cluster))
      ax1.set_xlim(np.min(B_Cluster-V_Cluster), np.max(B_Cluster-V_Cluster))
   elif var[0:3] == "ext":
      if "auto" in var:
         auto_ext = True
      else:
         auto_ext = False
         if var[3] == "=":
            val = float(var[4:])
            ext = val
         else:
            val = float(var[3:])
            ext += val
         ext = np.clip(ext, 0, 10)
         extp = ext
         extm = ext
   elif var[0:2] == "dd": 
      print("\033[33m"+"Warning: please use vd instead!"+ "\033[0m")
      if var[2] == "=":
         val = float(var[3:])
         dd_fit = val
      else:
         val = float(var[2:])
         dd_fit += val
      dd_fit = np.clip(dd_fit, 0, d_fit)
   elif var[0:2] == "dt": 
      print("\033[33m"+"Warning: please use vt instead!"+ "\033[0m")
      if var[2] == "=":
         val = float(var[3:])
         dt_fit = val
      else:
         val = float(var[2:])
         dt_fit += val
      dt_fit = np.clip(dt_fit, 0, t_fit)
   elif var[0] == "d": 
      if var[1] == "=":
         val = float(var[2:])
         d_fit = val
      else:
         val = float(var[1:])
         d_fit += val
      d_fit = np.clip(d_fit, 1, 1e12)
   elif var[0] == "t":
      if var[1] == "=":
         val = float(var[2:])
         t_fit = val
      else:
         val = float(var[1:])
         t_fit += val
      t_fit = np.clip(t_fit, t_min, t_max)
   elif var[0] == "v":
      if var[1] == "d":
         if var[2] == "=":
            val = float(var[3:])
            dd_fit = val
         else:
            val = float(var[2:])
            dd_fit += val
         dd_fit = np.clip(dd_fit, 0, d_fit)
      if var[1] == "t":
         if var[2] == "=":
            val = float(var[3:])
            dt_fit = val
         else:
            val = float(var[2:])
            dt_fit += val
         dt_fit = np.clip(dt_fit, 0, t_fit)
   else:
      print("\033[33m"+"Warning: not found!"+ "\033[0m")


   # * Extinction computation
   i = np.argmin(np.abs(Time_stars-t_fit))
   ip = np.argmin(np.abs(Time_stars-(t_fit+dt_fit)))
   im = np.argmin(np.abs(Time_stars-(t_fit-dt_fit)))
   X_BV = B_abs_stars[i] - V_abs_stars[i]
   Xp_BV = B_abs_stars[ip] - V_abs_stars[ip]
   Xm_BV = B_abs_stars[im] - V_abs_stars[im]
   
   if auto_ext:
      Coords_Cluster = SkyCoord(ra*u2.deg, de*u2.deg, distance=d_fit*u2.pc, frame='icrs')
      ext = ext_query(Coords_Cluster, mode='mean') 
      Coords_Clusterp = SkyCoord(ra*u2.deg, de*u2.deg, distance=(d_fit-dd_fit)*u2.pc, frame='icrs')
      extp = ext_query(Coords_Clusterp, mode='mean') 
      Coords_Clusterm = SkyCoord(ra*u2.deg, de*u2.deg, distance=(d_fit+dd_fit)*u2.pc, frame='icrs')
      extm = ext_query(Coords_Clusterm, mode='mean') 

   ext_BV = coeff_ext_BV * ext
   ext_V = coeff_ext_V * ext
   X_BV_corr = X_BV + ext_BV

   extp_BV = coeff_ext_BV * extp
   extp_V = coeff_ext_V * extp
   Xp_BV_corr = Xp_BV + extp_BV
   
   extm_BV = coeff_ext_BV * extm
   extm_V = coeff_ext_V * extm
   Xm_BV_corr = Xm_BV + extm_BV

   Y = V(t_fit, mu(d_fit)) + ext_V
   Yp = V(t_fit+dt_fit, mu(d_fit-dd_fit)) + extp_V
   Ym = V(t_fit-dt_fit, mu(d_fit+dd_fit)) + extm_V

   fit_BV[0].set_xdata(X_BV_corr)

   fit_BV[0].set_ydata(Y)

   err_BVp[0].set_xdata(Xp_BV_corr)
   err_BVm[0].set_xdata(Xm_BV_corr)

   err_BVp[0].set_ydata(Yp)
   err_BVm[0].set_ydata(Ym)

   sct[0].set_xdata([d_fit])
   sct[0].set_ydata([ext])
   sctp[0].set_xdata([d_fit-dd_fit])
   sctp[0].set_ydata([extp])
   sctm[0].set_xdata([d_fit+dd_fit])
   sctm[0].set_ydata([extm])


   ax1.legend([sct_BV, fit_BV[0]], ["Observational data", "$t =$ {:.2f} $\\pm$ {:.2f} Myr, $d =$ {:.2f} $\\pm$ {:.2f} kpc, $\\mathrm{{ext}} =$ {:.2f}".format(t_fit, dt_fit, d_fit/1000, dd_fit/1000, ext)], loc="upper right")

   fig1.canvas.draw_idle()
   fig1.canvas.flush_events()

   fig4.canvas.draw_idle()
   fig4.canvas.flush_events()


print_results(t_fit, dt_fit, d_fit, dd_fit, ext, ext_BV, ext_V)


save_fig = input("\033[32m"+"Save figure (y/n): "+"\033[0m").upper() in YES

if save_fig:
   fig1.savefig(fig_directory+source+"_HR_B-V_fit")
   fig4.savefig(fig_directory+source+"_HR_ext_fit")

plt.show(block=False)
input(("\033[33m"+"(Press enter to quit)"+"\033[0m"))
plt.close()

