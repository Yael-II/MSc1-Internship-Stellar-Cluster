"""
* COSMIC-VIS : Cluster Orbital SysteM Integration Code - Visualisation
* Version 2 - 2024-02-27
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (Intern)
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# * Parameters
if "YII" in plt.style.available: plt.style.use("YII")

if "YII_light" in plt.style.available: light = "YII_light"
else: light = "default"

# * Variables

in_directory = "./COSMIC_output/"
out_directory = "../Figs/COSMIC_figs/"

cluster_suffix = "_cluster.csv"
stars_suffix = "_stars.csv"

YES = ["Y", "YES", "1", "OUI"]
MARKER = ["o", "s", "D", "v", "^", "<", ">"]

Tmin = 2000  # K - Temperature truncation
Tmax = 10000 # K - Temperature truncation

# * Functions

def density_hist(X,Y,M,T,bins,i):
   x_cm = np.average(X[i], weights=M[i])
   y_cm = np.average(Y[i], weights=M[i])
   Radius = np.sqrt((X[i]-x_cm)**2 + (Y[i]-y_cm)**2)
   Rho = np.zeros(len(bins)-1)
   for k in range(len(bins)-1):
      r1 = bins[k]
      r2 = bins[k+1]
      w = np.where((Radius > r1) & (Radius < r2))
      Rho[k] = np.sum(M[i,w])/(4*np.pi*r2**3/3 - 4*np.pi*r1**3/3)
   return Rho

# * Listing files for selection 
filelist = np.sort(os.listdir(in_directory))
filelist = [f for f in filelist if f[0] != "."]
dates = []
for f in filelist:
   if len(f) == len("YYYY-MM-DD_HH:MM_cluster.csv") or len(f) == len("YYYY-MM-DD_HH:MM_stars.csv"):
      date = f[0:16]
      if not date in dates:
         dates.append(date)
   
# * File selection
print("\033[36m"+"Available dates:"+"\033[0m")
for i in range(len(dates)):
   print("\033[36m"+"\t{}: ".format(i+1) + "\033[0m" + "{}".format(dates[i].replace("h", ":").replace("_", " ")))
j = int(input("\033[32m"+"Select a date number: "+"\033[0m"))-1
print("")
date_prefix = dates[j]

# * Data collection
cluster = np.loadtxt(in_directory+date_prefix+cluster_suffix, skiprows=1, comments="#", delimiter=",")
stars = np.loadtxt(in_directory+date_prefix+stars_suffix, skiprows=1, comments="#", delimiter=",")

T = cluster[:,1]
E = cluster[:,2]
DE = cluster[:,3]
L = cluster[:,4]
DL = cluster[:,5]

M = np.zeros(np.int32(stars[-1,0:2])+1)
X = np.zeros(np.int32(stars[-1,0:2])+1)
Y = np.zeros(np.int32(stars[-1,0:2])+1)
Z = np.zeros(np.int32(stars[-1,0:2])+1)
VX = np.zeros(np.int32(stars[-1,0:2])+1)
VY = np.zeros(np.int32(stars[-1,0:2])+1)
VZ = np.zeros(np.int32(stars[-1,0:2])+1)
Lumi = np.zeros(np.int32(stars[-1,0:2])+1)
Radi = np.zeros(np.int32(stars[-1,0:2])+1)
Temp = np.zeros(np.int32(stars[-1,0:2])+1)

for star in stars:
   i = np.int32(star[0])
   j = np.int32(star[1])
   M[i,j] = star[2]
   X[i,j] = star[3]
   Y[i,j] = star[4]
   Z[i,j] = star[5]
   VX[i,j] = star[6]
   VY[i,j] = star[7]
   VZ[i,j] = star[8]
   Lumi[i,j] = star[9]
   Radi[i,j] = star[10]
   Temp[i,j] = star[11]
N_iter = i
N_part = j

x_cm = np.array([np.average(X[i], weights=M[i]) for i in range(N_iter)])
y_cm = np.array([np.average(Y[i], weights=M[i]) for i in range(N_iter)])
z_cm = np.array([np.average(Z[i], weights=M[i]) for i in range(N_iter)])

lim_min = np.min((X,Y))
lim_max = np.max((X,Y))
Alpha = np.clip((np.log(Lumi) - np.min(np.log(Lumi)))/(np.max(np.log(Lumi)) - np.min(np.log(Lumi))), 0.1, 1)

# * Output choice
print("\033[36m"+"Choose figures:"+"\033[0m")
figure_1 = input("\033[32m"+"\tFigure 1 (energy, angular momentum, time): "+"\033[0m").upper() in YES
figure_2 = input("\033[32m"+"\tFigure 2a and 2b (positions, speed): "+"\033[0m").upper() in YES
figure_3 = input("\033[32m"+"\tFigure 3 (mass distribution): "+"\033[0m").upper() in YES
figure_4 = input("\033[32m"+"\tFigure 4 (density distribution): "+"\033[0m").upper() in YES
figure_5 = input("\033[32m"+"\tFigure 5 (HR diagram): "+"\033[0m").upper() in YES
print("")

print("\033[36m"+"General settings:"+"\033[0m")
time =  np.int32(input("\033[32m"+"\tStatic images drawing instant ({} < int < {}): ".format(0,N_iter)+"\033[0m"))
if input("\033[32m"+"\tWindow size restriction (y/n): "+"\033[0m").upper() in YES:
      lim = np.float32(input("\033[32m"+"\t\tLimit [pc]: "+"\033[0m"))
      lim_min = -lim/2
      lim_max = +lim/2
print("")

draw_time = scatter_time = distrib_time = diagram_time = time

if figure_1:
   print("\033[36m"+"Figure 1 (energy, angular momentum, time):"+"\033[0m")
   plot_E = input("\033[32m"+"\tPlot energy vs. time (y/n): "+"\033[0m").upper() in YES
   plot_L = input("\033[32m"+"\tPlot angular momentum vs. time (y/n): "+"\033[0m").upper() in YES
   plot_save = input("\033[32m"+"\tSave figure (y/n): "+"\033[0m").upper() in YES
   print("")
else: plot_E = plot_L = plot_save = False

if figure_2:
   print("\033[36m"+"Figure 2a and 2b (positions, speed):"+"\033[0m")
   scatter_2D = input("\033[32m"+"\tScatter positions in 2D (y/n): "+"\033[0m").upper() in YES
   scatter_3D = input("\033[32m"+"\tScatter positions in 3D (y/n): "+"\033[0m").upper() in YES
   scatter_speed = input("\033[32m"+"\tShow speed vectors (y/n): "+"\033[0m").upper() in YES
   scatter_save = input("\033[32m"+"\tSave figure (y/n): "+"\033[0m").upper() in YES
   scatter_animate = input("\033[32m"+"\tAnimate over time (y/n): "+"\033[0m").upper() in YES
   print("")
else: scatter_2D = scatter_3D = scatter_speed = scatter_save = scatter_animate = False

if figure_3:
   print("\033[36m"+"Figure 3 (mass distribution):"+"\033[0m")
   distrib_mass = input("\033[32m"+"\tShow mass distribution (y/n): "+"\033[0m").upper() in YES
   distrib_bins = np.int32(input("\033[32m"+"\tNumber of bins (int): "+"\033[0m"))
   distrib_save = input("\033[32m"+"\tSave figure (y/n): "+"\033[0m").upper() in YES
   print("")
else: distrib_mass = distrib_bins = distrib_save = False

if figure_4:
   print("\033[36m"+"Figure 4 (density distribution):"+"\033[0m")
   draw_density = input("\033[32m"+"\tShow density distribution (y/n): "+"\033[0m").upper() in YES
   draw_bins = np.float32(input("\033[32m"+"\tWidth of bins [pc]: "+"\033[0m"))
   draw_max = np.float32(input("\033[32m"+"\tMaximum distance [pc]: "+"\033[0m"))
   draw_save = input("\033[32m"+"\tSave figure (y/n): "+"\033[0m").upper() in YES
   draw_animate = input("\033[32m"+"\tAnimate over time (y/n): "+"\033[0m").upper() in YES
   print("")
else: draw_density = draw_bins = draw_max = draw_save = draw_animate = False

if figure_5:
   print("\033[36m"+"Figure 5 (HR diagram):"+"\033[0m")
   diagram_HT = input("\033[32m"+"\tDraw HR diagram (y/n): "+"\033[0m").upper() in YES
   diagram_save = input("\033[32m"+"\tSave figure (y/n): "+"\033[0m").upper() in YES
   diagram_animate = input("\033[32m"+"\tAnimate over time (y/n): "+"\033[0m").upper() in YES
else : diagram_HT = diagram_save = diagram_animate = False

# * Plotting
if plot_E or plot_L:
   fig1, axs1 = plt.subplots(2, sharex=True)
   if plot_E and not plot_L: # Just the energy
      axs1[1].set_xlabel("$t\\ [\\mathrm{{Myr}}]$")
      axs1[0].set_ylabel("$E\\ [\\mathrm{{pc^2\\ M_\\odot\\ Myr^{{-2}}}}]$")
      axs1[1].set_ylabel("$\\Delta E/E_0$")
      axs1[0].plot(T, E, color="C3")
      axs1[1].plot(T, DE, color="C3")
   elif plot_L and not plot_E: # Just the angular momentum
      axs1[1].set_xlabel("$t\\ [\\mathrm{{Myr}}]$")
      axs1[0].set_ylabel("$L\\ [\\mathrm{{pc^2\\ M_\\odot\\ Myr^{{-1}}}}]$")
      axs1[1].set_ylabel("$\\Delta L/L_0$")
      axs1[0].plot(T, L, color="C4")
      axs1[1].plot(T, DL, color="C4")
   elif plot_L and plot_E: # ...both
      xas1 = [axs1[i].twinx() for i in range(len(axs1))]
      axs1[1].set_xlabel("$t\\ [\\mathrm{{Myr}}]$")
      axs1[0].set_ylabel("$E\\ [\\mathrm{{pc^2\\ M_\\odot\\ Myr^{{-2}}}}]$")
      axs1[1].set_ylabel("$\\Delta E/E_0$")
      xas1[0].set_ylabel("$L\\ [\\mathrm{{pc^2\\ M_\\odot\\ Myr^{{-1}}}}]$")
      xas1[1].set_ylabel("$\\Delta L/L_0$")
      axs1[0].plot(T, E, color="C3")
      axs1[1].plot(T, DE, color="C3")
      xas1[0].plot(T, L, color="C4")
      xas1[1].plot(T, DL, color="C4")
      axs1[0].tick_params(axis='y', which="both", colors='C3')
      axs1[0].yaxis.label.set_color('C3')
      xas1[0].tick_params(axis='y', which="both", colors='C4')
      xas1[0].yaxis.label.set_color('C4')
      axs1[1].tick_params(axis='y', which="both", colors='C3')
      axs1[1].yaxis.label.set_color('C3')
      xas1[1].tick_params(axis='y', which="both", colors='C4')
      xas1[1].yaxis.label.set_color('C4')
   if plot_save:
      fig1.tight_layout()
      fig1.savefig(out_directory+date_prefix+"_figure_1")
      print("\033[94m"+"Info: Figure saved as " + "\033[0m" + date_prefix+"_figure_1")

if scatter_2D or scatter_3D:
   if scatter_2D:
      fig2a, ax2a = plt.subplots(1)
      # ax2a.set_facecolor('k')
      i = scatter_time
      ax2a.set_title("$t =$ {:.2f} Myr".format(T[i]))
      if scatter_speed: 
         qvr2a = ax2a.quiver(X[i]-x_cm[i],Y[i]-y_cm[i],VX[i],VY[i], alpha=0.1, headwidth=2, headlength=2, headaxislength=2, color="k")
         ax2a.quiverkey(qvr2a, 0.95, 0.05, 1, "Velocity scale: $1\mathrm{{~km\\cdot s^{{-1}}}}$", labelpos='W', coordinates='figure', color="k", alpha=1)
      sct2a = ax2a.scatter(X[i]-x_cm[i],Y[i]-y_cm[i],s=np.clip(Radi[i]**2,0.5,30), c=Temp[i], alpha=Alpha[i], cmap="RdYlBu", vmin=Tmin, vmax=Tmax)
      ax2a.set_xlabel("$x$ [pc]")
      ax2a.set_ylabel("$y$ [pc]")
      ax2a.set_xlim(lim_min, lim_max)
      ax2a.set_ylim(lim_min, lim_max)
      ax2a.set_aspect("equal")
      fig2a.colorbar(sct2a, label="Stellar effective temperature $T$ [K]", extend='both')
      fig2a.tight_layout()
      if scatter_save:
         fig2a.savefig(out_directory+date_prefix+"_figure_2a")
         print("\033[94m"+"Info: Figure saved as " + "\033[0m" + date_prefix+"_figure_2a")
   if scatter_3D:
      fig2b, ax2b = plt.subplots(1, subplot_kw={"projection": "3d"})
      # ax2b.set_facecolor('k')
      i = scatter_time
      ax2b.set_title("$t =$ {:.2f} Myr".format(T[i]))
      if scatter_speed: 
         qvr2b = ax2b.quiver(X[i]-x_cm[i],Y[i]-y_cm[i],Z[i]-z_cm[i],VX[i],VY[i],VZ[i], alpha=0.1, color="k")
         # ax2b.quiverkey(qvr2b, 0.95, 0.05, 2, "Velocity scale: $2\mathrm{{~km\\cdot s^{{-1}}}}$", labelpos='W', coordinates='figure')
      sct2b = ax2b.scatter(X[i]-x_cm[i],Y[i]-y_cm[i],Z[i]-z_cm[i],s=np.clip(Radi[i]**2,0.5,30), c=Temp[i], alpha=Alpha[i], cmap="RdYlBu", vmin=Tmin, vmax=Tmax)
      ax2b.set_xlabel("$x$ [pc]")
      ax2b.set_ylabel("$y$ [pc]")
      ax2b.set_zlabel("$z$ [pc]")
      ax2b.set_xlim(lim_min, lim_max)
      ax2b.set_ylim(lim_min, lim_max)
      ax2b.set_zlim(lim_min, lim_max)
      ax2b.set_aspect("equal")
      fig2b.colorbar(sct2b, label="Stellar effective temperature $T$ [K]", location="left", extend='both')
      fig2b.tight_layout()
      if scatter_save:
         fig2b.savefig(out_directory+date_prefix+"_figure_2b")
         print("\033[94m"+"Info: Figure saved as " + "\033[0m" + date_prefix+"_figure_2b")

if distrib_mass:
   fig3, ax3 = plt.subplots(1)
   i = distrib_time
   ax3.set_title("$t =$ {:.2f} Myr".format(T[i]))
   nbins = distrib_bins
   hist, lbin = np.histogram(M[i], nbins)
   cbins = 0.5*(lbin[1:]+ lbin[: -1])
   ax3.scatter(cbins, hist, s=5)
   ax3.set_xscale ("log")
   ax3.set_yscale ("log")
   ax3.set_xlabel("Masses $M\\ [\\mathrm{{M_\\odot}}]$")
   ax3.set_ylabel("Star distribution")
   if distrib_save:
      fig3.savefig(out_directory+date_prefix+"_figure_3")
      print("\033[94m"+"Info: Figure saved as " + "\033[0m" + date_prefix+"_figure_3")

if draw_density:
   fig4, ax4 = plt.subplots(1)
   i = draw_time
   bin_width = draw_bins
   bin_max = draw_max
   bins = np.arange(0,bin_max,bin_width)
   Rho = density_hist(X,Y,M,T,bins,i)
   cbins = 0.5*(bins[1:]+ bins[: -1])
   ax4.set_title("$t =$ {:.2f} Myr".format(T[i]))
   ax4.set_xscale ("log")
   ax4.set_yscale ("log")
   ax4.set_xlabel("Distance from center $r$ [pc]")
   ax4.set_ylabel("Stellar masses density $\\rho\\ [\\mathrm{{M_\\odot\\ pc^{{-3}}}}]$")
   ax4.scatter(cbins, Rho, s=5)
   fig4.tight_layout()
   if draw_save:
      fig4.savefig(out_directory+date_prefix+"_figure_4")
      print("\033[94m"+"Info: Figure saved as " + "\033[0m" + date_prefix+"_figure_4")
   
if diagram_HT:
   fig5, ax5 = plt.subplots(1)
   ax5.invert_xaxis()
   i = diagram_time
   ax5.set_title("$t =$ {:.2f} Myr".format(T[i]))
   ax5.set_xlabel("Temperature $T$ [K]")
   ax5.set_ylabel("Luminosity $L\\ [\\mathrm{{L_\\odot}}]$")
   ax5.set_xscale ("log")
   ax5.set_yscale ("log")
   ax5.scatter(Temp[i], Lumi[i], s=5, alpha=Alpha[i], c=Temp[i], cmap="RdYlBu", vmin=Tmin, vmax=Tmax)
   fig5.tight_layout()
   if diagram_save:
      fig5.savefig(out_directory+date_prefix+"_figure_5")
      print("\033[94m"+"Info: Figure saved as " + "\033[0m" + date_prefix+"_figure_5")

if scatter_animate:
   figA2a,axA2a = plt.subplots(1)
   # ax2a.set_facecolor('k')
   axA2a.set_title("$t =$ {:.2f} Myr".format(T[0]))
   if scatter_speed:
      qvr = axA2a.quiver(X[0]-x_cm[0],Y[0]-y_cm[0],VX[0],VY[0], alpha=0.5, headwidth=2, headlength=2, headaxislength=2, color="k")
      axA2a.quiverkey(qvr, 0.95, 0.05, 1, "Velocity scale: $1\mathrm{{~km\\cdot s^{{-1}}}}$", labelpos='W', coordinates='figure', color="k", alpha=1)
   sctA2a = axA2a.scatter(X[0]-x_cm[0],Y[0]-y_cm[0],s=np.clip(Radi[0]**2,0.5,30), c=Temp[0], alpha=Alpha[0], cmap="RdYlBu", vmin=Tmin, vmax=Tmax)
   figA2a.colorbar(sctA2a, label="Stellar effective temperature $T$ [K]", extend='both')
   axA2a.set_xlabel("$x$ [pc]")
   axA2a.set_ylabel("$y$ [pc]")
   axA2a.set_aspect("equal")
   axA2a.set_xlim(lim_min, lim_max)
   axA2a.set_ylim(lim_min, lim_max)
   figA2a.tight_layout()
   def anim2a(i):
      if scatter_speed:
         qvr.set_offsets(np.array([X[i]-x_cm[i],Y[i]-y_cm[i]]).T)
         qvr.set_UVC(VX[i],VY[i])
      axA2a.set_title("$t =$ {:.2f} Myr".format(T[i]))
      sctA2a.set_offsets(np.array([X[i]-x_cm[i],Y[i]-y_cm[i]]).T)
      sctA2a.set_array(Temp[i])
      sctA2a.set_sizes(np.clip(Radi[i]**2,0.5,30))
      sctA2a.set_alpha(Alpha[i])
      return sctA2a,
      
   ani2a = animation.FuncAnimation(fig=figA2a, func=anim2a, frames=N_iter, interval=100)
   if scatter_save:
      ani2a.save(out_directory+date_prefix+"_animation_2a.png", dpi=300, fps=10)
      ani2a.save(out_directory+date_prefix+"_animation_2a.pdf", dpi=300, fps=10)
      print("\033[94m"+"Info: Animation saved as " + "\033[0m" + date_prefix+"_animation_2a")


if draw_animate:
   figA4, axA4 = plt.subplots(1)
   bin_width = draw_bins
   bin_max = draw_max
   bins = np.arange(0,bin_max,bin_width)
   cbins = 0.5*(bins[1:]+ bins[: -1])
   Rho = np.array([density_hist(X,Y,M,T,bins,k) for k in range(N_iter)])
   axA4.set_title("$t =$ {:.2f} Myr".format(T[0]))
   sctA4 = axA4.scatter(cbins, Rho[0], s=5)
   axA4.set_xscale ("log")
   axA4.set_yscale ("log")
   axA4.set_ylim(top=np.max(Rho))
   axA4.set_xlabel("Distance from center $r$ [pc]")
   axA4.set_ylabel("Stellar masses density $\\rho\\ [\\mathrm{{M_\\odot\\ pc^{{-3}}}}]$")
   figA4.tight_layout()
   def anim4(i):
      axA4.set_title("$t =$ {:.2f} Myr".format(T[i]))
      sctA4.set_offsets(np.array([cbins, Rho[i]]).T)
      return sctA4
   ani4 = animation.FuncAnimation(fig=figA4, func=anim4, frames=N_iter, interval=100)
   if scatter_save:
      ani4.save(out_directory+date_prefix+"_animation_4.png", dpi=300, fps=10)
      ani4.save(out_directory+date_prefix+"_animation_4.pdf", dpi=300, fps=10)
      print("\033[94m"+"Info: Animation saved as " + "\033[0m" + date_prefix+"_animation_4")

if diagram_animate:
   figA5, axA5 = plt.subplots(1)
   axA5.invert_xaxis()
   axA5.set_title("$t =$ {:.2f} Myr".format(T[0]))
   axA5.set_xlabel("Temperature $T$ [K]")
   axA5.set_ylabel("Luminosity $L\\ [\\mathrm{{L_\\odot}}]$")
   axA5.set_xscale ("log")
   axA5.set_yscale ("log")
   axA5.set_xlim(np.max(Temp), np.min(Temp))
   axA5.set_ylim(np.min(Lumi), np.max(Lumi))
   sctA5 = axA5.scatter(Temp[0], Lumi[0], s=5, alpha=Alpha[0], c=Temp[0], cmap="RdYlBu", vmin=Tmin, vmax=Tmax)
   fig5.tight_layout()
   def anim5(i):
      axA5.set_title("$t =$ {:.2f} Myr".format(T[i]))
      sctA5.set_offsets(np.array([Temp[i], Lumi[i]]).T)
      sctA5.set_array(Temp[i])
   ani5 = animation.FuncAnimation(fig=figA5, func=anim5, frames=N_iter, interval=100)
   if diagram_save:
      ani5.save(out_directory+date_prefix+"_animation_5.png", dpi=300, fps=10)
      ani5.save(out_directory+date_prefix+"_animation_5.pdf", dpi=300, fps=10)
      print("\033[94m"+"Info: Animation saved as " + "\033[0m" + date_prefix+"_animation_5")

# * SHOW !
plt.show(block=False)
input("(press enter to quit)")
