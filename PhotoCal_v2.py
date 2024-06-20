"""
* PhotoCal: Extract the instrumental magnitude coefficients and perform photometry calibration
* Version 2 - 2024-03-21
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (Intern)
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

if "YII_light" in plt.style.available: plt.style.use("YII_light")

Simbad_directory = "./Simbad_files/"
Joint_directory = "./Joint_files/"
PhotoCal_directory = "./PhotoCal_files/"
fig_directory = "../Figs/PhotoCal_figs/"
Joint_prefix = "Joint_"
Simbad_prefix = "Simbad_"
PhotoCal_prefix = "PhotoCal_"
BVR_suffix = "_BVR"
filtered_suffix = "_filtered"
file_ext = ".csv"

precision = 1/3600 # °

YES = ["Y", "YES", "1", "OUI"]

# * Functions

def calibration_function(I,c): return I + c 

# * Import

Joint_filelist = np.sort(os.listdir(Joint_directory))
Joint_filelist = [f for f in Joint_filelist if Joint_prefix in f]
filtered_filelist = np.sort(os.listdir(Joint_directory))
filtered_filelist = [f for f in filtered_filelist if (Joint_prefix in f and filtered_suffix in f)]
Simbad_filelist = np.sort(os.listdir(Simbad_directory))
Simbad_filelist = [f for f in Simbad_filelist if Simbad_prefix in f]
Joint_sources = []
filtered_sources = []
Simbad_sources = []

for i in range(len(Joint_filelist)):
   k = Joint_filelist[i][len(Joint_prefix):].find("_BVR")
   source = Joint_filelist[i][len(Joint_prefix):k+len(Joint_prefix)]
   if source not in Joint_sources: Joint_sources.append(source)
for i in range(len(filtered_filelist)):
   k = filtered_filelist[i][len(Joint_prefix):].find("_BVR")
   source = filtered_filelist[i][len(Joint_prefix):k+len(Joint_prefix)]
   if source not in filtered_sources: filtered_sources.append(source)
for i in range(len(Simbad_filelist)):
   k = Simbad_filelist[i][len(Simbad_prefix):].find("_BVR")
   source = Simbad_filelist[i][len(Simbad_prefix):k+len(Simbad_prefix)]
   if source not in Simbad_sources: Simbad_sources.append(source)
sources = np.intersect1d(Joint_sources, Simbad_sources)

print("\033[36m"+"Available source:"+"\033[0m")
for i in range(len(sources)):
   if sources[i] in filtered_sources:  print("\033[36m"+"\t{}: ".format(i+1) + "\033[0m" + "{} (filtered)".format(sources[i]))
   else: print("\033[36m"+"\t{}: ".format(i+1) + "\033[0m" + "{}".format(sources[i]))
k = int(input("\033[32m"+"Select source number: "+"\033[0m"))-1
print("")
source = sources[k]

i_ra = 0
i_de = 1
i_B_inst = 2
i_V_inst = 3
i_R_inst = 4
i_dB_inst = 5
i_dV_inst = 6
i_dR_inst = 7
i_B_ref = 37
i_V_ref = 38
i_R_ref = 39

cols = [i_ra, i_de, i_B_inst, i_V_inst, i_R_inst, i_dB_inst, i_dV_inst, i_dR_inst, i_B_ref, i_V_ref, i_R_ref]


Simbad_file = np.genfromtxt(Simbad_directory+Simbad_prefix+source+BVR_suffix+file_ext, skip_header=1, comments="#", delimiter=",", usecols=cols)
if source in filtered_sources: Joint_file = np.loadtxt(Joint_directory+Joint_prefix+source+BVR_suffix+filtered_suffix+file_ext, skiprows=1, comments="#", delimiter=",")
else: Joint_file = np.loadtxt(Joint_directory+Joint_prefix+source+BVR_suffix+file_ext, skiprows=1, comments="#", delimiter=",")

# * Calibration
w_B = np.intersect1d(np.argwhere(np.isfinite(Simbad_file[:,8])), np.argwhere(Simbad_file[:,2]!=99))
w_V = np.intersect1d(np.argwhere(np.isfinite(Simbad_file[:,9])), np.argwhere(Simbad_file[:,3]!=99))
w_R = np.intersect1d(np.argwhere(np.isfinite(Simbad_file[:,10])), np.argwhere(Simbad_file[:,4]!=99))

c_B, dc_B = curve_fit(calibration_function, Simbad_file[w_B,2].flatten(), Simbad_file[w_B,8].flatten())
c_V, dc_V = curve_fit(calibration_function, Simbad_file[w_V,3].flatten(), Simbad_file[w_V,9].flatten())
c_R, dc_R = curve_fit(calibration_function, Simbad_file[w_R,4].flatten(), Simbad_file[w_R,10].flatten())

dc_B = np.sqrt(np.diag(dc_B))
dc_V = np.sqrt(np.diag(dc_V))
dc_R = np.sqrt(np.diag(dc_R))

Simbad_file[:,2] += c_B[0]
Simbad_file[:,3] += c_V[0]
Simbad_file[:,4] += c_R[0]

Simbad_file[:,5] += dc_B[0]
Simbad_file[:,6] += dc_V[0]
Simbad_file[:,7] += dc_R[0]

Joint_file[:,2] += c_B[0]
Joint_file[:,3] += c_V[0]
Joint_file[:,4] += c_R[0]

Joint_file[:,5] += dc_B[0]
Joint_file[:,6] += dc_V[0]
Joint_file[:,7] += dc_R[0]

# * Output choices

plot_calibration = input("\033[32m"+"Calibration plot (y/n): "+"\033[0m").upper() in YES
plot_magnitude =  input("\033[32m"+"Show magnitude distribution (y/n): "+"\033[0m").upper() in YES
plot_save = input("\033[32m"+"Save figures (y/n): "+"\033[0m").upper() in YES
print("")
save_output = input("\033[32m"+"Save output file (y/n): "+"\033[0m").upper() in YES
print("")

# * Calibration plot
tab_mag = np.concatenate([Simbad_file[w_B,2].flatten(),Simbad_file[w_B,8].flatten(),Simbad_file[w_V,3].flatten(),Simbad_file[w_V,9].flatten(),Simbad_file[w_R,3].flatten(),Simbad_file[w_R,10].flatten()])

if plot_calibration:
   fig, ax = plt.subplots(1)
   ax.errorbar(Simbad_file[w_B,8].flatten(), Simbad_file[w_B,2].flatten(), Simbad_file[w_B,5].flatten(), markersize=2, color="#0000FF", linestyle="", marker="o", label='B filter', alpha=0.5)
   ax.errorbar(Simbad_file[w_V,9].flatten(), Simbad_file[w_V,3].flatten(), Simbad_file[w_V,6].flatten(), markersize=2, color="#00FF00", linestyle="", marker="s", label="V filter", alpha=0.5)
   ax.errorbar(Simbad_file[w_R,10].flatten(), Simbad_file[w_R,4].flatten(), Simbad_file[w_R,7].flatten(), markersize=2, color="#FF0000", linestyle="", marker="^", label="R filter", alpha=0.5)
   
   min_CAL = np.min(tab_mag)
   max_CAL = np.max(tab_mag)
   ax.plot([min_CAL-0.5,max_CAL+0.5], [min_CAL-0.5,max_CAL+0.5], "--", alpha=0.5)
   ax.set_xlabel("Reference magnitude")
   ax.set_ylabel("Calibrated magnitude")
   ax.legend()
   fig.tight_layout()
   if plot_save: 
      fig.savefig(fig_directory+source+"_calibration_plot")
      print("\033[94m"+"Info: Figure saved as " + "\033[0m" + source+"_calibration_plot")

# * Magnitude histograms

if plot_magnitude:
   fig, ax = plt.subplots(1)
   bins = np.linspace(np.quantile(tab_mag, 0.1), np.quantile(tab_mag,0.9), 20)
   hist_B, bins_B = np.histogram(Joint_file[:,2],bins)
   hist_V, bins_V = np.histogram(Joint_file[:,3],bins)
   hist_R, bins_R = np.histogram(Joint_file[:,4],bins)

   cbins_B = 0.5*(bins_B[1:] + bins_B[:-1])
   cbins_V = 0.5*(bins_V[1:] + bins_V[:-1])
   cbins_R = 0.5*(bins_R[1:] + bins_R[:-1])

   ax.scatter(cbins_B, hist_B, color="#0000FF", label="B filter", marker="o", alpha=0.5)
   ax.scatter(cbins_V, hist_V, color="#00FF00", label="V filter", marker="s", alpha=0.5)
   ax.scatter(cbins_R, hist_R, color="#FF0000", label="R filter", marker="^", alpha=0.5)

   ax.set_xlabel("Relative magnitude")
   ax.set_ylabel("Star distribution")
   ax.legend()
   fig.tight_layout()
   if plot_save: 
      fig.savefig(fig_directory+source+"_magnitude_hist")
      print("\033[94m"+"Info: Figure saved as " + "\033[0m" + source+"_magnitude_hist")

if save_output:
   file = open(PhotoCal_directory+PhotoCal_prefix+source+BVR_suffix+file_ext, "w")
   file.write("ra [deg],dec [deg],B_obs [mag],V_obs [mag],R_obs [mag],dB_obs [mag],dV_obs [mag],dR_obs [mag]\n")
   for i in range(len(Joint_file)):
      ra, dec, B, V, R, dB, dV, dR = Joint_file[i,:]
      file.write("{ra:e},{dec:e},{B:e},{V:e},{R:e},{dB:e},{dV:e},{dR:e}\n".format(ra=ra, dec=dec, B=B, V=V, R=R, dB=dB, dV=dV, dR=dR))
   print("\033[94m"+"Info: File saved as " + "\033[0m" + PhotoCal_prefix+source+BVR_suffix+file_ext)
   file.close()
plt.show()

