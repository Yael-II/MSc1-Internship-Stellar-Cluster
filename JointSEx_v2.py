"""
* JointSEx : Joint SExtractor files
* Version 2 - 2024-03-21
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (Intern)
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import FuncFormatter, MaxNLocator

# * Parameters

if "YII_light" in plt.style.available: plt.style.use("YII_light")

# * Variables
in_directory = "./SEx_files/"
out_directory = "./Joint_files/"
fig_directory = "../Figs/JointSEx_figs/"
SEx_prefix = "S-ex "
B_suffix = "-B_cal.cat"
V_suffix = "-G_cal.cat"
R_suffix = "-R_cal.cat"
BVR_suffix = "_BVR.csv"
out_prefix = "Tab_"
precision = 1/3600 # °

YES = ["Y", "YES", "1", "OUI"]

i_ra = 13
i_de = 14
i_mag = 1
i_dmag = 2
# MAG_ISO:  1, 2
# MAG_BEST: 9, 10
# TODO APERTURE

# * Functions
def calibration_function(I,c): return I + c 

# deg = FuncFormatter(lambda x, pos=None: f"{x:.0f}\N{DEGREE SIGN}{(x-np.floor(x))*60:02.0f}\'00\"")
# hrs = FuncFormatter(lambda x, pos=None: f"{x*24/360:.00f}:{(x*24/360-np.floor(x*24/360))*60:02.0f}:00")

def DMS_from_deg_format(x, pos=None):
   D = int(np.floor(x))
   M = int(np.floor((x-D)*60))
   S = int(np.floor(((x-D)*60-M)*60))
   return "{D:02d}°{M:02d}\'{S:02d}\"".format(D=D,M=M,S=S)
def HMS_from_deg_format(x, pos=None):
   y = x/15
   H = int(np.floor(y))
   M = int(np.floor((y-H)*60))
   S = int(np.floor(((y-H)*60-M)*60))
   return "{H:02d}:{M:02d}:{S:02d}".format(H=H,M=M,S=S)

DMS_from_deg = FuncFormatter(DMS_from_deg_format)
HMS_from_deg = FuncFormatter(HMS_from_deg_format)

def openfile(in_dir, SEx_prefix, source, suffix):
   tab = []
   with open(in_dir+SEx_prefix+source+suffix) as file:
      for i in range(20):
         file.readline()
      for line in file:
         columns = line.strip().split()
         tab.append([float(columns[i_ra]), float(columns[i_de]), float(columns[i_mag]), float(columns[i_dmag])])
      file.close()
   return np.array(tab)

# * Import 
filelist = np.sort(os.listdir(in_directory))
filelist = [f for f in filelist if SEx_prefix in f]
sources = []
for i in range(len(filelist)):
   k = filelist[i][5:].find("-")
   source = filelist[i][5:k+5]
   if source not in sources: sources.append(source)

print("\033[36m"+"Available source:"+"\033[0m")
for i in range(len(sources)):
   print("\033[36m"+"\t{}: ".format(i+1) + "\033[0m" + "{}".format(sources[i]))
k = int(input("\033[32m"+"Select source number: "+"\033[0m"))-1
print("")
source = sources[k]

tab_B = openfile(in_directory, SEx_prefix, source, B_suffix)
tab_V = openfile(in_directory, SEx_prefix, source, V_suffix)
tab_R = openfile(in_directory, SEx_prefix, source, R_suffix)

w1 = np.where(tab_B[:,0] -  tab_V[:,0] < precision)
w2 = np.where(tab_V[:,0] -  tab_R[:,0] < precision)
w3 = np.where(tab_B[:,0] -  tab_R[:,0] < precision)


if not (len(w1[0]) == len(w2[0]) == len(w3[0]) == len(tab_B) == len(tab_V) == len(tab_R)):
   print("\033[33m"+"Warning: not all tables have the same size, or a star may have moved !"+ "\033[0m")
   print("\033[33m"+"\tB table size:"+ "\033[0m"+"{:d}".format(len(tab_B)))
   print("\033[33m"+"\tV table size:"+ "\033[0m"+"{:d}".format(len(tab_V)))
   print("\033[33m"+"\tR table size:"+ "\033[0m"+"{:d}".format(len(tab_R)))
   print("\033[33m"+"\tSame position in B-V:"+ "\033[0m"+"{:d}".format(len(w1[0])))
   print("\033[33m"+"\tSame position in V-R:"+ "\033[0m"+"{:d}".format(len(w2[0])))
   print("\033[33m"+"\tSame position in B-R:"+ "\033[0m"+"{:d}".format(len(w3[0])))
else:
   N = len(tab_R)
   tab_BVR = [] # ra, dec, B, V, R, dB, dV, dR
   for i in range(N):
      ra = np.mean([tab_B[i,0], tab_V[i,0], tab_R[i,0]])
      de = np.mean([tab_B[i,1], tab_V[i,1], tab_R[i,1]])
      B = tab_B[i,2]
      V = tab_V[i,2]
      R = tab_R[i,2]
      dB = tab_B[i,3]
      dV = tab_V[i,3]
      dR = tab_R[i,3]
      tab_BVR.append([ra, de, B, V, R, dB, dV, dR])
   tab_BVR = np.array(tab_BVR)

# * Output choices

plot_map = input("\033[32m"+"Draw map (y/n): "+"\033[0m").upper() in YES
plot_save = input("\033[32m"+"Save figures (y/n): "+"\033[0m").upper() in YES
print("")
save_output = input("\033[32m"+"Save output file (y/n): "+"\033[0m").upper() in YES
print("")

# * Map
if plot_map:
   fig, ax = plt.subplots(1)
   ax.scatter(tab_BVR[:,0], tab_BVR[:,1], s=2, color='C3', alpha=1)
   ax.xaxis.set_major_formatter(HMS_from_deg)
   ax.xaxis.set_major_locator(MaxNLocator(5))
   ax.yaxis.set_major_formatter(DMS_from_deg)
   ax.set_aspect("equal")
   ax.set_xlabel("Right ascension")
   ax.set_ylabel("Declination")
   if plot_save: 
      fig.savefig(fig_directory+source+"_map")
      print("\033[94m"+"Info: Figure saved as " + "\033[0m" + source+"_map")

# * File output
if save_output:
   file = open(out_directory+out_prefix+source+BVR_suffix, "w")
   file.write("ra [deg],dec [deg],B_inst [mag],V_inst [mag],R_inst [mag],dB_inst [mag],dV_inst [mag],dR_inst [mag]\n")
   for i in range(len(tab_BVR)):
      ra, dec, B, V, R, dB, dV, dR = tab_BVR[i,:]
      file.write("{ra:e},{dec:e},{B:e},{V:e},{R:e},{dB:e},{dV:e},{dR:e}\n".format(ra=ra, dec=dec, B=B, V=V, R=R, dB=dB, dV=dV, dR=dR))
   print("\033[94m"+"Info: File saved as " + "\033[0m" + out_prefix+source+BVR_suffix)
   file.close()
plt.show()


