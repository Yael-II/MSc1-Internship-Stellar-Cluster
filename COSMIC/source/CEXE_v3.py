"""
* CEXE: COSMIC Extension for X-ray Emission
* Version 3
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (Intern)
"""

import numpy as np
import os 
import subprocess as sp
from astropy import units as u

directory = "./output/"

# SIMPUT (sources)
val_Elow = 0.1
val_Eup = 15
val_NBins = 1000
val_logEgrid = "yes"
val_Emin = 2
val_Emax = 10
val_nH = 0.2
val_flux_min = 2
val_flux_max = 10
# SIXTE (instrument)
val_ra = 0 
val_dec = 0
val_Exposure = 1000
XML_file = "TODO" # TODO !!

# * COSMIC
def open_file_terminal():
    filelist = np.sort([f for f in os.listdir(directory) if "stars" in f])
    print("\033[36m"+"Available files:"+"\033[0m")
    for i in range(len(filelist)):
        print("\033[36m"+"\t{}: ".format(i+1) + "\033[0m" + filelist[i])
    j = int(input("\033[32m"+"Select a file number: "+"\033[0m"))-1
    filename = filelist[j]
    date = filename[filename.find("COSMIC_")+len("COSMIC_") : filename.find("_stars")]
    time = filename[filename.find("stars_")+len("stars_") : filename.find("_Myr")]
    print("")
    file = np.loadtxt(directory+filename, delimiter=",", comments="#")
    with open(directory+filename) as f: head = f.readline().strip('\n')[1:]
    return file, head, date, time

# * XSpec
def xspec_build(date, time):
    with open(directory+"CEXE_XSPEC_{}_{}_Myr_0.6_0.0.sh".format(date, time), "w+") as xspec_file_0600:
        xspec_file_0600.write("model phabs*TODO\n")
        xspec_file_0600.write("phabs:nH>{}\n".format(val_nH))
        xspec_file_0600.write("TODO:...\n")
        xspec_file_0600.write("flux {} {}\n".format(val_flux_min, val_flux_max))
        xspec_file_0600.write("save model {}XSPEC_{}_{}_Myr_0.6_0.0.xcm\n".format(directory, date, time))
        xspec_file_0600.write("quit")
    with open(directory+"CEXE_XSPEC_{}_{}_Myr_0.5_0.0.sh".format(date, time), "w+") as xspec_file_0500:
        xspec_file_0500.write("model phabs*TODO\n")
        xspec_file_0500.write("phabs:nH>{}\n".format(val_nH))
        xspec_file_0500.write("TODO:...\n")
        xspec_file_0500.write("flux {} {}\n".format(val_flux_min, val_flux_max))
        xspec_file_0500.write("save model {}XSPEC_{}_{}_Myr_0.5_0.0.xcm\n".format(directory, date, time))
        xspec_file_0500.write("quit")
    with open(directory+"CEXE_XSPEC_{}_{}_Myr_0.4_1.0.sh".format(date, time), "w+") as xspec_file_0410:
        xspec_file_0410.write("model phabs*TODO\n")
        xspec_file_0410.write("phabs:nH>{}\n".format(val_nH))
        xspec_file_0410.write("TODO:...\n")
        xspec_file_0410.write("flux {} {}\n".format(val_flux_min, val_flux_max))
        xspec_file_0410.write("save model {}XSPEC_{}_{}_Myr_0.4_1.0.xcm\n".format(directory, date, time))
        xspec_file_0410.write("quit")
    with open(directory+"CEXE_XSPEC_{}_{}_Myr_0.2_0.8.sh".format(date, time), "w+") as xspec_file_0208:
        xspec_file_0208.write("model phabs*TODO\n")
        xspec_file_0208.write("phabs:nH>{}\n".format(val_nH))
        xspec_file_0208.write("TODO:...\n")
        xspec_file_0208.write("flux {} {}\n".format(val_flux_min, val_flux_max))
        xspec_file_0208.write("save model {}XSPEC_{}_{}_Myr_0.2_0.8.xcm\n".format(directory, date, time))
        xspec_file_0208.write("quit")
    return None

def xspec_exec(date, time):
    sp.run("echo \"xspec < {}CEXE_XSPEC_{}_{}_Myr_0.6_0.0.sh - TO EXECUTE\"".format(directory, date, time), shell=True)
    sp.run("echo \"xspec < {}CEXE_XSPEC_{}_{}_Myr_0.5_0.0.sh - TO EXECUTE\"".format(directory, date, time), shell=True)
    sp.run("echo \"xspec < {}CEXE_XSPEC_{}_{}_Myr_0.4_1.0.sh - TO EXECUTE\"".format(directory, date, time), shell=True)
    sp.run("echo \"xspec < {}CEXE_XSPEC_{}_{}_Myr_0.2_0.8.sh - TO EXECUTE\"".format(directory, date, time), shell=True)
    return None

# * SIMPUT
def simput_build(file, col, date, time):
    with open(directory+"CEXE_SIMPUT_{}_{}_Myr.sh".format(date, time), "w+") as simput_file:
        Infiles = ""
        for i in range(len(file)):
            star_flux = file[i,col["X_luminosity"]] / (4*np.pi*(file[i,col["dist"]] * u.pc).to(u.cm).value**2)
            star_ra = file[i,col["ra"]]
            star_dec = file[i,col["dec"]]
            star_T_0 = file[i,col["X_temperature_0"]]
            star_T_1 = file[i,col["X_temperature_1"]]

            Infiles += "{}SIMPUT_{}_{}_Myr_{}.fits,".format(directory, date, time, i) 
            
            simput_file.write("$SIXTE/bin/simputfile Simput={}SIMPUT_{}_{}_Myr_{}.fits ".format(directory, date, time, i))
            simput_file.write("Src_Name = star_{} ".format(i))
            simput_file.write("RA = {} Dec = {} ".format(star_ra, star_dec))
            simput_file.write("srcFlux = {} ".format(star_flux))
            simput_file.write("Elow={} Eup={} NBins={} logEgrid={} Emin={} Emax={} ".format(val_Elow, val_Eup, val_NBins, val_logEgrid, val_Emin, val_Emax))
            
            simput_file.write("XSPECFile={}XSPEC_{}_{}_Myr_{:2.01f}_{:2.01f}.xcm".format(directory, date, time, star_T_0, star_T_1))
            simput_file.write("\n")
        simput_file.write("$SIXTE/bin/simputmerge ")
        simput_file.write("Infiles={} ".format(Infiles[:-1]))
        simput_file.write("Outfile={}SIMPUT_{}_{}_Myr.fits".format(directory, date, time))
        simput_file.write("\n")
    return 0

def simput_exec(date, time):
    sp.run(["chmod", "a+x", "{}CEXE_SIMPUT_{}_{}_Myr.sh".format(directory, date, time)])
    sp.run(["echo" ,"{}CEXE_SIMPUT_{}_{}_Myr.sh - TO EXECUTE".format(directory, date, time)])
    return 0

def sixte_build(date, time):
    with open(directory+"CEXE_SIXTE_{}_{}_Myr.sh".format(date, time), "w+") as sixte_file:
        sixte_file.write("$SIXTE/bin/runsixt ")
        sixte_file.write("XMLFile = {} ".format(XML_file))
        sixte_file.write("RA = {} DEC = {} ".format(val_ra, val_dec))
        # sixte_file.write("Prefix = SIXTE_ ")
        sixte_file.write("Simput = {}SIMPUT_{}_{}_Myr.fits ".format(directory, date, time))
        sixte_file.write("EvtFile = {}SIXTE_Event_{}_{}_Myr.fits ".format(directory, date, time))
        sixte_file.write("Exposure = val_Exposure")
        sixte_file.write("\n")
    return 0

def sixte_exec(date, time):
    sp.run(["chmod", "a+x", "{}CEXE_SIXTE_{}_{}_Myr.sh".format(directory, date, time)])
    sp.run(["echo", "{}CEXE_SIXTE_{}_{}_Myr.sh - TO EXECUTE".format(directory, date, time)])
    return 0

file, head, date, time = open_file_terminal()
col = {head.split(",")[i]:i for i in range(len(head.split(",")))}

xspec_build(date, time)
xspec_exec(date, time)
simput_build(file, col, date, time)
simput_exec(date, time)
sixte_build(date, time)
sixte_exec(date, time)
