"""
* COSMIC - SAVE
* Version 3
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (Intern)
"""


import numpy as np

from amuse.io.base import write_set_to_file

def save_stars(stars, params, t=0):
    write_set_to_file(set=stars, filename=params["out_directory"]+params["filename"]+"_stars_{:04d}_Myr.{}".format(np.int32(t), params["format_type"]), format=params["format_type"])
    return 0

def save_binaries(binaries, params, t=0):
    write_set_to_file(set=binaries, filename=params["out_directory"]+params["filename"]+"_binaries_{:04d}_Myr.{}".format(np.int32(t), params["format_type"]), format=params["format_type"])
    return 0
