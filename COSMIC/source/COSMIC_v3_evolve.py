"""
* COSMIC - EVOLVE
* Version 3
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (Intern)
"""

import numpy as np
import amuse.units.units as u
import COSMIC_v3_output as output
import COSMIC_v3_emission as emission
import COSMIC_v3_coordinates as coords

def stellar_gravity(codes, stars, channels, params, save_stars = False):
    """
    Evolve a model with both gravity and stellar evolution.
    """
    for t in params["output_times"]:
        print("\033[36m"+"t: "+"\033[0m"+"{:02.2f} Myr/{:02.2f} Myr".format(t.value_in(u.Myr),np.max(params["output_times"].value_in(u.Myr))))
        codes["g"].evolve_model(t)
        channels.copy()
        if save_stars: output.save_stars(stars, params, t.value_in(u.Myr))
	 
    print("\033[36m"+"Simulation completed"+"\033[0m")
    return 0

def stellar_gravity_binaries(codes, stars, binaries, channels, params, save_stars = False, save_binaries = False, compute_X_emission = False):
    """
    Evolve a model with both gravity and stellar evolution.
    """
    for t in params["output_times"]:
        print("\033[36m"+"t: "+"\033[0m"+"{:02.2f} Myr/{:02.2f} Myr".format(t.value_in(u.Myr),np.max(params["output_times"].value_in(u.Myr))))
        codes["b"].evolve_model(t)
        channels.copy()
        stars = coords.xyz2radecdist(stars)
        if compute_X_emission: 
            stars = emission.rotation(stars, params)
            stars = emission.X_emission(stars, params)
            stars = emission.plasma_temperature(stars, params)
        if save_stars: output.save_stars(stars, params, t.value_in(u.Myr))
        if save_binaries: output.save_binaries(binaries, params, t.value_in(u.Myr))
	 
    print("\033[36m"+"Simulation completed"+"\033[0m")
    return codes, stars, channels


