"""
* COSMIC - EMISSION
* Version 3
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (Intern)
"""

import numpy as np
import amuse.units.units as u

def rotation(stars, params):
    N = params["N_stars"] + np.int32(np.round(params["N_stars"] * params["binary_fraction"]))
    stars.X_v = np.random.uniform(0, 10, N) | u.km * u.s**(-1)
    p = np.random.uniform(0, 1, N)
    stars.X_vsini = stars.X_v * np.sqrt(2*p-p**2)
    return stars

def X_emission(stars, params):
    for i in range(len(stars)):
        if stars[i].temperature.value_in(u.K) > 7300: # O B A
            stars[i].X_luminosity = params["X_coefficient_OBA"] * stars[i].luminosity.in_(u.erg * u.s**(-1))
        elif stars[i].temperature.value_in(u.K) < 6000: # G K M
            stars[i].X_luminosity = params["X_coefficient_GKM"] * (stars[i].X_vsini.value_in(u.km * u.s**(-1)))**(1.9) | u.erg * u.s**(-1)
        else: # F
            stars[i].X_luminosity = 0
    return stars
