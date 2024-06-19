"""
* COSMIC : Cluster Orbital SysteM Integration Code
* Version 3
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (Intern)
"""

import numpy as np
import amuse.units.units as u
# import matplotlib.pyplot as plt
import COSMIC_v3_config as config
import COSMIC_v3_init as init
import COSMIC_v3_evolve as evolve
# import COSMIC_v3_output as output
# import COSMIC_v3_emission as emission

params, params_list = config.default()
params = config.from_cfg_file(params)
stars, binaries, converter = init.king_salpeter_binaries(params)
# stars, converter = init.king_salpeter(params)
channels = init.create_channels()
codes = init.gravity_stellar_bridge(converter, params)
codes = init.add_stars(codes, stars)
codes, channels = init.add_binaries(codes, stars, binaries, channels)
codes, channels = init.commit(codes, stars, channels, params)
codes, stars, channels = evolve.stellar_gravity_binaries(codes, stars, binaries, channels, params, save_stars = True, save_binaries = True, compute_X_emission = True)
# evolve.stellar_gravity(codes, stars, channels, params, save_stars = True)
