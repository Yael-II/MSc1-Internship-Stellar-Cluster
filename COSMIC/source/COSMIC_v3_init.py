"""
* COSMIC - INIT
* Version 3 
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (Intern)
"""

import numpy as np
import amuse.units.units as u

import COSMIC_v3_galaxy as galaxy
import COSMIC_v3_coordinates as coords

from amuse.units import nbody_system
from amuse.ic.kingmodel import new_king_model
from amuse.ic.salpeter import new_salpeter_mass_distribution
# from amuse.community.huayno.interface import Huayno # FIXME DOES NOT WORK
from amuse.community.ph4.interface import Ph4
from amuse.community.seba.interface import SeBa
from amuse.couple.bridge import Bridge
from amuse.datamodel.particles import Channels
from amuse.datamodel import Particles

# def IMF_Besancon(M, alpha_1=1.6, alpha_2=3.0):
#     return (M.value_in(u.MSun))**(alpha_1-1) * ((M.value_in(u.MSun))<1) + (((M.value_in(u.MSun))-1)**(alpha_2-1) + 1) * ((M.value_in(u.MSun))>=1)

# def Besancon_disk(params):
#     """
#     Generate two IMF of stars and binaries
#     """
#     M_1 = []
#     M_2 = []
#     N_min = 0
#     N_max = IMF_Besancon(params["M_max_besancon"], params["a_besancon_1"], params["a_besancon_2"])
#     while len(M_1) < params["N_stars"]:
#         X = np.random.uniform(params["M_min_besancon"].value_in(u.MSun), params["M_max_besancon"].value_in(u.MSun))|u.MSun
#         Y = np.random.uniform(N_min, N_max)
#         if Y > IMF_Besancon(X, params["a_besancon_1"], params["a_besancon_2"]):
#             M_1.append(X)
#     while len(M_2) < params["N_stars"] * params["binary_fraction"]:
#         X = np.random.uniform(params["M_min_besancon"].value_in(u.MSun), params["M_max_besancon"].value_in(u.MSun))|u.MSun
#         Y = np.random.uniform(N_min, N_max)
#         if Y > IMF_Besancon(X, params["a_besancon_1"], params["a_besancon_2"]):
#             M_2.append(X)
#     return M_1, M_2

def orbital_parameters(params, m_tot):
    x_0 = params["mean_period"] * np.log(2)
    sigma_0 = params["std_period"] * np.log(2)
    T = np.random.lognormal(mean = x_0, sigma = sigma_0)
    a = (T**2 * u.constants.G.value_in(u.au**3*u.day**(-2)*u.MSun**(-1)) * m_tot.value_in(u.MSun) / (4*np.pi**2) )**(1/3) | u.au
    e = 1
    return a, e

def king_salpeter(params):
    """
    Create a cluster of N stars, with a King's law of parameter W0, with a (virial) radius R, stellar masses following a Salpeter law ranging between M_min and M_max with exponent a_salpeter.
    """
    M_stars = new_salpeter_mass_distribution(params["N_stars"], mass_min=params["M_min_salpeter"], mass_max=params["M_max_salpeter"], alpha=params["a_salpeter"])
    converter = nbody_system.nbody_to_si(M_stars.sum(), params["R0_king"])
    stars = new_king_model(params["N_stars"], params["W0_king"], converter)
    stars.mass = M_stars
    stars.zams_mass = M_stars
    return stars, converter

def king_salpeter_binaries(params):
    """
    
    """
    M_primary = new_salpeter_mass_distribution(params["N_stars"], mass_min=params["M_min_salpeter"], mass_max=params["M_max_salpeter"], alpha=params["a_salpeter"])
    M_secondary = new_salpeter_mass_distribution(np.int32(np.round(params["N_stars"]*params["binary_fraction"])), mass_min=params["M_min_salpeter"], mass_max=params["M_max_salpeter"], alpha=params["a_salpeter"])

    M_stars = np.concatenate([M_primary, M_secondary])
    N_binaries = len(M_secondary)
    converter = nbody_system.nbody_to_si(M_stars.sum(), params["R0_king"])
    stars_primary = new_king_model(params["N_stars"], params["W0_king"], converter)
    stars_primary.mass = M_primary
    stars_primary.zams_mass = M_primary
     
    stars_secondary = Particles(N_binaries)
    for i in range(N_binaries):
        stars_secondary[i].mass = M_secondary[i]
        stars_secondary[i].radius = 0|u.m
        stars_secondary[i].vx = stars_primary[i].vx
        stars_secondary[i].vy = stars_primary[i].vy
        stars_secondary[i].vz = stars_primary[i].vz
        stars_secondary[i].x = stars_primary[i].x
        stars_secondary[i].y = stars_primary[i].y
        stars_secondary[i].z = stars_primary[i].z
        stars_secondary[i].zams_mass = M_secondary[i]

    stars = Particles()
    binaries = Particles(N_binaries)

    stars.add_particles(stars_primary)
    stars.add_particles(stars_secondary)

    stars.ra = 0 | u.deg
    stars.dec = 0 | u.deg
    stars.dist = 0 | u.pc
    stars.X_v = 0 | u.km * u.s**(-1)
    stars.X_vsini = 0 | u.km * u.s**(-1)
    stars.X_temperature_0 = 0 | u.kilo(u.eV)
    stars.X_temperature_1 = 0 | u.kilo(u.eV)
    stars.X_luminosity = 0 | u.erg * u.s**(-1)

    if N_binaries == 0:
        binaries.child1 = [0]
        binaries.child2 = [0]
        binaries.mass1 = 0
        binaries.mass2 = 0
        binaries.mass_tot = 0
    else:
        binaries.child1 = list(stars_primary[:N_binaries])
        binaries.child2 = list(stars_secondary)
        binaries.mass1 = list(stars_primary[:N_binaries].mass)
        binaries.mass2 = list(stars_secondary.mass.in_(u.kg))
        binaries.mass_tot = binaries.mass1 + binaries.mass2

    for i in range(N_binaries):
        semi_major_axis, eccentricity = orbital_parameters(params, binaries[i].mass_tot)
        binaries[i].semi_major_axis = semi_major_axis
        binaries[i].eccentricity = eccentricity

    stars.move_to_center()
    return stars, binaries, converter


def gravity_stellar_bridge(converter, params):
    gravity = Ph4(converter, number_of_workers=params["workers_gravity"]) # TODO A CHANGER SI BESOIN
    stellar = SeBa(number_of_worker=params["workers_stellar"])
    bridge = Bridge()
    codes = {"g":gravity, "s":stellar, "b":bridge}
    return codes

def create_channels():
    channels = Channels()
    return channels

def generate_galaxy(stars, params):
    galaxy_model = galaxy.MilkyWay_AMUSE()
    stars.x += params["cluster_position_x"]
    stars.y += params["cluster_position_y"]
    stars.z += params["cluster_position_z"]
    vc = galaxy_model.vel_circ(stars.center_of_mass().length()).in_(u.km * u.s**(-1))
    phi = np.arctan2(stars.center_of_mass().y.value_in(u.pc), stars.center_of_mass().x.value_in(u.pc))
    stars.vx += - vc * np.sin(phi)
    stars.vy += vc * np.cos(phi)
    stars = coords.xyz2radecdist(stars)
    return stars, galaxy_model

def add_stars(codes, stars):
    codes["g"].particles.add_particles(stars)
    codes["s"].particles.add_particles(stars)
    return codes

def add_binaries(codes, stars, binaries, channels):
    codes["s"].binaries.add_particles(binaries)
    channels.add_channel(codes["s"].binaries.new_channel_to(stars))
    return codes, channels

def commit(codes, stars, channels, params):
    codes["s"].set_metallicity(params["metallicity_stellar"])
    codes["s"].commit_particles()
    codes["g"].commit_particles()
    # codes["g"].initialize_code() # NOT REQUIRED WITH Ph4

    codes["b"].add_system(codes["g"])
    codes["b"].add_system(codes["s"])

    codes["b"].synchronize_model()
    codes["b"].timestep = params["timestep"]

    codes["b"].channels.add_channel(codes["s"].particles.new_channel_to(codes["g"].particles, attributes=["mass", "radius"])) 

    channels.add_channel(codes["s"].particles.new_channel_to(stars))
    channels.add_channel(codes["g"].particles.new_channel_to(stars))
    return codes, channels

def commit_with_potential(codes, stars, channels, galaxy_model, params):
    codes["s"].set_metallicity(params["metallicity_stellar"])
    codes["s"].commit_particles()
    codes["g"].commit_particles()
    # codes["g"].initialize_code() # NOT REQUIRED WITH Ph4

    codes["b"].add_system(codes["g"], (galaxy_model,))
    codes["b"].add_system(codes["s"])

    codes["b"].synchronize_model()
    codes["b"].timestep = params["timestep"]

    codes["b"].channels.add_channel(codes["s"].particles.new_channel_to(codes["g"].particles, attributes=["mass", "radius"])) 

    channels.add_channel(codes["s"].particles.new_channel_to(stars))
    channels.add_channel(codes["g"].particles.new_channel_to(stars))
    return codes, channels

