"""
* COSMIC - COORDINATES
* Version 3
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (intern)
"""

import numpy as np
import amuse.units.units as u_amuse
import astropy.units as u_astropy

import astropy.coordinates as c

def xyz2radecdist(stars):
    X = stars.x.value_in(u_amuse.pc)
    Y = stars.y.value_in(u_amuse.pc)    
    Z = stars.z.value_in(u_amuse.pc)
    Galactocentric = c.SkyCoord(x = X * u_astropy.pc, y = Y * u_astropy.pc, z = Z * u_astropy.pc, frame=c.Galactocentric)
    ICRS = Galactocentric.transform_to(c.ICRS)
    RA = ICRS.ra.deg
    DEC = ICRS.dec.deg
    DIST = ICRS.distance.pc
    stars.ra = RA | u_amuse.deg
    stars.dec = DEC | u_amuse.deg
    stars.dist = DIST | u_amuse.pc
    return stars

