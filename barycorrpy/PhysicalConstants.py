import astropy.constants as ac
import numpy as np


# Parsing constants #
AU = ac.au.value # [m]
c = ac.c.value # Speed of light [m/s]
pctoau = 3600*180/np.pi # No. of AU in one parsec
SECS_PER_DAY = 3600 * 24
year = 365.25*SECS_PER_DAY # [s]
kmstoauyr = 1000 * year/(AU)

M_moon = 73476730924573500000000 # Mass of the Moon [kg]

# Mass and normalised mass of solar system bodies
ss_bodies = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Moon']
M = dict(zip(ss_bodies, [ac.M_sun.value, 0.3301e24, 4.867e24, ac.M_earth.value, 0.6417e24, ac.M_jup.value, 568.5e24, 86.82e24, 102.4e24, M_moon])) # [kg]

#Reduced number of bodies, to increase speed. Does not affect precision at 1 cm/s level
ss_bodies = ['Sun','Earth','Jupiter','Saturn']
M = dict(zip(ss_bodies, [ac.M_sun.value, ac.M_earth.value,ac.M_jup.value,568.5e24])) # [kg]

GM = {k:ac.G.value*M[k] for k in ss_bodies}
