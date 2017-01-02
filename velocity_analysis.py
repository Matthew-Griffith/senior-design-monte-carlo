import numpy as np
from math import sin, cos, pi, sqrt
from scipy.optimize import newton

# here we create a var for how many times we want to in the model
samples = 1032

# we need to keep track of time for each the cubesats
timeStep = 60               # secs
timeFinal = 60*60*24*14     # this is two weeks in seconds
timeArray = np.arange(0, timeFinal, timeStep)

# this is just an array of random values for the magnitude of velocity
stdDevMag = (1.2 - 0.8)/4      # this just finding the standard deviation using the range rule.
velMag = (np.random.normal(1.0, stdDevMag, samples))/1000    # this is in km/s
# for this project I have been told that the separation velocity vector has an pointing error of +/- 5 degrees
# based on this the velocity vector will be randomly generated in spherical coordinates and converted into cartesian
stdDevAng = (10*(pi/180))/4
# this might seem a bit odd but I could find a function that just created a half normal or more generally a fold normal
# random distribution in numpy or scipy so I wrote this while loop to create one
phi = np.array([])
while phi.size < samples:
    tempPhi = np.random.normal(0, stdDevAng)
    if tempPhi >= 0:
        phi = np.append(phi, tempPhi)
theta = np.random.uniform(0, 2*pi, samples)
# here is where we convert to cartesian coordinates
# note that the separation is happening in the z direction or normal to the orbit plane.
vx_sep = velMag*np.cos(theta)*np.sin(phi)
vy_sep = velMag*np.sin(theta)*np.sin(phi)
vz_sep = velMag*np.cos(phi)

# setting initial orbit parameters
stdGravPar = 398600             # (km^3)/(s^2)
radiusInitial = 6378 + 500      # km
angMoCombine = sqrt(radiusInitial*stdGravPar)   # angular momentum for initial orbit km^3/s^2
velTanInitial = angMoCombine/radiusInitial      # tangential velocity km/s
posCombine = np.array([[radiusInitial], [0], [0]])  # combine position vector at t = 0
velCombine = np.array([[0], [velTanInitial], [0]])  # combine velocity vector at t = 0
