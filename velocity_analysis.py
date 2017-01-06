import numpy as np
from math import sin, cos, cosh, sinh, pi, sqrt
from scipy.optimize import newton

def stumpC(z):
    '''
    :param z: z = alpha * X^2, where alpha is the reciprocal of the semimajor axis and X is the universal anomaly.
    :return: the value of the stumpff function C(z) at a given Z, note first you need to find the universal anomally.
    '''
    if z > 0:
        return (1- cos(sqrt(z)))/z
    elif z < 0:
        return (cosh(sqrt(-z))-1)/(-z)
    else:
        return 1/2

def stumpS(z):
    '''
    :param z: z = alpha * X^2, where alpha is the reciprocal of the semimajor axis and X is the universal anomaly.
    :return: the value of the stumpff function S(z) at a given Z, note first you need to find the universal anomally.
    '''
    if z > 0:
        return (sqrt(z) - sin(sqrt(z)))/((sqrt(z))**3)
    elif z < 0:
        return (sinh(sqrt(z)) - sqrt(-z))/((sqrt(-z))**3)
    else:
        return 1/6

def kepler_eqn(time, postion, velocity, semiMajorAxis, mu):
    '''
    :param time: time at which you are try to find X the universal anomaly at.
    :param postion: this the magnitude of the initial position vector
    :param velocity: the magnitude of the initial velocity vector
    :param semiMajorAxis: the semi-major axis of your orbits eclipse
    :param mu standard gravity parameter
    :return: X the universal anomaly at the given time.

    this function uses newton's method to solve for universal anomaly from kepler's equation. a iterative method is used
    here, the other options would be to use a graphical option or an approximation. But given that I want to make this
    analysis as general as possible this method seems the best and this can be reused to later analysis because of that.
    '''

    error = 1 * 10**(-8)
    maxIter = 1000
    x = sqrt(mu) * abs(semiMajorAxis) * time    # starting point
    alpha = 1/semiMajorAxis
    n = 0
    ratio = 1
    while abs(ratio) > error and n <= maxIter:
        n += 1
        C = stumpC(alpha*x**2)
        S = stumpS(alpha*x**2)
        F = postion * velocity / sqrt(mu) * x ** 2 * C + (1 - alpha * postion) * x ** 3 * S + postion * x - sqrt(mu) * time
        dFdx = postion * velocity / sqrt(mu) * x * (1 - alpha * x ** 2 * S) + (1 - alpha * postion) * x ** 2 * C + postion
        ratio = F/dFdx
        x = x - ratio

    return x

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
velTan = angMoCombine/radiusInitial      # tangential velocity km/s
posCombine = np.array([[radiusInitial], [0], [0]])  # combine position vector at t = 0
velCombine = np.array([[0], [velTan], [0]])  # combine velocity vector at t = 0

testing = kepler_eqn(250, radiusInitial, velTan, radiusInitial, stdGravPar)
print(testing)