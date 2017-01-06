import numpy as np
from math import sin, cos, cosh, sinh, pi, sqrt

def stumpC(z):
    '''
    :param z: z = alpha * X^2, where alpha is the reciprocal of the semimajor axis and X is the universal anomaly.
    :return: the value of the stumpff function C(z) at a given Z, note first you need to find the universal anomally.
    '''
    if z > 0:
        return (1 - cos(sqrt(z)))/z
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

def kepler_eqn(time, position, radialVector, semiMajorAxis, mu):
    '''
    :param time: time at which you are try to find X the universal anomaly at.
    :param position: this the magnitude of the initial position vector
    :param radialVector: the magnitude of the initial radial velocity
    :param semiMajorAxis: the semi-major axis of your orbits eclipse
    :param mu: standard gravity parameter
    :return: X the universal anomaly at the given time.

    this function uses newton's method to solve for universal anomaly from kepler's equation. a iterative method is used
    here, the other options would be to use a graphical option or an approximation. But given that I want to make this
    analysis as general as possible this method seems the best and this can be reused to later analysis because of that.
    '''

    error = 1 * 10**(-8)
    maxiter = 1000
    x = sqrt(mu) * abs(semiMajorAxis) * time    # starting point
    alpha = 1/semiMajorAxis
    n = 0
    ratio = 1
    while abs(ratio) > error and n <= maxiter:
        n += 1
        C = stumpC(alpha*x**2)
        S = stumpS(alpha*x**2)
        F = position * radialVector / sqrt(mu) * x ** 2 * C + (1 - alpha * position) * x ** 3 * S + position * x - sqrt(mu) * time
        dFdx = position * radialVector / sqrt(mu) * x * (1 - alpha * x ** 2 * S) + (1 - alpha * position) * x ** 2 * C + position
        ratio = F/dFdx
        x = x - ratio

    return x

def lagrange_coeff(x, time, positionMag, alpha, mu):
    '''
    :param x: universal anomaly at some time
    :param time: time you are trying to find the coeff for
    :param positionMag: magnitude of your position vector
    :param alpha: is the reciprocal of the semimajor axis
    :param mu: is the standard gravitational parameter
    :return: f and g the lagrange coeff
    '''

    z = alpha * x**2
    f = 1 - x**2 / positionMag * stumpC(z)
    g = time - 1 / sqrt(mu) * x**3 * stumpS(z)
    return f, g

def position_from_state_vector(positionVector, velocityVector, time, mu):
    '''
    :param positionVector: this is the initial position vector
    :param velocityVector: this is the initial velocity vector
    :param time: time you wanted to know the position for
    :param mu: the standard gravitational parameter
    :return: position_at_t
    '''

    magPosition = sqrt(positionVector[0]**2 + positionVector[1]**2 + positionVector[2]**2)
    magVelocity = sqrt(velocityVector[0]**2 + velocityVector[1]**2 + velocityVector[2]**2)
    radialVelocity = np.dot(positionVector, velocityVector)/magPosition
    alpha = 2 / magPosition - magVelocity**2 / mu
    a = 1 / alpha
    x = kepler_eqn(time, magPosition, radialVelocity, a, mu)
    f, g = lagrange_coeff(x, time, magPosition, alpha, mu)
    finalPositionVector = f * positionVector + g * velocityVector
    return finalPositionVector

# here we create a var for how many times we want to in the model
samples = 1032

# we need to keep track of time for each the cubesats
timeStep = 15               # secs
timeFinal = 60*60*24        # this is one day in seconds
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
stdGravPar = 398600                                 # (km^3)/(s^2)
radiusInitial = 6378 + 500                          # km
angMoCombine = sqrt(radiusInitial*stdGravPar)       # angular momentum for initial orbit km^3/s^2
velTan = angMoCombine/radiusInitial                 # tangential velocity km/s
posCombine = np.array([radiusInitial, 0, 0])  # combine position vector at t = 0
velCombine = np.array([0, velTan, 0])         # combine velocity vector at t = 0

time1km = np.array([])
for i in range(samples):
    vx_sc1 = velCombine[0] - vx_sep[i]
    vy_sc1 = velCombine[1] - vy_sep[i]
    vz_sc1 = velCombine[2] - vz_sep[i]
    vx_sc2 = velCombine[0] + vx_sep[i]
    vy_sc2 = velCombine[1] + vy_sep[i]
    vz_sc2 = velCombine[2] + vz_sep[i]
    vel_sc1 = np.array([vx_sc1, vy_sc1, vz_sc1])
    vel_sc2 = np.array([vx_sc2, vy_sc2, vz_sc2])

    timeCounter = 0
    for j in range(timeArray.size):
        pos_sc1 = position_from_state_vector(posCombine, vel_sc1, timeArray[j], stdGravPar)
        pos_sc2 = position_from_state_vector(posCombine, vel_sc2, timeArray[j], stdGravPar)

        distance = sqrt((pos_sc1[0] - pos_sc2[0])**2 + (pos_sc1[1] - pos_sc2[1])**2 + (pos_sc1[2] - pos_sc2[2])**2)
        if distance >= 1:
            timeCounter += 1

    time1km = np.append(timeArray, timeCounter * timeStep)

np.save('data/timeResults.npy', time1km)