import numpy as np
from math import sin, cos, pi, sqrt

# here we create a var for how many times we want to in the model
samples = 10

# we need to keep track of time for each the cubesats
timeStep = 60               # secs
timeFinal = 60*60*24*14     # this is two weeks in seconds
timeArray = np.arange(0, timeFinal, timeStep)

# this is just an array of random values for the magnitude of velocity
velMag = np.random.uniform(0.8, 1.2, size=samples)    # this is in m/s

# setting initial orbit parameters
initialEcc = 0
stdGravPar = 3.98574405 * 10**14            # (m^3)/(s^2)
initialRadius = 500000 + 6371000            # m

# here we will create functions for some of the eqns that we use multiple times

def angular_momentum_cross(position, velocity):
    '''this function take the position and velocity a point and find the magnitude of the angular momentum'''
    return np.cross(position, velocity)

def orbit_period(angularMomentum, eccentricity):
    '''simple eqn to find the period of an orbit.'''
    return ((angularMomentum**2)/stdGravPar)*(1/(1 - eccentricity**2))

def orbit_eccentricity(velocity, angularMomentum):
    '''this function finds the eccentricity for an orbit'''
    return velocity[1] * (angularMomentum/stdGravPar) - 1

def position_function(angularMomentum, eccentricity, period):
    ''' this function takes in angular momentum, eccentricity and period and returns
        an array for the x and y values in a numpy array.'''
    pos = np.array([])

    # the two following for loops find the position of the two cubesats at for the
    # timeArray and place them into their own arrays
    for i in range(timeArray.size):
        tempPos = ((angularMomentum ** 2) / stdGravPar) * (
            1 / (1 + eccentricity * cos((timeArray[i] / periodSC1) * (2 * pi)))) * np.array(
            [cos((timeArray[i] / periodSC1) * (2 * pi)), sin((timeArray[i] / period) * (2 * pi))])
        pos = np.append(pos, tempPos)

    # appending the array just makes it 1D, here with are separating x and y
    # into two different columns
    return np.reshape(pos, (timeArray.size, 2))

initialAngMo = sqrt(initialRadius*(stdGravPar*(1 + initialEcc * cos(0))))
initialPeriod = orbit_period(initialAngMo, initialEcc)
initialVel = np.array([0, (initialEcc + cos(0))]) * (stdGravPar/initialAngMo)
initialDis = np.array([initialRadius, 0])

# outside of the main for loop we need to create an emtpy array that will be used
# to store the results, it need to be created here because the values will be appended
# at the end of the for loop.
result = np.array([])

# this loop will cycle through our random array of values and compute the time the
# cubesats are within 1 km
for n in range(samples):
    # now we have to take these initial values and find them for the new orbits
    initialVelSC1 = initialVel + np.array([0, velMag[n]])
    initialVelSC2 = initialVel - np.array([0, velMag[n]])

    # now we will solve for the e, h and T for the new orbits
    angMoSC1 = angular_momentum_cross(initialDis, initialVelSC1)
    angMoSC2 = angular_momentum_cross(initialDis, initialVelSC2)
    eccSC1 = orbit_eccentricity(initialVelSC1, angMoSC1)
    eccSC2 = orbit_eccentricity(initialVelSC2, angMoSC2)
    periodSC1 = orbit_period(angMoSC1, eccSC1)
    periodSC2 = orbit_period(angMoSC2, eccSC2)

    # find the position of the satellites for each time step.
    posSC1 = position_function(angMoSC1, eccSC1, periodSC1)
    posSC2 = position_function(angMoSC2, eccSC2, periodSC2)

    # now we measure the distance between the two cubesats at each point in time
    # and then see if they are within 1 km. if they are a counter is incremented
    # and that will be multiplied by the timeStep to give us a time that they are
    # within 1 km of each other
    timeCounter = 0
    for i in range(timeArray.size):
        distance = sqrt((posSC2[i, 0] - posSC1[i, 0])**2 + (posSC2[i, 1] - posSC1[i, 1])**2)

        if distance <= 1000:
            timeCounter += 1

    time1km = timeCounter * timeStep
    result = np.append(result, [velMag[n], time1km])

# This program will take a very long time to run, so the results data is saved to a
# file for further analysis
result = np.reshape(result, (samples, 2))
np.save('testing', result)

