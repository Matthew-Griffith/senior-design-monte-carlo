import numpy as np
from math import sin, cos, pi, sqrt

# here we create a var for how many times we want to in the model
samples = 1000

# we need to keep track of time for each the cubesats
timeStep = 60       # secs
timeFinal = 60*60*24*14 # this is two weeks in seconds
timeArray = np.arange(0, timeFinal, timeStep)

# this is just an array of random values for the magnitude of velocity
velMag = np.random.uniform(0.8, 1.2, size=samples)    # this is in m/s

# setting initial orbit parameters
initialEcc = 0
stdGravPar = 3.98574405 * 10**14        # (m^3)/(s^2)
initialRadius = 500000 + 6371000            # m
initialAngMo = sqrt(initialRadius*(stdGravPar*(1 + initialEcc * cos(0))))
initialPeriod = ((initialAngMo**2)/stdGravPar)*(1/(1 - initialEcc**2))
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
    angMoSC1 = np.cross(initialDis, initialVelSC1)
    angMoSC2 = np.cross(initialDis, initialVelSC2)
    eccSC1 = initialVelSC1[1] * (angMoSC1/stdGravPar) - 1
    eccSC2 = initialVelSC2[1] * (angMoSC2/stdGravPar) - 1
    periodSC1 = ((angMoSC1**2)/stdGravPar)*(1/(1 - eccSC1**2))
    periodSC2 = ((angMoSC2**2)/stdGravPar)*(1/(1 - eccSC2**2))

    # similar to result array above these position arrays need to be created empty
    # outside of this loop so that they can appended inside
    disSC1 = np.array([])
    disSC2 = np.array([])

    # the two following for loops find the position of the two cubesats at for the
    # timeArray and place them into their own arrays
    # todo: make this into a function
    for i in range(timeArray.size):
        tempDis1 = ((angMoSC1**2)/stdGravPar) * (1/(1 + eccSC1*cos((timeArray[i]/periodSC1)*(2*pi)))) * np.array([cos((timeArray[i]/periodSC1)*(2*pi)), sin((timeArray[i]/periodSC1)*(2*pi))])
        disSC1 = np.append(disSC1, tempDis1)
        
        tempDis2 = ((angMoSC2**2)/stdGravPar) * (1/(1 + eccSC2*cos((timeArray[i]/periodSC2)*(2*pi)))) * np.array([cos((timeArray[i]/periodSC2)*(2*pi)), sin((timeArray[i]/periodSC2)*(2*pi))])
        disSC2 = np.append(disSC2, tempDis2)
    # appending the array just makes it 1D, here with are separating x and y
    # into two different columns
    disSC1 = np.reshape(disSC1, (timeArray.size, 2))
    disSC2 = np.reshape(disSC2, (timeArray.size, 2))

    # now we measure the distance between the two cubesats at each point in time
    # and then see if they are within 1 km. if they are a counter is incremented
    # and that will be multiplied by the timeStep to give us a time that they are
    # within 1 km of each other
    timeCounter = 0
    for i in range(timeArray.size):
        distance = sqrt((disSC2[i, 0] - disSC1[i, 0])**2 + (disSC2[i, 1] - disSC1[i, 1])**2)

        if distance <= 1000:
            timeCounter += 1

    time1km = timeCounter * timeStep
    result = np.append(result, [velMag[n], time1km])

# This program will take a very long time to run, so the results data is saved to a
# file for further analysis
result = np.reshape(result, (samples, 2))
np.save('testing', result)

