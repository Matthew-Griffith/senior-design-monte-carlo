import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.mlab as mlab

results = np.load('testing.npy')
time1km = results[:, 1]
cases = np.arange(0, 1000)
# stdDeviation = np.std(time1km)
# mean = np.mean(time1km)

# plt.plot(cases, time1km, label='Raw Data', marker='.', linestyle='none')
plt.hist(time1km, bins=45)

# y = mlab.normpdf(60, mean, stdDeviation)
# plt.plot(y, 'r--', linewidth=1)

plt.title('Monte Carlo Analysis for Velocity Magnitude')
plt.xlabel('Time within 1 km [secs]')
plt.ylabel('cases')
plt.legend(loc='upper right')
plt.grid()
plt.savefig('histogram_velocity_magnitude.png')
plt.show()
