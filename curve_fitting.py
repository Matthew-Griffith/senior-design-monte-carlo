import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

results = np.load('testing.npy')
time1km = results[:, 1]
cases = np.arange(0, 1000)
stdDeviation = np.std(time1km)
mean = np.mean(time1km)

plt.hist(time1km, bins=45, facecolor='blue')
x = np.sort(time1km)
y = (mlab.normpdf(x, mean, stdDeviation))*(65/(1.30081002398*10**(-5)))
plt.plot(x, y, 'r+', linewidth=1, label='Normal Distribution')

plt.title('Monte Carlo Analysis for Velocity Magnitude')
plt.xlabel('Time within 1 km [secs]')
plt.ylabel('Cases')
plt.legend(loc='upper right')
plt.grid()
plt.savefig('histogram_velocity_magnitude.png')
plt.show()

print(max(time1km))
print(min(time1km))
print(mean)
print(stdDeviation)
