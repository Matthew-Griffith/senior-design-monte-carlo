import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

results = np.load('testing.npy')

# Linear regression on the data set
slope, intercept, r_value, p_value, std_err = stats.linregress(results[:, 0], results[:, 1])
velRange = np.arange(0.8, 1.2, .001)
linTime = velRange * slope + intercept
# print(slope, intercept)

plt.plot(results[:, 0], results[:, 1], label='Raw Data', marker='.', linestyle='none')
plt.plot(velRange, linTime, label='Linear regression', color='red')
plt.title('Monte Carlo Analysis for Velocity Magnitude')
plt.xlabel('Velocity [m/s]')
plt.ylabel('Time within 1 km [secs]')
plt.legend(loc='upper right')
plt.grid()
plt.savefig('1000_samples_lin_regression.png')
plt.show()
