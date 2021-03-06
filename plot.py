import numpy             as np
import matplotlib.pyplot as plt

hist_data  = np.loadtxt("histogramdata.csv",dtype=float,delimiter=",",usecols=(1))
bin_labels = np.loadtxt("histogramdata.csv",dtype=float,delimiter=",",usecols=(0))
bin_labels = [round(x,2) for x in bin_labels]

#plt.bar(bin_labels,hist_data,width=1.0,align='edge')
plt.plot(bin_labels,hist_data)
plt.plot(bin_labels,[0 for x in bin_labels])
plt.plot(bin_labels,[4*x**2*2**3*np.exp(-2*2*x) for x in bin_labels])
plt.xlabel('radius/a.u.')
plt.ylabel('recovered probability (normalised frequency)')
plt.title('Radial Probability Density Recovered from a Hydrogenic Sampler')
plt.show()
