import numpy             as np
import matplotlib.pyplot as plt

hist_data  = np.loadtxt("histogramdata.csv",dtype=float,delimiter=",",usecols=(1))
bin_labels = np.loadtxt("histogramdata.csv",dtype=float,delimiter=",",usecols=(0))
bin_labels = [round(x,2) for x in bin_labels]

plt.bar(bin_labels,hist_data,width=1.0,align='edge')
plt.xlabel('radius, Hartree units')
plt.ylabel('recovered probability (normalised frequency)')
plt.show()