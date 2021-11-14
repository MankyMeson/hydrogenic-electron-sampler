import numpy             as np
import matplotlib.pyplot as plt

hist_data  = np.loadtxt("histogramdata.csv",dtype=float,delimiter=",",usecols=(1))
bin_labels = np.loadtxt("histogramdata.csv",dtype=float,delimiter=",",usecols=(0))

for i in range (len(bin_labels)):
  bin_labels[i] = round(bin_labels[i], 2)

plt.bar(bin_labels,hist_data,width=1.0,align='edge')
plt.xlabel('radius, Hartree units')
plt.ylabel('recovered probability (normalised frequency)')
plt.show()