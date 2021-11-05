import numpy as np
import matplotlib.pyplot as plt

MW= []
C = []
mort_child_min= []
mort_child_ave= []
mort_child_max= []

file = open('values.txt', 'r')
data = file.readlines()
for line in data:
    line = line.strip()
    columns = line.split()
    MW.append(np.float64(columns[1]))
    C.append(np.float64(columns[0]))
file.close()
plt.scatter(C, MW, alpha=0.5)
plt.xlabel(r'$log_{10}(C_0) (\mu g/m^{-3})$')
plt.ylabel('Molar mass (gr/mol)')
plt.show()
