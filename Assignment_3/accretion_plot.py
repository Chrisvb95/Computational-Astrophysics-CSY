import numpy as np
from matplotlib import pyplot as plt


fname = 'sim_output.csv'
data = np.genfromtxt(fname, delimiter=',')

# '(t.value_in(units.yr), a_Jup, e_Jup, disk_size, accreted_mass'

#plt.scatter(data[:,0]

time = data[:,0]
a_Jup = data[:,1]
e_Jup = data[:,2]
disk_size = data[:,3]
accreted_mass = data[:,4]

fig, ax1 = plt.subplots()
ax1.set_xlabel('Time (year)')
ax1.set_ylabel('Semi-major axis (AU)',color='tab:blue')
ax1.plot(time, a_Jup,color='tab:blue')
ax2 = ax1.twinx()
ax2.set_ylabel('Eccentricity',color='tab:red')
ax2.plot(time, e_Jup, color='tab:red')
fig.tight_layout()
plt.savefig('Jupiter_properties.png')
plt.show()
plt.close()

fig, axs = plt.subplots(2)
axs[0].plot(time, disk_size)
axs[0].set_ylabel('Disk size (AU)')
axs[0].axvline(x=101,color='k',linestyle='--')
axs[0].text(90, 93,'closest approach',rotation=90,fontsize=6)
axs[1].plot(time, accreted_mass)
axs[1].set_ylabel('Accreted mass ($M_{Jupiter}$)')
axs[1].set_xlabel('Time (years)')
axs[1].axvline(x=101,color='k',linestyle='--')
axs[1].text(90, 0.5,'closest approach',rotation=90,fontsize=6)
plt.savefig('accretion_plot.png')
plt.show()
plt.close()


