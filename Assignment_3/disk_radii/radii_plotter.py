import pickle as pkl
import numpy as np
from matplotlib import pyplot as plt
from amuse.units import units


if __name__ == "__main__":

    radii_dir = 'disk_radii/'   

    N500 = pkl.load(open(radii_dir+'lr9_N500.p','rb'))
    N1000 = pkl.load(open(radii_dir+'lr9_N1000.p','rb'))
    N2000 = pkl.load(open(radii_dir+'lr9_N2000.p','rb'))
    N4000 = pkl.load(open(radii_dir+'lr9_N4000.p','rb'))

    t_end = 1000 | units.yr
    n_steps = 1000
    dt = t_end/float(n_steps)
    t = np.arange(0,(t_end+dt).value_in(units.yr),dt.value_in(units.yr))    

    plt.plot(t,N500,label='N = 500')
    plt.plot(t,N1000,label='N = 1000')
    plt.plot(t,N2000,label='N = 2000')
    plt.plot(t,N4000,label='N = 4000')
    plt.xlabel('Time (years)')
    plt.ylabel('Disk size (AU)')
    plt.title('Disk radius over time')
    plt.legend()
    plt.savefig(radii_dir+'radius_time_plot.png')
    plt.show()

