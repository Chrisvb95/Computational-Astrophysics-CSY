from amuse.lab import *
from matplotlib import pyplot as plt
import numpy as np


def return_L9_radius(disk,Mstar,Rmin): 
    converter=nbody_system.nbody_to_si(Mstar, Rmin)     
    lr,mr = disk.LagrangianRadii(converter)
    return lr[7].value_in(units.AU)

if __name__ in '__main__':
     
    Mstar = 1|units.MSun
    Rmin = 1.0|units.AU
    t_end=1000|units.yr
    n_steps=10

    lr9 = []
    for i in range(1,12):
        filename = 'hydro_disk_i{0:04}.amuse'.format(i)
        disk = read_set_from_file(filename,format='amuse')
        lr9.append(return_L9_radius(disk, Mstar, Rmin))

    dt = t_end/float(n_steps)
    print "dt=", dt.in_(units.day)
    time = 0 | units.day

    t = np.arange(0,(t_end+dt).value_in(units.yr),dt.value_in(units.yr))
    
    plt.scatter(t,lr9)  
    plt.xlabel('Time (years)')
    plt.ylabel('Disk size (AU)')
    plt.savefig('disk_radius_time_plot.png')    
    plt.show()

