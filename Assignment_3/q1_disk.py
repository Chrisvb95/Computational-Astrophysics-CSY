# -*- coding: ascii -*-
"""
Creates a protoplanetary disk around a sun-like star
"""
from __future__ import print_function
import numpy as np
from matplotlib import pyplot as plt
from amuse.lab import *
from amuse.community.fi.interface import Fi
from amuse.units import units
from amuse.units import nbody_system
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.datamodel import Particles
import pickle as pkl

def return_L9_radius(disk,Mstar,Rmin): 
    converter=nbody_system.nbody_to_si(Mstar, Rmin)     
    lr,mr = disk.LagrangianRadii(converter)
    return lr[7].value_in(units.AU)

def plot_map(sph, title, N=100, L=1, show=True):

    print('Plotting', title)

    L = 200
   
    x, y = np.indices((N + 1, N + 1))

    x = L * (x.flatten() - N / 2.) / N
    y = L * (y.flatten() - N / 2.) / N
    z = x * 0.
    vx = 0. * x
    vy = 0. * x
    vz = 0. * x

    x = units.AU(x)
    y = units.AU(y)
    z = units.AU(z)
    vx = units.kms(vx)
    vy = units.kms(vy)
    vz = units.kms(vz)

    rho, rhovx, rhovy, rhovz, rhoe = sph.get_hydro_state_at_point(
        x, y, z, vx, vy, vz)
    rho = rho.reshape((N + 1, N + 1))
    
    plt.figure(figsize=(8, 8))
    plt.imshow(np.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)),
                  extent=[-L / 2, L / 2, -L / 2, L / 2], vmin=10, vmax=15)
    plt.title(title)
    plt.xlabel('AU')
    plt.savefig(title)

    if show:
        plt.show()
    
    plt.close()
        
if __name__ == "__main__":

    # Variables
    N = 4000
    t_end = 1000 | units.yr
    n_steps = 1000
    integration_dt = 0.125 | units.yr
    Mstar = 1. | units.MSun
    Rstar = 695510000 | units.m
    Mdisk = 0.01 | units.MSun
    Rmin = 1
    Rmax=100
    overwrite = True
    d_plot_dir = 'distribution_plots_4k/'
    
    # Setting a seed
    np.random.seed(1)

    # Initialising disk
    convert = nbody_system.nbody_to_si(Mstar , 1. | units.AU)
    proto = ProtoPlanetaryDisk(
        N, convert_nbody=convert, densitypower=1.5, Rmin=Rmin, Rmax=Rmax, q_out=1., discfraction=Mdisk/Mstar)
    disk_gas = proto.result
    
    # Initialising star
    sun = Particles(1)
    sun.mass = Mstar
    sun.radius = 0.1 | units.AU
    sun.x = 0. | units.AU
    sun.y = 0. | units.AU
    sun.z = 0. | units.AU
    sun.vx = 0. | units.kms
    sun.vy = 0. | units.kms
    sun.vz = 0. | units.kms

    # Initialising SPH code
    sph = Fi(convert)
    
    # Adding particles to code
    sph.gas_particles.add_particles(disk_gas)
    sph.particles.add_particles(sun)

    # Makking channels from SPH code to sun and disk
    sph_to_star = sph.particles.new_channel_to(sun)
    sph_to_disk = sph.particles.new_channel_to(disk_gas)

    lr9 = []
    lr9.append(return_L9_radius(disk_gas, Mstar, Rmin|units.AU))

    if overwrite: 

        dt = t_end/float(n_steps)
        time = 0 | units.yr

        # Evolving the system
        while time < t_end:
            time += dt
            print('Time:', time)
            sph.evolve_model(time)
            #plot_map(sph,d_plot_dir+'test{0}.png'.format(time.value_in(units.yr)),show=False)
            sph_to_star.copy()
            sph_to_disk.copy()
            lr9.append(return_L9_radius(disk_gas, Mstar, Rmin|units.AU))
            
        sph.stop()
    
    # Plotting the radius size
    pkl.dump(lr9,open('disk_radii/lr9_N4000.p','wb'))    
    t = np.arange(0,(t_end+dt).value_in(units.yr),dt.value_in(units.yr))    
    plt.scatter(t,lr9)  
    plt.xlabel('Time (years)')
    plt.ylabel('Disk size (AU)')
    plt.savefig('disk_radius_time_plot.png')    
    plt.show()

