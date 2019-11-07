# -*- coding: ascii -*-
"""
Task 3: Sun + Jupiter(sink) + disk
"""
#from __future__ import print_function
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle

from amuse.lab import *
from amuse.units import units, nbody_system, quantities
from amuse.couple import bridge
from amuse.community.fi.interface import Fi
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ext.orbital_elements import orbital_elements_from_binary

def plot_map(sph, Sun_and_Jupiter, title, N=100, L=1, show=True):

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

    Sx = Sun_and_Jupiter[0].x.value_in(units.AU)
    Sy = Sun_and_Jupiter[0].y.value_in(units.AU)
    Sz = Sun_and_Jupiter[0].z.value_in(units.AU)
    Sr = 1    
    sun = Circle((Sx,Sy),Sr,color='y')

    Jx = Sun_and_Jupiter[1].x.value_in(units.AU)
    Jy = Sun_and_Jupiter[1].y.value_in(units.AU)
    Jz = Sun_and_Jupiter[1].z.value_in(units.AU)
    Jr = 0.5     
    jup = Circle((Jx,Jy),Jr,color='orange')

    rho, rhovx, rhovy, rhovz, rhoe = sph.get_hydro_state_at_point(
        x, y, z, vx, vy, vz)
    rho = rho.reshape((N + 1, N + 1))

    fig,ax = plt.subplots(1,figsize=(8, 8))    
    ax.imshow(np.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)),
                  extent=[-L / 2, L / 2, -L / 2, L / 2], vmin=10, vmax=15,origin='lower')
    ax.add_patch(sun)
    ax.add_patch(jup)
    plt.title(title)
    plt.xlabel('AU')
    plt.savefig(title)

    if show:
        plt.show()
    
    plt.close()


def return_L9_radius(disk,Mstar,Rmin): 
    converter=nbody_system.nbody_to_si(Mstar, Rmin)     
    lr,mr = disk.LagrangianRadii(converter)
    return lr[7].value_in(units.AU)

def Hill_radius(a, e, Sun_Jupiter):
    return a * (1-e) * ((Sun_Jupiter[1].mass / Sun_Jupiter[0].mass) / 3.) ** (1. / 3.)

def hydro_sink_particles(sink, gas):

    removed_particles = Particles()
    xs, ys, zs = sink.x, sink.y, sink.z
    radius_squared = sink.radius**2
    insink = gas.select_array(lambda x,y,z: (x-xs)**2+(y-ys)**2+(z-zs)**2 < radius_squared,['x','y','z'])  
    if len(insink)==0:
        return insink
    
    cm = sink.position*sink.mass
    p = sink.velocity*sink.mass
    sink.mass += insink.total_mass()
    sink.position = (cm+insink.center_of_mass()*insink.total_mass())/sink.mass
    sink.velocity = (p+insink.total_momentum())/sink.mass
    removed_particles.add_particles(insink)
    gas.remove_particles(insink)
    
    return removed_particles


def init_sun_jupiter():
    particles = Particles(2)

    sun = particles[0]
    sun.mass = 1.0 | units.MSun
    sun.radius = 1.0 | units.RSun
    sun.position = (855251, -804836, -3186) |units.km
    sun.velocity = (7.893, 11.894, 0.20642) |(units.m/units.s)
       
    jupiter = particles[1]
    jupiter.mass = 1 | units.MJupiter
    jupiter.radius = 1.0 | units.RJupiter
    jupiter.position = (-4.9829, 2.062, -0.10990) | units.AU
    jupiter.velocity = (-5.158, -11.454, 0.13558) | units.kms

    particles.move_to_center()

    return particles


def init_sink_particle(m, r, position, velocity):
    sink = Particles(1)
    sink.mass = 0. | units.MSun
    sink.radius = r
    sink.x, sink.y, sink.z = position
    sink.vx, sink.vy, sink.vz = velocity
    return sink


def gravity_hydro_bridge(gravity, hydro, sink, local_particles, Rmin, t_end=1000.|units.yr, dt=10.|units.yr):

    Sun_and_Jupiter, disk_gas = local_particles
    Mstar = 1.0 | units.MSun

    print 'Bridging...'
    # Build up the bridge between gravity and hydrodynamics
    grav_hydro = bridge.Bridge(use_threading=False)
    grav_hydro.add_system(gravity, (hydro,))
    grav_hydro.add_system(hydro, (gravity,))
    grav_hydro.timestep = dt

    # Set up channels for updating the particles
    channel_from_grav = gravity.particles.new_channel_to(Sun_and_Jupiter)
    channel_from_hydro = hydro.gas_particles.new_channel_to(disk_gas)
    channel_to_grav = Sun_and_Jupiter.new_channel_to(gravity.particles)
    channel_to_hydro = disk_gas.new_channel_to(hydro.gas_particles)

    # Preparing lists for data-recording
    a_Jup = []
    e_Jup = []
    disk_size = []
    accreted_mass = []
    accreted_mass.append((sink.mass).value_in(units.MJupiter)[0])
    sink0_mass = 0 | units.MJupiter

    # Start evolution
    print 'Start evolving...'
    times = quantities.arange(0.|units.yr, t_end+1*dt, dt)
    model_time = 0.0 | units.yr
    while model_time <= t_end:
        
        # Save the data for plots
        orbit = orbital_elements_from_binary(Sun_and_Jupiter, G=constants.G)
        a = orbit[2].value_in(units.AU)
        e = orbit[3]
        lr9 = return_L9_radius(disk_gas, Mstar, Rmin|units.AU)
        a_Jup.append(a)
        e_Jup.append(e)
        disk_size.append(lr9)
       
        # Plotting system
        print 'Time = %.1f yr:'%model_time.value_in(units.yr), \
                  'a = %.2f au, e = %.2f,'%(a, e), \
                  'disk size = %.2f au'%lr9
        plot_map(hydro,Sun_and_Jupiter,'q3_distribution_plot_4K/{0}.png'.format(int(model_time.value_in(units.yr))),show=False)
        
        # Evolve the bridge system for one step
        model_time += dt
        grav_hydro.evolve_model(model_time)
        channel_from_grav.copy()
        channel_from_hydro.copy()

        # Calculating accreted mass in new position
        Jupiter = gravity.particles[0]   
        sink.position = Jupiter.position
        sink.radius = Hill_radius(a, e, Sun_and_Jupiter) | units.AU       
        removed_particles = hydro_sink_particles(sink, disk_gas)     
        Jupiter.mass += sink.mass - sink0_mass
        sink0_mass = sink.mass.copy()
        channel_to_grav.copy()
        channel_to_hydro.copy()

    gravity.stop()
    hydro.stop()

    return a_Jup, e_Jup, disk_size, times, accreted_mass


def main(t_end=1000.|units.yr, dt=10.|units.yr):

    print 'Initializing...'
    converter = nbody_system.nbody_to_si(1.|units.MSun, 1.|units.AU)

    # Initialize the gravity system
    Sun_and_Jupiter = init_sun_jupiter()
    Sun = Sun_and_Jupiter[0]
    Jupiter = Sun_and_Jupiter[1]
    orbit0 = orbital_elements_from_binary(Sun_and_Jupiter, G=constants.G)
    a0 = orbit0[2].in_(units.AU)
    e0 = orbit0[3]
    Hill_radius0 = a0 * (1-e0) * ((1.0|units.MJupiter).value_in(units.MSun)/3.)**(1./3.)

    # Initialising the direct N-body integrator
    gravity = ph4(converter)
    gravity.particles.add_particle(Jupiter)
    gravity.timestep = dt
    
    # Setting proto-disk parameters
    N = 4000
    Mstar = 1. | units.MSun
    Mdisk = 0.01 * Mstar
    Rmin = 2.
    Rmax = 100.

    # Initialising the proto-planetary disk    
    np.random.seed(1)
    disk = ProtoPlanetaryDisk(N, convert_nbody=converter, densitypower=1.5, Rmin=Rmin, Rmax=Rmax, q_out=1., discfraction=Mdisk/Mstar)
    disk_gas = disk.result
    
    # Initialize the sink particle  
    sink = init_sink_particle(0.|units.MSun, Hill_radius0, Jupiter.position, Jupiter.velocity)
   
    # Initialising the SPH code
    sph = Fi(converter, mode="openmp")
    sph.gas_particles.add_particles(disk_gas)
    sph.dm_particles.add_particle(Sun)
    sph.parameters.timestep = dt

    # bridge and evolve
    a_Jup, e_Jup, disk_size, times, accreted_mass = gravity_hydro_bridge(gravity, sph, sink,
                                     [Sun_and_Jupiter, disk_gas], Rmin, t_end, dt)

    return a_Jup, e_Jup, disk_size, times, accreted_mass


if __name__ == "__main__":
    t_end = 1000. | units.yr
    dt = 1. |units. yr
    a_Jup, e_Jup, disk_size, times, accreted_mass = main(t_end, dt)
    
    sim_output = np.column_stack((times.value_in(units.yr), a_Jup, e_Jup, disk_size, accreted_mass))
    np.savetxt('sim_output_q3.csv',sim_output,delimiter=',')
















