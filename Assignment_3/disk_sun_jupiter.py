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

def plot_map(sph, particles, title, N=100, L=200, show=True):

    #print 'Plotting', title
   
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

    Sun = particles[0]
    Sx = Sun.x.value_in(units.AU)
    Sy = Sun.y.value_in(units.AU)
    Sz = Sun.z.value_in(units.AU)
    #Sr = Sun.radius.value_in(units.AU)
    Sr = 1    
    sun = Circle((Sx,Sy),Sr,color='y')

    Jupiter = particles[1]
    Jx = Jupiter.x.value_in(units.AU)
    Jy = Jupiter.y.value_in(units.AU)
    Jz = Jupiter.z.value_in(units.AU)
    #Jr = Jupiter.radius.value_in(units.AU)
    Jr = 0.5     
    jup = Circle((Jx,Jy),Jr,color='orange')

    rho, rhovx, rhovy, rhovz, rhoe = sph.get_hydro_state_at_point(
        x, y, z, vx, vy, vz)
    rho = rho.reshape((N + 1, N + 1))

    fig,ax = plt.subplots(1,figsize=(8, 8))    
    ax.imshow(np.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)),
                  extent=[-L / 2, L / 2, -L / 2, L / 2], vmin=10, vmax=15)
    ax.add_patch(sun)
    ax.add_patch(jup)

    plt.title(title)
    plt.xlabel('AU')
    plt.savefig('task3/1000_1/new_'+title+'.png', dpi=150)

    if show:
        plt.show()
    
    plt.close()


def return_L9_radius(disk, Mstar, Rmin): 
    converter=nbody_system.nbody_to_si(Mstar, Rmin)     
    lr,mr = disk.LagrangianRadii(converter)
    return lr[7].value_in(units.AU)


def hydro_sink_particles(sink, gas):
    removed_particles = Particles()
    #for s in sinks:
    #    xs, ys, zs = [s.x, s.y, s.z]
    #    radius_squared = s.radius**2
    #    insink = gas.select_array(lambda x,y,z: (x-xs)**2+(y-ys)**2+(z-zs)**2 < radius_squared,['x','y','z'])  
    #    if len(insink)==0:
    #        return insink
    
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


def init_star_planet():
    particles = Particles(2)

    # Iniialize the Sun
    Sun = particles[0]
    Sun.mass = 1.0 | units.MSun
    Sun.radius = 1.0 | units.RSun
    Sun.position = (855251, -804836, -3186) |units.km
    Sun.velocity = (7.893, 11.894, 0.20642) |(units.m/units.s)
       
    # Initialize the Jupiter
    Jupiter = particles[1]
    Jupiter.mass = 1 | units.MJupiter
    Jupiter.radius = 1.0 | units.RJupiter
    Jupiter.position = (-4.9829, 2.062, -0.10990) | units.AU
    Jupiter.velocity = (-5.158, -11.454, 0.13558) | units.kms

    particles.move_to_center()

    return particles


def init_sink_particle(m, r, position, velocity):
    sink = Particles(1)
    sink.mass = 0. | units.MSun
    sink.radius = r
    sink.x, sink.y, sink.z = position
    sink.vx, sink.vy, sink.vz = velocity
    return sink


def gravity_hydro_bridge(gravity, hydro, sink, local_particles, Rmin, t_end=1000.|units.yr, dt=1.|units.yr):

    # Set the local particles in each system
    gravity_particles, disk_star, disk_gas = local_particles
    Sun = disk_star[0]
    Jupiter = gravity_particles[0]
    Mstar = Sun.mass.in_(units.MSun)

    print 'Bridging...'
    # Build up the bridge between gravity and hydrodynamics
    grav_hydro = bridge.Bridge(use_threading=False)
    grav_hydro.add_system(gravity, (hydro,))
    grav_hydro.add_system(hydro, (gravity,))
    grav_hydro.timestep = dt

    # Set up channels for updating the particles
    channel_to_grav = gravity.particles.new_channel_to(gravity_particles)
    channel_from_grav = gravity_particles.new_channel_to(gravity.particles)
    channel_to_hydro_gas = hydro.gas_particles.new_channel_to(disk_gas)
    channel_from_hydro_gas = disk_gas.new_channel_to(hydro.gas_particles)
    channel_to_hydro_star = hydro.dm_particles.new_channel_to(disk_star)
    channel_from_hydro_star = disk_star.new_channel_to(hydro.dm_particles)

    # Sanity checks:
    print 'Sanity checks:'
    print 'Sun\'s position', (Sun.position).in_(units.AU)
    print 'Jupiter\'s position', (Jupiter.position).in_(units.AU)
    print 'Disk particle map saved to: initial_check_disk.png'
    plot_map(hydro, [Sun, Jupiter], N=400, L=200,\
             title='initial_check_disk', show=False)

    times = list(np.arange(0.0, t_end.value_in(units.yr)+dt.value_in(units.yr), dt.value_in(units.yr)))
    a_Jup = list()
    e_Jup = list()
    disk_size = list()
    accreted_mass = list()

    # start evolotuion
    print 'Start evolving...'
    model_time = 0.0 | units.yr
    while model_time <= t_end:

        # Evolve the bridge system for one step
        grav_hydro.evolve_model(model_time)
        channel_to_grav.copy()
        channel_to_hydro_gas.copy()
        channel_to_hydro_star.copy()

        Sun = disk_star[0]
        Jupiter = gravity_particles[0]

        # Save the data for plots
        Sun_and_Jupiter = Particles()
        Sun_and_Jupiter.add_particle(Sun)
        Sun_and_Jupiter.add_particle(Jupiter)
        orbit = orbital_elements_from_binary(Sun_and_Jupiter, G=constants.G)
        a = orbit[2].value_in(units.AU)
        a_Jup.append(a)
        e = orbit[3]
        e_Jup.append(e)
        lr9 = return_L9_radius(disk_gas, Mstar, Rmin)
        disk_size.append(lr9)
        sink_mass = (sink.mass).value_in(units.MJupiter)[0]
        accreted_mass.append(sink_mass)
        with open('task3/1000_1/task3_data_new.txt', 'a') as output_file:
            output_file.write('%f,%f,%f,%f,%f\n'%(model_time.value_in(units.yr), a, e, lr9, sink_mass))
        
        print 'Time = %.1f yr:'%model_time.value_in(units.yr), \
              'a = %.2f au, e = %.3f,'%(a, e), \
              'disk size = %.2f au,'%lr9, \
              'accreted mass = %.4f M_Jupiter.'%sink.mass.value_in(units.MJupiter)
        plot_map(hydro, [Sun, Jupiter], N=400, L=200, title='%.1f_year'%model_time.value_in(units.yr), show=False)
        

        # Add the 'sinked' mass to Jupiter & keep the sink particle along with Jupiter
        sink.position = Jupiter.position
        sink.radius = a * (1-e) * ((Jupiter.mass).value_in(units.MSun)/3.)**(1./3.) | units.AU
        removed_particles = hydro_sink_particles(sink, disk_gas)
        Jupiter.mass += removed_particles.mass.sum()
        channel_from_grav.copy()
        channel_from_hydro_gas.copy()

        model_time += dt

    gravity.stop()
    hydro.stop()

    return times, a_Jup, e_Jup, disk_size, accreted_mass


def main(t_end=1000.|units.yr, dt=10.|units.yr):

    print 'Initializing...'
    #converter = nbody_system.nbody_to_si(1.|units.MSun, 1.|units.AU)

    # Initialize the Sun and Jupiter system
    Sun, Jupiter = init_star_planet()
    Sun_and_Jupiter0 = Particles()
    Sun_and_Jupiter0.add_particle(Sun)
    Sun_and_Jupiter0.add_particle(Jupiter)
    orbit0 = orbital_elements_from_binary(Sun_and_Jupiter0, G=constants.G)
    a0 = orbit0[2].in_(units.AU)
    e0 = orbit0[3]
    Hill_radius0 = a0 * (1-e0) * ((1.0|units.MJupiter).value_in(units.MSun)/3.)**(1./3.)


    # Initialize the direct N-body integrator
    gravity_particles = Particles()
    gravity_particles.add_particle(Jupiter)
    converter_gravity = nbody_system.nbody_to_si(gravity_particles.mass.sum(), gravity_particles.position.length())
    gravity = ph4(converter_gravity)
    gravity.particles.add_particles(gravity_particles)
    gravity.timestep = dt


    # Set proto-disk parameters
    N = 4000
    Mstar = 1. | units.MSun
    Mdisk = 0.01 * Mstar
    Rmin = 1. | units.AU
    Rmax = 100. | units.AU

    # Initialize the proto-planetary disk    
    np.random.seed(1)
    converter_disk = nbody_system.nbody_to_si(Mstar, Rmin)
    disk = ProtoPlanetaryDisk(N, convert_nbody=converter_disk, densitypower=1.5, Rmin=1., Rmax=Rmax/Rmin, q_out=1., discfraction=Mdisk/Mstar)
    disk_gas = disk.result
    disk_star = Particles()
    disk_star.add_particle(Sun)
    
    # Initialize the sink particle  
    sink = init_sink_particle(0.|units.MSun, Hill_radius0, Jupiter.position, Jupiter.velocity)
   
    # Initialising the SPH code
    #converter_hydro = nbody_system.nbody_to_si(Mdisk, Rmax)
    sph = Fi(converter_disk, mode="openmp")
    sph.gas_particles.add_particles(disk_gas)
    sph.dm_particles.add_particles(disk_star)
    sph.parameters.timestep = dt


    # bridge and evolve
    times, a_Jup, e_Jup, disk_size, accreted_mass = gravity_hydro_bridge(gravity, sph, sink,\
                                     [gravity_particles, disk_star, disk_gas], Rmin, t_end, dt)

    return times, a_Jup, e_Jup, disk_size, accreted_mass



if __name__ == "__main__":
    t_end = 1000. | units.yr
    dt = 1. |units. yr
    times, a_Jup, e_Jup, disk_size, accreted_mass = main(t_end, dt)

    # save data
    #output_data = np.column_stack((times, a_Jup, e_Jup, disk_size, accreted_mass))
    #np.savetxt('task3/100_0.5/task3_data.txt', output_data, delimiter=',')
    #print 'Saved the data of timestamp, semi-major-axis, eccentricity, disk size, and accreted mass in task3_data.txt'

    # make plots
    fig = plt.figure(figsize=(15, 10))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    ax1.plot(times, a_Jup)
    ax2.plot(times, e_Jup)
    ax3.plot(times, disk_size)
    ax4.plot(times, accreted_mass)

    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xlabel('Time [yr]')
    
    ax1.set_ylabel('a [au]')
    ax2.set_ylabel('e')
    ax3.set_ylabel('disk size [au]')
    ax4.set_ylabel('accreted mass [MJupiter]')

    fig.savefig('task3/1000_1/task3_data_new.png', bbox_inches='tight')
    #plt.show()

