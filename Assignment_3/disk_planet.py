# -*- coding: ascii -*-
"""
Task 3: Sun + Jupiter(sink) + disk
"""
#from __future__ import print_function
import numpy as np
from matplotlib import pyplot as plt

from amuse.lab import *
from amuse.units import units, nbody_system, quantities
from amuse.couple import bridge
from amuse.community.fi.interface import Fi
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ext.orbital_elements import orbital_elements_from_binary


def make_map(sph, N=100, L=1):

    x, y = numpy.indices((N + 1, N + 1))
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

    return numpy.transpose(rho)


def return_L9_radius(disk,Mstar,Rmin): 
    converter=nbody_system.nbody_to_si(Mstar, Rmin)     
    lr,mr = disk.LagrangianRadii(converter)
    return lr[7].value_in(units.AU)


def hydro_sink_particles(sinks, gas):
    removed_particles = Particles()
    for s in sinks:
        xs, ys, zs = [s.x, s.y, s.z]
        radius_squared = s.radius**2
        insink = gas.select_array(lambda x,y,z: (x-xs)**2+(y-ys)**2+(z-zs)**2 < radius_squared,['x','y','z'])  
        if len(insink)==0:
            return insink

        cm = s.position*s.mass
        p = s.velocity*s.mass
        s.mass += insink.total_mass()
        s.position = (cm+insink.center_of_mass()*insink.total_mass())/s.mass
        s.velocity = (p+insink.total_momentum())/s.mass
        removed_particles.add_particles(insink)
    return removed_particles


def init_sun_jupiter():
    particles = Particles(2)

    sun = particles[0]
    sun.mass = 1.0 | units.MSun
    sun.radius = 1.0 | units.RSun
    sun.position = (855251, -804836, -3186) |units.km
    sun.velocity = (7.893, 11.894, 0.20642) |(units.m/units.s)
    
    jupiter = particles[1]
    jupiter.mass = 1.0 | units.MJupiter
    jupiter.radius = 1.0 | units.RJupiter
    jupiter.position = (-4.9829, 2.062, -0.10990) | units.AU
    jupiter.velocity = (-5.158, -11.454, 0.13558) | units.kms
    
    particles.move_to_center()
	
    return particles


def init_sink_particle(m, r, position, velocity):
    sink = Particles(1)
    sink.mass = m
    sink.radius = r
    sink.x, sink.y, sink.z = position
    sink.vx, sink.vy, sink.vz = velocity
    return sink


def gravity_hydro_bridge(gravity, hydro, local_particles, Rmin, t_end=1000.|units.yr, dt=10.|units.yr):

    Sun_and_Jupiter, disk_gas = local_particles
    Mstar = gravity.particles[1].mass.in_(units.MSun)

    print 'Bridging...'
    # build up the bridge between gravity and hydrodynamics
    grav_hydro = bridge.Bridge(use_threading=False)
    grav_hydro.add_system(gravity, (hydro,))
    grav_hydro.add_system(hydro, (gravity,))
    #grav_hydro.timestep = dt

    # set up channels for updating the particles
    channel_from_grav = gravity.particles.new_channel_to(Sun_and_Jupiter)
    channel_from_hydro = hydro.gas_particles.new_channel_to(disk_gas)

    a_Jup = list()
    e_Jup = list()
    disk_size = list()

    # start evolotuion
    print 'Start evolving...'
    times = quantities.arange(0.|units.yr, t_end+1*dt, 2*dt)
    model_time = 0.0 | units.yr
    while model_time <= t_end:
        # save the data for plots
        orbit = orbital_elements_from_binary(Sun_and_Jupiter, G=constants.G)
        a = orbit[2].value_in(units.AU)
        e = orbit[3]
        lr9 = return_L9_radius(disk_gas, Mstar, Rmin|units.AU)
        a_Jup.append(a)
        e_Jup.append(e)
        disk_size.append(lr9)

        if model_time.value_in(units.yr) % 50 == 0:
            print 'Time = %.1f yr:'%model_time.value_in(units.yr), \
                  'a = %.2f au, e = %.2f,'%(a, e), \
                  'disk size = %.2f au'%lr9
        # evolve the bridge system for one step
        model_time += dt
        grav_hydro.evolve_model(model_time, timestep=2*dt)
        channel_from_grav.copy()
        channel_from_hydro.copy()

        # add the 'sinked' mass to Jupiter & keep the sink particle along with Jupiter
        sink = hydro.particles[0]
        _ = hydro_sink_particles([sink], disk_gas)
        Jupiter = gravity.particles[1]
        Jupiter.mass += sink.mass
        sink.position = Jupiter.position
        sink.radius = a * (1-e) * ((1.0|units.MJupiter).value_in(units.MSun)/3.)**(1./3.) | units.au
        channel_from_grav.copy()
        channel_from_hydro.copy()

    gravity.stop()
    hydro.stop()

    return a_Jup, e_Jup, disk_size, times


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

    gravity = ph4(converter)
    gravity.particles.add_particles(Sun_and_Jupiter)
    gravity.timestep = dt

    # initialize the hydrodynamics(SPH) system
    # initialize the protoplanetary disk
    N = 1000
    Mstar = 1. | units.MSun
    Mdisk = 0.01 * Mstar
    Rmin = 2.
    Rmax = 100.

    np.random.seed(1)
    disk = ProtoPlanetaryDisk(N, convert_nbody=converter, densitypower=1.5, Rmin=Rmin, Rmax=Rmax, q_out=1.)
    disk_gas = disk.result
    disk_gas.h_smooth = 0.06 | units.AU

    # initialize the sink particle  
    sink = init_sink_particle(0.|units.MSun, Hill_radius0, Jupiter.position, Jupiter.velocity)
   
    sph = Fi(converter)
    sph.particles.add_particles(sink)
    sph.gas_particles.add_particles(disk_gas)
    sph.parameters.timestep = dt

    # bridge and evolve
    a_Jup, e_Jup, disk_size, times = gravity_hydro_bridge(gravity, sph, \
                                     [Sun_and_Jupiter, disk_gas], Rmin, t_end, dt)

    return a_Jup, e_Jup, disk_size, times


if __name__ == "__main__":
    t_end = 1000. | units.yr
    dt = 10. |units. yr
    a_Jup, e_Jup, disk_size, times = main(t_end, dt)
    

    #N_rho = 1000
    #L = 2 * Rmax
    #rho = make_map(sph, N=N_rho, L=L)
    
    fig = plt.figure(figsize=(8, 20))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    ax1.plot(times, a_Jup)
    ax2.plot(times, e_Jup)
    ax3.plot(times, disk_size)

    for ax in [ax1, ax2, ax3]:
        ax.set_xlabel('Time [yr]')
    
    ax1.set_ylabel('a [au]')
    ax2.set_ylabel('e')
    ax3.set_ylabel('disk size [au]')

    plt.show()
    
    #plt.imshow(numpy.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)),
    #              extent=[-L / 2, L / 2, -L / 2, L / 2], vmin=10, vmax=15)
    #plt.title(t_end)
    #plt.savefig('test.png')
















