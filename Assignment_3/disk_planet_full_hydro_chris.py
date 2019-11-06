# -*- coding: ascii -*-
"""
Task 3: Sun + Jupiter(sink) + disk
"""
#from __future__ import print_function
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
import pickle as pkl

from amuse.lab import *
from amuse.units import units, nbody_system, quantities
from amuse.couple import bridge
from amuse.community.fi.interface import Fi
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ext.orbital_elements import orbital_elements_from_binary

def plot_map(sph, Sun_and_Jupiter, Pstar, title, N=100, L=1, show=True):

    print('Plotting', title)

    L = 400
   
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
    #Sr = Sun_and_Jupiter[0].radius.value_in(units.AU)
    Sr = 1    
    sun = Circle((Sx,Sy),Sr,color='y')

    Jx = Sun_and_Jupiter[1].x.value_in(units.AU)
    Jy = Sun_and_Jupiter[1].y.value_in(units.AU)
    Jz = Sun_and_Jupiter[1].z.value_in(units.AU)
    #Jr = Sun_and_Jupiter[1].radius.value_in(units.AU)
    Jr = 0.5     
    jup = Circle((Jx,Jy),Jr,color='orange')

    Starx = Pstar.x.value_in(units.AU)
    Stary = Pstar.y.value_in(units.AU)
    Starz = Pstar.z.value_in(units.AU)
    Starr = 1.0
    star = Circle((Starx,Stary),Starr,color='r')

    rho, rhovx, rhovy, rhovz, rhoe = sph.get_hydro_state_at_point(
        x, y, z, vx, vy, vz)
    rho = rho.reshape((N + 1, N + 1))

    fig,ax = plt.subplots(1,figsize=(8, 8))    
    ax.imshow(np.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)),
                  extent=[-L / 2, L / 2, -L / 2, L / 2], vmin=10, vmax=15)
    ax.add_patch(sun)
    ax.add_patch(jup)
    ax.add_patch(star)
    #plt.title(title)
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
    sink.mass = 0. | units.MJupiter
    sink.radius = r
    sink.x, sink.y, sink.z = position
    sink.vx, sink.vy, sink.vz = velocity

    return sink

def init_passing_star(pericenter, e, M):
  
    Pstar = Particles(1)
    Pstar.mass = M
    Pstar.radius = 1.0 | units.RSun
    Pstar.position = (200, -400, 0) | units.AU
    Pstar.eccentricity = e
    Pstar.argument_of_pericenter = 200 | units.AU
    ds = 400/100
    Pstar.velocity = (0, ds, 0) | (units.AU/units.yr)  
  
    return Pstar
    
def evolve(Sun_Jupiter, disk_gas, sink, Pstar, dt_gravity, dt_sph, dt_diagnostic, dt_bridge, tend):

    Sun = Sun_Jupiter[0]
    Jupiter = Sun_Jupiter[1]

    
    # Initialising the SPH code
    sph_converter = nbody_system.nbody_to_si(1|units.MSun, 1|units.AU)    
    sph = Fi(sph_converter, mode="openmp")
    sph.gas_particles.add_particles(disk_gas)
    sph.dm_particles.add_particle(Sun)
    sph.dm_particles.add_particle(Jupiter)
    sph.dm_particles.add_particle(Pstar)
 
    # Set up channels for updating the particles
    sph_to_disk = sph.gas_particles.new_channel_to(disk_gas)
    sph_to_Sun_Jupiter = sph.dm_particles.new_channel_to(Sun_Jupiter)
    sph_to_Pstar = sph.dm_particles.new_channel_to(Pstar)
    
    Sun_Jupiter_to_sph = Sun_Jupiter.new_channel_to(sph.dm_particles)
    disk_to_sph = disk_gas.new_channel_to(sph.gas_particles)
    Pstar_to_sph = Pstar.new_channel_to(sph.dm_particles)
 
    # Preparing lists for data-recording
    a_Jup = []
    e_Jup = []
    disk_size = []
    accreted_mass = []
    accreted_mass.append((sink.mass).value_in(units.MJupiter)[0])
    sink0_mass = 0 | units.MJupiter
    
    # Start evolution
    print 'Start evolving...'
    times = quantities.arange(0.|units.yr, tend, dt_diagnostic)
    for i,t in enumerate(times):
	    
        # Save the data for plots
        orbit = orbital_elements_from_binary(Sun_Jupiter, G=constants.G)
        a = orbit[2].value_in(units.AU)
        e = orbit[3]
        lr9 = return_L9_radius(disk_gas, Sun_Jupiter[0].mass, Rmin)
        a_Jup.append(a)
        e_Jup.append(e)
        disk_size.append(lr9)
 
        # Plotting system        
        print 'Time = %.1f yr:'%t.value_in(units.yr), \
                  'a = %.2f au, e = %.2f,'%(a, e), \
                  'disk size = %.2f au'%lr9
        plot_map(sph, Sun_Jupiter, Pstar, 
		'distribution_plot_passing_star_full_hydro_4k/{0}.png'.format(int(t.value_in(units.yr))),show=False)
    
        # Evolve the bridge system for one step
        sph.evolve_model(t, dt_diagnostic)
        sph_to_disk.copy()
        sph_to_Sun_Jupiter.copy()
        sph_to_Pstar.copy()

        # Calculating accreted mass in new position
        sink.position = Jupiter.position
        sink.radius = Hill_radius(a, e, Sun_Jupiter) | units.AU       
        removed_particles = hydro_sink_particles(sink, disk_gas) 
        Jupiter.mass += sink.mass - sink0_mass
        accreted_mass.append((sink.mass).value_in(units.MJupiter)[0])            
        sink0_mass = sink.mass.copy()
        Sun_Jupiter_to_sph.copy()
        disk_to_sph.copy()

	Sun_Jupiter_to_sph.copy()
	Pstar_to_sph.copy()
    sph.stop()

    return times, a_Jup, e_Jup, disk_size, accreted_mass


def init_body_solar_disk_planetary(Ndisk, Mdisk, Rmin, Rmax):

    print 'Initializing...'
    #converter = nbody_system.nbody_to_si(1.|units.MJupiter, 1.|units.AU)

    # Initialize the solar system
    Sun_Jupiter = init_sun_jupiter()
    Sun = Sun_Jupiter[0]
    Jupiter = Sun_Jupiter[1]

    orbit = orbital_elements_from_binary(Sun_Jupiter, G=constants.G)
    a = orbit[2].in_(units.AU)
    e = orbit[3]
    hr = Hill_radius(a, e, Sun_Jupiter)

    # Initialize the proto-planetary disk  
    np.random.seed(1)
    disk_converter = nbody_system.nbody_to_si(Sun.mass, Rmin)
    disk = ProtoPlanetaryDisk(Ndisk, convert_nbody=disk_converter, 
	densitypower=1.5, Rmin=1.0, Rmax=Rmax/Rmin, q_out=1.0,
	discfraction=Mdisk/Sun.mass)
    disk_gas = disk.result

    # Initialize the sink particle  
    sink = init_sink_particle(0.|units.MSun, hr, Jupiter.position, Jupiter.velocity)
    
    # Initizalize the passing star
    e = 1.2
    Mstar = 2.0 | units.MSun
    pericenter = 200.0 | units.AU
    Pstar = init_passing_star(pericenter, e, Mstar)

    return Sun_Jupiter, disk_gas, sink, Pstar


if __name__ == "__main__":
    
    # Setting parameters
    tend = 500. | units.yr
    dt_gravity = 0.1 | units.yr
    dt_sph = 0.01 |units.yr
    dt_diagnostic = 1.0 | units.yr
    dt_bridge = 0.1 | units.yr
    # Setting proto-disk parameters
    Ndisk = 4000
    Mdisk = 0.01 | units.MSun
    Rmin = 2. | units.AU
    Rmax = 100. | units.AU

    Sun_Jupiter, disk_gas, sink, Pstar = init_body_solar_disk_planetary(Ndisk, Mdisk, Rmin, Rmax)
    t, a_Jup, e_Jup, disk_size, accreted_mass= evolve(Sun_Jupiter, 
			disk_gas, sink, Pstar, dt_gravity, dt_sph, dt_diagnostic, dt_bridge,tend)

    sim_output = np.column_stack((t.value_in(units.yr), a_Jup, e_Jup, disk_size, accreted_mass))
    np.savetxt('sim_output.csv',sim_output,delimiter=',')
















