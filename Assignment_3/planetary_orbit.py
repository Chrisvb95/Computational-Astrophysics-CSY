from amuse.lab import Particles, units, nbody_system, constants
from amuse.units import quantities
import numpy as np
import matplotlib.pyplot as plt
from amuse.community.ph4.interface import ph4
from amuse.ext.orbital_elements import orbital_elements_from_binary

def body_init_sun_jupiter():
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

def evolve_sun_jupiter(particles, tend, dt):
    
    SunJupiter = Particles()
    SunJupiter.add_particle(particles[0])
    SunJupiter.add_particle(particles[1])

    converter=nbody_system.nbody_to_si(particles.mass.sum(), particles[1].position.length())
    gravity = ph4(converter)
    gravity.particles.add_particles(particles)

    channel_from_to_SunJupiter = gravity.particles.new_channel_to(SunJupiter)
    
    semi_major_axis = []
    eccentricity = []
    times = quantities.arange(0|units.yr, tend, dt)
    for i,t in enumerate(times):
        print "Time=", t.in_(units.yr)
	channel_from_to_SunJupiter.copy()
	orbital_elements = orbital_elements_from_binary(SunJupiter, G=constants.G)
	a = orbital_elements[2]
	e = orbital_elements[3]
	semi_major_axis.append(a.value_in(units.AU))
	eccentricity.append(e)
        gravity.evolve_model(t, timestep=dt)
    
    # Plot
    
    fig = plt.figure(figsize = (8, 6))
    ax1 = fig.add_subplot(211, xlabel = 'time(yr)', ylabel = 'semi major axis (AU)')
    ax1.plot(times.value_in(units.yr), semi_major_axis)
    ax2 = fig.add_subplot(212, xlabel = 'time(yr)', ylabel = 'eccentriity')
    ax2.plot(times.value_in(units.yr), eccentricity)
    plt.show()
    
    gravity.stop()

if __name__ in ('__main__','__plot__'):

    dt = 1.0 | units.yr
    tend = 1000.0 | units.yr

    particles = body_init_sun_jupiter()
    evolve_sun_jupiter(particles, tend, dt) 

