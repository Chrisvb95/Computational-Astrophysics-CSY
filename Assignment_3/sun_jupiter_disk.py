import numpy as np
from amuse.lab import *
from amuse.units import quantities
from amuse.couple import bridge
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.orbital_elements import orbital_elements_from_binary

###BOOKLISTSTART1###
class BaseCode:
    def __init__(self, code, particles, eps=0|units.RSun):

        self.local_particles = particles
        m = self.local_particles.mass.sum()
        l = self.local_particles.position.length()
        self.converter = nbody_system.nbody_to_si(m, l)
        #self.converter = nbody_system.nbody_to_si(1.|units.MSun, 1.|units.AU)
        self.code = code(self.converter)
        self.code.parameters.epsilon_squared = eps**2

    def evolve_model(self, time):
        self.code.evolve_model(time)
    def copy_to_framework(self):
        self.channel_to_framework.copy()
    def get_gravity_at_point(self, r, x, y, z):
        return self.code.get_gravity_at_point(r, x, y, z)
    def get_potential_at_point(self, r, x, y, z):
        return self.code.get_potential_at_point(r, x, y, z)
    def get_timestep(self):
        return self.code.parameters.timestep
    @property
    def model_time(self):            
        return self.code.model_time
    @property
    def particles(self):
        return self.code.particles
    @property
    def total_energy(self):
        return self.code.kinetic_energy + self.code.potential_energy
    @property
    def stop(self):
        return self.code.stop
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
class Gravity(BaseCode):
    def __init__(self, code, particles, eps=0|units.RSun):
        BaseCode.__init__(self, code, particles, eps)
        self.code.particles.add_particles(self.local_particles)
        self.channel_to_framework \
            = self.code.particles.new_channel_to(self.local_particles)
        self.channel_from_framework \
            = self.local_particles.new_channel_to(self.code.particles)
        self.initial_total_energy = self.total_energy
###BOOKLISTSTOP2###

###BOOKLISTSTART3###
class Hydro(BaseCode):
    def __init__(self, code, particles, eps=0|units.RSun,
                 dt=None, Rbound=None):
        BaseCode.__init__(self, code, particles, eps)
        self.channel_to_framework \
            = self.code.gas_particles.new_channel_to(self.local_particles)
        self.channel_from_framework \
            = self.local_particles.new_channel_to(self.code.gas_particles)
        self.code.gas_particles.add_particles(particles)
        m = self.local_particles.mass.sum()
        l = self.code.gas_particles.position.length()
        if Rbound is None:
            Rbound = 10*l
        self.code.parameters.periodic_box_size = Rbound
        if dt is None:
            dt = 0.01*numpy.sqrt(l**3/(constants.G*m))
        self.code.parameters.timestep = dt/5.
        self.initial_total_energy = self.total_energy
    @property
    def total_energy(self):
        return self.code.kinetic_energy \
            + self.code.potential_energy \
            + self.code.thermal_energy


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


def init_sink_particle(m, r, position, velocity):
    sink = Particles(1)
    sink.mass = m
    sink.radius = r
    sink.x, sink.y, sink.z = position
    sink.vx, sink.vy, sink.vz = velocity
    return sink

    
def gravity_hydro_bridge(gravity, hydro, t_end, dt):
    
    model_time = 0 | units.yr
    filename = "gravhydro.hdf5"
    #write_set_to_file(stars.savepoint(model_time), filename, 'amuse')
    #write_set_to_file(gas, filename, 'amuse')

    gravhydro = bridge.Bridge(use_threading=False)
    gravhydro.add_system(gravity, (hydro,))
    gravhydro.add_system(hydro, (gravity,))
    gravhydro.timestep = 2*hydro.get_timestep()


    a_Jup = list()
    e_Jup = list()
    disk_size = list()
    
    # start evolotuion
    print 'Start evolving...'
    times = quantities.arange(0.|units.yr, t_end+dt, dt)
    while model_time < t_end:
        stars = gravity.local_particles
        orbit = orbital_elements_from_binary(stars, G=constants.G)
        a = orbit[2].value_in(units.AU)
        e = orbit[3]
        if model_time.value_in(units.yr)%50 == 0:
            print "Time:", model_time.in_(units.yr), \
                  "ae=", a ,e
        a_Jup.append(a)
        e_Jup.append(e)

        model_time += gravhydro.timestep
        gravhydro.evolve_model(model_time)
        gravity.copy_to_framework()
        hydro.copy_to_framework()

        sink = hydro.local_particles[0]
        _ = hydro_sink_particles([sink], hydro.local_particles[1:])
        Jupiter = gravity.local_particles[1]
        Jupiter.mass += sink.mass
        sink.position = Jupiter.position
        sink.radius = a * (1-e) * ((1.0|units.MJupiter).value_in(units.MSun)/3.)**(1./3.) | units.au
        gravity.copy_to_framework()
        hydro.copy_to_framework()
        write_set_to_file(stars.savepoint(model_time), filename, 'amuse')
        write_set_to_file(ism, filename, 'amuse')
        print "P=", model_time.in_(units.yr), gravity.particles.x.in_(units.au)

    gravity.stop()
    hydro.stop()

    return a_Jup, e_Jup, times


def main(t_end, dt):
    stars = init_sun_jupiter()
    Sun = stars[0]
    Jupiter = stars[1]
    orbit0 = orbital_elements_from_binary(stars, G=constants.G)
    a0 = orbit0[2].in_(units.AU)
    e0 = orbit0[3]
    Hill_radius0 = a0 * (1-e0) * ((1.0|units.MJupiter).value_in(units.MSun)/3.)**(1./3.)
    eps = 1 | units.RSun
    gravity = Gravity(ph4, stars, eps)

    N = 500
    Rmin = 1.0
    Rmax = 100.0
    converter = nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.AU)
    np.random.seed(1)
    disk = ProtoPlanetaryDisk(N, convert_nbody=converter, densitypower=1.5,\
                              Rmin=Rmin, Rmax=Rmax, q_out=1., discfraction=0.01)
    gas = disk.result
    gas.h_smooth = 0.06 | units.AU
    hydro = Hydro(Fi, gas, eps, dt)

    sink = init_sink_particle(0.|units.MSun, Hill_radius0, Jupiter.position, Jupiter.velocity)
    hydro.particles.add_particles(sink)

    a_Jup, e_Jup, times = gravity_hydro_bridge(gravity, hydro, t_end, dt)

    return a_Jup, e_Jup, times


if __name__ in ('__main__', '__plot__'):
    #o, arguments  = new_option_parser().parse_args()
    #gravity_hydro_bridge(**o.__dict__)
    a_Jup, e_Jup, times = main(t_end=1000.|units.yr, dt=10.|units.yr)



def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="int", default = 1000,
                      help="number of diagnostics time steps [%default]")
    result.add_option("-N", dest="Ngas", type="int", default = 1024,
                      help="number of gas particles [%default]")
    result.add_option("--Mprim", unit=units.MSun,
                      dest="Mprim", type="float", default = 3.2|units.MSun,
                      help="Mass of the primary star [%default]")
    result.add_option("--Msec", unit=units.MSun,
                      dest="Msec", type="float", default = 3.1|units.MSun,
                      help="Mass of the secondary star [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mgas", type="float", default = 1|units.MSun,
                      help="Mass of the gas [%default]")
    result.add_option("-R", unit=units.AU,
                      dest="Rgas", type="float", default = 10|units.AU,
                      help="Size of the gas distribution [%default]")
    result.add_option("-a", unit=units.AU,
                      dest="a", type="float", default = 1|units.AU,
                      help="initial orbital separation [%default]")
    result.add_option("-e", dest="ecc", type="float", default = 0.6,
                      help="initial orbital eccentricity [%default]")
    result.add_option("-t", unit=units.yr, 
                      dest="t_end", type="float", default = 10000|units.yr,
                      help="end time of the simulation [%default]")
    return result

