from amuse.lab import *
from amuse.community.simplex.interface import SimpleXInterface, SimpleX, \
     SimpleXSplitSet
from amuse.ext.protodisk import ProtoPlanetaryDisk
from matplotlib import pyplot as plt
import numpy as np

def mu(X = None, Y = 0.25, Z = 0.02, x_ion = 0.1):
    """
    Compute the mean molecular weight in kg (the average weight of
    particles in a gas) X, Y, and Z are the mass fractions of
    Hydrogen, of Helium, and of metals, respectively.  x_ion is the
    ionisation fraction (0 < x_ion < 1), 1 means fully ionised.
    """
    if X is None:
        X = 1.0 - Y - Z
    elif abs(X + Y + Z - 1.0) > 1e-6:
        raise Exception("Error in calculating mu: mass "
                         + "fractions do not sum to 1.0")
    return constants.proton_mass / (X*(1.0+x_ion) + Y*(1.0+2.0*x_ion)/4.0
                                     + Z*x_ion/2.0)

class RadHydro:
    def __init__(self, rad, hydro, star, disk):

        self.time = 0|units.day
        self.star = star
        self.disk = disk

        disk.r2 = disk.x**2 + disk.y**2
        disk = disk.sorted_by_attributes("r2")
        Rmax = disk.r2.max().sqrt()
        print "MaxR=", Rmax.in_(units.AU)

        self.hydro = hydro(nbody_system.nbody_to_si(self.disk.mass.sum(), Rmax),
                           number_of_workers=4)
        self.hydro.parameters.epsilon_squared = (10|units.AU)**2
        
        self.hydro.gas_particles.add_particles(self.disk)
        self.hydro.gas_particles.new_channel_to(self.disk)
        self.hydro.dm_particles.add_particles(self.star)
        self.hydro.dm_particles.new_channel_to(self.star)

        self.hydro_to_star = self.hydro.dm_particles.new_channel_to(self.star)
        self.hydro_to_disk = self.hydro.gas_particles.new_channel_to(self.disk)
        self.star_to_hydro = self.star.new_channel_to(self.hydro.dm_particles)
        self.disk_to_hydro = self.disk.new_channel_to(self.hydro.gas_particles)
            
        self.hydro.evolve_model(1|units.s)
        self.hydro_to_star.copy()
        self.hydro_to_disk.copy()

        self.rad = rad()
        for si in self.star:
            if si.mass>=5|units.MSun:
                self.rad.src_particles.add_particle(si)
        self.rad.gas_particles.add_particles(self.disk)
        self.rad.parameters.box_size=2.01*Rmax
        self.rad.parameters.timestep=1|units.day
        self.rad.set_source_Teff(star.temperature)

        self.rad_to_disk = self.rad.gas_particles.new_channel_to(self.disk,
                                                    attributes=["xion", "u"])
        self.star_to_rad = self.star.new_channel_to(self.rad.src_particles,
                                                    attributes=["x", "y", "z"])
        self.disk_to_rad = self.disk.new_channel_to(self.rad.gas_particles,
                                                    attributes=["x", "y", "z"])

        self.rad.stop()
        self.index = 0

    def write_file(self):
        self.index += 1
        filename = "hydro_disk_i{0:04}.amuse".format(self.index)
        write_set_to_file(self.star, filename, "amuse", append_to_file=False)
        write_set_to_file(self.disk, filename, "amuse")
        
    def evolve_model(self, model_time):
        dt = model_time - self.time
        self.old_time = self.time
        self.time += dt/2.
      
        self.time += dt/2.

        self.disk_to_hydro.copy()
        self.star_to_hydro.copy()
        self.hydro.evolve_model(self.time)
        self.hydro_to_disk.copy()
        self.hydro_to_star.copy()

        print "RT done at time:", self.time.in_(units.day)

    def print_diagnostics(self):
        umin = self.disk.u.min()
        umean = self.disk.u.mean()
        umax = self.disk.u.max()
        Tmin =  mu() / constants.kB * umax
        Tmean =  mu() / constants.kB * umean
        Tmax =  mu() / constants.kB * umin

        print "Time=", self.time.in_(units.day)
        #print "Ionization:", self.disk.xion.min(), self.disk.xion.mean(), \
        #      self.disk.xion.max()
        print "Intenal energy:", umin, umean, umax
        print "Temperature:", Tmin, Tmean, Tmax
        print "Density:", self.disk.density.min().in_(units.amu/units.cm**3), \
              self.disk.density.mean().in_(units.amu/units.cm**3), \
              self.disk.density.max().in_(units.amu/units.cm**3)
        print "scaleheight:", abs(self.disk.z.value_in(units.AU)).mean()

    
    def return_L9_radius(self,Mstar,Rmin): 
        converter=nbody_system.nbody_to_si(Mstar, Rmin)     
        self.lr,self.mr = self.disk.LagrangianRadii(converter)
        return self.lr[7].value_in(units.AU)

    def stop(self):
        self.hydro.stop()


def new_disk(Mstar=1|units.MSun, Ndisk=500,
             Rmin=1.0|units.AU, Rmax=100.0|units.AU):

    Mdisk = 0.01*Mstar #setting the mass of the disk to 1% of star

    # Initialise disk
    converter=nbody_system.nbody_to_si(Mstar, Rmin)
    disk = ProtoPlanetaryDisk(Ndisk, convert_nbody=converter,
                              densitypower=1.5, Rmin=1, Rmax=Rmax/Rmin,
                              q_out=1.0, discfraction=Mdisk/Mstar).result
    com = disk.center_of_mass()

    return disk

def evolve_star(Mstar, tstar):

    print 'Particles:', Particles(2)
        
    star = Particles(1)
    star[0].mass = Mstar
    stellar = SeBa()
    stellar.particles.add_particle(star)
    stellar.evolve_model(tstar)
    star.mass = stellar.particles.mass
    star.position = (0, 0, 0) |units.AU
    star.velocity = (0, 0, 0) |units.kms
    star.luminosity = stellar.particles.luminosity/(20. | units.eV)
    star.temperature = stellar.particles.temperature
    star.flux = star.luminosity
    star.rho = 1.0|(units.g/units.cm**3)
    star.xion = 0.0 #ionization_fraction
    star.u = (9. |units.kms)**2 #internal_energy

    stellar.stop()
    print 'stars', star
    return star
 
def hydro_disk(Mstar=1|units.MSun, Ndisk=500,
               Rmin=1.0|units.AU, Rmax=100.0|units.AU,
               t_end=1000|units.yr, n_steps=10):

    star = evolve_star(Mstar, t_end)
    disk = new_disk()
    lr9 = []

    radhydro = RadHydro(SimpleXSplitSet, Gadget2, star, disk)
    radhydro.write_file()
    lr9.append(radhydro.return_L9_radius(Mstar,Rmin))

    dt = t_end/float(n_steps)
    print "dt=", dt.in_(units.day)
    time = 0 | units.day

    while time<t_end:
        time += dt
        radhydro.evolve_model(time)
        radhydro.print_diagnostics()
        radhydro.write_file()
        lr9.append(radhydro.return_L9_radius(Mstar,Rmin))
    radhydro.stop()

    print dt.value_in(units.yr)
    print t_end.value_in(units.yr)
    print lr9

    
    t = np.arange(0,(t_end+dt).value_in(units.yr),dt.value_in(units.yr))
    print(t)
    
    plt.scatter(t,lr9)
    plt.xlabel('Time (years)')
    plt.ylabel('Disk size (AU)')
    plt.savefig('disk_radius_time_plot.png')    
    plt.show()

if __name__ in '__main__':
    hydro_disk()



