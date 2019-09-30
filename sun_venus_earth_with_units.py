###BOOKLISTSTART1###
from amuse.lab import Particles, units
from amuse.lab import SmallN, Huayno, nbody_system
    

v1 = 0.2009656237
v2 = 0.2431076328
T = 20 | nbody_system.time

def sun_venus_and_earth():
    particles = Particles(3)
    sun = particles[0]
    sun.mass = 1.0 | nbody_system.mass
    #sun.radius = 1.0 | units.RSun
    sun.position = (-1, 0., 0.) | nbody_system.length
    sun.velocity = (v1, v2, 0) | nbody_system.speed
    
    venus = particles[1]
    venus.mass = 1.0 | nbody_system.mass
    #venus.radius = 1.0 | units.RSun
    venus.position = (1., 0., 0.) | nbody_system.length
    venus.velocity = (v1, v2, 0.) | nbody_system.speed
    
    earth = particles[2]
    earth.mass = 0.5 | nbody_system.mass
    #earth.radius = 1.0 | units.RSun
    earth.position = (0., 0., 0.) | nbody_system.length
    earth.velocity = (-4 * v1, -4 * v2, 0.) | nbody_system.speed
    
    particles.move_to_center()
    return particles
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def integrate_solar_system(particles, end_time):
    from amuse.units import constants, units
    #print(particles)
    #convert_body = nbody_system.nbody_to_si(particles.mass.sum(),
    #                                        particles[1].position.length())

    gravity = Huayno()
    gravity.particles.add_particles(particles)
    sun = gravity.particles[0]
    venus = gravity.particles[1]
    earth = gravity.particles[2]
    
    
    x_earth = [] | nbody_system.length
    y_earth = [] | nbody_system.length
    x_venus = [] | nbody_system.length
    y_venus = [] | nbody_system.length
    x_sun = [] | nbody_system.length
    y_sun = [] | nbody_system.length
    
    while gravity.model_time < end_time:
        
        gravity.evolve_model(gravity.model_time + (0.001 | nbody_system.time))
        x_sun.append(sun.x)
        y_sun.append(sun.y)
        x_earth.append(earth.x)
        y_earth.append(earth.y)
        x_venus.append(venus.x)
        y_venus.append(venus.y)
    gravity.stop()
    return x_earth, y_earth, x_venus, y_venus, x_sun, y_sun
###BOOKLISTSTOP2###

###BOOKLISTSTART3###
def plot_track(xs,ys,xe,ye,xv,yv, output_filename):

    from matplotlib import pyplot
    figure = pyplot.figure(figsize=(10, 10))
    pyplot.rcParams.update({'font.size': 30})
    plot = figure.add_subplot(1,1,1)
    ax = pyplot.gca()
    ax.minorticks_on() 
    ax.locator_params(nbins=3)

    x_label = 'x'
    y_label = 'y'
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)

    
    plot.plot(xe.value_in(nbody_system.length), ye.value_in(nbody_system.length), color = 'b')
    plot.plot(xv.value_in(nbody_system.length), yv.value_in(nbody_system.length), color = 'r')
    plot.plot(xs.value_in(nbody_system.length), ys.value_in(nbody_system.length), color = 'k')
        
    plot.set_xlim(-1.3, 1.3)
    plot.set_ylim(-1.3, 1.3)

    save_file = 'sun_venus_earth_test.png'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    pyplot.show()
###BOOKLISTSTOP3###

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-o", 
                      dest="output_filename", default ="SunVenusEarth",
                      help="output filename [%default]")
    return result
    
if __name__ in ('__main__','__plot__'):
    o, arguments  = new_option_parser().parse_args()

    particles = sun_venus_and_earth()
    xs,ys, xe,ye, xv,yv = integrate_solar_system(particles, 2 * T)
    plot_track(xs, ys, xe, ye, xv, yv, o.output_filename)
    
