###BOOKLISTSTART1###
from amuse.lab import Particles, units 

print(type(units.RSun))

v1 = 0.201
v2 = 0.243

def sun_venus_and_earth():
    particles = Particles(3)
    sun = particles[0]
    sun.mass = 1.0 | units.MSun
    sun.radius = 1e-5 | units.m
    sun.position = (-1, 0, 0) | units.AU
    sun.velocity = (v1, v2, 0) | units.kms

    venus = particles[1]
    venus.mass = 1.0 | units.MSun
    venus.radius = 1e-5 | units.m
    venus.position = (1, 0, 0) | units.AU
    venus.velocity = (v1, v2, 0) | units.kms

    earth = particles[2]
    earth.mass = 0.5 | units.MSun
    earth.radius = 1e-5 | units.m
    earth.position = (0, 0, 0) | units.AU
    earth.velocity = (-4 * v1, -4 * v2, 0) | units.kms

    particles.move_to_center()
    return particles
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def integrate_solar_system(particles, end_time):
    from amuse.lab import Huayno, nbody_system
    #print(particles[1].position())
    convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun,
                                             1.0 | units.AU)
    print(convert_nbody)

    gravity = Huayno(convert_nbody)
    gravity.particles.add_particles(particles)
    venus = gravity.particles[1]
    earth = gravity.particles[2]
    sun = gravity.particles[0]
    
    
    x_earth = [] | units.AU
    y_earth = [] | units.AU
    x_venus = [] | units.AU
    y_venus = [] | units.AU
    x_sun = [] | units.AU
    y_sun = [] | units.AU
    
    #print('Error is here')
    #print(gravity.model_time + 1.|units.s)

    while gravity.model_time < end_time:
        
        gravity.evolve_model(gravity.model_time + (1 | units.day))
        x_earth.append(earth.x)
        y_earth.append(earth.y)
        x_venus.append(venus.x)
        y_venus.append(venus.y)
        x_sun.append(sun.x)
        y_sun.append(sun.y)
    gravity.stop()
    return x_earth, y_earth, x_venus, y_venus, x_sun, y_sun
###BOOKLISTSTOP2###

###BOOKLISTSTART3###
def plot_track(xe,ye,xv,yv,xs,ys, output_filename):

    from matplotlib import pyplot
    figure = pyplot.figure(figsize=(10, 10))
    pyplot.rcParams.update({'font.size': 30})
    plot = figure.add_subplot(1,1,1)
    ax = pyplot.gca()
    ax.minorticks_on() 
    ax.locator_params(nbins=3)

    x_label = 'x [au]'
    y_label = 'y [au]'
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)

    #plot.scatter([0.0], [0.0], color='y', lw=8)
    plot.plot(xe.value_in(units.m), ye.value_in(units.m), color = 'b')
    plot.plot(xv.value_in(units.m), yv.value_in(units.m), color = 'r')
    plot.plot(xs.value_in(units.m), ys.value_in(units.m), color = 'k')
    
    plot.set_xlim(-1.3, 1.3)
    plot.set_ylim(-1.3, 1.3)

    save_file = 'sun_venus_earth.png'
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
    xe,ye, xv,yv, xs,ys = integrate_solar_system(particles, 2 | units.yr)
    plot_track(xe, ye, xv, yv, xs, ys, o.output_filename)
    
