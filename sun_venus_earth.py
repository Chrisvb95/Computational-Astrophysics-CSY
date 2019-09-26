###BOOKLISTSTART1###
from amuse.lab import Particles, units

v1 = 0.2009656237
v2 = 0.2431076328

def sun_venus_and_earth():
    particles = Particles(3)
    sun = particles[0]
    sun.mass = 1.0 | units.kg
    #sun.radius = 1e-2 | units.m
    sun.position = (-1, 0, 0) | units.m
    sun.velocity = (v1, v2, 0) | (units.m/units.s)

    venus = particles[1]
    venus.mass = 1.0 | units.kg
    #venus.radius = 1e-2 | units.m
    venus.position = (1, 0, 0) | units.m
    venus.velocity = (v1, v2, 0) | (units.m/units.s)

    earth = particles[2]
    earth.mass = 0.5 | units.kg
    #earth.radius = 1e-2 | units.m
    earth.position = (0, 0, 0) | units.m
    earth.velocity = (-4 * v1, -4 * v2, 0) | (units.m/units.s)

    particles.move_to_center()
    return particles
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def integrate_solar_system(particles, end_time):
    from amuse.lab import Huayno, nbody_system
    from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
    from amuse.units import constants, units
    #print(particles)
    #converter_body = nbody_system.nbody_to_si(particles.mass.sum(),
    #                                       	particles[1].position.length())
    converter_body = nbody_system.nbody_to_generic(1|units.kg, 1|units.m)
    #print(convert_nbody)
    converter = ConvertBetweenGenericAndSiUnits(particles)#????????????????????
    #gravity = Huayno(converter.to_si(particles))
    gravity = Huayno(converr_body)
    gravity.particles.add_particles(particles)
    sun = gravity.particles[0]
    venus = gravity.particles[1]
    earth = gravity.particles[2]
    
    
    x_earth = [] | units.m
    y_earth = [] | units.m
    x_venus = [] | units.m
    y_venus = [] | units.m
    x_sun = [] | units.m
    y_sun = [] | units.m


    while gravity.model_time < end_time:
        
        gravity.evolve_model(gravity.model_time + (1 | units.s))
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

    x_label = 'x [m]'
    y_label = 'y [m]'
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)

    #plot.scatter([0.0], [0.0], color='y', lw=8)
    plot.plot(xe.value_in(units.m), ye.value_in(units.m), color = 'b')
    #print(xe.value_in(units.m), ye.value_in(units.m))
    plot.plot(xv.value_in(units.m), yv.value_in(units.m), color = 'r')
    plot.plot(xs.value_in(units.m), ys.value_in(units.m), color = 'k')
    
    #plot.set_xlim(-1.3, 1.3)
    #plot.set_ylim(-1.3, 1.3)

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
    xs,ys, xe,ye, xv,yv = integrate_solar_system(particles, 50 | units.s)
    plot_track(xs, ys, xe, ye, xv, yv, o.output_filename)
    
