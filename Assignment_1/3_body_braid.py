###BOOKLISTSTART1###
from amuse.lab import Particles, units
from amuse.lab import SmallN, Huayno, Hermite, nbody_system
import numpy as np

v1_init = 0.2009656237
v2_init = 0.2431076328
init = np.array([-1,0,v1_init,v2_init])   
P = 19.0134164290
T = P*1 | nbody_system.time

def sigma(x1,x2,v1,v2):
    current = np.array([x1, x2, v1, v2])
    return np.sqrt(np.sum((current - init) ** 2))     
    

def body_init():
    particles = Particles(3)
    body0 = particles[0]
    body0.mass = 1.0 | nbody_system.mass
    body0.position = (-1, 0., 0.) | nbody_system.length
    body0.velocity = (v1_init, v2_init, 0) | nbody_system.speed
    
    body1 = particles[1]
    body1.mass = 1.0 | nbody_system.mass
    body1.position = (1., 0., 0.) | nbody_system.length
    body1.velocity = (v1_init, v2_init, 0.) | nbody_system.speed
    
    body2 = particles[2]
    body2.mass = 0.5 | nbody_system.mass
    body2.position = (0., 0., 0.) | nbody_system.length
    body2.velocity = (-4 * v1_init, -4 * v2_init, 0.) | nbody_system.speed
    
    particles.move_to_center()
    return particles
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def integrate_solar_system(particles, end_time):
    from amuse.units import constants, units
   
    gravity = Hermite()
    gravity.particles.add_particles(particles)
    body0 = gravity.particles[0]
    body1 = gravity.particles[1]
    body2 = gravity.particles[2]
    
    
    x_body0 = [] | nbody_system.length
    y_body0 = [] | nbody_system.length
    x_body1 = [] | nbody_system.length
    y_body1 = [] | nbody_system.length
    x_body2 = [] | nbody_system.length
    y_body2 = [] | nbody_system.length

    vx_body0 = [] | nbody_system.speed
    vy_body0 = [] | nbody_system.speed

    x_body0.append(body0.x)
    y_body0.append(body0.y)
    x_body1.append(body1.x)
    y_body1.append(body1.y)
    x_body2.append(body2.x)
    y_body2.append(body2.y)    

    vx_body0.append(body0.vx)
    vy_body0.append(body0.vy)

    while gravity.model_time < end_time:
        
        gravity.evolve_model(gravity.model_time + (0.005 | nbody_system.time))
        x_body0.append(body0.x)
        y_body0.append(body0.y)
        x_body1.append(body1.x)
        y_body1.append(body1.y)
        x_body2.append(body2.x)
        y_body2.append(body2.y)

        vx_body0.append(body0.vx)
        vy_body0.append(body0.vy)

    gravity.stop()
    return x_body0, y_body0, x_body1, y_body1, x_body2, y_body2, vx_body0, vy_body0
###BOOKLISTSTOP2###

###BOOKLISTSTART3###
def plot_track(x0,y0,x1,y1,x2,y2,vx0,vy0, output_filename):

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
    
    print('Printing x-values')
    print(np.array(x0.value_in(nbody_system.length)))
    print(len(x0.value_in(nbody_system.length)))
    
    sig = sigma(np.array(x0.value_in(nbody_system.length))[-2],
                np.array(y0.value_in(nbody_system.length))[-2],
                np.array(vx0.value_in(nbody_system.speed))[-2],
                np.array(vy0.value_in(nbody_system.speed))[-2])
    
    print('Sigma:',sig) 
    
    plot.plot(x0.value_in(nbody_system.length), y0.value_in(nbody_system.length), color = 'b')
    plot.plot(x1.value_in(nbody_system.length), y1.value_in(nbody_system.length), color = 'r')
    plot.plot(x2.value_in(nbody_system.length), y2.value_in(nbody_system.length), color = 'k')
        
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
                      dest="output_filename", default ="Assignment_1",
                      help="output filename [%default]")
    return result
    
if __name__ in ('__main__','__plot__'):
    o, arguments  = new_option_parser().parse_args()

    particles = body_init()
    x0,y0, x1,y1, x2,y2, vx0,vy0 = integrate_solar_system(particles, T)
    plot_track(x0, y0, x1, y1, x2, y2, vx0, vy0, o.output_filename)
    
