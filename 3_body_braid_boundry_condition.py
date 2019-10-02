from amuse.lab import Particles, units
from amuse.lab import Hermite, nbody_system
import numpy as np

#Initial Condition
v1_init = 0.9911981217
v2_init = 0.7119472124
m3 = 4
Cycle = 1
lower_boundry = 1e-1
upper_boundry = 25
#P = 17.6507807837 # uncomment, if you want to use time boundry
#T = P * 1 | nbody_system.time # uncomment, if you want to use time boundry


def sigma(x1,x2, v1,v2): #Calculates the return proximity function at any point
    init = np.array([-1, 0, v1_init, v2_init])
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
    body2.mass = m3 | nbody_system.mass
    body2.position = (0., 0., 0.) | nbody_system.length
    body2.velocity = (-2 * v1_init / m3, -2 * v2_init / m3, 0.) | nbody_system.speed
    
    particles.move_to_center()
    return particles

def integrate_3body_system(particles):#, end_time):
    from amuse.units import constants, units
    import matplotlib.pyplot as plt
   
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
    
    # Plotting the return proximity funcion
    fig = plt.figure(figsize=(8,16))
    ax1 = fig.add_subplot(211, xlabel = 'Integraion time step', ylabel = 'Return proximity function', 
				title='Deviation from the initial point after each period')
    ax2 = fig.add_subplot(212, xlabel = 'Intergration time step', ylabel = 'Return proximity function',
				title='Return proximity function per time step')
	
    sig = sigma(np.array(body0.x.value_in(nbody_system.length)),
	        np.array(body0.y.value_in(nbody_system.length)),
	        np.array(body0.vx.value_in(nbody_system.speed)),
	        np.array(body0.vy.value_in(nbody_system.speed)))
    c = 0 # cycle
    sig_array = []
    time_step = 0
    sig_O_x = [0]
    sig_O_y = [sig]
    while (c < Cycle and sig < upper_boundry): # gravity.model_time < end_time:
	
	gravity.evolve_model(gravity.model_time + (0.005 | nbody_system.time))
	x_body0.append(body0.x)
	y_body0.append(body0.y)
	x_body1.append(body1.x)
	y_body1.append(body1.y)
	x_body2.append(body2.x)
	y_body2.append(body2.y)

	vx_body0.append(body0.vx)
	vy_body0.append(body0.vy)
	
	sig = sigma(np.array(body0.x.value_in(nbody_system.length)),
	        np.array(body0.y.value_in(nbody_system.length)),
	        np.array(body0.vx.value_in(nbody_system.speed)),
	        np.array(body0.vy.value_in(nbody_system.speed)))
	#print(sig)
	sig_array.append(sig)
	if time_step > 2: 
		slope = (sig_array[-2] - sig_array[-3]) * (sig_array[-1] - sig_array[-2])
	else: slope = 1
		
	if (sig < lower_boundry) and (slope < 0):
		c += 1
		print(x_body0.value_in(nbody_system.length)[-1], y_body0.value_in(nbody_system.length)[-1],
		vx_body0.value_in(nbody_system.speed)[-1], vy_body0.value_in(nbody_system.speed)[-1])
		print(sig)
		sig_O_x.append(time_step)
		sig_O_y.append(sig)
		
	time_step += 1
    gravity.stop()
    ax1.plot(sig_O_x, sig_O_y, color = 'k')

    ax2.plot(sig_array)
    ax2.plot(range(time_step), np.zeros(time_step))
    ax2.scatter(sig_O_x, sig_O_y, color = 'k')

    plt.show()
    return x_body0, y_body0, x_body1, y_body1, x_body2, y_body2
    
def plot_track(x0,y0,x1,y1,x2,y2, output_filename):

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
    
    plot.plot(x0.value_in(nbody_system.length), y0.value_in(nbody_system.length), color = 'b')
    plot.plot(x1.value_in(nbody_system.length), y1.value_in(nbody_system.length), color = 'r')
    plot.plot(x2.value_in(nbody_system.length), y2.value_in(nbody_system.length), color = 'k')
        
    plot.set_xlim(-1.3, 1.3)
    plot.set_ylim(-1.3, 1.3)

    save_file = 'sun_venus_earth_test.png'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    pyplot.show()

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
    x0,y0, x1,y1, x2,y2= integrate_3body_system(particles)#, T)
    plot_track(x0, y0, x1, y1, x2, y2, o.output_filename)
    
