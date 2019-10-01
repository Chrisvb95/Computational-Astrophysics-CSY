###BOOKLISTSTART1###
from amuse.lab import Particles, units
from amuse.lab import Hermite, Huayno, nbody_system


pattern = ['I.A_4(0.5)', 'I.B_68(0.5)', 'I.B_59(0.75)', 'I.A_1(2)', 'II.D_2(2)', 'I.A_2(4)']
m3_list = [0.5, 0.5, 0.75, 2., 2., 4.]
v1_list = [0.2009656237, 0.2138410831, 0.4101378717, 0.6649107583, 0.3057224330, 0.9911981217]
v2_list = [0.2431076328, 0.0542938396, 0.1341894173, 0.8324167864, 0.5215124257, 0.7119472124]
T_list = [19.01341164290, 83.8472907647, 121.0976361440, 12.6489061509, 8.8237067653, 17.6507807837] | nbody_system.time

yrange = [1.5, 0.3, 0.6, 1.5, 0.8, 0.8]

t_step = 0.1
n_T = 100
    

def body_init():
    particles = Particles(3)
    body0 = particles[0]
    body0.mass = 1.0 | nbody_system.mass
    body0.position = (-1., 0., 0.) | nbody_system.length
    body0.velocity = (v1, v2, 0.) | nbody_system.speed
    
    body1 = particles[1]
    body1.mass = 1.0 | nbody_system.mass
    body1.position = (1., 0., 0.) | nbody_system.length
    body1.velocity = (v1, v2, 0.) | nbody_system.speed
    
    body2 = particles[2]
    body2.mass = m3 | nbody_system.mass
    body2.position = (0., 0., 0.) | nbody_system.length
    body2.velocity = (-2/m3 * v1, -2/m3 * v2, 0.) | nbody_system.speed
    
    particles.move_to_center()
    return particles
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def integrate_solar_system(particles, end_time):
   
    gravity = Hermite(number_of_workers=2)
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

    x_body0.append(body0.x)
    y_body0.append(body0.y)
    x_body1.append(body1.x)
    y_body1.append(body1.y)
    x_body2.append(body2.x)
    y_body2.append(body2.y)    


    while gravity.model_time < end_time:
        
        gravity.evolve_model(gravity.model_time + (t_step | nbody_system.time))
        x_body0.append(body0.x)
        y_body0.append(body0.y)
        x_body1.append(body1.x)
        y_body1.append(body1.y)
        x_body2.append(body2.x)
        y_body2.append(body2.y)

    gravity.stop()
    return x_body0, y_body0, x_body1, y_body1, x_body2, y_body2
###BOOKLISTSTOP2###

###BOOKLISTSTART3###
def plot_track(x0,y0,x1,y1,x2,y2, output_filename):
    import matplotlib as mpl
    from matplotlib import pyplot
    mpl.rc('xtick', direction='in', top=True)
    mpl.rc('xtick.major', size=6, width=2.5)
    mpl.rc('xtick.minor', size=3.5, width=2)
    mpl.rc('ytick', direction='in', right=True)
    mpl.rc('ytick.major', size=6, width=2.5)
    mpl.rc('ytick.minor', size=3.5, width=2)
    mpl.rc('axes', linewidth=3)

    figure = pyplot.figure(figsize=(10, 10))
    pyplot.rcParams.update({'font.size': 25})
    plot = figure.add_subplot(1,1,1)
    ax = pyplot.gca()
    ax.minorticks_on()
    ax.tick_params(labelsize=22, pad=6)
    ax.locator_params(nbins=9)

    x_label = 'x'
    y_label = 'y'
    pyplot.xlabel(x_label)
    pyplot.ylabel(y_label)

    plot.plot(x0.value_in(nbody_system.length), y0.value_in(nbody_system.length), color = 'b', label='body1')
    plot.plot(x1.value_in(nbody_system.length), y1.value_in(nbody_system.length), color = 'r', label='body2')
    plot.plot(x2.value_in(nbody_system.length), y2.value_in(nbody_system.length), color = 'k', label='body3')
        
    plot.set_xlim(-1.5, 1.5)
    plot.set_ylim(-yr, yr)
    lg = plot.legend(fontsize=15, title='step=%.1f, n_T=%d'%(t_step,n_T))
    lg.get_title().set_weight('bold')
    lg.get_title().set_fontsize(15)
    plot.set_title(pattern[i], fontsize=25, weight='bold', pad=10)

    save_file = '%s_%.1f_%d.png'%(pattern[i], t_step, n_T)
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

    for i in [3,5]:
	v1 = v1_list[i]
	v2 = v2_list[i]
	m3 = m3_list[i]
	T = T_list[i]
	yr = yrange[i]
        
	particles = body_init()
	x0,y0, x1,y1, x2,y2 = integrate_solar_system(particles, n_T * T)
	plot_track(x0, y0, x1, y1, x2, y2, o.output_filename)
    

