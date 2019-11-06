import numpy as np
import sys
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from amuse.units import units
from amuse.units import quantities
from amuse.units import constants
from amuse.units import nbody_system
from amuse.ext.bridge import bridge
from amuse.community.ph4.interface import ph4
from amuse.community.bhtree.interface import BHTree
from amuse.lab import new_plummer_model
from amuse.lab import new_powerlaw_mass_distribution


def cluster_init(N,mMin,mMax,alpha,mCluster,rCluster):

    converter = nbody_system.nbody_to_si(mCluster,rCluster)
    cluster = new_plummer_model(N,convert_nbody=converter)
    mZAMS = new_powerlaw_mass_distribution(N,mMin,mMax,alpha)
    cluster.mass = mZAMS
    cluster.scale_to_standard(converter)
    return cluster


def evolve_cluster(cluster,mCut,dt,tend,mCluster,rCluster):
    
    converter=nbody_system.nbody_to_si(mCluster,rCluster)
    # Splitting particle according to mass_cut
    low_mass_stars = cluster.select(lambda m: m < mCut,["mass"])
    high_mass_stars = cluster.select(lambda m: m >= mCut,["mass"])

    # Sanity checks
    print('dynamical timescale:', cluster.dynamical_timescale().value_in(units.Myr))
    print('Number of low-mass stars:',low_mass_stars.__len__())
    print('Number of high-mass stars:',high_mass_stars.__len__())
    #plot_cluster(low_mass_stars, high_mass_stars)
    
    # Making models and assigning particles    
    code_tree = BHTree(converter)
    code_direct = ph4(converter)
    code_tree.particles.add_particles(low_mass_stars)
    code_direct.particles.add_particles(high_mass_stars)    
    channel_from_code_tree = code_tree.particles.new_channel_to(low_mass_stars)
    channel_from_code_direct = code_direct.particles.new_channel_to(high_mass_stars)    
    
    # Making bridge
    combined_gravity = bridge()
    combined_gravity.add_system(code_tree,(code_direct,))
    combined_gravity.add_system(code_direct,(code_tree,))
    combined_gravity.timestep = dt

    # Evolving the model
    time_series = []
    half_mass_radii = []
    core_radii = [] 
    times = quantities.arange(0.|units.Myr, tend, dt)
    for i,t in enumerate(times):
        print "Time =", t.in_(units.Myr)
        channel_from_code_tree.copy()
        channel_from_code_direct.copy()
        time_series.append(t.value_in(units.Myr))
        pos,coreradius,coredens=cluster.densitycentre_coreradius_coredens(converter)
        lr,mf=cluster.LagrangianRadii(converter) # outputs are radii, mass fractions
        half_mass_radii.append(lr[5].value_in(units.parsec)) # 5th argument attributes to half-mass radius
        core_radii.append(coreradius.value_in(units.parsec))
        combined_gravity.evolve_model(t, timestep=dt)

    code_tree.stop()
    code_direct.stop()    
    
    # Plotting results
    #plot_cluster(low_mass_stars, high_mass_stars)
    
    # Plotting radii
    plot_radii(time_series, half_mass_radii, x_label='time(Myr)', y_label='half-mass radius (pc)')
    plot_radii(time_series, core_radii, x_label='time(Myr)', y_label='core radius (pc)')

def plot_radii(time, radii, x_label, y_label):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel=x_label,  ylabel=y_label)
    ax.plot(time, radii) 
    plt.show()



def plot_cluster(low_mass, high_mass):
    
    m_low = low_mass.mass.value_in(units.MSun)
    x_low = low_mass.x.value_in(units.parsec)
    y_low = low_mass.y.value_in(units.parsec)
    z_low = low_mass.z.value_in(units.parsec)
    
    m_high = high_mass.mass.value_in(units.MSun)
    x_high = high_mass.x.value_in(units.parsec)
    y_high = high_mass.y.value_in(units.parsec)
    z_high = high_mass.z.value_in(units.parsec)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x_low,y_low, c='b',zorder=0)
    ax.scatter(x_high,y_high, c='r',zorder=1)    
    plt.show()
    
    masses = np.append(m_high,m_low)
    plt.hist(masses,range=[0.1,100],bins=100)
    plt.yscale('log')    
    plt.show()


if __name__ == "__main__":
    
    # Determining number of particles for which the code must be run
    if len(sys.argv) > 1:
        N = int(sys.argv[1])
    else:
        print('N-particles not given, setting N = 10000...')
        N = 10000

    # Setting parameters for cluster
    mMin = 0.1|units.MSun
    mMax = 100|units.MSun
    alpha = -2.35
    mCluster = N|units.MSun
    rCluster = 3|units.parsec    

    # Setting parameters for evolving the model
    mCut = 1|units.MSun
    dt = 0.1 |units.Myr
    tend = 10 |units.Myr

    cluster = cluster_init(N,mMin,mMax,alpha,mCluster,rCluster)
    evolve_cluster(cluster,mCut,dt,tend,mCluster,rCluster)    

    #plot_cluster(cluster)

    
