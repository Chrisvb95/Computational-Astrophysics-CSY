import sys
import numpy as np
import matplotlib as mpl
from time import time
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
   
    print 'Initializing the cluster...'
    converter = nbody_system.nbody_to_si(mCluster,rCluster)
    cluster = new_plummer_model(N,convert_nbody=converter)
    mZAMS = new_powerlaw_mass_distribution(N,mMin,mMax,alpha)
    cluster.mass = mZAMS
    cluster.scale_to_standard(converter)
    return cluster


def evolve_cluster(cluster,mCut,dt,tend,mCluster,rCluster):

    print 'Start evolving the cluster with mCut =', mCut
    # timing
    t1 = time()

    converter=nbody_system.nbody_to_si(mCluster,rCluster)

    # split particle according to mass_cut
    low_mass_stars = cluster.select(lambda m: m < mCut,["mass"])
    high_mass_stars = cluster.select(lambda m: m >= mCut,["mass"])
    print 'Number of low-mass stars:', low_mass_stars.__len__()
    print 'Number of high-mass stars:', high_mass_stars.__len__()

    
    # make models and assigning particles    
    code_tree = BHTree(converter)
    code_direct = ph4(converter)
    code_tree.particles.add_particles(low_mass_stars)
    code_direct.particles.add_particles(high_mass_stars)
    channel_from_code_tree = code_tree.particles.new_channel_to(low_mass_stars)
    channel_from_code_direct = code_direct.particles.new_channel_to(high_mass_stars)    
    
    # make bridge
    combined_gravity = bridge()
    combined_gravity.add_system(code_tree,(code_direct,))
    combined_gravity.add_system(code_direct,(code_tree,))
    combined_gravity.timestep = dt

    # the initial total energy
    E0 = combined_gravity.potential_energy + combined_gravity.kinetic_energy
                     
    # evolve the model
    times = quantities.arange(0.|units.Myr, tend+dt, dt)
    n_step = len(times)
    dE = np.empty(n_step)
    Qs = np.empty(n_step)

    for i,t in enumerate(times):
        if t.value_in(units.Myr)%1 == 0:
            print "Time =", t.in_(units.Myr)
        channel_from_code_tree.copy()
        channel_from_code_direct.copy()

        combined_gravity.evolve_model(t, timestep=dt)

        U = cluster.potential_energy()
        T = cluster.kinetic_energy()
        Et = U + T
        dE[i] = np.abs((E0-Et)/E0)
        Qs[i] = -T / U

    code_tree.stop()
    code_direct.stop()

    t2 = time()
    print 'Finished. Time cost: %.2f s'%(t2-t1)
    output_filename = 'ee_q_%.1f.txt'%mCut.value_in(units.MSun)
    with open(output_filename, 'w'):
        np.savetxt(output_filename, np.append(dE,Qs).reshape(2,n_step).transpose())
    print 'Saved the data of energy error and virial ratios in', output_filename

    # plot energy
    plot_energy_error(times.value_in(units.Myr), dE, Qs, mCut)
    return dE[-1], Qs[-1], t2-t1


def plot_energy_error(times, dE, Qs, mCut):
    fig_ee = plt.figure(figsize=(8,6))
    ax = fig_ee.add_subplot(111)
    ax.plot(times, dE, '-o', label='dE/E0')
    ax.plot(times, Qs, '-o', label='Q = -T/U')
    ax.set_xlabel('Time [Myr]')
    ax.set_ylabel('Relative Energy Error')
    ax.legend()
    #plt.show()
    fig_ee.savefig('energy_error_%.1f.jpg'%mCut.value_in(units.MSun))
   

if __name__ == "__main__":

    
    # determine the number of particles for which the code must be run
    if len(sys.argv) > 1:
        N = int(sys.argv[1])
    else:
        print('N-particles not given, setting N = 10000...')
        N = 10000

    # set parameters for cluster
    mMin = 0.1|units.MSun
    mMax = 100|units.MSun
    alpha = -2.35
    mCluster = N|units.MSun
    rCluster = 3|units.parsec    

    # Setting parameters for evolving the model
    mCut_list = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]) |units.MSun
    #mCut_low = np.array([0.3, 0.4, 0.5])
    #mCut_mid = np.array([0.6, 0.7, 0.8, 0.9, 1.0, 2.0]) |units.MSun
    #mCut_high = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]) |units.MSun
    #mCut_list = [10|units.MSun]
    dt = 0.1 |units.Myr
    tend = 10 |units.Myr

    # Initialize empty arrays for Energy Error and Virial Ratio
    n_mCut = len(mCut_list)
    EE = np.empty(n_mCut)
    Q = np.empty(n_mCut)
    wc_time = np.empty(n_mCut)

    # set a random seed
    np.random.seed(1)
    cluster0 = cluster_init(N,mMin,mMax,alpha,mCluster,rCluster)
    
    for i, mCut in enumerate(mCut_list):
        # set a random seed
        np.random.seed(1)
        cluster = cluster0.copy()
        EE[i], Q[i], wc_time[i] = evolve_cluster(cluster,mCut,dt,tend,mCluster,rCluster)

    output_filename = 'ee_q_t.txt'
    with open(output_filename, 'w'):
        np.savetxt(output_filename, np.append(np.append(EE,Q),wc_time).reshape(3, n_mCut).transpose())
    print 'Saved the data of energy error, virial ratios, and wall-clock time in', output_filename
            


    


    
