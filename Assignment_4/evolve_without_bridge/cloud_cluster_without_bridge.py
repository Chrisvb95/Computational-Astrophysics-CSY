'''
Task 4: cloud + cluster
'''
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
mpl.rc('xtick', direction='in', top=True, labelsize=14)
mpl.rc('xtick.major', size=5, width=2)
mpl.rc('xtick.minor', size=3, width=2)
mpl.rc('ytick', direction='in', right=True, labelsize=14)
mpl.rc('ytick.major', size=5, width=2)
mpl.rc('ytick.minor', size=3, width=2)
mpl.rc('axes', linewidth=2, labelsize=18, labelpad=10,\
       titlesize=20, titlepad=10)
mpl.rc('image', origin='lower')

from amuse.lab import * 
from amuse.units import units, nbody_system, quantities
from amuse.couple import bridge
from amuse.community.fi.interface import Fi
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube  #what's this?


def make_map(sph, L=200):

    N = L
   
    x, y = np.indices((N + 1, N + 1))

    x = L * (x.flatten() - N / 2.) / N
    y = L * (y.flatten() - N / 2.) / N
    z = x * 0.
    vx = 0. * x
    vy = 0. * x
    vz = 0. * x

    x = units.parsec(x)
    y = units.parsec(y)
    z = units.parsec(z)
    vx = units.kms(vx)
    vy = units.kms(vy)
    vz = units.kms(vz)

    rho, rhovx, rhovy, rhovz, rhoe = sph.get_hydro_state_at_point(
        x, y, z, vx, vy, vz)
    rho = rho.reshape((N+1,N+1))

    return rho

def plot_cloud_cluster(cluster, sph, title, vrange=[0,5], L=200):

    rho = make_map(sph, L)
    #rhomin = int(np.log10(rho.value_in(units.amu/units.cm**3)).min())
    #rhomax = np.ceil(np.log10(rho.value_in(units.amu/units.cm**3)).max())

    fig,ax = plt.subplots(1,figsize=(8,8))
    vmin, vmax = vrange   
    cax = ax.imshow(np.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)).transpose(),
                    extent=[-L/2, L/2, -L/2, L/2], vmin=vmin, vmax=vmax, cmap='jet')
    
    ax.scatter(cluster.x.value_in(units.parsec), cluster.y.value_in(units.parsec), s=5, c='C3')   

    #ticks = np.arange(rhomin, rhomax+1, 1, dtype=int)
    ticks = np.arange(vmin, vmax+1, 1, dtype=int)
    ticklabels = [str(i) for i in ticks]
    cbar = fig.colorbar(cax, ticks=ticks, orientation='vertical', fraction=0.045)
    cbar.ax.set_yticklabels(ticklabels)  # horizontal colorbar
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('log projected density [$amu/cm^3$]', fontsize=15, rotation=270, labelpad=25)

    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.set_title(title+' Myr', weight='bold')
    fig.savefig(title+'.png', dpi=200)

    #plt.close()

    #plt.show()

    return None
    


def cluster_init(N, Mtot, Rvir, init_v, init_position, alpha=-2.35, seed=0):
    '''Initialises a star cluster that follows a powerlaw mass distribution.
       `N` is the number of star members, `Mtot` is the total mass of the cluster,
       and `Rvir` is the virial radius of the model.'''
    converter = nbody_system.nbody_to_si(Mtot,Rvir)

    np.random.seed(seed)
    cluster = new_plummer_model(N,convert_nbody=converter)
    m_average = Mtot/N

    # new_powerlaw_mass_distribution(N, Mmin, Mmax, alpha), alpha is the powerlaw index
    mZAMS = new_powerlaw_mass_distribution(N,m_average*0.3, m_average*80, alpha)
    cluster.mass = mZAMS
    cluster.scale_to_standard(converter)
    cluster.velocity += init_v
    cluster.position += init_position
    return cluster, converter


def cloud_init(N=100, M=100.|units.MSun, R=1.|units.parsec, seed=0):
    '''Initialises the cloud'''
    converter = nbody_system.nbody_to_si(M,R)
    cloud = molecular_cloud(targetN=N, convert_nbody=converter,\
                            base_grid=body_centered_grid_unit_cube,\
                            seed=seed).result
    return cloud, converter


def evolve(cluster,cloud, converter_grav,converter_hydro, t_end,dt_bridge):

    # Initialising the direct N-body integrator
    print 'Setting up the gravity code...'
    gravity = ph4(converter_grav)
    gravity.particles.add_particles(cluster)

    # Initialising the hydro code
    print 'Setting up the hydrodynamics code...'
    sph = Fi(converter_hydro, mode="openmp")
    sph.gas_particles.add_particles(cloud)

    sph.parameters.use_hydro_flag = True
    sph.parameters.radiation_flag = False
    sph.parameters.gamma = 1
    sph.parameters.isothermal_flag = True
    sph.parameters.integrate_entropy_flag = False
    sph.parameters.timestep = dt  
    sph.parameters.verbosity = 0
    sph.parameters.eps_is_h_flag = False    # h_smooth is constant
    eps = 0.1 | units.parsec
    sph.parameters.gas_epsilon = eps
    sph.parameters.sph_h_const = eps
    cloud.h_smooth= eps


    # Build a bridge between hydro and grav
    print 'Bridging...'
    grav_hydro = bridge.Bridge(use_threading=False)
    grav_hydro.add_system(gravity, (sph,))
    grav_hydro.add_system(sph, (gravity,))
    grav_hydro.timestep = dt_bridge

    # Set up channels
    channel_from_grav_to_cluster = gravity.particles.new_channel_to(cluster)
    channel_from_sph_to_cloud = sph.gas_particles.new_channel_to(cloud)
    
    # Start evolution
    print 'Start evolving...'
    model_time = 0.|units.yr
    while model_time <= t_end:

        # Evolve for one step        
        model_time += dt
        grav_hydro.evolve_model(model_time)
        channel_from_grav_to_cluster.copy()
        channel_from_sph_to_cloud.copy()

    gravity.stop()
    sph.stop()

    return cluster, cloud 


if __name__ == "__main__":

    # Set initial parameters
    print 'Initializing the star cluster and the molecular cloud...'
    # Molecular cloud (typically, mass = 1e3-1e7 MSun, diameter = 15-600 pc ?)
    N_cloud = 2000
    M_cloud = 1e4 | units.MSun
    R_cloud = 10 | units.parsec

    # Cluster (typically, number of stars = 1e5-1e6, diameter = 3-100 pc ?)
    N_cluster = 1000
    Mtot_cluster = N_cluster | units.MSun
    Rvir_cluster = 5. | units.parsec
    # iniial velocity and position of the cluster's COM
    v_cluster = (0,0,0) | units.km/units.s
    p_cluster = (-(R_cloud*2+Rvir_cluster*5).value_in(units.parsec),0,0)|units.parsec

    # Initialising the two systems
    seed = 50
    cluster, converter_cluster = cluster_init(N_cluster, Mtot_cluster, Rvir_cluster,\
                                              v_cluster, p_cluster, seed=seed)
    cloud, converter_cloud = cloud_init(N_cloud, M_cloud, R_cloud, seed=seed)

    t_end = 10. | units.Myr #??
    dt_bridge = 0.1 | units.Myr
    dt_grav = 0.1 | units.Myr
    dt_sph = 0.01 | units.Myr

    #evolve(cluster,cloud, converter_cluster,converter_cloud, t_end,dt_bridge)
    # Initialising the hydro code
    #converter = nbody_system.nbody_to_si(1|units.MSun,1|units.parsec)
    #hydro = Fi(converter, mode="openmp")
    #hydro.gas_particles.add_particles(cloud)
    #hydro.stop()

    #print('CLUSTER')
    #print(cluster)
    #print('CLOUD')
    #print(cloud)

    print 'Setting up the gravity code and the hydro code...'
    gravity = ph4(converter_cluster)
    gravity.particles.add_particles(cluster)

    sph = Fi(converter_cloud, mode="openmp")
    sph.gas_particles.add_particles(cloud)
    sph.parameters.use_hydro_flag = True
    sph.parameters.radiation_flag = False
    sph.parameters.gamma = 1
    sph.parameters.isothermal_flag = True
    sph.parameters.integrate_entropy_flag = False
    sph.parameters.timestep = dt_sph
    sph.parameters.verbosity = 0
    sph.parameters.eps_is_h_flag = False    # h_smooth is constant
    eps = 0.1 | units.parsec
    sph.parameters.gas_epsilon = eps
    sph.parameters.sph_h_const = eps
    cloud.h_smooth= eps

    channel_from_grav_to_cluster = gravity.particles.new_channel_to(cluster)
    channel_from_sph_to_cloud = sph.gas_particles.new_channel_to(cloud)


    print 'Start evolving the molecular cloud...'
    model_time = 0.|units.yr
    while model_time <= t_end:

        # Evolve for one step        
        print model_time.in_(units.Myr)
        gravity.evolve_model(model_time, dt=dt_grav)
        sph.evolve_model(model_time, dt=dt_sph)
        channel_from_grav_to_cluster.copy()
        channel_from_sph_to_cloud.copy()
        plot_cloud_cluster(cluster, sph, title='%.1f'%model_time.value_in(units.Myr),\
                           L=160, vrange=[0,5])
        model_time += dt_bridge

    gravity.stop()
    sph.stop()


    # Evolve model
    #evolve(cluster,cloud)



    
    
    
    
    
    
    
