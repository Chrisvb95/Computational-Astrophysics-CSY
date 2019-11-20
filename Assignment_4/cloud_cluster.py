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

def plot_cloud_cluster(cluster, sph, title, vrange=[-5,2], L=400):

    rho = make_map(sph, L)
    #rhomin = int(np.log10(rho.value_in(units.amu/units.cm**3)).min())
    #rhomax = np.ceil(np.log10(rho.value_in(units.amu/units.cm**3)).max())

    fig,ax = plt.subplots(1,figsize=(8,8))
    #vmin, vmax = vrange   
    #cax = ax.imshow(np.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)).transpose(),
    #                extent=[-L/2, L/2, -L/2, L/2], vmin=vmin, vmax=vmax, cmap='jet')
    
    #cax = ax.imshow((rho.value_in(units.amu / units.cm**3)).transpose(),
    #                extent=[-L/2, L/2, -L/2, L/2], vmin=vmin, vmax=vmax, cmap='jet')

    ax.scatter(cluster.x.value_in(units.parsec), cluster.y.value_in(units.parsec), s=5, c='C3')   
    ax.scatter(cloud.x.value_int(units.parsec), cluster.y.value_in(units.parsec))
    
    #ticks = np.arange(rhomin, rhomax+1, 1, dtype=int)
    #ticks = np.arange(vmin, vmax+1, 1, dtype=int)
    #ticklabels = [str(i) for i in ticks]
    #cbar = fig.colorbar(cax, ticks=ticks, orientation='vertical', fraction=0.045)
    #cbar.ax.set_yticklabels(ticklabels)  # horizontal colorbar
    #cbar.ax.tick_params(labelsize=13)
    #cbar.set_label('log projected density [$amu/cm^3$]', fontsize=15, rotation=270, labelpad=25)

    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.set_title(title+' Myr', weight='bold')
    plt.savefig('cloud_cluster/'+title+'.png', dpi=200)

    plt.close()

    return None
    

def return_L9_radius(system,M,R): 
    converter=nbody_system.nbody_to_si(M, R)     
    lr,mr = system.LagrangianRadii(converter)
    return lr[7].value_in(units.parsec)

def cluster_init(N, Mtot, Rvir, init_v, init_position, alpha=-2.35):
    '''Initialises a star cluster that follows a powerlaw mass distribution.
       `N` is the number of star members, `Mtot` is the total mass of the cluster,
       and `Rvir` is the virial radius of the model.'''
    converter = nbody_system.nbody_to_si(Mtot,Rvir)

    cluster = new_plummer_model(N,convert_nbody=converter)
    m_average = Mtot/N

    # new_powerlaw_mass_distribution(N, Mmin, Mmax, alpha), alpha is the powerlaw index
    mZAMS = new_powerlaw_mass_distribution(N,m_average*0.3, m_average*80, alpha)
    cluster.mass = mZAMS
    cluster.scale_to_standard(converter)
    cluster.velocity += init_v
    cluster.position += init_position
    return cluster, converter


def cloud_init(N=100, M=100.|units.MSun, R=1.|units.parsec):
    '''Initialises the cloud'''
    converter = nbody_system.nbody_to_si(M,R)
    cloud = molecular_cloud(targetN=N, convert_nbody=converter,\
                            base_grid=body_centered_grid_unit_cube,\
                            ).result
    #cloud = new_plummer_model(N, convert_nbody=converter)
    return cloud, converter


def evolve(cluster,cloud, converter_grav,converter_hydro, t_end, dt_bridge, dt_diag):

    converter = nbody_system.nbody_to_si(1. | units.MSun,1. | units.parsec)
    # Initialising the direct N-body integrator
    print 'Setting up the gravity code and the hydro code...'
    gravity = ph4(converter)
    gravity.particles.add_particles(cluster)

    # Initialising the hydro code
    sph = Fi(converter, mode="openmp")
    sph.gas_particles.add_particles(cloud)
    sph.parameters.use_hydro_flag = True
    sph.parameters.radiation_flag = False
    sph.parameters.gamma = 1
    sph.parameters.isothermal_flag = True
    sph.parameters.integrate_entropy_flag = False
    sph.parameters.timestep = dt_diag
    sph.parameters.verbosity = 0
    sph.parameters.eps_is_h_flag = False    # h_smooth is constant
    eps = 0.1 | units.parsec
    sph.parameters.gas_epsilon = eps
    sph.parameters.sph_h_const = eps
    cloud.h_smooth= eps

    print "cloud:", sph.gas_particles
    #from matplolib import pyplot
    plt.scatter(np.log10(sph.gas_particles.density.value_in(units.g/units.cm**3)), sph.gas_particles.pressure.value_in(units.kg/units.m/units.s**2), s=10)
    plt.show()
    # Building a bridge between hydro and grav
    print 'Bridging...'
    grav_sph = bridge.Bridge(use_threading=False)
    grav_sph.add_system(gravity, (sph,))
    grav_sph.add_system(sph, (gravity,))
    grav_sph.timestep = dt_bridge

    # Setting up channels from code to cloud and cluster
    channel_from_grav_to_cluster = gravity.particles.new_channel_to(cluster)
    channel_from_sph_to_cloud = sph.gas_particles.new_channel_to(cloud)

    # Open file 
    f = open('l9_radius.csv','wb')

    # Initializing 90 percent lagrangian radius 
    lr9 = []
    lr9.append(return_L9_radius(cloud, 1|units.MSun, 1|units.AU))
    f.write(str(lr9[-1])+"\n")

    # Evolving
    print 'Start evolving the molecular cloud...'
    times = quantities.arange(0.|units.Myr, t_end, dt_diag)
    for i,t in enumerate(times):
      
        print t.in_(units.Myr)
        grav_sph.evolve_model(t, dt_diag)        
        channel_from_grav_to_cluster.copy()
        channel_from_sph_to_cloud.copy()
        #plot_cloud_cluster(cluster, sph, title='{0}'.format(float(t.value_in(units.Myr))),\
        #                   L=400, vrange=[-5,1])
        lr9.append(return_L9_radius(cloud, 1|units.MSun, 1|units.AU))
        f.write(str(lr9[-1])+"\n")        
        print('Lagrangian radius:', lr9[-1])
        print "cloud:", sph.gas_particles        
        plt.scatter(np.log10(sph.gas_particles.density.value_in(units.g/units.cm**3)), sph.gas_particles.pressure.value_in(units.kg/units.m/units.s**2), s=10)
        plt.show()
    gravity.stop()
    sph.stop()

    return cluster, cloud 


if __name__ == "__main__":

    # Set initial parameters
    print 'Initializing the star cluster and the molecular cloud...'
    # Molecular cloud (typically, mass = 1e3-1e7 MSun, diameter = 5-200 pc ?)
    N_cloud = 4000
    M_cloud = 1e4 | units.MSun
    R_cloud = 100 | units.parsec

    # Cluster (typically, number of stars = 1e5-1e6, diameter = 3-100 pc ?)
    N_cluster = 1000
    Mtot_cluster = N_cluster | units.MSun
    Rvir_cluster = 5. | units.parsec
    # iniial velocity and position of the cluster's COM
    v_cluster = (100,0,0) | units.km/units.s
    #p_cluster = (-(R_cloud*2+Rvir_cluster*5).value_in(units.parsec),0,0)|units.parsec
    p_cluster = (-200, -100, 0) | units.parsec

    # Setting a seed
    np.random.seed(1)

    # Initialising the two systems
    cluster, converter_cluster = cluster_init(N_cluster, Mtot_cluster, Rvir_cluster,\
                                              v_cluster, p_cluster)
    cloud, converter_cloud = cloud_init(N_cloud, M_cloud, R_cloud)

    t_end = 10 | units.Myr
    dt_bridge = 0.1 | units.Myr
    dt_diag = 0.1 | units.Myr

    evolve(cluster,cloud, converter_cluster,converter_cloud, t_end,dt_bridge, dt_diag)


    
    
    
    
    
    
    
