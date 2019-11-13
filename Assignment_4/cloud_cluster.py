'''
Task 4: cloud + cluster
'''
import numpy as np
from amuse.lab import * 
from amuse.units import units, nbody_system, quantities
from amuse.couple import bridge
from amuse.community.fi.interface import Fi
from matplotlib import pyplot as plt

def plot_cloud_cluster(cloud,cluster,N=200):

    #print('Plotting', title)
    '''
    L = 400
   
    x, y = np.indices((N + 1, N + 1))

    x = L * (x.flatten() - N / 2.) / N
    y = L * (y.flatten() - N / 2.) / N
    z = x * 0.
    vx = 0. * x
    vy = 0. * x
    vz = 0. * x

    x = units.AU(x)
    y = units.AU(y)
    z = units.AU(z)
    vx = units.kms(vx)
    vy = units.kms(vy)
    vz = units.kms(vz)

    rho, rhovx, rhovy, rhovz, rhoe = cloud.get_hydro_state_at_point(
        x, y, z, vx, vy, vz)
    rho = rho.reshape((N + 1, N + 1))

    fig,ax = plt.subplots(1,figsize=(8, 8))    
    ax.imshow(np.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)).transpose(),
                  extent=[-L / 2, L / 2, -L / 2, L / 2], vmin=10, vmax=15, origin='lower')
    '''
    fig,ax = plt.subplots(1,figsize=(8, 8))   
    ax.scatter(cloud.x.value_in(units.parsec),cloud.y.value_in(units.parsec))
    ax.scatter(cluster.x.value_in(units.parsec),cluster.y.value_in(units.parsec))    
    #plt.title(title+' yr', fontsize=15, weight='bold', pad=10)
    plt.xlabel('AU')
    #plt.savefig('distribution_plot_passing_star_full_hydro_4k/'+title+'.png', dpi=200)
    plt.show()
    


def cluster_init(N,m,hmr,v,init_position,alpha=-2.35):
    # Initialises a star cluster that follows a powerlaw mass distrib.
    converter = nbody_system.nbody_to_si(m,hmr)
    cluster = new_plummer_model(N,convert_nbody=converter)
    mZAMS = new_powerlaw_mass_distribution(N,m*0.3,m*80,alpha)
    cluster.mass = mZAMS
    cluster.scale_to_standard(converter)
    cluster.vy += v
    cluster.position += init_position
    return cluster

def cloud_init(N,m,hmr):
    # Initialises the cloud
    converter = nbody_system.nbody_to_si(m,hmr)
    cloud = new_plummer_model(N,convert_nbody=converter)
    return cloud

def evolve(cluster,cloud,dt_bridge):
    converter = nbody_system.nbody_to_si(1|units.MSun,1|units.parsec)

    # Initialising the direct N-body integrator
    gravity = ph4(converter)
    gravity.particles.add_particles(cluster)

    # Initialising the hydro code
    hydro = Fi(converter, mode="openmp")
    hydro.gas_particles.add_particles(cloud)

    # Build a bridge between hydro and grav
    grav_hydro = bridge.Bridge(use_threading=False)
    grav_hydro.add_system(gravity, (hydro,))
    grav_hydro.add_system(hydro, (gravity,))
    grav_hydro.timestep = dt_bridge

    # Set up channels
    channel_from_grav = gravity.particles.new_channel_to(cluster)
    channel_from_hydro = hydro.gas_particles.new_channel_to(cloud)
    
    # Start evolution
    while model <= t_end:

        # Evolve for one step        
        model_time += dt
        grav_hydro.evolve_model(model_time)
        channel_from_grav.copy()
        channel_from_hydro.copy()

    gravity.stop()
    hydro.stop()

    return cluster, cloud 


if __name__ == "__main__":

    # Set initial parameters
    # Cluster
    N_cluster = 1000
    m_cluster = N_cluster | units.MSun
    hmr_cluster = 5 | units.parsec   
    v_cluster = 100 | units.km/units.s
    init_position_cluster = (100,-100,0)|units.parsec
    # Cloud    
    N_cloud = 2000
    m_cloud = 100 | units.MSun
    hmr_cloud = 20 | units.parsec

    # Initialising the two systems
    cluster = cluster_init(N_cluster,m_cluster,hmr_cluster,v_cluster,init_position_cluster)
    cloud = cloud_init(N_cloud,m_cloud,hmr_cloud)

    # Initialising the hydro code
    #converter = nbody_system.nbody_to_si(1|units.MSun,1|units.parsec)
    #hydro = Fi(converter, mode="openmp")
    #hydro.gas_particles.add_particles(cloud)
    #hydro.stop()
    #print('CLUSTER')
    #print(cluster)
    #print('CLOUD')
    #print(cloud)

    plot_cloud_cluster(cloud,cluster)

    # Evolve model
    #evolve(cluster,cloud)



    
    
    
    
    
    
    
