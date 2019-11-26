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
from amuse.ic.gasplummer import new_plummer_gas_model



def make_map(sph, L=200):

    N = 500

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

def plot_cloud_cluster(cluster, sph, title, vrange=[-6,2]):

    gas_position = (sph.gas_particles.position).value_in(units.parsec)
    L = 2 * (np.ceil(np.max(np.abs(gas_position))/100.) * 100.)
    rho = make_map(sph, L)
 
    rhomin = int(np.log10(1.e-5+rho.value_in(units.amu/units.cm**3)).min())
    rhomax = np.ceil(np.log10(1.e-5+rho.value_in(units.amu/units.cm**3)).max())

    fig,ax = plt.subplots(1,figsize=(8,8))
    vmin, vmax = vrange   
    cax = ax.imshow(np.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)).transpose(),
                    extent=[-L/2, L/2, -L/2, L/2], vmin=vmin, vmax=vmax, cmap='jet')
    
    #cax = ax.imshow((rho.value_in(units.amu / units.cm**3)).transpose(),
    #                extent=[-L/2, L/2, -L/2, L/2], vmin=vmin, vmax=vmax, cmap='jet')
   
    ax.scatter(cloud.x.value_in(units.parsec), cloud.y.value_in(units.parsec), s=1, c='magenta',\
               label='cloud', alpha=0.3)
    ax.scatter(cluster.x.value_in(units.parsec), cluster.y.value_in(units.parsec), s=2, c='C0',\
               label='cluster')
    if len(sph.dm_particles) > 0:
        ax.scatter(sph.dm_particles.x.value_in(units.parsec), sph.dm_particles.y.value_in(units.parsec), s=5,\
                   marker='*', c='k', label='formed star', alpha=0.6)
    
    ticks = np.arange(rhomin, rhomax+1, 1, dtype=int)
    ticks = np.arange(vmin, vmax+1, 1, dtype=int)
    ticklabels = [str(i) for i in ticks]
    cbar = fig.colorbar(cax, ticks=ticks, orientation='vertical', fraction=0.045)
    cbar.ax.set_yticklabels(ticklabels)  # horizontal colorbar
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('log projected density [$amu/cm^3$]', fontsize=15, rotation=270, labelpad=25)

    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.set_title(title+' Myr', weight='bold')
    ax.legend(fontsize=14)
    plt.savefig(title+'.png', dpi=200)

    plt.close()

    return None
    

def L9_radius(model, converter):      
    lr,mr = model.LagrangianRadii(converter)
    return lr[7].value_in(units.parsec)


def Jeans_density(M, m=2*1.66e-24|units.g, T=15|units.K):
    rho_J = 81*constants.kB**3 / (32*np.pi*constants.G**3) * T**3 / (m**3 * M**2)
    return rho_J.in_(units.kg/units.m**3)


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


def cloud_init(Ngas, Mgas, Rgas):
    '''Initialises the cloud'''
    converter = nbody_system.nbody_to_si(Mgas,Rgas)
    #cloud = molecular_cloud(targetN=Ngas, convert_nbody=converter,\
    #                         base_grid=body_centered_grid_unit_cube,\
    #                        ).result
    cloud = new_plummer_gas_model(
            Ngas, convert_nbody=converter)
    return cloud, converter


def evolve(cluster,cloud, converter_grav,converter_sph, t_end, dt_bridge, dt_diag,\
           sink=False):

    with open('print_out.txt', 'a') as pf:
        pf.write('Setting up the gravity code and the hydro code...\n')
    #converter = nbody_system.nbody_to_si(1. | units.MSun, 1. | units.parsec)
    # Initialising the direct N-body integrator
    gravity = ph4(converter_grav)
    gravity.particles.add_particles(cluster)

    # Initialising the hydro code
    sph = Fi(converter_sph, mode="openmp")
    sph.gas_particles.add_particles(cloud)
    #sph.parameters.use_hydro_flag = True
    sph.parameters.radiation_flag = False
    #sph.parameters.isothermal_flag = True
    #sph.parameters.integrate_entropy_flag = False
    sph.parameters.timestep = dt_bridge
    #sph.parameters.verbosity = 0
    #sph.parameters.eps_is_h_flag = False    # h_smooth is constant
    #eps = 0.1 | units.parsec
    #sph.parameters.gas_epsilon = eps
    #sph.parameters.sph_h_const = eps
    #cloud.h_smooth= eps

    #print("cloud:", sph.gas_particles)
    #plt.scatter(np.log10(sph.gas_particles.density.value_in(units.g/units.cm**3)), 
	#	sph.gas_particles.pressure.value_in(units.kg/units.m/units.s**2), s=10)
    #plt.show()

    # Building a bridge between hydro and grav
    with open('print_out.txt', 'a') as pf:
        pf.write('Bridging...\n')
    grav_sph = bridge.Bridge(use_threading=False)
    grav_sph.add_system(gravity, (sph,))
    grav_sph.add_system(sph, (gravity,))
    grav_sph.timestep = dt_bridge

    # Setting up channels from code to cloud and cluster
    channel_from_grav_to_cluster = gravity.particles.new_channel_to(cluster)
    channel_from_sph_to_cloud = sph.gas_particles.new_channel_to(cloud)

    if sink == True:
        with open('print_out.txt','a') as pf:
            pf.write('star formation is considered\n')
        stars = Particles(0)
        sph.dm_particles.add_particles(stars)
        channel_from_sph_to_star = sph.dm_particles.new_channel_to(stars)

        density_threshold = Jeans_density(M=sph.gas_particles.mass.max())
        sph.parameters.stopping_condition_maximum_density = density_threshold
        density_limit_detection = sph.stopping_conditions.density_limit_detection
        density_limit_detection.enable()

        merge_radius = 1e-2 | units.parsec  # around 20 AU

    # Initializing 90 percent lagrangian radius 
    lr9 = []
    lr9.append(L9_radius(cloud, converter_sph))

    # Evolving
    with open('print_out.txt', 'a') as pf:
        pf.write('Start evolving the molecular cloud...\n')


    times = quantities.arange(0.|units.Myr, t_end+dt_diag, dt_diag)
    for i,t in enumerate(times):
        with open('print_out.txt', 'a') as pf:
            pf.write(str(t.value_in(units.Myr))+' Myr\n')

        if sink == True:
            # make sure at this time 'stars' and 'sph.gas_particles' are the same
            resolve_sinks(sph, stars, cloud, density_threshold, t)

        grav_sph.evolve_model(t, timestep=dt_bridge)
        channel_from_grav_to_cluster.copy()
        channel_from_sph_to_cloud.copy()

        if sink == True:
            merge_stars(sph, stars, merge_radius)

        #  make plots
        plot_cloud_cluster(cluster, sph, title='{0}'.format(float(t.value_in(units.Myr))),\
                           vrange=[-5,3])

        # save data (energy will be added afterwards)
        lr9.append(L9_radius(cloud, converter_sph))
        #print("cloud:", sph.gas_particles)
	#xx  
        #plt.scatter(np.log10(sph.gas_particles.density.value_in(units.g/units.cm**3)),
	#		 sph.gas_particles.pressure.value_in(units.kg/units.m/units.s**2), s=10)
        #plt.show()

    gravity.stop()
    sph.stop()

    return lr9


def resolve_sinks(sph, stars, cloud, density_thr, model_time):
    high_dens = sph.gas_particles.select_array(lambda rho: rho > density_thr, ["rho"])
    if len(high_dens) <= 0:
        with open('print_out.txt','a') as pf:
            pf.write('(no high density particle has been identified in this step)\n')
    else:
        candidate_stars = high_dens.copy()
        sph.gas_particles.remove_particles(high_dens)
        sph.gas_particles.synchronize_to(cloud)
    
    #if len(candidate_stars) > 0:
        with open('print_out.txt','a') as pf:
            pf.write('identified %d candidates; '%len(candidate_stars))
        newstars_in_code = sph.dm_particles.add_particles(candidate_stars)
        newstars = Particles()
        for nsi in newstars_in_code:
            if nsi not in stars:
                newstars.add_particle(nsi)
        with open('print_out.txt','a') as pf:
            pf.write('confirmed %d new stars\n'%len(candidate_stars))

        newstars.birth_age = model_time
        newstars.Lx = 0 | (units.g * units.m**2)/units.s
        newstars.Ly = 0 | (units.g * units.m**2)/units.s
        newstars.Lz = 0 | (units.g * units.m**2)/units.s

        with open('print_out.txt','a') as pf:
            pf.write('previous N: %d in stars, %d in sph.dm_particles;\n'%(len(stars), len(sph.dm_particles)))
        stars.add_particles(newstars)
        with open('print_out.txt','a') as pf:
            pf.write('afterwards N: %d in stars, %d in sph.dm_particles;\n'%(len(stars), len(sph.dm_particles)))


def merge_stars(sph, stars, merge_radius):
    if len(stars) <= 0 or stars.radius.max() <= (0.|units.AU):
        with open('print_out.txt','a') as pf:
           pf.write('(no merge happens in this step)\n')
        return

    with open('print_out.txt','a') as pf:
        pf.write('identifying merging groups...\n')
    # select out star pairs in which the distance between two members is smaller than `merge_radius` 
    ccs = stars.copy().connected_components(threshold=merge_radius)
    if len(ccs) <= 0:
        with open('print_out.txt','a') as pf:
            pf.write('(no merge happens in this step)\n')
        return

    n_merge = 0
    newstars = Particles()
    for cc in ccs:
        if len(cc) > 1:
            n_merge += 1
            merge_two_stars(stars, cc.copy())
            stars.synchronize_to(sph.dm_particles)
            with open('print_out.txt','a') as pf:
                pf.write('two stars just merged\n')
    

def merge_two_stars(stars, particles_in_encounter):
    com_pos = particles_in_encounter.center_of_mass()
    com_vel = particles_in_encounter.center_of_mass_velocity()
    new_particle = Particles(1)
    mu = particles_in_encounter[0].mass / particles_in_encounter.mass.sum()

    new_particle.birth_age = particles_in_encounter.birth_age.min()
    new_particle.mass = particles_in_encounter.total_mass()
    new_particle.position = com_pos
    new_particle.velocity = com_vel
    new_particle.name = 'Star'
    new_particle.radius = particles_in_encounter.radius.max()
    stars.add_particles(new_particle)
    stars.remove_particles(particles_in_encounter)



if __name__ == '__main__':

    with open('print_out.txt', 'w') as pf:
        pf.write('START RUNNING\n')

    # Set initial parameters
    with open('print_out.txt', 'a') as pf:
        pf.write('Initializing the star cluster and the molecular cloud...\n')
    # Molecular cloud (typically, mass = 1e3-1e7 MSun, diameter = 5-200 pc ?)
    N_cloud = 1000
    Mtot_cloud = 1e4 | units.MSun
    Rvir_cloud = 10 | units.parsec

    # Cluster (typically, number of stars = 1e5-1e6, diameter = 3-100 pc ?)
    N_cluster = 100
    Mtot_cluster = 1e3 | units.MSun
    Rvir_cluster = 5. | units.parsec
    # initial velocity and position of the cluster's COM
    v_cluster = (10,10,0) | units.km/units.s
    p_cluster = (-20,-20,0) | units.parsec

    # Setting a seed
    np.random.seed(1)

    # Initialising the two systems
    # don't want to include cluster yet
    #Mtot_cluster = 1e-3|units.MSun
    #Rvir_cluster = 1.|units.parsec
    cluster, conv_grav = cluster_init(N_cluster, Mtot_cluster, Rvir_cluster,\
                                              v_cluster, p_cluster)
    cloud, conv_sph = cloud_init(N_cloud, Mtot_cloud, Rvir_cloud)

    t_end = 40. | units.Myr
    dt_bridge = 0.1 | units.Myr
    dt_diag = 1.0 | units.Myr
    

    lr9 = evolve(cluster,cloud, conv_grav,conv_sph, t_end,dt_bridge,dt_diag, sink=True)

    with open('print_out.txt','a') as pf:
        pf.write('END RUNNING\n')

    
    
    
    
    
    
    
