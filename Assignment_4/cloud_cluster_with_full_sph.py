'''after removing the bridge, we should use only one-directional channel from `stars` to `sph.dm_particles`,
and always remove the `stars` from the `cloud` before calculating the energy, Lagrange radii and so on.
Also remember to plot the `stars` on top of `sph.dm_particles`.
Anyway, try everything to distinguish between the cluster stars and newly formed stars!'''

import numpy as np
import matplotlib as mpl

from amuse.lab import * 
from amuse.units import units, nbody_system, quantities
from amuse.couple import bridge
from amuse.community.fi.interface import Fi
from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube  #what's this?
from amuse.ic.gasplummer import new_plummer_gas_model

# plot customization
from matplotlib import pyplot as plt
mpl.rc('xtick', direction='in', top=True, labelsize=14)
mpl.rc('xtick.major', size=5, width=2)
mpl.rc('xtick.minor', size=3, width=2)
mpl.rc('ytick', direction='in', right=True, labelsize=14)
mpl.rc('ytick.major', size=5, width=2)
mpl.rc('ytick.minor', size=3, width=2)
mpl.rc('axes', linewidth=2, labelsize=20, labelpad=10,\
       titlesize=25, titlepad=15)
mpl.rc('image', origin='lower')


printout_file = 'print_out.txt'


def make_map(sph, L=200):
    '''Makes 2D-projected density map of the gas particles in the sph code'''

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


def plot_cloud_cluster(cluster, sph, stars, title, vrange=[-6,2]):
    '''Makes snapshot of each evolutionary stage'''

    #gas_position = (sph.gas_particles.position).value_in(units.parsec)
    #L = 2 * (np.ceil(np.max(np.abs(gas_position))/100.) * 100.)
    L = 600

    rho = make_map(sph, L) 
    rhomin = int(np.log10(1.e-5+rho.value_in(units.amu/units.cm**3)).min())
    rhomax = np.ceil(np.log10(1.e-5+rho.value_in(units.amu/units.cm**3)).max())

    fig = plt.figure(figsize=(17,7))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    # plot scatters of gas particles only in the first subplot
    ax1.scatter(cloud.x.value_in(units.parsec), cloud.y.value_in(units.parsec), s=1, c='r',\
                label='cloud')

    # make colormap for the projected cloud density in the second subplot
    vmin, vmax = vrange
    cax = ax2.imshow(np.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)).transpose(),
                    extent=[-L/2, L/2, -L/2, L/2], vmin=vmin, vmax=vmax, cmap='jet', aspect="auto")

    # plot scatters of cluster stars and formed stars in both subplots
    for ax in [ax1,ax2]:
        if ax == ax1:
            cc = 'b'
        else:
            cc = 'm'
        ax.scatter(cluster.x.value_in(units.parsec), cluster.y.value_in(units.parsec), s=1, c=cc,\
                   label='cluster')
        if len(stars) > 0:
            ax.scatter(stars.x.value_in(units.parsec), stars.y.value_in(units.parsec),\
                       s=3, marker='*', c='k', label='formed star')

        ax.set_xlim(-300, 300)
        ax.set_ylim(-300, 300)


    # customize colorbar ticks for the colormap in the second subplot
    ticks = np.arange(rhomin, rhomax+1, 1, dtype=int)
    ticks = np.arange(vmin, vmax+1, 1, dtype=int)
    ticklabels = [str(i) for i in ticks]
    cbar = fig.colorbar(cax, ax=ax2, ticks=ticks, orientation='vertical', fraction=0.045)
    cbar.ax.set_yticklabels(ticklabels)  # horizontal colorbar
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('log(density) [$amu/cm^3$]', fontsize=15, rotation=270, labelpad=25)

    ax1.set_title(title+' Myr', weight='bold')
    ax2.set_title('projected density map of cloud', fontsize=15, pad=10)
    ax1.set_xlabel('x [pc]', weight='bold')
    ax1.set_ylabel('y [pc]', weight='bold')
    ax1.legend(fontsize=14)

    fig.subplots_adjust(wspace=0.2)
    plt.savefig(title+'.png', dpi=200, bbox_inches='tight')
    plt.close()

    return None
    

def Lagrange_radii(model, converter):
    ''' Returns the 0.2, 0.5 and 0.9 Lagrangia radii'''      
    lr,mr = model.LagrangianRadii(converter)
    return np.take(lr.value_in(units.parsec),[4,5,7])


def get_energy(particles, norm=1.):
    '''Calculates the kinetic/potential/total energies of a particle set'''
    Ek = particles.kinetic_energy()/norm
    Ep = particles.potential_energy()/norm
    Etot = Ek+Ep
    return np.array([Ek, Ep, Etot])


def Jeans_density(M, m=2*1.66e-24|units.g, T=15|units.K):
    '''calculates the Jeans density of a molecular cloud:
       `M` is the total mass of the cloud,
       `m` is the average mass of gas particles, here estimated as the mass of a molecular hydrogen,
       `T` is the equilibrium temperature of the cloud, here estimated as 15K.'''

    rho_J = 81*constants.kB**3 / (32*np.pi*constants.G**3) * T**3 / (m**3 * M**2)
    return rho_J.in_(units.kg/units.m**3)


def MRR(mass):
    '''Empirical mass-radius relation for a zero-age mein-sequence star (ZAMS)'''
    if mass > 1.0 | units.MSun:
        R = (mass.value_in(units.MSun))**0.57
    else:
        R = (mass.value_in(units.MSun))**0.80
    return (R | units.RSun)


def timescale(particles, unit=units.Myr, G=constants.G):
    '''Calculates the dynamical/relaxation/free-fall timescales for a particles set:
       N = number of particles, M = total mass, R = virial radius'''
    # dynamical timescale
    N = len(particles)
    M = particles.total_mass()
    R = particles.virial_radius()
    t_dyn = np.sqrt(R**3/(G*M)).value_in(unit)
    # relaxation timescale
    t_rh = 0.138*N/np.log(0.4*N) * t_dyn
    # free-fall timescale
    t_ff = particles.dynamical_timescale().value_in(unit)

    return t_dyn, t_rh, t_ff


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
    cluster.birth_age = -1. | units.Myr  # distinguish from the newly-formed stars in the cloud
    return cluster, converter


def cloud_init(Ngas, Mgas, Rgas):
    '''Initialises the cloud'''
    converter = nbody_system.nbody_to_si(Mgas,Rgas)
    #cloud = molecular_cloud(targetN=Ngas, convert_nbody=converter,\
    #                        base_grid=body_centered_grid_unit_cube).result
    cloud = new_plummer_gas_model(Ngas, convert_nbody=converter)
    return cloud, converter



def resolve_sinks(sph, stars, cloud, density_thr, model_time):
    '''Deals with potential star formation events in the cloud'''

    # select out the gas particles of which the density exceeds Jeans density
    high_dens = sph.gas_particles.select_array(lambda rho: rho > density_thr, ["rho"])

    if len(high_dens) <= 0:
        with open(printout_file,'a') as pf:
            pf.write('(no high density particle has been identified in this step)\n')
    # if there is any gas particles that has been selected out, consider a star formation event
    else:
        candidate_stars = high_dens.copy()
        # remove these high-density particles out of the sph code, and synchronize the change to the cloud
        sph.gas_particles.remove_particles(high_dens)
        sph.gas_particles.synchronize_to(cloud)

        # prevent stars forming duplicately
        with open(printout_file,'a') as pf:
            pf.write('Identified %d candidates; '%len(candidate_stars))
        newstars = Particles()
        for cs in candidate_stars:
            # if a candidate star has not been considered as a newly formed star before, add it
            if cs not in stars:
                # update the radius of the new star, otherwise it will inherit the same radius 
                # from the gas particle, which is unresonably large for a ZAMS star
                cs.radius = MRR(cs.mass)
                newstars.add_particle(cs)
        with open(printout_file,'a') as pf:
            pf.write('Confirmed %d new stars.\n'%len(candidate_stars))

        newstars.birth_age = model_time
        newstars.Lx = 0 | (units.g * units.m**2)/units.s
        newstars.Ly = 0 | (units.g * units.m**2)/units.s
        newstars.Lz = 0 | (units.g * units.m**2)/units.s

        N_star_pre = len(stars)
        stars.add_particles(newstars)
        sph.dm_particles.add_particles(newstars)
        N_star_after = len(stars)
        with open(printout_file,'a') as pf:
            pf.write('previous: %d stars, afterwards: %d stars\n'%(N_star_pre, N_star_after))


def merge_stars(sph, stars, merge_radius):
    if len(stars) <= 0 or stars.radius.max() <= (0.|units.AU):
        with open(printout_file,'a') as pf:
           pf.write('(no merge happens in this step)\n')
        return

    with open(printout_file,'a') as pf:
        pf.write('identifying merging groups...\n')
    # select out star pairs in which the distance between two members is smaller than `merge_radius` 
    ccs = stars.copy().connected_components(threshold=merge_radius)
    if len(ccs) <= 0:
        with open(printout_file,'a') as pf:
            pf.write('(no merge happens in this step)\n')
        return

    n_merge = 0
    for cc in ccs:
        if len(cc) > 1:
            n_merge += 1
            new_star = merge_two_stars(stars, cc.copy())
            sph.dm_particles.remove_particles(cc)
            sph.dm_particles.add_particles(new_star)
    

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
    with open(printout_file,'a') as pf:
        pf.write(f'Two stars with {np.round(particles_in_encounter.mass.value_in(units.MSun),1)} MSun '+\
                 f'collided at {np.round(com_pos.value_in(units.parsec),2)} pc;\n'+\
                 f'    new radius = {np.round(particles_in_encounter.radius.value_in(units.RSun),2)} RSun, '+\
                 f'new radius = {np.round(new_particle.radius.value_in(units.RSun),2)} RSun.\n')
    stars.remove_particles(particles_in_encounter)
    return new_particle




def evolve(cluster,cloud, converter_grav,converter_sph, t_end, dt_sph, dt_diag,\
           sink=True, merge=False):

    with open(printout_file, 'a') as pf:
        pf.write('Setting up the gravity code and the hydro code...\n')

    '''# Initialising the direct N-body integrator
    gravity = ph4(converter_grav)
    gravity.particles.add_particles(cluster)'''

    # Initialising the hydro code
    sph = Fi(converter_sph, mode="openmp")
    sph.gas_particles.add_particles(cloud)
    sph.dm_particles.add_particles(cluster)

    sph.parameters.radiation_flag = False
    #sph.parameters.isothermal_flag = True
    #sph.parameters.integrate_entropy_flag = False  #should we consider isothermal process or adiabatic one?
    sph.parameters.timestep = dt_sph
    #sph.parameters.eps_is_h_flag = False    # h_smooth is constant
    #eps = 0.1 | units.parsec
    #sph.parameters.gas_epsilon = eps
    #sph.parameters.sph_h_const = eps
    #cloud.h_smooth= eps


    '''# Building a bridge between hydro and grav
    with open(printout_file, 'a') as pf:
        pf.write('Bridging...\n')
    grav_sph = bridge.Bridge(use_threading=False)
    grav_sph.add_system(gravity, (sph,))
    grav_sph.add_system(sph, (gravity,))
    grav_sph.timestep = dt_bridge'''

    # Setting up channels from code to cloud and cluster
    #channel_from_grav_to_cluster = gravity.particles.new_channel_to(cluster)
    channel_from_sph_to_cloud = sph.gas_particles.new_channel_to(cloud)
    channel_from_sph_to_cluster = sph.dm_particles.new_channel_to(cluster)

    # consider star formation
    if sink == True:
        with open(printout_file,'a') as pf:
            pf.write('[Star formation is considered.]\n')
        stars = Particles(0)
        sph.dm_particles.add_particles(stars)

        density_threshold = Jeans_density(M=sph.gas_particles.mass.max())
        sph.parameters.stopping_condition_maximum_density = density_threshold
        density_limit_detection = sph.stopping_conditions.density_limit_detection
        density_limit_detection.enable()


    if merge == True:
        merge_radius = 1e-4 | units.parsec  # 1 pc = 206265 AU
        with open(printout_file,'a') as pf:
            pf.write('[Star merging is considered.]\n')

    # Initialize data lists (Lagrangian radii and relative total energy)
    lr_cloud_list = [Lagrange_radii(cloud, converter_sph)]
    lr_cluster_list = [Lagrange_radii(cluster, converter_grav)]

    E0_cloud = get_energy(cloud)
    Etot0_cloud = E0_cloud[-1]
    E_cloud_list = [E0_cloud/Etot0_cloud]

    E0_cluster = get_energy(cluster)
    Etot0_cluster = E0_cluster[-1]
    E_cluster_list = [E0_cluster/Etot0_cluster]

    n_formed_star = [0]


    # Start Evolving!
    # unify the unit of times
    unit_time = units.Myr
    t_end = t_end.in_(unit_time)
    dt_diag = dt_diag.in_(unit_time)
    dt_sph = dt_sph.in_(unit_time)
    times = quantities.arange(0.|unit_time, t_end+dt_diag, dt_diag)

    with open(printout_file, 'a') as pf:
        pf.write('\nStart evolving the molecular cloud...\n')
        pf.write(f'End time = {t_end}; Diagnostic timestep = {dt_diag}; sph code timestep = {dt_sph}.\n')

    for i,t in enumerate(times):
        with open(printout_file, 'a') as pf:
            pf.write(f'---------- Time = {t.in_(units.Myr)} ----------\n')

        # calculate the dynamical(dyn), half-mass relaxation(rh), and free-fall(ff) timescales
        # of both the cloud and the cluster systems
        t_dyn_cloud, t_rh_cloud, t_ff_cloud = timescale(cloud, unit_time)
        t_dyn_cluster, t_rh_cluster, t_ff_cluster = timescale(cluster, unit_time)
        all_timescales = [t_dyn_cloud, t_rh_cloud, t_ff_cloud, t_dyn_cluster, t_rh_cluster, t_ff_cluster] | unit_time
        # check if the bridge timestep should be reduced
        if all_timescales.amax() < 10 * dt_sph:
            dt_sph_old = dt_sph
            dt_sph = all_timescales.amax() / 10.
            with open(printout_file, 'a') as pf:
                pf.write(f'Change the bridge timestep from {dt_sph_old} to {dt_sph}.\n')

        if sink == True:
            # make sure at this time, elements in 'stars' and 'sph.gas_particles' are the same
            resolve_sinks(sph, stars, cloud, density_threshold, t)

        # evolve for one diagnostic timestep
        #grav_sph.evolve_model(t, timestep=dt_bridge)
        sph.evolve_model(t, timestep=dt_sph)

        #channel_from_grav_to_cluster.copy()
        channel_from_sph_to_cloud.copy()
        channel_from_sph_to_cluster.copy()
        if len(stars) > 0:
            newstars = cluster.select_array(lambda birth_age: birth_age >= 0.|unit_time , ["birth_age"])
            cluster.remove_particles(newstars)

        # merge two stars if they are too close to each other (distance < merge_radius)
        # typically we don't consider merging event since it is very rare given a reasonable merge radius
        if merge == True:
            merge_stars(sph, stars, merge_radius)

        with open(printout_file,'a') as pf:
            pf.write('Number of stars in `cluster`= %d; in `stars` = %d; in `sph.dm_particles`= %d.\n'\
                     %(len(cluster), len(stars), len(sph.dm_particles)))

        # make snapshots
        plot_cloud_cluster(cluster, sph, stars, title='{0}'.format(float(t.value_in(units.Myr))),\
                           vrange=[-5,3])

        # save data (energy will be added afterwards)
        lr_cloud = Lagrange_radii(cloud, converter_sph)
        lr_cloud_list.append(lr_cloud)
        lr_cluster = Lagrange_radii(cluster, converter_grav)
        lr_cluster_list.append(lr_cluster)

        E_cloud = get_energy(cloud, norm=Etot0_cloud)
        E_cloud_list.append(E_cloud)
        E_cluster = get_energy(cluster, norm=Etot0_cluster)
        E_cluster_list.append(E_cluster)

        n_formed_star.append(len(stars))

        # save data instantaneously
        lr_data = np.concatenate(([t.value_in(unit_time)], lr_cloud, lr_cluster))
        E_data = np.concatenate(([t.value_in(unit_time)], E_cloud, E_cluster))
        with open('lr_data.txt', 'a') as f_lr_data:
            f_lr_data.write(','.join(str(x) for x in lr_data)+'\n')
        with open('E_data.txt', 'a') as f_E_data:
            f_E_data.write(','.join(str(x) for x in E_data)+'\n')
        with open('formed_star_data.txt', 'a') as f_fs_data:
            f_fs_data.write('%f,%d\n'%(t.value_in(unit_time), len(stars)))

        with open(printout_file, 'a') as pf:
            pf.write(f'Finished.\n')

    #gravity.stop()
    sph.stop()

    return times.value_in(unit_time), n_formed_star, lr_cloud_list, lr_cluster_list, E_cloud_list, E_cluster_list





if __name__ == '__main__':

    with open(printout_file, 'w') as pf:
        pf.write('START RUNNING\n')
    with open('lr_data.txt', 'w') as f_lr_data:
        f_lr_data.write('# Time, lr2_cd, lr5_cd, lr9_cd, lr2_cr, lr5_cr, lr9_cr\n')
    with open('E_data.txt', 'w') as f_E_data:
        f_E_data.write('# Time, Ek_cd, Ep_cd, Etot_cd; Ek_cr, Ep_cr, Etot_cr\n')
    with open('formed_star_data.txt', 'w') as f_fs_data:
        f_fs_data.write('# Time, Number of newly formed star\n')
        

    # Set initial parameters
    with open(printout_file, 'a') as pf:
        pf.write('Initializing the star cluster and the molecular cloud...\n')
    # Molecular cloud (typically, mass = 1e3-1e7 MSun, diameter = 5-200 pc ?)
    N_cloud = int(5e4)
    Mtot_cloud = 2. | units.MSun * N_cloud
    Rvir_cloud = 10. | units.parsec

    # Cluster (typically, number of stars = 1e4-1e6, diameter = 3-100 pc ?)
    N_cluster = int(5e3)
    Mtot_cluster = 2. | units.MSun * N_cluster
    Rvir_cluster = 3. | units.parsec
    # initial velocity and position of the cluster's COM
    v_cluster = (3,3,0) | units.km/units.s
    p_cluster = (-200,-200,0) | units.parsec

    # Setting a seed
    np.random.seed(1)

    # Initialising the two systems
    cluster, conv_grav = cluster_init(N_cluster, Mtot_cluster, Rvir_cluster,\
                                              v_cluster, p_cluster)
    cloud, conv_sph = cloud_init(N_cloud, Mtot_cloud, Rvir_cloud)

    t_end = 100. | units.Myr
    dt_sph = 0.1 | units.Myr
    dt_diag = 1.0 | units.Myr
    
    for filename in [printout_file, 'lr_data.txt', 'E_data.txt', 'formed_star_data.txt']:
        with open(filename,'a') as f:
            f.write('# Cloud: number of gas particles: %d, total mass: %d MSun, virial radius: %d pc;\n'\
                     %(N_cloud, Mtot_cloud.value_in(units.MSun), Rvir_cloud.value_in(units.parsec)))
            f.write('# Cluster: number of stars: %d, total mass: %d MSun, virial radius: %d pc;\n'\
                     %(N_cluster, Mtot_cluster.value_in(units.MSun), Rvir_cluster.value_in(units.parsec)))
            f.write(f'# Initial position of cluster: {p_cluster.in_(units.parsec)}; '+\
                    f'Initial velocity of cluster: {v_cluster.in_(units.kms)}.\n')

    times, n_formed_star, lr_cloud, lr_cluster, E_cloud, E_cluster = evolve(cluster,cloud, conv_grav,conv_sph,\
                                                      t_end,dt_sph,dt_diag, sink=True, merge=False)

    #with open(printout_file,'a') as pf:
    #    pf.write('Saving data...\n')
    #lr_data = np.hstack((lr_cloud, lr_cluster))
    #np.savetxt('lr.csv', lr_data, delimiter=',')
    #E_data = np.column_stack((E_cloud, E_cluster))
    #np.savetxt('Etot.csv', Etot_data, delimiter=',')

    with open(printout_file,'a') as pf:
        pf.write('END RUNNING\n')

    
    
    
    
    
    
    
