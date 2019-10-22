import numpy as np
import sys
import pickle as pkl
from matplotlib import pyplot as plt
from amuse.units import units
from amuse.units import quantities
from amuse.units import constants
from amuse.units import nbody_system

def read_dump(dump_dir, mCut, star_group):
    print('Reading data stored data...')
    data = []
    time_steps = np.arange(0,10,0.1)        
    for i in time_steps:
        fname = 't'+str(i)+'_m'+mCut+'_'+star_group+'.p'
        data.append(pkl.load(open(dump_dir+'/'+fname,'rb')))
        
    return time_steps, data        

def plot_cluster(cluster, title, fig_dir, save=False):
   
    m = cluster.mass.value_in(units.MSun)
    x = cluster.x.value_in(units.parsec)
    y = cluster.y.value_in(units.parsec)
    z = cluster.z.value_in(units.parsec)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x,y, c='b',zorder=0,marker='.')
    #ax.scatter(x_high,y_high, c='r',zorder=1,marker='.') 

    plt.xlim(left=-50,right=50)  
    plt.ylim(top=50,bottom=-50)
    
    if save:
        plt.savefig(fig_dir+title+'.png')
    else:    
        plt.show()
    plt.close()    

def plot_energy_error(times, dE_direct, dE_tree):
    fig_ee = plt.figure()
    ax = fig_ee.add_subplot(111, xlabel = 'time(Myr)', ylabel = 'relative energy error')
    ax.plot(times, dE_direct, linestyle = '-', label='all stars in direct code', color='k')
    ax.plot(times, dE_tree, linestyle = ':', label='all stars in tree code', color='k')
    #ax.set_yscale('log')
    ax.legend()
    plt.show()
    #fig_ee.savefig('energy_error_%.1f.jpg'%mCut.value_in(units.MSun))

def plot_radii(times, radii_direct, radii_tree, x_label, y_label):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel=x_label,  ylabel=y_label)
    ax.plot(times, radii_direct, linestyle = '-', label = 'all stars in direct code', color='k')
    ax.plot(times, radii_tree, linestyle = ':', label = 'all stars in tree code', color='k') 
    ax.legend()
    plt.show()

if __name__ == "__main__":
    '''
    dump_dir = 'data_dump_direct'
    fig_dir = 'figs_ndirect/'
    mCut = '0'
    star_group = 'high_mass' 

    data_direct = read_dump(dump_dir, mCut, star_group)
    for i in range(len(data_direct)):
        plot_cluster(data_direct[i],str(i),fig_dir,save=True)
    '''
    
    N = 10000
    mCluster = N|units.MSun
    rCluster = 3|units.parsec  
    converter=nbody_system.nbody_to_si(mCluster,rCluster)  
    
    # All stars in direct code
    dump_dir_direct = 'data_dump_direct'
    mCut = '0'
    star_group = 'high_mass'
    time_series, data_direct = read_dump(dump_dir_direct, mCut, star_group)
    
    # All stars in tree code
    dump_dir_tree = 'data_dump_tree'
    mCut = '100'
    star_group = 'low_mass'
    time_series, data_tree = read_dump(dump_dir_tree, mCut, star_group)
    
    print('Plotting energy error and radii...')
    # Saving the initial total energy
    print(data_direct[0].potential_energy)
    E0_direct = data_direct[0].potential_energy() + data_direct[0].kinetic_energy()    
    dE_direct = []
    E0_tree = data_tree[0].potential_energy() + data_tree[0].kinetic_energy()
    dE_tree = []
    
    print(E0_direct)
    
    direct_code_half_mass_radii = []
    direct_code_core_radii = [] 
    tree_code_half_mass_radii = []
    tree_code_core_radii = []

    for i in range(len(data_direct)):        

        print "time step : " + str(i * 0.1) + "Myr"      
        
        # direct code
        
        ## energy calculation
        
        U = data_direct[i].potential_energy()
        T = data_direct[i].kinetic_energy()
        Et = U + T
        dE_direct.append(np.abs((E0_direct-Et)/E0_direct))
        
        ## radii calculation
        pos,coreradius,coredens=data_direct[i].densitycentre_coreradius_coredens(converter)
        lr,mf=data_direct[i].LagrangianRadii(converter) # outputs are radii, mass fractions
        direct_code_half_mass_radii.append(lr[5]
                            .value_in(units.parsec)) # 5th argument attributes to half-mass radius
        direct_code_core_radii.append(coreradius.value_in(units.parsec))

        # tree code

        U = data_tree[i].potential_energy()
        T = data_tree[i].kinetic_energy()
        Et = U + T
        dE_tree.append(np.abs((E0_tree-Et)/E0_tree))

        ## radii calculation
        pos,coreradius,coredens=data_tree[i].densitycentre_coreradius_coredens(converter)
        lr,mf=data_tree[i].LagrangianRadii(converter) # outputs are radii, mass fractions
        tree_code_half_mass_radii.append(lr[5]
                            .value_in(units.parsec)) # 5th argument attributes to half-mass radius
        tree_code_core_radii.append(coreradius.value_in(units.parsec))       

        #plot_cluster(data_direct[i],str(i),fig_dir,save=False)
    # Plotting Energy Error
    plot_energy_error(time_series, dE_direct, dE_tree)
    # Plotting radii
    fig = plt.figure(figsize=(12, 8))
    x_label = 'time(Myr)'
    y_label = 'radius(pc)'
    ax = fig.add_subplot(111, xlabel=x_label,  ylabel=y_label)
    ax.plot(time_series, direct_code_half_mass_radii, linestyle = '-', 
        label = 'half-mass radius direct code', color='b')
    ax.plot(time_series, tree_code_half_mass_radii, linestyle = ':', 
        label = 'half-mass radius tree code', color='b')
    ax.plot(time_series, direct_code_core_radii, linestyle = '-', 
        label = 'core radius direct code', color='r')
    ax.plot(time_series, tree_code_core_radii, linestyle = ':', 
        label = 'core radius tree code', color='r') 
    ax.legend(bbox_to_anchor=(1.05, 0.5), loc='center left') 
    plt.tight_layout()
    plt.show()
    #plot_radii(time_series, direct_code_half_mass_radii,
    #        tree_code_half_mass_radii, x_label='time(Myr)', y_label='half-mass radius(pc)')
    #plot_radii(time_series, direct_code_core_radii,
    #        tree_code_core_radii, x_label='time(Myr)', y_label='core radius(pc)')
    



