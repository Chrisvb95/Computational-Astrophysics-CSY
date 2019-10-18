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
    for i in np.arange(0,10,0.1):
        fname = 't'+str(i)+'_m'+mCut+'_'+star_group+'.p'
        data.append(pkl.load(open(dump_dir+'/'+fname,'rb')))
    return data        

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
    dump_dir = 'data_dump_tree'
    fig_dir = 'figs_tree/'
    mCut = '100'
    star_group = 'low_mass' 

    data_direct = read_dump(dump_dir, mCut, star_group)
    print('Plotting data...')
    for i in range(len(data_direct)):
        plot_cluster(data_direct[i],str(i),fig_dir,save=True)
