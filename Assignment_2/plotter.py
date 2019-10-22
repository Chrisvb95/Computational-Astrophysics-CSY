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

def plot_coordinates(cluster, title, fig_dir, save=False):
   
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
    
    dump_dir = 'data_dump_direct'
    fig_dir = 'figs_ndirect/'
    mCut = '0'
    star_group = 'high_mass' 

    data_direct = read_dump(dump_dir, mCut, star_group)
    E0 = data_direct[0].potential_energy() + data_direct[0].kinetic_energy()
    Eer = np.zeros(len(data_direct))
    for i in range(len(data_direct)):
        Et = data_direct[i].potential_energy() + data_direct[i].kinetic_energy()     
        Eer[i] = np.abs((Et-E0)/E0)
    t = np.arange(0,10,0.1)
    
    plt.title('Relative energy error over time \n for direct integration')
    plt.xlabel('Time (Myr)')
    plt.ylabel('Relative energy error')
    plt.plot(t,Eer)
    plt.show() 
    
