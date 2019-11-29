import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

if __name__ == "__main__":

    data = np.loadtxt("lr.csv", delimiter=",")

    t_end = len(data[:,[0]])
    dt = 1
    t = np.arange(0, t_end, dt)    

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, xlabel= 'Time (Myr)', ylabel='Radius (pc)', title = 'Cluster radius over time')
    ax.plot(t,data[:, [0]], label='0.2 lr')
    ax.plot(t,data[:, [1]], label='0.5 lr')
    ax.plot(t,data[:, [2]], label='0.9 lr')
    ax.legend()
    plt.savefig('Cluster_radius_time.png')
    plt.close()
'''
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, xlabel= 'Time (Myr)', ylabel='Radius (pc)', title = 'Cluster radius over time')
    ax.plot(t,data[:, [3]], label='0.2 lr')
    ax.plot(t,data[:, [4]], label='0.5 lr')
    ax.plot(t,data[:, [5]], label='0.9 lr')
    ax.legend()
    plt.savefig('Cluster_radius_time.png')
    plt.close()
'''
