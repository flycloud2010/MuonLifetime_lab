import os
import argparse
import itertools
from tqdm import tqdm

import numpy as np
import pandas as pd
import scipy as sp
import scipy.optimize as op
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator, MultipleLocator



# define ticks
def set_ticks(ax, xMaj, yMaj):
    ax.xaxis.set_major_locator(MultipleLocator(xMaj))
    ax.yaxis.set_major_locator(MultipleLocator(yMaj))
    ax.minorticks_on()
    ax.tick_params(which='major', width=1.0, length=8, direction='in', labelsize=14)
    ax.tick_params(which='minor', width=1.0, length=4, direction='in', labelsize=14)
    ax.yaxis.get_offset_text().set_fontsize(14)
    ax.xaxis.get_offset_text().set_fontsize(14)

def ang_sin(theta):
    '''
        Angular distribution formula 
    '''
    return(np.cos(theta) * np.cos(theta) * np.sin(theta))


def generate_angular(dist, limits, N):
    
    # extremes in x direction 
    u1 = limits[0] # minimum of x sampling 
    u2 = limits[1] # maximum of x sampling 

    u = np.random.uniform(u1,u2,5*N)

    lower = 0                                       # accept-reject minimum for each point
    upper = np.max(dist(u)) # accept-reject maximum for each point 

    v = np.random.uniform(lower, upper, 5*N)

    # accepted points 
    points = u[v < dist(u)]

    # keep only N events if we have more 
    if len(points) >= N:
        points = points[:N]

    return points


def generate_xy(theta, x0,y0,phi,N, h=8.):
    
    x_ = x0 - (h/np.cos(theta)) * np.sin(theta)*np.cos(phi)
    y_ = y0 - (h/np.cos(theta)) * np.sin(theta)*np.sin(phi)


    y_select =  np.logical_and(y_ >= 0, y_ <= 20)
    x_select = np.logical_and(x_ >= 0, x_ <= 182)
    mask = np.logical_and(x_select,y_select)

    x = x_[mask]
    y = y_[mask]

    n_coincidences = N * (len(x)/N) / 100

    return x_, y_, x, y, n_coincidences



def propagate_muons(seed):

    N = 6067

    np.random.seed((seed))

    x0 = np.random.uniform(0,182, size = N)
    y0 = np.random.uniform(0,20, size = N)
    phi = np.random.uniform(0,2*np.pi, size = N)

    thetas = generate_angular(ang_sin,(0.,0.5*np.pi), N)

    coincidences = []

    for i in [8,16,24,32]:
        _, _, _, _, n_coinc = generate_xy(thetas, x0, y0, phi, N, i)
        coincidences.append(n_coinc)

    return coincidences



def main(N_sim):

    data = np.zeros((N_sim,4))

    for i in tqdm(range(0,N_sim)):
        seed = i+1
        rate = propagate_muons(seed)
        data[i,:] = rate
             
    print("Number of simulations: ", N_sim)


    for i in range(0,4):
        mean = np.mean(data[:,i])

        mean_error = np.std(data[:,i])/np.sqrt(len(data[:,i]))

        print("slab %1.0f rate : %1.4f ± %1.4f"% (i+2 ,mean, mean_error))



def parse_args():
    '''Parse the args.'''
    parser = argparse.ArgumentParser()

    parser.add_argument('-N', '--simulations', type = int, required = True,
                        help = 'Number of simulations',dest='N_sim')

    return parser.parse_args()

    
if __name__ == '__main__':
    args = parse_args()
    main(N_sim = args.N_sim)
