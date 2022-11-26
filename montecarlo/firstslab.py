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
import seaborn as sns 
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator, MultipleLocator


## Constants 

L = 182 # cm
W = 20  # cm
H = 2.68


PI = np.pi

sns.set(style = 'white')
mpl.rc('xtick.minor', visible = True)
mpl.rc('ytick.minor', visible = True)
mpl.rc('xtick', direction='in', top=True, bottom = True)
mpl.rc('ytick', direction='in', right=True, left = True)

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



def coincidences(X_tilde):

    L2 = L+2*X_tilde
    W2 = W+2*X_tilde 
    D = 15

    N_events = int((L2*W2) * (100/60))

    #print('N events: ',N_events)

    counter, entries = 0,0
    
    x0, y0, z0, theta, phi = np.zeros(N_events), np.zeros(N_events), np.zeros(N_events), np.zeros(N_events), np.zeros(N_events)

    i = 0

    counter = 0
    xf, yf, zf = np.zeros(N_events), np.zeros(N_events), np.zeros(N_events)
    x, y, z = np.zeros(N_events), np.zeros(N_events), np.zeros(N_events)


    i = 0

    while(i<N_events):

        if counter >= N_events: break

        x0[i] = np.random.uniform(-X_tilde,L + X_tilde, size = 1)
        y0[i]= np.random.uniform(-X_tilde,W + X_tilde, size = 1)
        phi[i] = np.random.uniform(0,2*np.pi, size = 1)
        theta[i] = generate_angular(ang_sin,(0.,0.5*np.pi), 1)

        entries = entries + 1
        if x0[i] >= 0 and x0[i] <= L and y0[i] >= 0 and y0[i] <= W:
            counter = counter + 1 
            z[i] = 0  
            x[i] = x0[i]
            y[i] = y0[i]
            i = i+1
        elif x0[i] < 0 and phi[i] > PI / 2 and phi[i] < PI and (y0[i] - x0[i] * np.tan(phi[i])) > 0 and (y0[i] - x0[i] * np.tan(phi[i])) < W and theta[i] > np.arctan(x0[i] / (np.cos(phi[i]) * H)):
            counter = counter + 1 
            z[i] = -x0[i] / (np.cos(phi[i]) * np.tan(theta[i]))
            x[i] = 0
            y[i] = (y0[i] - x0[i] * np.tan(phi[i]))
            i = i+1
        elif x0[i] < 0 and phi[i] > PI / 2 and phi[i] < PI and (y0[i] + x0[i] * np.tan(PI-phi[i])) > W and x0[i]-(D-y0[i])/ np.tan(PI-phi[i]) > 0 and x0[i] - (W - y0[i]) / np.tan(PI-phi[i]) <L and theta[i] > np.arctan((y0[i]-W) / (np.sin(PI-phi[i]) * H)):
            counter = counter + 1 
            z[i] = -(y0[i] - W) / (np.sin(PI - phi[i]) * np.tan(theta[i]))
            x[i] = x0[i] - (W - y0[i]) / np.tan(PI-phi[i])
            y[i] = W
            i = i+1
        elif x0[i] < 0 and phi[i] >PI and phi[i] < 3*PI/2 and (y0[i] - x0[i] * np.tan(phi[i])) > 0 and (y0[i] - x0[i] * np.tan(phi[i])) < D and theta[i] > np.arctan(-x0[i] / (np.cos(PI-phi[i]) * H)):
            counter = counter + 1 
            z[i] = -x0[i] / (np.cos(phi[i]) * np.tan(theta[i]))
            x[i] = 0 
            y[i] = y0[i] - x0[i] * np.tan(phi[i])
            i = i+1
        elif x0[i] < 0 and phi[i] >PI and phi[i] < 3 * PI / 2 and x0[i] - y0[i] / np.tan(phi[i]) > 0 and x0[i] - y0[i] / np.tan(phi[i]) < L and (y0[i] - x0[i] * np.tan(phi[i])) < 0 and theta[i] > np.arctan(-y0[i] / (np.sin(phi[i]-PI) * H)):
            counter = counter + 1 

            z[i] = y0[i] / (np.sin(phi[i] - PI) * np.tan(theta[i]))
            y[i] = 0
            x[i] = x0[i] - y0[i] / np.tan(phi[i])
            i = i+1

        elif x0[i] > L and phi[i] > 0 and phi[i] < PI / 2 and (y0[i] - (x0[i] - L) * np.tan(phi[i])) > 0 and (y0[i] - (x0[i] - L) * np.tan(phi[i])) < W and theta[i] > np.arctan((x0[i] - L) / (np.cos(phi[i]) * H)):
            counter = counter + 1 

            z[i] = - (x0[i] - L) / (np.cos(phi[i]) * np.tan(theta[i]))
            x[i] = L
            y[i] = y0[i] - ( x0[i] - L ) * np.tan(phi[i])
            i = i+1
        elif x0[i] > L and phi[i] > 0 and phi[i] < PI / 2 and (y0[i] - (x0[i] - L) * np.tan(phi[i])) > W and x0[i] - (y0[i] - W) / np.tan(phi[i]) > 0 and x0[i] - (y0[i] - W) / np.tan(phi[i]) < L and theta[i] > np.arctan((y0[i]-W) / (np.sin(phi[i]) * H)):
            counter = counter + 1 

            z[i] = -(y0[i] - W) / (np.sin(phi[i]) * np.tan(theta[i]))
            y[i] = W
            x[i] = x0[i] - (y0[i] - W) / np.tan(phi[i])
            i = i+1
        elif x0[i] > L and phi[i] > 3 * PI / 2 and (y0[i] + (x0[i]-L) / np.tan(phi[i] - 3 * PI / 2)) > 0  and (y0[i] + (x0[i]-L) / np.tan(phi[i] - 3 * PI / 2)) < W and theta[i] > np.arctan((x0[i]-L) / (np.sin(phi[i] - 3 * PI / 2) * H)):
            counter = counter + 1 
    
            z[i]  = -(x0[i] - L) / (np.sin(phi[i] - 3 * PI / 2) * np.tan(theta[i]))
            x[i] = L
            y[i] = y0[i] + (x0[i] - L) / np.tan(phi[i] - 3 * PI / 2)
            i = i+1
        elif x0[i] > L and phi[i] > 3 * PI / 2 and x0[i] + y0[i] * np.tan(phi[i] - 3 * PI / 2) > 0 and x0[i] + y0[i] * np.tan(phi[i] - 3 * PI / 2) < L and (y0[i] + (x0[i] - L) / np.tan(phi[i] - 3 * PI / 2)) < 0 and theta[i] > np.arctan(-y0[i] / (np.sin(2*PI-phi[i]) * H)):
            counter = counter + 1 
    
            z[i]  = y0[i] / (np.sin(2 * PI - phi[i]) * np.tan(theta[i]))
            x[i] = x0[i] + y0[i] * np.tan(phi[i]-3*PI/2); 
            y[i] = 0
            i = i+1
        elif y0[i] < 0 and phi[i] > PI and phi[i] < 3 * PI/2 and (x0[i] - y0[i] / np.tan(phi[i]-PI)) > 0 and (x0[i] - y0[i] / np.tan(phi[i]-PI)) < L and theta[i] > np.arctan(-y0[i] / (np.sin(phi[i]-PI) * H)):
            counter = counter + 1 
    
            z[i]  = y0[i] / (np.sin(phi[i]-PI) * np.tan(theta[i]))
            x[i] = x0[i] - y0[i] / np.tan(phi[i]-PI); 
            y[i] = 0
            i = i+1
        elif y0[i] < 0 and phi[i] > 3 * PI / 2 and phi[i] < 2 * PI and (x0[i] + y0[i] / np.tan(2 * PI - phi[i])) > 0 and (x0[i] + y0[i] / np.tan(2*PI-phi[i])) < L and theta[i] > np.arctan(-y0[i] / (np.sin(2*PI-phi[i]) * H)):
            counter = counter + 1 
            z[i]  = y0[i] / (np.sin(2 * PI - phi[i]) * np.tan(theta[i]))
            x[i] = x0[i] + y0[i] / np.tan(2 * PI - phi[i])
            y[i] = 0
            i = i+1
        elif y0[i] > W and phi[i] > 0 and phi[i]  < PI / 2 and (x0[i] - (y0[i] - W) / np.tan(phi[i])) > 0 and (x0[i] - (y0[i] - W) / np.tan(phi[i])) < L and theta[i] > np.arctan((y0[i] - W) / (np.sin(phi[i]) * H)):
            counter = counter + 1 
            z[i] = -(y0[i] - W) / (np.sin(phi[i]) * np.tan(theta[i]))
            x[i] = x0[i] - (y0[i] - W) / np.tan(phi[i])
            y[i] = W
            i = i+1
        elif y0[i] > W and phi[i] > PI / 2 and phi[i] < PI and (x0[i] + (y0[i] - W) / np.tan(PI - phi[i])) > 0 and (x0[i] + (y0[i] - W) / np.tan(PI - phi[i])) < L  and theta[i] > np.arctan ( (y0[i] - W) / (np.sin(PI-phi[i]) * H)):
            counter = counter + 1 
            z[i] = -(y0[i] - W) / (np.sin(PI-phi[i]) * np.tan(theta[i]))
            x[i] = x0[i] + (y0[i] - W) / np.tan(PI - phi[i])
            y[i] = W
            i = i+1
        else:
            i = i+1
            continue

        
    x_f,y_f,z_f = [],[],[]


    x_f.append(x0)
    y_f.append(y0)
    z_f.append(z0)


    
    for i in range(0, N_events):
        xf[i]= x0[i] - (D + z[i]) * np.tan(theta[i]) * np.cos(phi[i])
        yf[i]= y0[i] - (D + z[i]) * np.tan(theta[i]) * np.sin(phi[i]) 
        zf[i]= -D


    x_f.append(xf)
    y_f.append(yf)
    z_f.append(zf)


    
    counterz = 0
    counterzcoin = 0
    countercoinc = 0

    for i in range(0,N_events):

        if 0 <= xf[i] and xf[i] <= L and yf[i] >= 0 and yf[i] <= W:
            countercoinc+=1
        if z[i] != 0: 
            counterz+=1
        if 0 <= xf[i] and xf[i] <= L and yf[i] >= 0 and yf[i] <= W and z[i] != 0:
            counterzcoin+=1

           

    #print('Extension: ', X_tilde)        
    #print('entries: ', entries)
    print('Coincidences: ',countercoinc )
    print('Perc. coinc: ', countercoinc / N_events * 100)
    print('Lateral events: ', counterz)
    print('Lateral events in coincidence: ', counterzcoin)

  
    return N_events, countercoinc
 

def main():

    n, ntot = [],[]

    for i in tqdm(range(5,60,5)):
        for j in range(0,100):
            Ntot, rate = coincidences(i)
            n.append(rate)
        ntot.append(Ntot)


    ext = range(5,60,5)

    plt.plot(ext, n)
    plt.show()


   

    
def parse_args():
    '''Parse the args.'''
    parser = argparse.ArgumentParser()

    parser.add_argument('-ex', '--extension', type = int, required = True,
                        help = 'Extension size',dest='X_tilde')


    return parser.parse_args()

    
if __name__ == '__main__':
    #args = parse_args()
    main()

		


