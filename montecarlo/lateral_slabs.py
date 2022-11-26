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

#N_events = 10*int((L*W + 2*L*H + 2*W*H) * (100/60)) # events in 100s in whole surface 
N_events = 6016
#N_events = 6183 # events in 100s in whole surface 

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



def lateral(X_tilde):

    L2 = L+2*X_tilde
    W2 = W+2*X_tilde 

    #N_events = int((L*W + 2*L*H + 2*W*H) * (100/60)) # events in 100s in whole surface 
    N_events = 6016

    #print("N events: ",N_events)

    counter, entries = 0,0
    
    x0, y0, z0, theta, phi = np.zeros(N_events), np.zeros(N_events), np.zeros(N_events), np.zeros(N_events), np.zeros(N_events)

    i = 0

    for D in [8,16,24,32]:
        globals()[f'counter{D}'] = 0
        globals()[f'xf{D}'], globals()[f'yf{D}'], globals()[f'zf{D}'] = np.zeros(N_events), np.zeros(N_events), np.zeros(N_events)
        globals()[f'x{D}'], globals()[f'y{D}'], globals()[f'z{D}'] = np.zeros(N_events), np.zeros(N_events), np.zeros(N_events)


    for D in [8,16,24,32]:
        i = 0

        while(i>=-1):

            if globals()[f'counter{D}'] >= N_events: break

            if D == 8:
                x0[i] = np.random.uniform(-X_tilde,L + X_tilde, size = 1)
                y0[i]= np.random.uniform(-X_tilde,W + X_tilde, size = 1)
                phi[i] = np.random.uniform(0,2*np.pi, size = 1)
                theta[i] = generate_angular(ang_sin,(0.,0.5*np.pi), 1)


            entries = entries + 1
            #print(D,i, globals()[f'counter{D}'],entries)
            if x0[i] >= 0 and x0[i] <= L and y0[i] >= 0 and y0[i] <= W:
                globals()[f'counter{D}'] = globals()[f'counter{D}'] + 1 
                globals()[f'z{D}'][i] = 0  
                globals()[f'x{D}'][i] = x0[i]
                globals()[f'y{D}'][i] = y0[i]
                i = i+1
            elif x0[i] < 0 and phi[i] > PI / 2 and phi[i] < PI and (y0[i] - x0[i] * np.tan(phi[i])) > 0 and (y0[i] - x0[i] * np.tan(phi[i])) < W and theta[i] > np.arctan(x0[i] / (np.cos(phi[i]) * H)):
                globals()[f'counter{D}'] = globals()[f'counter{D}'] + 1 
                globals()[f'z{D}'][i] = -x0[i] / (np.cos(phi[i]) * np.tan(theta[i]))
                globals()[f'x{D}'][i] = 0
                globals()[f'y{D}'][i] = (y0[i] - x0[i] * np.tan(phi[i]))
                i = i+1
            elif x0[i] < 0 and phi[i] > PI / 2 and phi[i] < PI and (y0[i] + x0[i] * np.tan(PI-phi[i])) > W and x0[i]-(D-y0[i])/ np.tan(PI-phi[i]) > 0 and x0[i] - (W - y0[i]) / np.tan(PI-phi[i]) <L and theta[i] > np.arctan((y0[i]-W) / (np.sin(PI-phi[i]) * H)):
                globals()[f'counter{D}'] = globals()[f'counter{D}'] + 1 
                globals()[f'z{D}'][i] = -(y0[i] - W) / (np.sin(PI - phi[i]) * np.tan(theta[i]))
                globals()[f'x{D}'][i] = x0[i] - (W - y0[i]) / np.tan(PI-phi[i])
                globals()[f'y{D}'][i] = W
                i = i+1
            elif x0[i] < 0 and phi[i] >PI and phi[i] < 3*PI/2 and (y0[i] - x0[i] * np.tan(phi[i])) > 0 and (y0[i] - x0[i] * np.tan(phi[i])) < D and theta[i] > np.arctan(-x0[i] / (np.cos(PI-phi[i]) * H)):
                globals()[f'counter{D}'] = globals()[f'counter{D}'] + 1 
                globals()[f'z{D}'][i] = -x0[i] / (np.cos(phi[i]) * np.tan(theta[i]))
                globals()[f'x{D}'][i] = 0 
                globals()[f'y{D}'][i] = y0[i] - x0[i] * np.tan(phi[i])
                i = i+1
            elif x0[i] < 0 and phi[i] >PI and phi[i] < 3 * PI / 2 and x0[i] - y0[i] / np.tan(phi[i]) > 0 and x0[i] - y0[i] / np.tan(phi[i]) < L and (y0[i] - x0[i] * np.tan(phi[i])) < 0 and theta[i] > np.arctan(-y0[i] / (np.sin(phi[i]-PI) * H)):
                globals()[f'counter{D}'] = globals()[f'counter{D}'] + 1 

                globals()[f'z{D}'][i] = y0[i] / (np.sin(phi[i] - PI) * np.tan(theta[i]))
                globals()[f'y{D}'][i] = 0
                globals()[f'x{D}'][i] = x0[i] - y0[i] / np.tan(phi[i])
                i = i+1

            elif x0[i] > L and phi[i] > 0 and phi[i] < PI / 2 and (y0[i] - (x0[i] - L) * np.tan(phi[i])) > 0 and (y0[i] - (x0[i] - L) * np.tan(phi[i])) < W and theta[i] > np.arctan((x0[i] - L) / (np.cos(phi[i]) * H)):
                globals()[f'counter{D}'] = globals()[f'counter{D}'] + 1 

                globals()[f'z{D}'][i] = - (x0[i] - L) / (np.cos(phi[i]) * np.tan(theta[i]))
                globals()[f'x{D}'][i] = L
                globals()[f'y{D}'][i] = y0[i] - ( x0[i] - L ) * np.tan(phi[i])
                i = i+1
            elif x0[i] > L and phi[i] > 0 and phi[i] < PI / 2 and (y0[i] - (x0[i] - L) * np.tan(phi[i])) > W and x0[i] - (y0[i] - W) / np.tan(phi[i]) > 0 and x0[i] - (y0[i] - W) / np.tan(phi[i]) < L and theta[i] > np.arctan((y0[i]-W) / (np.sin(phi[i]) * H)):
                globals()[f'counter{D}'] = globals()[f'counter{D}'] + 1 

                globals()[f'z{D}'][i] = -(y0[i] - W) / (np.sin(phi[i]) * np.tan(theta[i]))
                globals()[f'y{D}'][i] = W
                globals()[f'x{D}'][i] = x0[i] - (y0[i] - W) / np.tan(phi[i])
                i = i+1
            elif x0[i] > L and phi[i] > 3 * PI / 2 and (y0[i] + (x0[i]-L) / np.tan(phi[i] - 3 * PI / 2)) > 0  and (y0[i] + (x0[i]-L) / np.tan(phi[i] - 3 * PI / 2)) < W and theta[i] > np.arctan((x0[i]-L) / (np.sin(phi[i] - 3 * PI / 2) * H)):
                globals()[f'counter{D}'] = globals()[f'counter{D}'] + 1 
    
                globals()[f'z{D}'][i]  = -(x0[i] - L) / (np.sin(phi[i] - 3 * PI / 2) * np.tan(theta[i]))
                globals()[f'x{D}'][i] = L
                globals()[f'y{D}'][i] = y0[i] + (x0[i] - L) / np.tan(phi[i] - 3 * PI / 2)
                i = i+1
            elif x0[i] > L and phi[i] > 3 * PI / 2 and x0[i] + y0[i] * np.tan(phi[i] - 3 * PI / 2) > 0 and x0[i] + y0[i] * np.tan(phi[i] - 3 * PI / 2) < L and (y0[i] + (x0[i] - L) / np.tan(phi[i] - 3 * PI / 2)) < 0 and theta[i] > np.arctan(-y0[i] / (np.sin(2*PI-phi[i]) * H)):
                globals()[f'counter{D}'] = globals()[f'counter{D}'] + 1 
    
                globals()[f'z{D}'][i]  = y0[i] / (np.sin(2 * PI - phi[i]) * np.tan(theta[i]))
                globals()[f'x{D}'][i] = x0[i] + y0[i] * np.tan(phi[i]-3*PI/2); 
                globals()[f'y{D}'][i] = 0
                i = i+1
            elif y0[i] < 0 and phi[i] > PI and phi[i] < 3 * PI/2 and (x0[i] - y0[i] / np.tan(phi[i]-PI)) > 0 and (x0[i] - y0[i] / np.tan(phi[i]-PI)) < L and theta[i] > np.arctan(-y0[i] / (np.sin(phi[i]-PI) * H)):
                globals()[f'counter{D}'] = globals()[f'counter{D}'] + 1 
    
                globals()[f'z{D}'][i]  = y0[i] / (np.sin(phi[i]-PI) * np.tan(theta[i]))
                globals()[f'x{D}'][i] = x0[i] - y0[i] / np.tan(phi[i]-PI); 
                globals()[f'y{D}'][i] = 0
                i = i+1
            elif y0[i] < 0 and phi[i] > 3 * PI / 2 and phi[i] < 2 * PI and (x0[i] + y0[i] / np.tan(2 * PI - phi[i])) > 0 and (x0[i] + y0[i] / np.tan(2*PI-phi[i])) < L and theta[i] > np.arctan(-y0[i] / (np.sin(2*PI-phi[i]) * H)):
                globals()[f'counter{D}'] = globals()[f'counter{D}'] + 1 
                globals()[f'z{D}'][i]  = y0[i] / (np.sin(2 * PI - phi[i]) * np.tan(theta[i]))
                globals()[f'x{D}'][i] = x0[i] + y0[i] / np.tan(2 * PI - phi[i])
                globals()[f'y{D}'][i] = 0
                i = i+1
            elif y0[i] > W and phi[i] > 0 and phi[i]  < PI / 2 and (x0[i] - (y0[i] - W) / np.tan(phi[i])) > 0 and (x0[i] - (y0[i] - W) / np.tan(phi[i])) < L and theta[i] > np.arctan((y0[i] - W) / (np.sin(phi[i]) * H)):
                globals()[f'counter{D}'] = globals()[f'counter{D}'] + 1 
                globals()[f'z{D}'][i] = -(y0[i] - W) / (np.sin(phi[i]) * np.tan(theta[i]))
                globals()[f'x{D}'][i] = x0[i] - (y0[i] - W) / np.tan(phi[i])
                globals()[f'y{D}'][i] = W
                i = i+1
            elif y0[i] > W and phi[i] > PI / 2 and phi[i] < PI and (x0[i] + (y0[i] - W) / np.tan(PI - phi[i])) > 0 and (x0[i] + (y0[i] - W) / np.tan(PI - phi[i])) < L  and theta[i] > np.arctan ( (y0[i] - W) / (np.sin(PI-phi[i]) * H)):
                globals()[f'counter{D}'] = globals()[f'counter{D}'] + 1 
                globals()[f'z{D}'][i] = -(y0[i] - W) / (np.sin(PI-phi[i]) * np.tan(theta[i]))
                globals()[f'x{D}'][i] = x0[i] + (y0[i] - W) / np.tan(PI - phi[i])
                globals()[f'y{D}'][i] = W
                i = i+1
            else: continue

        
    x_f,y_f,z_f = [],[],[]
    
    x_f.append(x0)
    y_f.append(y0)
    z_f.append(z8)



    for D in [8,16,24,32]:
        for i in range(0, N_events):
            globals()[f'xf{D}'][i]= x0[i] - (D + globals()[f'z{D}'][i]) * np.tan(theta[i]) * np.cos(phi[i])
            globals()[f'yf{D}'][i]= y0[i] - (D + globals()[f'z{D}'][i]) * np.tan(theta[i]) * np.sin(phi[i]) 
            globals()[f'zf{D}'][i]= -D


    for i in [8,16,24,32]:
        x_f.append(globals()[f'xf{i}'])
        y_f.append(globals()[f'yf{i}'])
        z_f.append(globals()[f'zf{i}'])


    index = []

    for D in [8,16,24,32]:

        globals()[f'countercoinc{D}'] = 0
        counterz = 0
        counterzcoin = 0
        
        for i in range(0,N_events):
    
            if 0 <= globals()[f'xf{D}'][i] and globals()[f'xf{D}'][i] <= L and globals()[f'yf{D}'][i] >= 0 and globals()[f'yf{D}'][i] <= W:
                globals()[f'countercoinc{D}']+=1

            if globals()[f'z{D}'][i] != 0: 
                counterz+=1

            if 0 <= globals()[f'xf{D}'][i] and globals()[f'xf{D}'][i] <= L and globals()[f'yf{D}'][i] >= 0 and globals()[f'yf{D}'][i] <= W and globals()[f'z{D}'][i] != 0:
                counterzcoin+=1
                if D == 32:
                    index.append(i)
               
                
        ##print("entries: ", entries)
        print("Coincidences between first and considered slab: ", globals()[f'countercoinc{D}'] )
        ##print("Perc. coinc: ", globals()[f'countercoinc{D}'] / N_events * 100)
        print("Lateral events: ", counterz)
        ##print("Lateral events coincidence: ", counterzcoin)

    counterzcoin = []
    for D in [8,16,24,32]:
        counterzcoin.append(globals()[f'countercoinc{D}'])

    return x_f, y_f, z_f, index, counterzcoin

def main(X_tilde):
    
    xf, yf, zf, idx = lateral(X_tilde)

    x, y, z = [],[], []

    for i in range(0,5):
        for j in range(0,N_events):
            x.append(xf[i][j])
            y.append(yf[i][j])
            z.append(zf[i][j])


    fig = plt.figure(figsize=(10,10))

    ax = fig.add_subplot(111, projection='3d')
    cmap = mpl.cm.get_cmap('YlGnBu')

    sctt = ax.scatter3D(x, y, z,
                    alpha = 0.8,
                    c = z,
                    cmap = cmap,
                    marker = '.', s = 1)


    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    ax.set_zlabel("z [cm]")
    ax.yaxis.set_major_locator(MultipleLocator(200))

    fig.colorbar(sctt, ax = ax,pad = 0.027, shrink=0.47)

    x = [0,182]
    y = [20,20]
    z = [-32,-32]
    ax.plot(x,y,z, c = 'red', zorder = 1000)
    x = [0,182]
    y = [0,0]
    z = [-32,-32]
    ax.plot(x,y,z, c = 'red', zorder = 1000)
    x = [0,0]
    y = [0,20]
    z = [-32,-32]
    ax.plot(x,y,z, c = 'red', zorder = 1000)
    x = [182,182]
    y = [0,20]
    z = [-32,-32]
    ax.plot(x,y,z, c = 'red', zorder = 1000)


    x = [0,182]
    y = [20,20]
    z = [0,0]
    ax.plot(x,y,z, c = 'red', zorder = 1000)
    x = [0,182]
    y = [0,0]
    z = [0,0]
    ax.plot(x,y,z, c = 'red', zorder = 1000)
    x = [0,0]
    y = [0,20]
    z = [0,0]
    ax.plot(x,y,z, c = 'red', zorder = 1000)
    x = [182,182]
    y = [0,20]
    z = [0,0]
    ax.plot(x,y,z, c = 'red', zorder = 1000)


    x = [-X_tilde,182+X_tilde]
    y = [20+X_tilde,20+X_tilde]
    z = [0,0]
    ax.plot(x,y,z, c = 'black', zorder = 1000, ls = '--')
    x = [-X_tilde,182+X_tilde]
    y = [-X_tilde,-X_tilde]
    z = [0,0]
    ax.plot(x,y,z, c = 'black', zorder = 1000, ls = '--')
    x = [-X_tilde,-X_tilde]
    y = [-X_tilde,20+X_tilde]
    z = [0,0]
    ax.plot(x,y,z, c = 'black', zorder = 1000, ls = '--')
    x = [182+X_tilde,182+X_tilde]
    y = [-X_tilde,20+X_tilde]
    z = [0,0]
    ax.plot(x,y,z, c = 'black', zorder = 1000, ls = '--')


    ax.set_xlim(-300,400)
    ax.set_ylim(-400,400)


    index = idx[0]

    xline = [xf[0][index], xf[4][index]]
    yline = [yf[0][index], yf[4][index]]
    zline = [zf[0][index], zf[4][index]]
    
    ax.plot(xline, yline, zline, color = 'b')

    xline = [xf[0][index], xf[4][index]]
    yline = [yf[0][index], yf[4][index]]
    zline = [15, zf[4][index]]
    
    ax.plot(xline, yline, zline, color = 'r')

    ax.view_init(20, 120)
    fig.tight_layout()

    fig.savefig('./plots/3dslabs.png', dpi=400)
    plt.show()

    return
 

    
    
def parse_args():
    '''Parse the args.'''
    parser = argparse.ArgumentParser()

    parser.add_argument('-ex', '--extension', type = int, required = True,
                        help = 'Extension size',dest='X_tilde')


    return parser.parse_args()

    
if __name__ == '__main__':
    args = parse_args()
    main(X_tilde = args.X_tilde)#, D = args.D)

		


