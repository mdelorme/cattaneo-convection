import pyPLUTO as pp
import numpy as np
import matplotlib.pyplot as plt
import os, shutil, sys
from multiprocessing import Pool

def create_directory(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.mkdir(path)

def count_pluto_snaps(path='.', ext='.dbl'):
    snap_count = 0
    for f in os.listdir(path):
        if f.startswith('data') and f.endswith(ext):
            snap_count += 1
    return snap_count

def extract_Tbar(sid):
    d = pp.pload(sid)
    T = d.prs / d.rho
    Tbar = np.average(T, axis=0)
    return d.time, Tbar

def plot_2D(sid):
    d = pp.pload(sid)
    fig, ax = plt.subplots(3, 2, figsize = (15, 7))

    rho = d.rho
    T   = d.prs / d.rho

    rho_bar = np.average(rho, axis=0)
    T_bar   = np.average(T,   axis=0)

    rho_dev = rho-rho_bar
    T_dev   = T-T_bar

    T = T[:,::-1].T
    rho = rho[:,::-1].T
    T_dev = T_dev[:,::-1].T
    rho_dev = rho_dev[:,::-1].T
    uz      = d.vx2[:,::-1].T
    uh      = d.vx1[:,::-1].T

    ext = [d.x1[0], d.x1[-1], d.x2[0], d.x2[-1]]

    ax[0,0].imshow(rho, origin='lower', extent=ext, cmap='viridis')
    ax[0,1].imshow(rho_dev, origin='lower', extent=ext, cmap='bwr', clim=(-2.0, 2.0))
    ax[1,0].imshow(T,   origin='lower', extent=ext, cmap='viridis')
    ax[1,1].imshow(T_dev, origin='lower', extent=ext, cmap='bwr', clim=(-2.0, 2.0))
    ax[2,0].imshow(uz, origin='lower', extent=ext, cmap='bwr', clim=(-0.5, 0.5))
    ax[2,1].imshow(uh, origin='lower', extent=ext, cmap='bwr', clim=(-0.5, 0.5))

    for i in range(3):
        for j in range(2):
            ax[i,j].set_xlabel('X')
            ax[i,j].set_ylabel('Y')

    ax[0,0].set_title(r'$\rho$')
    ax[1,0].set_title(r'$T$')
    ax[0,1].set_title(r'$\rho - \bar{\rho}$')
    ax[1,1].set_title(r'$T-\bar{T}$')
    ax[2,0].set_title(r'$u_z$')
    ax[2,1].set_title(r'$u_h$')
    plt.tight_layout()
    plt.savefig('render/img.{:04}.png'.format(sid))
    plt.close()

def extract_quantities_2D(d, gamma=5.0/3.0):
    T  = d.time
    dV = d.dx1[0] * d.dx2[0]
    mass = (dV * d.rho).sum()

    Ek = 0.5 * d.rho * (d.vx1**2.0 + d.vx2**2.0)
    e  = d.prs / (d.rho * (gamma-1.0))
    E  = Ek + d.rho * e

    Ek = Ek.sum()
    e  = e.sum()
    E  = E.sum()

    return T, mass, Ek, e, E

create_directory('render')
snap_count = count_pluto_snaps()

if not '--no-render' in sys.argv:
    p = Pool()
    p.map(plot_2D, range(snap_count))

T    = []
mass = []
Ek   = []
e    = []
E    = []

for sid in range(snap_count):
    d = pp.pload(sid)
    ndims = len(d.rho.shape)

    if ndims == 2:
        T_, mass_, Ek_, e_, E_ = extract_quantities_2D(d)
        T.append(T_)
        mass.append(mass_)
        Ek.append(Ek_)
        e.append(e_)
        E.append(E_)

print('Plotting')

fig, ax = plt.subplots(2, 2, figsize=(15, 15))
ax[0,0].plot(T, mass, '-k')
ax[0,0].axhline(mass[0], linestyle='--')
ax[0,0].set_xlabel('T')
ax[0,0].set_ylabel('Mass')

ax[0,1].plot(T, Ek, '-k')
ax[0,1].set_xlabel('T')
ax[0,1].set_ylabel('Kinetic energy')

ax[1,0].plot(T, e, '-k')
ax[1,0].set_xlabel('T')
ax[1,0].set_ylabel('Internal energy')

ax[1,1].plot(T, E, '-k')
ax[1,1].set_xlabel('T')
ax[1,1].set_ylabel('Total energy')

plt.show()

T_snaps = range(0, snap_count, 50)
T_bar = []
z = []
for sid in T_snaps:
    if sid == 0:
        d = pp.pload(sid)
        z = d.x2
    t, Tbar = extract_Tbar(sid)
    plt.plot(z, Tbar, label='t={:.3f}'.format(t))
    plt.legend()
plt.show()


        

    

