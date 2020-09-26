# Sean Keenan 5th Year MPhys
# Nanophotnics Assesment 2

import numpy as np
import matplotlib.pyplot as mp

# define variables: wavelength, refractive index, angle of incidence, material thickness
freq = np.linspace(start=170e12, stop=600e12, num=100, endpoint=True)
c = 2.998e8
n = [3, 1]
d = [100e-9, 300e-9]
theta_rad = 0.0

# pre-assign lists
s = []
AD_abs = []
Keig = []
# loop to generate S matrix
for i in range(len(freq)):
    theta = [theta_rad]
    k = (2 * np.pi * freq[i]) / c

    # calculate transmission angles
    for index in range(len(n) - 1):
        theta.append(np.arcsin(n[index] * np.sin(theta[index]) / n[index + 1]))

    # calculate psi
    ncos = [n[i] * np.cos(theta[i]) for i in range(len(theta))]
    psi = [k * d[i] * ncos[i] for i in range(len(d))]

    # build matrices
    d_1 = np.array([[1, 1], [n[index] * np.cos(theta[index]), -n[index] * np.cos(theta[index])]])
    d_2 = np.array([[1, 1], [n[index + 1] * np.cos(theta[index + 1]), -n[index + 1] * np.cos(theta[index + 1])]])
    s_12 = (np.linalg.inv(d_1) @ d_2)
    s_21 = (np.linalg.inv(d_2) @ d_1)
    g_1 = np.array([[np.exp(1j * psi[0]), 0], [0, np.exp(-1j * psi[0])]])
    g_2 = np.array([[np.exp(1j * psi[1]), 0], [0, np.exp(-1j * psi[1])]])

    # calculate transfer matrix
    s.append(np.linalg.multi_dot([s_12, g_2, s_21, g_1]))
    AD = (s[i][0,0] + s[i][1,1])/2
    AD_abs.append(abs(AD))
    Keig.append((1/(d[0] + d[1])) * np.arccos(AD))

# set global fontsize
fsizes = {'axes.labelsize' : 12, 'axes.titlesize' : 14}
mp.rcParams.update(fsizes)

# lpot of bandgap
fig_1, ax_1 = mp.subplots()
ax_1.plot(freq*1e-12, AD_abs, label='$\\frac{|A+D|}{2}$')
ax_1.fill_between(x=[170, 600], y1=1.0, y2=2.0, facecolor='none', hatch='x', edgecolor='r', label='Forbidden Solutions')
ax_1.set_title('Bandgap of Si/SiO$_{2}$ Unit Cell')
ax_1.set(xlabel='Frequency (THz)', ylabel='$^\\frac{1}{2}|A+D|$')
ax_1.set_xlim(left=np.amin(freq*1e-12), right=np.amax(freq*1e-12))
ax_1.set_ylim(top=np.amax(AD_abs)*1.1 , bottom=0)
ax_1.legend(loc='lower left')
ax_1.grid(True)

# plot of dispersion relation
fig_2, ax_2 = mp.subplots()
ax_2.plot(Keig, freq*1e-12)
ax_2.set_title('Dispersion Relation for Si/SiO$_{2}$ Unit Cell')
ax_2.set(xlabel='K(f)', ylabel='f (THz)')
ax_2.set_xlim(left=np.min(Keig), right=np.amax(Keig))
ax_2.set_ylim(top=np.amax(freq*1e-12), bottom=np.amin(freq*1e-12))
ax_2.grid(True)
mp.show()

# save figs
fig_1.savefig(fname='/Users/Message/Documents/' + 'AD.pdf', format='pdf')
fig_2.savefig(fname='/Users/Message/Documents/' + 'Disp.pdf', format='pdf')
