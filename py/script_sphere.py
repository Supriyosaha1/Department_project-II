#!/usr/bin/env python3.10

# imports 
import argparse
import sys
from matplotlib import pyplot as plt
import numpy as np
import jphot as pp
import lya_utils as lya


# command-line arguments
parser = argparse.ArgumentParser()

parser.add_argument("-fres","--fileResults", help="Results file to read (default is ./photons_done.dat)",default="./photons_done.dat")

parser.add_argument("-fics","--fileICs", help="ICs file to read (default is ./ppic.dat)",default="./ppic.dat")

parser.add_argument("-temp","--temperature", type=float, help="temparature of the gas (default 1.e4 K)",default="1.e4")

parser.add_argument("-tau","--tau", type=float, help="line centre optical depth (default 1.e5)",default="1.e5")

parser.add_argument("-nbin","--nbin",type=int,help="number of bins for histogram (default is nbin = 100)",default="100")

#if len(sys.argv)==1:
#    parser.print_help()
#    sys.exit(1)
args = parser.parse_args()

fileRes = args.fileResults
fileICs = args.fileICs
nbin = args.nbin

p = pp.photonlist(fileICs,fileRes)

# Dijkstra 2006, formulae C-17 
# analytic solution from Dijkstra+06

temperature = args.temperature
tau_0       = args.tau

print("nu_0 [Hz] =",lya.nu0)
print("c [cm/s]  =",lya.clight)
print("T         =",temperature)
print("tau_0     =",tau_0)
print("nbin      =",nbin)

xlim = 100.
x = np.arange(-xlim,xlim,0.01)
    
vth_kms = 12.9 * np.sqrt(temperature / 1.e4) # km/s
vth_cgs = vth_kms * 1.e5 # cm/s
a       = 4.71e-4 / np.sqrt(temperature / 1.e4)
J_x     = x**2 * np.sqrt(np.pi) / (np.sqrt(24.)*a*tau_0) / (1+np.cosh(np.sqrt(2.*np.pi**3/27.)*np.abs(x**3)/a/tau_0)) / 2.


# normalisation of this is 2*PI
# in Dijkstra2014 norm = 4.*PI => factor 2

delta_nu = lya.nu0 * vth_cgs / lya.clight
xx = (p.nu - lya.nu0) / delta_nu

hh, xb = np.histogram(xx,bins=nbin,density=True,range=(-xlim,xlim))

plt.step(xb[:-1],hh/4/np.pi,where='post')
##plt.plot(xb[:-1],hh/4/np.pi)
##plt.bar(xb[:-1],hh/4/np.pi,width=0.99)

plt.plot(x,J_x,color='red',lw=1)#2)


plt.show()

