"""Routines for processing shear-layer stability runs with no magnetic field
but with different omega1s and omega2s.

There are two types of runs: kscan, in which a number of runs are made for
fixed shear layer width and fixed deltaomega=omega1-omega2, but with varying
k_z and omega2. The files are named with the pattern k123W456.nc, where the
number after k indicates the value of k_z, and the number after W indicates
omega2. These files are arranged in directories for each shear layer width
and value of deltaomega.

In the kcritscan runs, a minimization routine finds the most unstable k_z and
largest unstable k_z. These modes are saved to files that look like
w123W456_kmax.nc and w123W456_kpeak.nc, where the number after w indicates
the inner cylinder speed, and the number after W indicates the outer cylinder
speed.
"""

import scipy.io.netcdf as netcdf
import matplotlib.pyplot as pyplot
import numpy
import glob
import re
import os.path

def grab_data(rule):
    files = glob.glob(rule)

    filebases = []
    for i in range(0, len(files)):
        #Grab the part after the / of any path, and before the _k*.nc
        filebase = re.split('_', 
                            re.split('/', files[i])[-1]
                            )[0]
        #If we haven't seen this base before, keep it around.
        if not(filebases.count(filebase)):
            filebases.append(filebase)

    numbases = len(filebases)

    data =  numpy.ones(numbases, dtype = [('omega1', float), ('omega2', float),
                                          ('deltaomega', float),
                                          ('kmin', float), ('kmax', float),
                                          ('kpeak', float),
                                          ('peakgr', float)])

    path = ""
    for dir in (re.split('/', files[0])[0:-1]):
        path = path + dir + '/'

    for i in range(0,numbases):
        data[i]['omega1'] = float(re.split('W',
                                           re.split('w', filebases[i])[1])[0])
        data[i]['omega2'] = float(re.split('W', filebases[i])[1])
        data[i]['deltaomega'] = data[i]['omega1'] - data[i]['omega2']
        
        kminname = path + filebases[i] + "_kmin.nc"
        kmaxname = path + filebases[i] + "_kmax.nc"
        kpeakname = path + filebases[i] + "_kpeak.nc"
        
        if(os.path.exists(kminname)):
            ncfile = netcdf.netcdf_file(kminname, 'r')
            data[i]['kmin'] = ncfile.kz
            ncfile.close()
        else:
            data[i]['kmin'] = numpy.nan

        if(os.path.exists(kmaxname)):
            ncfile = netcdf.netcdf_file(kmaxname, 'r')
            data[i]['kmax'] = ncfile.kz
            ncfile.close()
        else:
            data[i]['kmax'] = numpy.nan

        if(os.path.exists(kpeakname)):
            ncfile = netcdf.netcdf_file(kpeakname, 'r')
            data[i]['kpeak'] = ncfile.kz
            data[i]['peakgr'] = ncfile.variables['lambda'][0,0]
            ncfile.close()
        else:
            data[i]['kpeak'] = numpy.nan
            data[i]['peakgr'] = 0

    #Sort the array so things are nicer
    data = numpy.sort(data, order=['deltaomega', 'omega2'])
        
    return data


def plot_quantities_const_deltaomega(rule, *deltaomegas):
    if len(deltaomegas) <= 0:
        print "Need to specify at least one omega."
        return

    data = grab_data(rule)
    idxs = []
    axs = []

    #We want to be sure that we plot a growth rate for all omega2s,
    #even if the growthrate there was 0 so that there's no output
    #file. So, first find all the omega2s that we use.

    omega2s = []
    for point in data:
        if(omega2s.count(point['omega2']) == 0):
            omega2s.append(point['omega2'])
    omega2s.sort()
    
    numdeltaomegas = len(deltaomegas)
    fig = pyplot.figure(figsize=(8, 4*numdeltaomegas + 2))
    fig.subplots_adjust(hspace=0.5)

    for i in range(0, len(deltaomegas)):
        tempax = []
        idxs.append(numpy.equal(data[:]['deltaomega'],
                                deltaomegas[i]*numpy.ones(data.size)))
        if(i == 0):
            tempax.append(fig.add_subplot(2*numdeltaomegas, 1, i*2 + 1))
            tempax.append(fig.add_subplot(2*numdeltaomegas, 1, i*2 + 2,
                                          sharex = tempax[0]))
        else:
            tempax.append(fig.add_subplot(2*numdeltaomegas, 1, i*2 + 1,
                                          sharex = axs[0][0]))
            tempax.append(fig.add_subplot(2*numdeltaomegas, 1, i*2 + 2,
                                          sharex = axs[0][0]))
                
    
        axs.append(tempax)

    for i in range(0, len(deltaomegas)):
        axs[i][0].plot(data[idxs[i]]['omega2'], data[idxs[i]]['kmax'],
                       'r.-', label=r"$k_{max}$")
        axs[i][0].plot(data[idxs[i]]['omega2'], data[idxs[i]]['kpeak'],
                       'g.-', label=r"$k_{peak}$")
        #Special treatment from the growth rate, where we want to make
        #sure that we plot all zeros.
        zerogromega2s = []
        for omega2 in omega2s:
            if (list(data[idxs[i]]['omega2']).count(omega2) == 0):
                zerogromega2s.append(omega2)
        axs[i][1].plot(numpy.concatenate((data[idxs[i]]['omega2'],
                                          numpy.array(zerogromega2s))),
                       (numpy.concatenate((data[idxs[i]]['peakgr'],
                                           numpy.zeros(len(zerogromega2s)))) /
                        (2*numpy.pi*deltaomegas[i]/60)),
                       'b.-', label=r"Peak Re{$\gamma$}")
        axs[i][0].loglog()
        axs[i][1].set_xscale('log')
        axs[i][0].axvline(deltaomegas[i]/2.35, color='k')
        axs[i][1].axvline(deltaomegas[i]/2.35, color='k')
        axs[i][0].text(0.1, 0.1,
                       r"$\Omega_1-\Omega_2$= %g rpm" % deltaomegas[i],
                       transform = axs[i][0].transAxes)
        axs[i][1].text(0.1, 0.1,
                       r"$\Omega_1-\Omega_2$= %g rpm" % deltaomegas[i],
                       transform = axs[i][1].transAxes)
        axs[i][0].set_ylabel(r"$k_z$ [1/cm]")
        axs[i][1].set_ylabel(r"Re{$\gamma$}/$\Delta\Omega$")
        
    axs[0][0].legend(loc='upper right')
    axs[0][1].legend(loc='upper right')
    
    axs[-1][1].set_xlabel(r"$\Omega_2$ [rpm]")


def plot_quantities_const_omega2(rule, *omega2s):
    if len(omega2s) <= 0:
        print "Need to specify at least one omega."
        return

    data = grab_data(rule)
    idxs = []
    axs = []
    
    numomega2s = len(omega2s)
    fig = pyplot.figure(figsize=(8, 4*numomega2s + 2))
    fig.subplots_adjust(hspace=0.5)

    for i in range(0, len(omega2s)):
        tempax = []
        idxs.append(numpy.equal(data[:]['omega2'],
                                omega2s[i]*numpy.ones(data.size)))
        if(i == 0):
            tempax.append(fig.add_subplot(2*numomega2s, 1, i*2 + 1))
            tempax.append(fig.add_subplot(2*numomega2s, 1, i*2 + 2,
                                          sharex = tempax[0]))
        else:
            tempax.append(fig.add_subplot(2*numomega2s, 1, i*2 + 1,
                                          sharex = axs[0][0]))
            tempax.append(fig.add_subplot(2*numomega2s, 1, i*2 + 2,
                                          sharex = axs[0][0]))
                
    
        axs.append(tempax)

    for i in range(0, len(omega2s)):
        axs[i][0].plot(data[idxs[i]]['deltaomega'], data[idxs[i]]['kmax'],
                       'r.-', label=r"$k_{max}$")
        axs[i][0].plot(data[idxs[i]]['deltaomega'], data[idxs[i]]['kpeak'],
                       'g.-', label=r"$k_{peak}$")
        axs[i][1].plot(data[idxs[i]]['deltaomega'],
                       data[idxs[i]]['peakgr']/(2*numpy.pi * 
                                                data[idxs[i]]['deltaomega'] /
                                                60),
                       'b.-', label=r"Peak Re{$\gamma$}")
        axs[i][0].loglog()
        axs[i][1].set_xscale('log')
        axs[i][0].text(0.1, 0.1, r"$\Omega_2$= %g rpm" % omega2s[i],
                       transform = axs[i][0].transAxes)
        axs[i][1].text(0.1, 0.1, r"$\Omega_2$= %g rpm" % omega2s[i],
                       transform = axs[i][1].transAxes)
        axs[i][0].set_ylabel(r"$k_z$ [1/cm]")
        axs[i][1].set_ylabel(r"Re{$\gamma$}/$\Delta\Omega$")
        
    axs[0][0].legend(loc='upper right')
    axs[0][1].legend(loc='upper right')
    
    axs[-1][1].set_xlabel(r"$\Omega_1-\Omega_2$ [rpm]")


def grab_data_kscan(rule):
    files = glob.glob(rule)

    numfiles = len(files)

    data =  numpy.ones(numfiles, dtype = [('kz', float), ('omega2', float),
                                          ('gr', float)])

    for i in range(0,numfiles):        
        #Make sure there's not a path stuck on the front of this.
        fname = re.split('/', files[i])[-1]
        data[i]['kz'] = float(re.split('W',
                                       re.split('k', fname)[1])[0])
        data[i]['omega2'] = float(re.split('.nc',
                                      re.split('W', fname)[1])[0])

        ncfile = netcdf.netcdf_file(files[i], 'r')

        if(ncfile.arpack_modes_converged > 0):
            data[i]['gr'] = ncfile.variables['lambda'][0,0]
        else:
            data[i]['gr'] = numpy.nan

        ncfile.close()

    #Sort the array so things are nicer
    data = numpy.sort(data, order=['omega2', 'kz'])
    
    return data

def plot_kscan_curves(rule):
    data = grab_data_kscan(rule)
    
    omega2s = []

    #Find the unique omega2s.
    for i in range (0, data.size):
        if not(omega2s.count(data[i]['omega2'])):
            omega2s.append(data[i]['omega2'])
    
    #Now for each omega, find the indexes that correspond to it.
    omega2idx = []
    for i in range(0, len(omega2s)):
        omega2idx.append(numpy.equal(data[:]['omega2'],
                                     omega2s[i]*numpy.ones(data.size)))
    
    #Now plot each of these.
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    for i in range(0, len(omega2s)):
        ax.plot(data[omega2idx[i]]['kz'], data[omega2idx[i]]['gr'],
                '.-', label=r"$\Omega_2$=%g rpm" % omega2s[i])
    
    ax.set_xlabel(r"$k_z$ [1/cm]")
    ax.set_ylabel(r"Re[$\gamma$] [1/s]")
    ax.autoscale(axis='y', tight=True)
    ytop = ax.set_ylim()[1]
    if (ytop > 0):
        ax.set_ylim(bottom=-0.2*ytop)
    ax.set_xscale('log')
    ax.legend(loc='best')


def plot_gr_contour(rule, maxgr=0):
    data = grab_data(rule)

    #Find all of the unique Bs and deltaomegas
    omega2s = []
    deltaomegas = []

    for point in data:
        if (omega2s.count(point['omega2']) == 0):
            omega2s.append(point['omega2'])
        if (deltaomegas.count(point['deltaomega']) == 0):
            deltaomegas.append(point['deltaomega'])
    omega2s.sort()
    deltaomegas.sort()

    #Now make an array to hold the growthrates.
    #We'll have deltaomega along the x-axis, and omega2 along the y-axis,
    #keeping in mind that this gets transposed in the call to contourf.
    grs = numpy.zeros((len(omega2s), len(deltaomegas)))
    normgrs = numpy.zeros((len(omega2s), len(deltaomegas)))
    
    #Now for each datapoint, stuff the data into the appropriate array.
    for point in data:
        i = deltaomegas.index(point['deltaomega'])
        j = omega2s.index(point['omega2'])        
        grs[j,i] = point['peakgr']
        normgrs[j,i] = point['peakgr']/(point['deltaomega']*2*numpy.pi/60)
    
    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)

    if (maxgr != 0):
        normgrs = normgrs.clip(-maxgr, maxgr)

    cf = ax.contourf(deltaomegas, omega2s, normgrs, 50)
    ax.scatter(data[:]['deltaomega'], data[:]['omega2'], s=1, c='k')
    cb = fig.colorbar(cf, ax=ax)
    cb.set_label(r"Re{$\gamma$}/($\Omega_1-\Omega_2$)", rotation=270)
    ax.set_xlabel(r"$\Omega_1-\Omega_2$ [rpm]")
    ax.set_ylabel(r"$\Omega_2$ [rpm]")


def elsasser(omega, B):
    rho = 6.36 #Density in g/cm^3
    eta = 2.57e3 #Resistivity in cm^2/s
    return B**2/(4*numpy.pi*eta*rho*(omega*2*numpy.pi/60))

def Bcrit(omega, elsasser=1):
    rho = 6.36 #Density in g/cm^3
    eta = 2.57e3 #Resistivity in cm^2/s
    return numpy.sqrt(4*numpy.pi*eta*rho*(omega*2*numpy.pi/60)*elsasser)

def plot_omega2crit_scaling():
    #Note that these measurements were all made by hand from the plots!
    deltaomegas = numpy.array((5, 10, 50, 100, 200, 400, 600, 800))
    omega_75 = numpy.array((2.9, 6.0, 32, 64, 128, 264, 396, 528))
    omega_50 = numpy.array((3.3, 7.1, 34, 73, 149, 299, 450, 600))
    omega_25 = numpy.array((3.6, 7.6, 36, 79, 154, 310, 467, 620))
    omega_10 = numpy.array((3.8, 7.8, 38, 81, 160, 315, 474, 630))

    rossby235 = deltaomegas/2.35

    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    ax.loglog()
    ax.plot(deltaomegas, rossby235, 'k-',
            label=r"$(\Omega_1-\Omega_2)/\Omega=2.35$", lw=2)
    ax.plot(deltaomegas, omega_75, 'bs-',
            label=r"Re{$\gamma$}$=0.75(\Omega_1-\Omega_2)$")
    ax.plot(deltaomegas, omega_50, 'go-',
            label=r"Re{$\gamma$}$=0.50(\Omega_1-\Omega_2)$")
    ax.plot(deltaomegas, omega_25, 'r*-',
            label=r"Re{$\gamma$}$=0.25(\Omega_1-\Omega_2)$")
    ax.plot(deltaomegas, omega_10, 'cp-',
            label=r"Re{$\gamma$}$=0.10(\Omega_1-\Omega_2)$")
    
    ax.set_xlabel(r"$\Omega_1-\Omega_2$ [rpm]")
    ax.set_ylabel(r"$\Omega_2$ [rpm]")
    ax.legend(loc='upper left')
    ax.set_xlim(3,1000)
