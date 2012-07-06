"""Routines for processing shear-layer stability runs with magnetic
field and a fixed omega2, and varying deltaomega=omega1-omega2.

In the kcritscan runs, a minimization routine finds the most unstable
k_z and largest unstable k_z for a fixed shear layer width and fixed
omega2. These modes are saved to files that look like
dw123B456_kmax.nc and w123W456_kpeak.nc, where the number after dw
indicates the differential speed between the inner and outer
cylinders, and the number after B indicates the magnetic field.
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

    data =  numpy.ones(numbases, dtype = [('B', float),
                                          ('deltaomega', float),
                                          ('kmin', float), ('kmax', float),
                                          ('kpeak', float),
                                          ('peakgr', float)])

    path = ""
    for dir in (re.split('/', files[0])[0:-1]):
        path = path + dir + '/'

    for i in range(0,numbases):
        data[i]['deltaomega'] = float(re.split('B',
                                           re.split('dw', filebases[i])[1])[0])
        data[i]['B'] = float(re.split('B', filebases[i])[1])

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
            data[i]['peakgr'] = numpy.nan

    #Sort the array so things are nicer
    data = numpy.sort(data, order=['deltaomega', 'B'])
        
    return data

def plot_quantities_const_deltaomega(rule, *deltaomegas):
    if len(omegas) <= 0:
        print "Need to specify at least one deltaomega."
        return

    data = grab_data(rule)
    idxs = []
    axs = []
    
    numdeltaomegas = len(deltaomegas)
    fig = pyplot.figure(figsize=(8, 4*numdeltaomegas))
    fig.subplots_adjust(hspace=0.5)

    for i in range(0, len(omegas)):
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
        axs[i][0].plot(data[idxs[i]]['B'], data[idxs[i]]['kmax'],
                       '.-', label=r"$k_{max}$")
        axs[i][0].plot(data[idxs[i]]['B'], data[idxs[i]]['kpeak'],
                       '.-', label=r"$k_{peak}$")
        axs[i][1].plot(data[idxs[i]]['B'], data[idxs[i]]['peakgr'],
                       '.-', label="peak gr")
        axs[i][0].loglog()
        axs[i][1].set_xscale('log')
        axs[i][0].axvline(Bcrit(deltaomegas[i]), color='k')
        axs[i][1].axvline(Bcrit(deltaomegas[i]), color='k')
        axs[i][0].text(0.1, 0.1, r"$\Delta\Omega$= %g rpm" % deltaomegas[i],
                       transform = axs[i][0].transAxes)
        axs[i][1].text(0.1, 0.1, r"$\Delta\Omega$= %g rpm" % deltaomegas[i],
                       transform = axs[i][1].transAxes)
        axs[i][0].set_ylabel(r"$k_z$ [1/cm]")
        axs[i][1].set_ylabel(r"Re{$\gamma$} [1/s]")
        axs[i][0].axhline(2*numpy.pi, color='k')
        
    axs[0][0].legend(loc='upper right')
    axs[0][1].legend(loc='upper right')

    axs[-1][1].set_xlabel("B [gauss]")


def plot_quantities_const_B(rule, *Bs):
    if len(Bs) <= 0:
        print "Need to specify at least one B."
        return

    data = grab_data(rule)
    idxs = []
    axs = []
    
    numBs = len(Bs)
    fig = pyplot.figure(figsize=(8, 4*numBs))
    fig.subplots_adjust(hspace=0.5)

    for i in range(0, len(Bs)):
        tempax = []
        idxs.append(numpy.equal(data[:]['B'],
                                Bs[i]*numpy.ones(data.size)))
        if(i == 0):
            tempax.append(fig.add_subplot(2*numBs, 1, i*2 + 1))
            tempax.append(fig.add_subplot(2*numBs, 1, i*2 + 2,
                                          sharex = tempax[0]))
        else:
            tempax.append(fig.add_subplot(2*numBs, 1, i*2 + 1,
                                          sharex = axs[0][0]))
            tempax.append(fig.add_subplot(2*numBs, 1, i*2 + 2,
                                          sharex = axs[0][0]))
                
    
        axs.append(tempax)

    for i in range(0, len(Bs)):
        axs[i][0].plot(data[idxs[i]]['deltaomega'], data[idxs[i]]['kmax'],
                       '.-', label=r"$k_{max}$")
        axs[i][0].plot(data[idxs[i]]['deltaomega'], data[idxs[i]]['kpeak'],
                       '.-', label=r"$k_{peak}$")
        axs[i][1].plot(data[idxs[i]]['deltaomega'], data[idxs[i]]['peakgr'],
                       '.-', label="peak gr")
        axs[i][0].loglog()
        axs[i][1].set_xscale('log')
        axs[i][0].text(0.1, 0.1, r"$B$= %g gauss" % Bs[i],
                       transform = axs[i][0].transAxes)
        axs[i][1].text(0.1, 0.1, r"$B$= %g gauss" % Bs[i],
                       transform = axs[i][1].transAxes)
        axs[i][0].set_ylabel(r"$k_z$ [1/cm]")
        axs[i][1].set_ylabel(r"Re{$\gamma$} [1/s]")
        axs[i][0].axhline(2*numpy.pi, color='k')
        
    axs[0][0].legend(loc='upper right')
    axs[0][1].legend(loc='upper right')

    axs[-1][1].set_xlabel(r"$\Delta\Omega$ [RPM]")

def plot_gr_contour(rule):
    data = grab_data(rule)

    #Find all of the unique Bs and deltaomegas
    Bs = []
    deltaomegas = []

    for point in data:
        if (Bs.count(point['B']) == 0):
            Bs.append(point['B'])
        if (deltaomegas.count(point['deltaomega']) == 0):
            deltaomegas.append(point['deltaomega'])
    
    #Make sure we have the 0,0 point, even there are no saved files.
    if (Bs.count(0) == 0):
        Bs.append(0)
    if (deltaomegas.count(0) == 0):
        deltaomegas.append(0)
      
    Bs.sort()
    deltaomegas.sort()

    #Now make an array to hold the growthrates.
    #We'll have deltaomega along the x-axis, and B along the y-axis,
    #keeping in mind that this gets transposed in the call to contourf.
    grs = numpy.zeros((len(Bs), len(deltaomegas)))
    normgrs = numpy.zeros((len(Bs), len(deltaomegas)))
    
    #Now for each datapoint, stuff the data into the appropriate array.
    for point in data:
        i = deltaomegas.index(point['deltaomega'])
        j = Bs.index(point['B'])        
        grs[j,i] = point['peakgr']
        normgrs[j,i] = point['peakgr']/(point['deltaomega']*2*numpy.pi/60)
    
    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)

    cf = ax.contourf(deltaomegas, Bs, normgrs, 50)
    ax.scatter(data[:]['deltaomega'], data[:]['B'], s=5, c='k')
    cb = fig.colorbar(cf, ax=ax)
    cb.set_label(r"$\gamma/\Delta\Omega$", rotation=270)
    ax.set_xlabel(r"$\Delta\Omega$ [RPM]")
    ax.set_ylabel("B [gauss]")

    


def elsasser(omega, B):
    rho = 6.36 #Density in g/cm^3
    eta = 2.43e3 #Resistivity in cm^2/s
    return B**2/(4*numpy.pi*eta*rho*(omega*2*numpy.pi/60))

def Bcrit(omega, elsasser=1):
    rho = 6.36 #Density in g/cm^3
    eta = 2.43e3 #Resistivity in cm^2/s
    return numpy.sqrt(4*numpy.pi*eta*rho*(omega*2*numpy.pi/60)*elsasser)
