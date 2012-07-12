"""Routines for processing shear-layer stability runs with omega2=0
but with different omega1s and applied fields B.

There are two types of runs: kscan, in which a number of runs are made
for fixed shear layer width and fixed omega1, but with varying k_z and
B. The files are named with the pattern k123B456.nc, where the number
after k indicates the value of k_z, and the number after B indicates
B. These files are arranged in directories for each shear layer width
and value of omega1.
      
In the kcritscan runs, a minimization routine finds the most unstable
k_z and largest unstable k_z. These modes are saved to files that look
like w123B456_kmax.nc and w123B456_kpeak.nc, where the number after w
indicates the inner cylinder speed, and the number after W indicates
the magnetic field.
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

    data =  numpy.ones(numbases, dtype = [('omega1', float), ('B', float),
                                          ('kmin', float), ('kmax', float),
                                          ('kpeak', float),
                                          ('peakgr', float)])

    path = ""
    for dir in (re.split('/', files[0])[0:-1]):
        path = path + dir + '/'

    for i in range(0,numbases):
        data[i]['omega1'] = float(re.split('B',
                                           re.split('w', filebases[i])[1])[0])
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
    data = numpy.sort(data, order=['omega1', 'B'])
        
    return data

def plot_quantities_const_omega1(rule, *omegas):
    if len(omegas) <= 0:
        print "Need to specify at least one omega."
        return

    data = grab_data(rule)
    idxs = []
    axs = []
    
    numomegas = len(omegas)
    fig = pyplot.figure(figsize=(8, 4*numomegas))
    fig.subplots_adjust(hspace=0.5)

    for i in range(0, len(omegas)):
        tempax = []
        idxs.append(numpy.equal(data[:]['omega1'],
                                omegas[i]*numpy.ones(data.size)))
        if(i == 0):
            tempax.append(fig.add_subplot(2*numomegas, 1, i*2 + 1))
            tempax.append(fig.add_subplot(2*numomegas, 1, i*2 + 2,
                                          sharex = tempax[0]))
        else:
            tempax.append(fig.add_subplot(2*numomegas, 1, i*2 + 1,
                                          sharex = axs[0][0]))
            tempax.append(fig.add_subplot(2*numomegas, 1, i*2 + 2,
                                          sharex = axs[0][0]))
                
    
        axs.append(tempax)

    for i in range(0, len(omegas)):
        axs[i][0].plot(data[idxs[i]]['B'], data[idxs[i]]['kmax'],
                       'r.-', label=r"$k_{max}$")
        axs[i][0].plot(data[idxs[i]]['B'], data[idxs[i]]['kpeak'],
                       'g.-', label=r"$k_{peak}$")
        axs[i][1].plot(data[idxs[i]]['B'], data[idxs[i]]['peakgr'],
                       'b.-', label="peak gr")
        axs[i][0].loglog()
        axs[i][1].set_xscale('log')
        axs[i][0].axvline(Bcrit(omegas[i]), color='k')
        axs[i][1].axvline(Bcrit(omegas[i]), color='k')
        axs[i][0].text(0.1, 0.1, r"$\Omega_{ic}$= %g rpm" % omegas[i],
                       transform = axs[i][0].transAxes)
        axs[i][1].text(0.1, 0.1, r"$\Omega_{ic}$= %g rpm" % omegas[i],
                       transform = axs[i][1].transAxes)
        axs[i][0].set_ylabel(r"$k_z$ [1/cm]")
        axs[i][1].set_ylabel(r"Re{$\gamma$} [1/s]")
        
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
        axs[i][0].plot(data[idxs[i]]['omega1'], data[idxs[i]]['kmax'],
                       'r.-', label=r"$k_{max}$")
        axs[i][0].plot(data[idxs[i]]['omega1'], data[idxs[i]]['kpeak'],
                       'g.-', label=r"$k_{peak}$")
        axs[i][1].plot(data[idxs[i]]['omega1'], data[idxs[i]]['peakgr'],
                       'b.-', label="peak gr")
        axs[i][0].loglog()
        axs[i][1].set_xscale('log')
        axs[i][0].text(0.1, 0.1, r"$B$= %g gauss" % Bs[i],
                       transform = axs[i][0].transAxes)
        axs[i][1].text(0.1, 0.1, r"$B$= %g gauss" % Bs[i],
                       transform = axs[i][1].transAxes)
        axs[i][0].set_ylabel(r"$k_z$ [1/cm]")
        axs[i][1].set_ylabel(r"Re{$\gamma$} [1/s]")
        
    axs[0][0].legend(loc='upper right')
    axs[0][1].legend(loc='upper right')

    axs[-1][1].set_xlabel(r"$\Omega_{ic}$ [rpm]")


def grab_data_kscan(rule):
    files = glob.glob(rule)

    numfiles = len(files)

    data =  numpy.ones(numfiles, dtype = [('kz', float), ('B', float),
                                          ('gr', float)])

    for i in range(0,numfiles):        
        #Make sure there's not a path stuck on the front of this.
        fname = re.split('/', files[i])[-1]
        data[i]['kz'] = float(re.split('B',
                                       re.split('k', fname)[1])[0])
        data[i]['B'] = float(re.split('.nc',
                                      re.split('B', fname)[1])[0])

        ncfile = netcdf.netcdf_file(files[i], 'r')
        data[i]['gr'] = ncfile.variables['lambda'][0,0]
        ncfile.close()

    #Sort the array so things are nicer
    data = numpy.sort(data, order=['B', 'kz'])
    
    return data

def plot_kscan_curves(rule):
    data = grab_data_kscan(rule)
    
    Bs = []

    #Find the unique Bs.
    for i in range (0, data.size):
        if not(Bs.count(data[i]['B'])):
            Bs.append(data[i]['B'])
    
    #Now for each B, find the indexes that correspond to it.
    Bidx = []
    for i in range(0, len(Bs)):
        Bidx.append(numpy.equal(data[:]['B'],
                                Bs[i]*numpy.ones(data.size)))
    
    #Now plot each of these.
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    for i in range(0, len(Bs)):
        ax.plot(data[Bidx[i]]['kz'], data[Bidx[i]]['gr'],
                '.-', label="B=%g" % Bs[i])
    
    ax.set_xlabel(r"$k_z$ [1/cm]")
    ax.set_ylabel(r"Re[$\gamma$] [1/s]")
    ax.autoscale(axis='y', tight=True)
    ytop = ax.set_ylim()[1]
    if (ytop > 0):
        ax.set_ylim(bottom=-0.2*ytop)
    ax.set_xscale('log')
    ax.legend(loc='best')


def plot_gr_contour(rule):
    data = grab_data(rule)

    #Find all of the unique Bs and deltaomegas
    Bs = []
    omega1s = []

    for point in data:
        if (Bs.count(point['B']) == 0):
            Bs.append(point['B'])
        if (omega1s.count(point['omega1']) == 0):
            omega1s.append(point['omega1'])
    Bs.sort()
    omega1s.sort()

    #Now make an array to hold the growthrates.
    #We'll have Omega1 along the x-axis, and B along the y-axis,
    #keeping in mind that this gets transposed in the call to contourf.
    grs = numpy.zeros((len(Bs), len(omega1s)))
    normgrs = numpy.zeros((len(Bs), len(omega1s)))
    
    #Now for each datapoint, stuff the data into the appropriate array.
    for point in data:
        i = omega1s.index(point['omega1'])
        j = Bs.index(point['B'])
        grs[j,i] = point['peakgr']
        normgrs[j,i] = point['peakgr']/(point['omega1']*2*numpy.pi/60)
    
    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)

    cf = ax.contourf(omega1s, Bs, normgrs, 50)
    ax.scatter(data[:]['omega1'], data[:]['B'], s=2, c='k')
    cb = fig.colorbar(cf, ax=ax)
    cb.set_label(r"Re{$\gamma$}/($\Omega_{ic}-\Omega_{oc}$)", rotation=270)
    ax.set_xlabel(r"$\Omega_{ic}$ [rpm]")
    ax.set_ylabel("B [gauss]")


def elsasser(omega, B):
    rho = 6.36 #Density in g/cm^3
    eta = 2.57e3 #Resistivity in cm^2/s
    return B**2/(4*numpy.pi*eta*rho*(omega*2*numpy.pi/60))

def Bcrit(omega, elsasser=1):
    rho = 6.36 #Density in g/cm^3
    eta = 2.57e3 #Resistivity in cm^2/s
    return numpy.sqrt(4*numpy.pi*eta*rho*(omega*2*numpy.pi/60)*elsasser)

def plot_Bcrit_scaling():
    #Note that these measurements were all made by hand from the plots!
    omega1s = numpy.array((1, 5, 10, 50, 100, 200, 400, 600, 800))
    omega_75 = numpy.array((60.4, 169, 248, 610, 868,
                            1258, 1793, 2213, 2584))
    omega_50 = numpy.array((99.8, 240, 352, 805, 1156,
                            1641, 2330, 2861, 3342))
    omega_25 = numpy.array((131, 318, 455, 1073, 1566,
                            2237, 3211, 4032, 4642))
    omega_10 = numpy.array((169, 400, 630, 1580, 2426,
                            3500, 4975, 6304, 7400))
    elsasserone = numpy.array([Bcrit(omega) for omega in omega1s])

    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    ax.loglog()
    ax.plot(omega1s, elsasserone, 'k-',
            label=r"$\Lambda=1$", lw=2)
    ax.plot(omega1s, omega_75, 'bs-',
            label=r"Re{$\gamma$}$=0.75\Omega_{ic}$")
    ax.plot(omega1s, omega_50, 'go-',
            label=r"Re{$\gamma$}$=0.50\Omega_{ic}$")
    ax.plot(omega1s, omega_25, 'r*-',
            label=r"Re{$\gamma$}$=0.25\Omega_{ic}$")
    ax.plot(omega1s, omega_10, 'cp-',
            label=r"Re{$\gamma$}$=0.10\Omega_{ic}$")
    
    ax.set_xlabel(r"$\Omega_{ic}$ [rpm]")
    ax.set_ylabel(r"$B$ [gauss]")
    ax.legend(loc='upper left')

    ax.set_ylim(40, 10000)
