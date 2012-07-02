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

def plot_quantities_const_omega(rule, *omegas):
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
                       '.-', label=r"$k_{max}$")
        axs[i][0].plot(data[idxs[i]]['B'], data[idxs[i]]['kpeak'],
                       '.-', label=r"$k_{peak}$")
        axs[i][1].plot(data[idxs[i]]['B'], data[idxs[i]]['peakgr'],
                       '.-', label="peak gr")
        axs[i][0].loglog()
        axs[i][1].set_xscale('log')
        axs[i][0].axvline(Bcrit(omegas[i]), color='k')
        axs[i][1].axvline(Bcrit(omegas[i]), color='k')
        axs[i][0].text(0.1, 0.1, r"$\Omega$= %g rpm" % omegas[i],
                       transform = axs[i][0].transAxes)
        axs[i][1].text(0.1, 0.1, r"$\Omega$= %g rpm" % omegas[i],
                       transform = axs[i][1].transAxes)
        axs[i][0].set_ylabel(r"$k_z$ [1/cm]")
        axs[i][1].set_ylabel(r"Re{$\gamma$} [1/s]")
        axs[i][0].axhline(2*numpy.pi, color='k')
        
    axs[0][0].legend(loc='upper right')
    axs[0][1].legend(loc='upper right')

    axs[-1][1].set_xlabel("B [gauss]")


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

def elsasser(omega, B):
    rho = 6.36 #Density in g/cm^3
    eta = 2.43e3 #Resistivity in cm^2/s
    return B**2/(4*numpy.pi*eta*rho*(omega*2*numpy.pi/60))

def Bcrit(omega, elsasser=1):
    rho = 6.36 #Density in g/cm^3
    eta = 2.43e3 #Resistivity in cm^2/s
    return numpy.sqrt(4*numpy.pi*eta*rho*(omega*2*numpy.pi/60)*elsasser)