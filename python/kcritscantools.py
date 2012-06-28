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
        filebase = re.split('_', files[i])[0]
        if not(filebases.count(filebase)):
            filebases.append(filebase)
    
    numbases = len(filebases)

    data =  numpy.ones(numbases, dtype = [('omega1', float), ('B', float),
                                          ('kmin', float), ('kmax', float),
                                          ('kpeak', float),
                                          ('peakgr', float)])

    for i in range(0,numbases):
        data[i]['omega1'] = float(re.split('B',
                                           re.split('w', filebases[i])[1])[0])
        data[i]['B'] = float(re.split('B', filebases[i])[1])

        kminname = filebases[i] + "_kmin.nc"
        kmaxname = filebases[i] + "_kmax.nc"
        kpeakname = filebases[i] + "_kpeak.nc"
        
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

def plot_quantities_const_omega(rule, omega1):
    data = grab_data(rule)
    idxs = numpy.equal(data[:]['omega1'], omega1*numpy.ones(data.size))
    
    fig = pyplot.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)

    ax1.plot(data[idxs]['B'], data[idxs]['kmin'], '.-', label=r"$k_{min}$")
    ax1.plot(data[idxs]['B'], data[idxs]['kmax'], '.-', label=r"$k_{max}$")
    ax1.plot(data[idxs]['B'], data[idxs]['kpeak'], '.-', label=r"$k_{peak}$")

    ax2.plot(data[idxs]['B'], data[idxs]['peakgr'], '.-', label="peak gr")
    ax1.loglog()
    ax2.set_xscale('log')
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    ax1.axvline(Bcrit(omega1), color='k')
    ax2.axvline(Bcrit(omega1), color='k')
    ax2.set_xlabel("B [gauss]")
    ax1.set_ylabel(r"$k_z$ [1/cm]")
    ax2.set_ylabel(r"Re{$\gamma$} [1/s]")


def grab_data_kscan(rule):
    files = glob.glob(rule)

    numfiles = len(files)

    data =  numpy.ones(numfiles, dtype = [('kz', float), ('B', float),
                                          ('gr', float)])

    for i in range(0,numfiles):
        data[i]['kz'] = float(re.split('B',
                                       re.split('k', files[i])[1])[0])
        data[i]['B'] = float(re.split('.nc',
                                      re.split('B', files[i])[1])[0])

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
    ax.set_ylim(bottom=-5)
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
