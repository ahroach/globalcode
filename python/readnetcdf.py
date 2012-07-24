import scipy.io.netcdf as netcdf
import numpy
from pylab import *
import math
import cmath
import scipy
import scipy.interpolate
import scipy.stats
import glob
import re

def list_attributes(filename):
    ncfile = netcdf.netcdf_file(filename, 'r')
    print 'eta =', ncfile.eta
    print 'nu = ', ncfile.nu
    print 'rho = ', ncfile.rho
    print 'r1 = ', ncfile.r1
    print 'r2 = ', ncfile.r2
    print 'height = ', ncfile.height
    print 'B0 = ', ncfile.B0
    print 'kz = ', ncfile.kz
    print 'm = ', ncfile.m
    print 'va = ', ncfile.va
    print 'Pm = ', ncfile.Pm
    print 'Re = ', ncfile.Re
    print 'Rm = ', ncfile.Rm
    print 'Ha = ', ncfile.Ha
    print 'nmode = ', ncfile.nmode
    print 'numcells = ', ncfile.numcells
    print 'Magnetic BC = ' + ncfile.magnetic_bc
    print 'Arpack tolerance = ', ncfile.arpack_tol
    print 'Real(sigma) = ', ncfile.arpack_sigma_r
    print 'Imag(sigma) = ', ncfile.arpack_sigma_i
    print '# Modes Requested = ', ncfile.arpack_modes_requested
    print '# Modes Converged = ', ncfile.arpack_modes_converged
    print '# Iterations Used = ', ncfile.arpack_itersused
    print 'Max # iterations  = ', ncfile.arpack_maxiters
    print ncfile.dimensions
    ncfile.close()


def list_grs(filename):
    ncfile = netcdf.netcdf_file(filename, 'r')
    for i in range(ncfile.dimensions['index']):
        print i, ":", ncfile.variables['lambda'][i][0], "+", ncfile.variables['lambda'][i][1],"*I", "resid= ", ncfile.variables['residual'][i]
    ncfile.close()
    del ncfile

def find_eigenvalue_number(filename, gr, omega, tol=2.0):
    ncfile = netcdf.netcdf_file(filename, 'r')
    indices = list()
    distances = list()
    for i in range(ncfile.dimensions['index']):
        tempgr = ncfile.variables['lambda'][i][0]
        tempomega = -ncfile.variables['lambda'][i][1]
        distance = sqrt((gr - tempgr)**2 +
                        (omega - tempomega)**2)
        if (sqrt(distance) < tol):
            indices.append(i)
            distances.append(distance)
    
        
    sorted_indices = argsort(distances)
    print sorted_indices

    for index in sorted_indices:
        print "Eigenvalue", indices[index], ":", ncfile.variables['lambda'][indices[index]][0], "+", ncfile.variables['lambda'][indices[index]][1], "*I", ": distance=", sqrt(distances[index])
    
    ncfile.close()
    del ncfile


def plot_grs(filename):
    ncfile = netcdf.netcdf_file(filename, 'r')

    subplot(2,1,1)
    plot(ncfile.variables['lambda'][:,0])
    ylabel("Real growth rate [1/sec]")
    title("Growth rates for " + filename)    

    subplot(2,1,2)
    plot(-ncfile.variables['lambda'][:,1])
    xlabel("Mode number")
    ylabel("Oscillation frequency [rad/sec]")
    ncfile.close()
    del ncfile

def plot_eigenvalues(filename, logx=0):
    ncfile = netcdf.netcdf_file(filename, 'r')

    ax = subplot(111)

    if(logx):
        ax.set_xscale('log')
        plot(-ncfile.variables['lambda'][:,0], -ncfile.variables['lambda'][:,1], 'o')    
        xlabel("-Growth rate [1/sec]")
        ylabel("Oscillation frequency [rad/sec]")

    else:
        plot(ncfile.variables['lambda'][:,0], -ncfile.variables['lambda'][:,1], 'o')    
        xlabel(r"Re[$\gamma$] [1/s]")
        ylabel(r"-Im[$\gamma$] [1/s]")
    grid(b=1)
    ncfile.close()


def plot_eigenvalues_residual(filename):
    ncfile = netcdf.netcdf_file(filename, 'r')

    scatter(ncfile.variables['lambda'][:,0], -ncfile.variables['lambda'][:,1], c=ncfile.variables['residual'], cmap=cm.spectral)
    colorbar()
    xlabel("Growth rate [1/sec]")
    ylabel("Oscillation frequency [rad/sec]")
    grid(b=1)
    ncfile.close()

def plot_eigenvalues_number(filename):
    ncfile = netcdf.netcdf_file(filename, 'r')

    index = arange(0, ncfile.dimensions['index'])
    scatter(ncfile.variables['lambda'][:,0], -ncfile.variables['lambda'][:,1], c=index, cmap=cm.spectral)
    colorbar()
    xlabel("Growth rate [1/sec]")
    ylabel("Oscillation frequency [rad/sec]")
    grid(b=1)
    ncfile.close()


def plot_mode_spectrogram(filename, component, min_mode, max_mode):
    ncfile = netcdf.netcdf_file(filename, 'r')

    subplot(2,1,1)
    plottable_array = ncfile.variables[component][min_mode:max_mode+1,:,0]
    transpose(plottable_array)
    contourf(ncfile.variables['r'][:],
             range(min_mode, max_mode + 1),
             plottable_array)
    xlabel("r [cm]")
    ylabel("Abs(" + component + ")")
    title("Mode structures of " + component + " for " + filename)

    subplot(2,1,2)
    plottable_array = ncfile.variables[component][min_mode:max_mode+1,:,1]
    transpose(plottable_array)
    contourf(ncfile.variables['r'][:],
             range(min_mode, max_mode + 1),
             plottable_array)
    xlabel("r [cm]")
    ylabel("Phase(" + component + ")")
    ncfile.close()
    del ncfile

def plot_single_component(filename, component, mode_number, norm=1.0):
    ncfile = netcdf.netcdf_file(filename, 'r')

    subplot(2,1,1)
    plot(ncfile.variables['r'][:], 
         norm*ncfile.variables[component][mode_number,:,0], '.-')
    xlabel("r [cm]")
    ylabel("Abs(" + component + ")")
    titlestring = "%s for mode %d, GR = %.6g + %.6g*i 1/s, resid= %.6g" % (component, mode_number,  ncfile.variables['lambda'][mode_number,0], ncfile.variables['lambda'][mode_number,1], ncfile.variables['residual'][mode_number])
    title(titlestring)

    subplot(2,1,2)
    plot(ncfile.variables['r'][:], ncfile.variables[component][mode_number,:,1], '.-')
    xlabel("r [cm]")
    ylabel("Phase(" + component + ")")
    ncfile.close()
    del ncfile


def plot_vzderiv(filename, mode_number):
    ncfile = netcdf.netcdf_file(filename, 'r')
    vzderiv_mag, vzderiv_phase = find_vz_deriv(filename, mode_number)

    subplot(2,1,1)
    plot(ncfile.variables['r'][:], vzderiv_mag, '.-')
    xlabel("r [cm]")
    ylabel("Abs(dvz/dz)")
    titlestring = "vzderiv for mode %d, GR = %.6g + %.6g*i 1/s, resid= %.6g" % (mode_number,  ncfile.variables['lambda'][mode_number,0], ncfile.variables['lambda'][mode_number,1], ncfile.variables['residual'][mode_number])
    title(titlestring)

    subplot(2,1,2)
    plot(ncfile.variables['r'][:], vzderiv_phase, '.-')
    xlabel("r [cm]")
    ylabel("Phase(dvz/dz)")
    ncfile.close()
    del ncfile

def plot_velocity_components(filename, mode_number):
    ncfile = netcdf.netcdf_file(filename, 'r')
    r = ncfile.variables['r'][:]
    vr_mag = ncfile.variables['vr'][mode_number, :, 0]
    vr_arg = ncfile.variables['vr'][mode_number, :, 1]

    vt_mag = ncfile.variables['vt'][mode_number, :, 0]
    vt_arg = ncfile.variables['vt'][mode_number, :, 1]

    subplot(2,1,1)
    plot(r, vr_mag, '-', label=r"$v_r$")
    plot(r, vt_mag, '-', label=r"$v_{\theta}$")
    legend()
    ylabel("Magnitude")

    subplot(2,1,2)
    plot(r, vr_arg, '-', label=r"Phase{$v_r$)")
    plot(r, vt_arg, '-', label=r"Phase{$v_{\theta}$)")
    xlabel("r [cm]")
    ylabel("Phase")
    nfile.close()
    del ncfile


def plot_all_components(filename, mode_number):
    ncfile = netcdf.netcdf_file(filename, 'r')
    r = ncfile.variables['r'][:]
    vr_mag = ncfile.variables['vr'][mode_number, :, 0]
    vr_arg = ncfile.variables['vr'][mode_number, :, 1]

    vt_mag = ncfile.variables['vt'][mode_number, :, 0]
    vt_arg = ncfile.variables['vt'][mode_number, :, 1]

    br_mag = ncfile.variables['br'][mode_number, :, 0]
    br_arg = ncfile.variables['br'][mode_number, :, 1]

    bt_mag = ncfile.variables['bt'][mode_number, :, 0]
    bt_arg = ncfile.variables['bt'][mode_number, :, 1]

    pi_mag = ncfile.variables['pi'][mode_number, :, 0]
    pi_arg = ncfile.variables['pi'][mode_number, :, 1]

    vr_real = zeros(vr_mag.size)
    vt_real = zeros(vr_mag.size)
    br_real = zeros(vr_mag.size)
    bt_real = zeros(vr_mag.size)
    pi_real = zeros(vr_mag.size)
    for i in range(0, vr_real.size):
        vr_real[i] = real(vr_mag[i]*exp(1j*vr_arg[i]))
        vt_real[i] = real(vt_mag[i]*exp(1j*vt_arg[i]))
        br_real[i] = real(br_mag[i]*exp(1j*br_arg[i]))
        bt_real[i] = real(bt_mag[i]*exp(1j*bt_arg[i]))
        pi_real[i] = real(pi_mag[i]*exp(1j*pi_arg[i]))

    subplot(5,1,1)
    plot(r, vr_real, '.-')
    ylabel(r"$v_r$")
    titlestring = "mode %d, GR = %.6g + %.6g*i 1/s" % (mode_number,  ncfile.variables['lambda'][mode_number,0], ncfile.variables['lambda'][mode_number,1])
    title(titlestring)

    subplot(5,1,2)
    plot(r, vt_real, '.-')
    ylabel(r"$v_\theta$")

    subplot(5,1,3)
    plot(r, br_real, '.-')
    ylabel(r"$b_r$")

    subplot(5,1,4)
    plot(r, bt_real, '.-')
    ylabel(r"$b_\theta$")


    subplot(5,1,5)
    plot(r, pi_real, '.-')
    ylabel(r"$\pi$")
    xlabel("r [cm]")

    ncfile.close()
    del ncfile


def plot_shear_dependence(rule):
    files = glob.glob(rule)
    numfiles = len(files)

    data = zeros(numfiles, dtype = [('width', float), ('m', int),
                                    ('growthrate', float)])
    for i in range(0, numfiles):
        ncfile = netcdf.netcdf_file(files[i], 'r')
        filename = re.split('/', files[i])[-1]
        #Now grab width out of the title string
        data[i]['width'] = float(re.split('d', re.split('m',
                                                        filename)[0])[1])
        #Now grab m out of the file name
        data[i]['m'] = int(re.split('m', re.split('.nc',
                                                  filename)[0])[1])
        data[i]['growthrate'] = ncfile.variables['lambda'][0,0]
        ncfile.close()

    #Sort the data array so that things plot correctly later.
    data = sort(data, order=['width', 'm'])

    #Find the list of unique ms
    ms = []
    for point in data:
        if(ms.count(point['m']) == 0):
            ms.append(point['m'])
    ms.sort()

    invwidth = 1/data[:]['width']

    fig = figure()
    ax = fig.add_subplot(111)
    #Now plot all of the azimuthal mode numbers as separate lines
    #Now plot all of the azimuthal mode numbers as separate lines
    for m in ms:
        indices = equal(data[:]['m'], m*ones(numfiles))
        ax.plot(invwidth[indices], data[indices]['growthrate'],
                '.-', label="m=%i" % m)

    ax.set_xscale('log')
    ax.set_ylabel("Growth rate [1/sec]")
    ax.set_xlabel("Steepness (1/layer width) [1/cm]")
    ax.legend(loc='upper left')
    ax.axhline(0, color='k')

def plot_width_dependence(rule):
    files = glob.glob(rule)
    numfiles = len(files)

    data = zeros(numfiles, dtype = [('width', float), ('m', int),
                                    ('growthrate', float)])
    for i in range(0, numfiles):
        ncfile = netcdf.netcdf_file(files[i], 'r')
        filename = re.split('/', files[i])[-1]
        #Now grab width out of the title string
        data[i]['width'] = float(re.split('d', re.split('m',
                                                        filename)[0])[1])
        #Now grab m out of the file name
        data[i]['m'] = int(re.split('m', re.split('.nc',
                                                  filename)[0])[1])
        data[i]['growthrate'] = ncfile.variables['lambda'][0,0]
        ncfile.close()

    #Sort the data array so that things plot correctly later.
    data = sort(data, order=['width', 'm'])


    #Find the list of unique ms
    ms = []
    for point in data:
        if(ms.count(point['m']) == 0):
            ms.append(point['m'])
    ms.sort()

    fig = figure()
    ax = fig.add_subplot(111)
    #Now plot all of the azimuthal mode numbers as separate lines
    for m in ms:
        indices = equal(data[:]['m'], m*ones(numfiles))
        ax.plot(data[indices]['width'], data[indices]['growthrate'],
                '.-', label="m=%i" % m)

    ax.set_xscale('log')
    ax.set_ylabel("Growth rate [1/sec]")
    ax.set_xlabel("Shear layer half-width [cm]")
    ax.legend(loc='best')
    ax.axhline(0, color='k')

def find_vz_deriv(filename, mode_number):
    ncfile = netcdf.netcdf_file(filename, 'r')
    vr_mag = ncfile.variables['vr'][mode_number,:,0]
    vr_phase = ncfile.variables['vr'][mode_number,:,1]
    vt_mag = ncfile.variables['vt'][mode_number,:,0]
    vt_phase = ncfile.variables['vt'][mode_number,:,1]

    vzderiv=zeros(vr_mag.size, dtype=complex)
    m = ncfile.m
    r = ncfile.variables['r'][:]

    #Now calculate from incompressibility condition:
    for i in range(0, r.size-1):
        vrcontrib = (1/r[i])*((r[i+1]*vr_mag[i+1]*e**(1j*vr_phase[i+1]) -
                               r[i]*vr_mag[i]*e**(1j*vr_phase[i]))/
                              (r[i+1]-r[i]))
        vtcontrib = 1j*m*vt_mag[i]*e**(1j*vt_phase[i])/r[i]
        vzderiv[i] = -vrcontrib - vtcontrib

    vzderiv_mag = abs(vzderiv)
    vzderiv_phase = angle(vzderiv)
    ncfile.close()
    del ncfile

    return vzderiv_mag, vzderiv_phase

def plot_velocity_contours(filename, mode_number):
    figure(figsize=(12,4), dpi=120, facecolor='w')
    axes([0,0,1,1]) #Use the entire canvas
    subplot(1,3,1)
    axis('off')
    axis('equal')
    plot_component_contour(filename, 'vt', mode_number, showcolorbar=1,
                           axisequalize=0)
    subplot(1,3,2)
    axis('off')
    axis('equal')
    plot_component_contour(filename, 'vr', mode_number, showcolorbar=1,
                           axisequalize=0)
    subplot(1,3,3)
    axis('off')
    axis('equal')
    vzderiv_mag, vzderiv_phase = find_vz_deriv(filename, mode_number)
    plot_derived_quantity_contour(filename, vzderiv_mag, vzderiv_phase,
                                  showcolorbar=1, axisequalize=0)

def plot_avg_torque(filename, mode_number):
    ncfile = netcdf.netcdf_file(filename, 'r')

    vr_mag = ncfile.variables['vr'][mode_number,:,0]
    vr_phase = ncfile.variables['vr'][mode_number,:,1]
    vt_mag = ncfile.variables['vt'][mode_number,:,0]
    vt_phase = ncfile.variables['vt'][mode_number,:,1]
    ncells = ncfile.numcells
    
    var_r = ncfile.variables['r'][:]
    m = ncfile.m
    ncfile.close()
    del ncfile

    #First let's renormalize vr and vt so that the expected value of
    #any given element is 1.
    #scalefac^2*(sum(vr^2) + sum(vt^2)) = 2*numcells
    x = sqrt((vr_mag*vr_mag).sum()
             + (vt_mag*vt_mag).sum())
    scalefac = sqrt(2*ncells)/x
    vr_mag = vr_mag*scalefac
    vt_mag = vt_mag*scalefac

    #Now calculate the average torque
    avg_torque = zeros(var_r.size)

    for i in range(0, var_r.size):
        avg_torque[i] = 0.5*(vr_mag[i]*vt_mag[i]*
                             (cos(vr_phase[i]-vt_phase[i])))
    
    ax1=gca()
    ax1.plot(var_r, vr_mag, 'r-', label=r"$|v_r|$")
    ax1.plot(var_r, vt_mag, 'g-', label=r"$|v_{\theta}|$")
    ax1.plot(var_r, avg_torque, 'b-',
             label=r"$<\tilde{v_{r}}\tilde{v_{\theta}}>$")

    ax1.set_ylabel(r"$|v_r|$, $|v_{\theta}|$, $<\tilde{v_{r}}\tilde{v_{\theta}}>$ [arb.]")

    ax2 = ax1.twinx()
    phasedifference = vt_phase-vr_phase
    for i in range(0,phasedifference.size):
        if phasedifference[i] > pi:
            phasedifference[i] = phasedifference[i] - 2.0*pi
        if phasedifference[i] < -pi:
            phasedifference[i] = phasedifference[i] + 2.0*pi
    
    ax2.plot(var_r, phasedifference, 'c-', label=r"$\Delta\phi$")
    ax2.set_ylabel(r"$\Delta\phi=$Arg($v_{\theta}$) - Arg($v_r$) [rad]", rotation=270)
    ax2.set_ylim(-pi,pi)
    ax2.axhline(pi/2.0, color='black')
    ax2.axhline(-pi/2.0, color='black')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    legend(lines1 + lines2, labels1 + labels2,
           loc='upper right', ncol=2)
    ax1.set_xlabel("r [cm]")


def plot_avg_torque_m0(filename, mode_number):
    ncfile = netcdf.netcdf_file(filename, 'r')

    vr_mag = ncfile.variables['vr'][mode_number,:,0]
    vr_phase = ncfile.variables['vr'][mode_number,:,1]
    vt_mag = ncfile.variables['vt'][mode_number,:,0]
    vt_phase = ncfile.variables['vt'][mode_number,:,1]
    ncells = ncfile.numcells
    
    var_r = ncfile.variables['r'][:]
    m = ncfile.m
    ncfile.close()
    del ncfile

    #First let's renormalize vr and vt so that the expected value of
    #any given element is 1.
    #scalefac^2*(sum(vr^2) + sum(vt^2)) = 2*numcells
    x = sqrt((vr_mag*vr_mag).sum()
             + (vt_mag*vt_mag).sum())
    scalefac = sqrt(2*ncells)/x
    vr_mag = vr_mag*scalefac
    vt_mag = vt_mag*scalefac

    #Now calculate the average torque

    vr = zeros(var_r.size)
    vt = zeros(var_r.size)
    avg_torque = zeros(var_r.size)

    for i in range(0, var_r.size):
        vr[i] = real(vr_mag[i]*exp(1j*vr_phase[i]))
        vt[i] = real(vt_mag[i]*exp(1j*vt_phase[i]))
        avg_torque[i] = 0.5*(vr[i]*vt[i])
    

    ax1=gca()
    ax1.plot(var_r, vr, 'r-', label=r"$v_r$")
    ax1.plot(var_r, vt, 'g-', label=r"$v_{\theta}$")
    ax1.plot(var_r, avg_torque, 'b-',
             label=r"$<\tilde{v_{r}}\tilde{v_{\theta}}>$")

    ax1.set_ylabel(r"$v_r$, $v_{\theta}$, $<\tilde{v_{r}}\tilde{v_{\theta}}>$ [arb.]")
    ax1.set_xlabel("r [cm]")
    ax1.legend(loc='upper right')


def plot_torque(filename, mode_number):
    ncfile = netcdf.netcdf_file(filename, 'r')
    vr_mag = ncfile.variables['vr'][mode_number,:,0]
    vr_phase = ncfile.variables['vr'][mode_number,:,1]
    vr_unwrapped_phase = numpy.unwrap(vr_phase)

    vt_mag = ncfile.variables['vt'][mode_number,:,0]
    vt_phase = ncfile.variables['vt'][mode_number,:,1]
    vt_unwrapped_phase = numpy.unwrap(vt_phase)

    var_r = ncfile.variables['r'][:]
    m = ncfile.m

    #You might be tempted to just multiply the complex vr and vt to get the
    #torque and then take the real part, but that would not be correct!
    #The torque T=Re{v_r}*Re{v_t}

    #Set up the arrays to store the x and y coordinates, and the values there
    desiredcells = 100
    xi = linspace(-ncfile.r2, ncfile.r2, desiredcells)
    yi = linspace(-ncfile.r2, ncfile.r2, desiredcells)
    torque = zeros((desiredcells, desiredcells))

    interp_vr_mag = scipy.interpolate.interp1d(var_r, vr_mag, kind='nearest')
    interp_vr_phase = scipy.interpolate.interp1d(var_r, vr_unwrapped_phase,
                                                 kind='nearest')
    interp_vt_mag = scipy.interpolate.interp1d(var_r, vt_mag, kind='nearest')
    interp_vt_phase = scipy.interpolate.interp1d(var_r, vt_unwrapped_phase,
                                                 kind='nearest')

    #Avoid lookups of these functions every iteration.
    lsqrt = math.sqrt
    latan2 = math.atan2

    i=0
    for x in xi:
        k=0
        for y in yi:
            r = lsqrt(x**2 + y**2)
            if ((r > ncfile.r1) and (r < ncfile.r2)):
                angle = latan2(y, x)
                torque[i][k] = 0.5*(interp_vr_mag(r)*interp_vt_mag(r)*
                              (cos(interp_vr_phase(r)-interp_vt_phase(r)) -
                               cos(2*m*angle + interp_vr_phase(r) +
                                   interp_vt_phase(r))))
            else:
                torque[i][k] = numpy.nan
                
            k = k + 1
        i = i + 1

    contourf(xi, yi, torque, 50)

    xlim(-ncfile.r2, ncfile.r2)
    ylim(-ncfile.r2, ncfile.r2)
    axes().set_aspect('equal')

    xlabel("x [cm]")
    ylabel("y [cm]")
    titlestring = "Torque for mode %d, GR = %.6g + %.6g*i 1/s" % (mode_number, ncfile.variables['lambda'][mode_number,0], ncfile.variables['lambda'][mode_number,1])
    title(titlestring)
    colorbar()

    ncfile.close()
    del ncfile

def calculate_stream(filename, mode_number):
    ncfile = netcdf.netcdf_file(filename, 'r')
    vt_mag = ncfile.variables['vt'][mode_number,:,0]
    vt_phase = ncfile.variables['vt'][mode_number,:,1]
    
    r = ncfile.variables['r'][:]
    
    vt = vt_mag*exp(1j*vt_phase)
    
    stream = zeros(r.size, dtype=complex)
    for i in range(1, stream.size):
        stream[i] = stream[i-1] - (r[i] - r[i-1])*vt[i]
    
    return abs(stream), angle(stream)

def calculate_potential(filename, mode_number, B=1000):
    ncfile = netcdf.netcdf_file(filename, 'r')
    vt_mag = ncfile.variables['vt'][mode_number,:,0]
    vt_phase = ncfile.variables['vt'][mode_number,:,1]

    var_r = ncfile.variables['r'][:]
    
    potential = zeros(var_r.size, dtype=complex)

    potential[0] = 0
    for i in range(1, potential.size):
        potential[i] = (potential[i-1] + 
                        ((var_r[i] - var_r[i-1]) *
                         vt_mag[i]*(cos(vt_phase[i]) +
                                    1j*sin(vt_phase[i]))))
    
    c = 2.998e10
    potential = potential*B/c
    potential_mag = abs(potential)
    potential_phase = angle(potential)
    ncfile.close()

    return potential_mag, potential_phase                         
    
def calculate_currents_kz0(filename, mode_number, B=1000):
    ncfile = netcdf.netcdf_file(filename, 'r')
    vt_mag = ncfile.variables['vt'][mode_number,:,0]
    vt_phase = ncfile.variables['vt'][mode_number,:,1]
    vr_mag = ncfile.variables['vr'][mode_number,:,0]
    vr_phase = ncfile.variables['vr'][mode_number,:,1]
    
    eta = ncfile.eta
    m = ncfile.m
    ncfile.close()

    c = 2.998e10

    r = ncfile.variables['r'][:]
    jr_mag = vt_mag*c*B/(4*pi*eta)
    jr_phase = vt_phase
    jt_mag = vr_mag*c*B/(4*pi*eta)
    jt_phase = vr_phase + pi

    jr = zeros(jr_mag.size, dtype=complex)
    jt = zeros(jr_mag.size, dtype=complex)
    for i in range(jr.size):
        jr[i] = jr_mag[i]*(cos(jr_phase[i]) + 1j*sin(jr_phase[i]))
        jt[i] = jt_mag[i]*(cos(jt_phase[i]) + 1j*sin(jt_phase[i]))

    delj_inplane = zeros(jr.size, dtype=complex)

    for i in range(1,delj_inplane.size):
        delj_inplane[i] = ((1/r[i])*jr[i] +
                           (jr[i] - jr[i-1])/(r[i]-r[i-1]) +
                           1j*m*jt[i]/r[i])
    
    delj_inplane_mag = abs(delj_inplane)
    delj_inplane_phase = angle(delj_inplane)

    return jr_mag, jr_phase, jt_mag, jt_phase, delj_inplane_mag, delj_inplane_phase


def plot_derived_quantity_contour(filename, quantity_mag, quantity_phase,
                                  levels=50, showcolorbar=0, axisequalize=1,
                                  lines=0, phase_offset=0, autorotate=0,
                                  mode_number=0):
    ncfile = netcdf.netcdf_file(filename, 'r')

    if(autorotate):
        #Rotate things so that the peak vr is at the right of the plot.
        ridx = abs(ncfile.variables['r'][:]-14).argmin()
        phase_offset = -ncfile.variables['vr'][mode_number,ridx,1]

    (xi,yi,zi)=regrid_component_two(ncfile.variables['r'][:],
                                    quantity_mag,
                                    quantity_phase + phase_offset,
                                    ncfile.m, ncfile.r1, ncfile.r2)
    print "Regrid complete"


    cf=gca().contourf(xi, yi, zi.T, levels)
    if(lines):
        ct=gca().contour(xi, yi, zi.T, lines, colors='k')

    gca().set_xlim(-ncfile.r2, ncfile.r2)
    gca().set_ylim(-ncfile.r2, ncfile.r2)
    if(axisequalize):
        gca().set_aspect('equal')

    gca().set_xlabel("x [cm]")
    gca().set_ylabel("y [cm]")
    if(showcolorbar):
        gca().get_figure().colorbar(mappable=cf, ax=gca())
    ncfile.close()
    del ncfile


    
def plot_component_contour(filename, component, mode_number, showcolorbar=0,
                           axisequalize=1, showtitle=1, desiredcells=400,
                           phase_offset=0.0, autorotate=0, levels=50,
                           showlabels=1):
    ncfile = netcdf.netcdf_file(filename, 'r')

    if(autorotate):
        #Rotate things so that the peak vr is at the right of the plot.
        ridx = abs(ncfile.variables['r'][:]-14).argmin()
        phase_offset = -ncfile.variables['vr'][mode_number,ridx,1]
    
    (xi,yi,zi)=regrid_component_two(ncfile.variables['r'][:],
                                    ncfile.variables[component][mode_number,:,0],
                                    ncfile.variables[component][mode_number,:,1] + phase_offset,
                                    ncfile.m, ncfile.r1, ncfile.r2,
                                    desiredcells=desiredcells)
    print "Regrid complete"


    cf=gca().contourf(xi, yi, zi.T, levels)

    gca().set_xlim(-ncfile.r2, ncfile.r2)
    gca().set_ylim(-ncfile.r2, ncfile.r2)
    if(axisequalize):
        gca().set_aspect('equal')

    if(showlabels):
        gca().set_xlabel("x [cm]")
        gca().set_ylabel("y [cm]")
    
    titlestring = "%s for mode %d, GR = %.6g + %.6g*i 1/s" % (component, mode_number, ncfile.variables['lambda'][mode_number,0], ncfile.variables['lambda'][mode_number,1])
    if(showtitle):
        title(titlestring)
    if(showcolorbar):
        gca().get_figure().colorbar(mappable=cf,  ax=gca())
    ncfile.close()
    del ncfile
    return cf

def plot_all_component_contours(filename, mode_number):
    subplot(2,3,1)
    plot_component_contour(filename, 'vr', mode_number, showcolorbar=1)
    subplot(2,3,2)
    plot_component_contour(filename, 'vt', mode_number, showcolorbar=1)
    subplot(2,3,3)
    plot_component_contour(filename, 'vz', mode_number, showcolorbar=1)
    subplot(2,3,4)
    plot_component_contour(filename, 'br', mode_number, showcolorbar=1)
    subplot(2,3,5)
    plot_component_contour(filename, 'bt', mode_number, showcolorbar=1)
    subplot(2,3,6)
    plot_component_contour(filename, 'bz', mode_number, showcolorbar=1)


def regrid_component(filename, component, mode_number):
    ncfile = netcdf.netcdf_file(filename, 'r')
    print "Hello!"

    #First, figure out how densely we should calculate values in theta

    m = ncfile.m

    #Now figure out how many x and y values we want to eventually return
    #desiredcells = int(2*ncfile.numcells*(ncfile.r2/(ncfile.r2-ncfile.r1)))
    desiredcells=100
    x = linspace(-ncfile.r2, ncfile.r2, desiredcells)
    y = linspace(-ncfile.r2, ncfile.r2, desiredcells)
    z = zeros([desiredcells, desiredcells])


    tck_mag = scipy.interpolate.splrep(ncfile.dimensions['r'],
                                       ncfile.variables[component][mode_number,:,0],
                                       s=0)
    temp_phase = unwrap(ncfile.variables[component][mode_number,:,1])
    tck_phase = scipy.interpolate.splrep(ncfile.dimensions['r'], temp_phase, s=0)
    
    
    #Iterate over all the rs and thetas
    rs = ncfile.variables['r'][:]
    for i in range(0, x.size):
        print "i = " + str(i)
        for j in range(0, y.size):
            r = sqrt(x[i]^2 + y[j]^2)
            if (r < ncfile.r1) | (r > ncfile.r2):
                z[i,j] = numpy.nan
            else:
                mag = scipy.interpolate.splev(r, tck_mag, der=0)
                phase = scipy.interpolate.splev(r, tck_phase, der=0)
                theta = arctan(y[j]/x[i])
                value = mag*e**(1j*(phase+m*theta))
                z[i,j] = value.real

    ncfile.close()
    del ncfile
    return (x, y, z)



def regrid_component_two(var_r, var_mag, var_phase, m, r1, r2,
                         desiredcells=400):

    unwrapped_phase = numpy.unwrap(var_phase)

    #Set up the arrays to store the x and y coordinates, and the values there
    xi = linspace(-r2, r2, desiredcells)
    yi = linspace(-r2, r2, desiredcells)
    zi = zeros((desiredcells, desiredcells))

    interp_mag = scipy.interpolate.interp1d(var_r, var_mag, kind='nearest')
    interp_phase = scipy.interpolate.interp1d(var_r, unwrapped_phase,
                                              kind='nearest')

    #Avoid lookups of these functions every iteration.
    lsqrt = math.sqrt
    latan2 = math.atan2

    i=0
    for x in xi:
        k=0
        for y in yi:
            r = lsqrt(x**2 + y**2)
            if ((r > r1) and (r < r2)):
                angle = latan2(y, x)
                total_mag = interp_mag(r)*e**(1j*(interp_phase(r) + m*angle))
                zi[i][k] = total_mag.real
            else:
                zi[i][k] = numpy.nan
                
            k = k + 1
        i = i + 1

    return (xi, yi, zi)


    
def plot_k_dependence(filename):
    ncfile = netcdf.netcdf_file(filename, 'r')
    n = size(ncfile.variables['lambda'][:,0])
    modes = range(n)
    k = zeros(n)
    r = zeros(n)
    omega_at_peak = zeros(n)
    if((ncfile.m == 0) or (ncfile.variables['omega'][-1] == ncfile.variables['omega'][0])):
        for mode in modes:
            mode_attrs = get_mode_attributes_mequalzero(filename, mode)
            k[mode] = mode_attrs['k']
            r[mode] = mode_attrs['peak_radius']
            omega_at_peak[mode] = ncfile.variables['omega'][0]
    else:
        #Expect the mode to be sheared
        for mode in modes:
            mode_attrs = get_mode_attributes_mnotzero(filename, mode)
            k[mode] = mode_attrs['k']
            r[mode] = mode_attrs['peak_radius']
            omega_at_peak[mode] = mode_attrs['omega_at_peak']

    subplot(2,1,1)
    scatter(k, ncfile.variables['lambda'][:,0], c=modes, cmap=cm.spectral)
    xlabel("k_r [1/cm]")
    ylabel("Growth rate [1/sec]")

    undoppler_frequency = -ncfile.variables['lambda'][:,1] - ncfile.m*omega_at_peak
    
    subplot(2,1,2)
    scatter(k, undoppler_frequency, c=modes, cmap=cm.spectral)
    xlabel("k_r [1/cm]")
    ylabel("Oscillation frequency [rad/sec] - m*Omega")

    figure(3)
    scatter(r, -ncfile.variables['lambda'][:,1], c=modes, cmap=cm.spectral)
    xlabel("peak radius [cm]")
    ylabel("Oscillation frequency [rad/sec]")

    figure(4)
    scatter(r, k, c=modes, cmap=cm.spectral)
    xlabel("peak radius [cm]")
    ylabel("k_r [1/cm]")

    ncfile.close()
    del ncfile


def plot_rotation_profile(filename):
    ncfile = netcdf.netcdf_file(filename, 'r')
    plot(ncfile.variables['r'][:], ncfile.variables['omega'][:])
    xlabel("r [cm]")
    ylabel("Omega [rads/sec]")
    title("Background rotation profile for " + filename)
    ncfile.close()
    del ncfile

def plot_mode_power_spectrum(filename, mode):
    ncfile = netcdf.netcdf_file(filename, 'r')
    k, fourier = get_mode_fft(filename, mode)
    power = fourier*fourier.conjugate()/(2.0*pi)
    titlestring = "Mode %d, GR = %.6g + %.6g*i 1/s" % (mode, ncfile.variables['lambda'][mode,0], ncfile.variables['lambda'][mode,1])
    loglog(k, power, '.-')
    title(titlestring)
    xlabel("k_r [1/cm]")
    ylabel("Power spectrum of vr")
    axes = axis()
    newaxes = [0.0, axes[1], axes[2], axes[3]]
    axis(newaxes)
    grid(b=1)
    
    ncfile.close()
    del ncfile

def get_mode_fft(filename, mode):
    ncfile = netcdf.netcdf_file(filename, 'r')
    
    vr_mag = ncfile.variables['vr'][mode,:,0]
    vr_arg = ncfile.variables['vr'][mode,:,1]

    #Calculate the real part along a chord
    vr_real = zeros(vr_mag.size)

    for i in range(0, vr_real.size):
        vr_real[i] = vr_mag[i]*exp(1j*vr_arg[i])

    #Take off the DC part
    vr_real = vr_real - vr_real.mean()

    fourier = fft(vr_real)
    spatialstep = (ncfile.r2 - ncfile.r1)/ncfile.numcells
    n = ncfile.numcells
    k = fftfreq(n, d=spatialstep)

    fourier = fftshift(fourier)
    k = fftshift(k)*2.0*pi
    ncfile.close()
    del ncfile
    return k, fourier

def get_mode_attributes_mequalzero(filename, mode):
    ncfile = netcdf.netcdf_file(filename, 'r')

    growthrate = ncfile.variables['lambda'][mode,0]
    frequency = -ncfile.variables['lambda'][mode,1]

    ncfile.close()
    del ncfile

    #Not finding a peak here
    peak_radius = numpy.nan
    peak_width = numpy.nan
    local_omega = numpy.nan

    #Get the mode power spectrum
    #k, fourier = get_mode_fft(filename, mode)
    #power = fourier*fourier.conjugate()/(2.0*pi)

    #Clear out the negative part of the spectrum
    #numpoints = power.size

    #trimmedpower = power[numpoints/2:]
    #trimmedk = k[numpoints/2:]

    #Find the peak in the power spectrum
    #max_element = argmax(trimmedpower)
    #max_k = trimmedk[max_element]

    #k = max_k
    k = find_kr(filename, mode, 'vr', localized=0)
    
    data = {'growthrate': growthrate, 'frequency': frequency, 'peak_radius': peak_radius, 'peak_width': peak_width, 'k': k, 'omega_at_peak': local_omega}
    return data

    
def get_mode_attributes_mnotzero(filename, mode):
    ncfile = netcdf.netcdf_file(filename, 'r')

    #First get growthrate & frequency
    growthrate = ncfile.variables['lambda'][mode,0]
    frequency = -ncfile.variables['lambda'][mode,1]

    #Next find approximate location of maximum vr amplitude
    r = ncfile.variables['r'][:]
    vr_mag = ncfile.variables['vr'][mode,:,0]
    vr_arg = ncfile.variables['vr'][mode,:,1]
    max_element = argmax(vr_mag)
    max_radius = r[max_element]

    #Find fluid rotation speed at radius of maximum vr
    local_omega = ncfile.variables['omega'][max_element]

    ncfile.close()
    del ncfile

    #Fit a Gaussian to the magnitude, giving the approximate
    #spatial-extent of the mode, and a better estimate of the maximum point

    gaussfunc = lambda p, x: p[0]*exp(-0.5*((x-p[1])/p[2])**2.0)
    errfunc = lambda p, x, y: gaussfunc(p, x) - y
    
    p_init = [vr_mag[max_element], max_radius, 1.0]
    out = scipy.optimize.leastsq(errfunc, p_init, args = (r, vr_mag), epsfcn=0.0005, xtol=0.005, full_output=1)

    p_fit = out[0]

    peak_radius = p_fit[1]
    peak_width = p_fit[2]

    #Now look at the phase change inside this region.

    #dr = r[max_element] - r[max_element - 1]
    
    #minroi = int(max_element -1.5*peak_width/dr)
    #maxroi = int(max_element + 1.5*peak_width/dr)

    #vr_freq = vr_arg[minroi:maxroi]
    #r_freq = r[minroi:maxroi]

    #Unwrap the phase.  Then fit a line to it to see what rads/cm we get
    #vr_freq = unwrap(vr_freq)

    #out = scipy.stats.linregress(r_freq, vr_freq)
 
    #k = out[0] #m from the linear regression fit => k in rads/cm
    #fit = out[0]*r_freq + out[1]

    k = find_kr(filename, mode, 'vr')
    
    data = {'growthrate': growthrate, 'frequency': frequency, 'peak_radius': peak_radius, 'peak_width': peak_width, 'k': k, 'omega_at_peak': local_omega}
    return data


def find_dispersion_relation_k(filename, kr):
    ncfile = netcdf.netcdf_file(filename, 'r')

    eta = ncfile.eta
    nu = ncfile.nu
    kz = ncfile.kz
    m = ncfile.m
    r1 = ncfile.r1
    r2 = ncfile.r2
    va = ncfile.va
    omega1 = ncfile.variables['omega'][0]
    omega2 = ncfile.variables['omega'][-1]

    ncfile.close()
    del ncfile

    k = sqrt(kr**2 + kz**2 + (m*2.0/(r2-r1))**2)

    a = (omega2*r2**2-omega1*r1**2)/(r2**2-r1**2)
    b = r1**2*r2**2*(omega1-omega2)/(r2**2-r1**2)
    
    zeta = 2-2/(1+(a/b)*(0.5*r1 + 0.5*r2)**2)
    if(omega1 == omega2):
        zeta=2.0
    omega = (omega1+omega2)/2

    p = zeros(5)

    p[0] = 1
    p[1] = 2*(eta + nu)*k**2
    p[2] = 2*eta*nu*k**4 + 2*(kz*va)**2 + ((nu + eta)*k**2)**2 + 2*zeta*(omega*kz/k)**2
    p[3] = 2*nu*eta*(eta + nu)*k**6 + 2*(nu + eta)*(k*kz*va)**2 + 4*zeta*eta*(omega*kz)**2
    p[4] = (nu*eta*k**4)**2 + 2*nu*eta*(kz*va*k**2)**2 + (kz*va)**4
    p[4] = p[4] + 2*zeta*(omega*eta*k*kz)**2 -2*(2-zeta)*(omega*kz**2*va/k)**2

    coeff = numpy.roots(p)
    return coeff

def find_dispersion_relation_range(filename, krmin, krmax):
    krmin=log10(krmin)
    krmax=log10(krmax)
    krs = logspace(krmin, krmax, 1000 )
    rootsarray = zeros([4,1000], complex)
    n = 0
    for kr in krs:
        roots = find_dispersion_relation_k(filename, kr)
        #Remove ridiculous roots
        for i in range(0,4):
            if roots[i].real < -1e8:
                roots[i] = numpy.nan

        rootsarray[:,n] = roots
        n = n + 1

    data = {'krs': krs, 'roots': rootsarray}
    return data
        
def plot_dispersion_relation_eigenvalues_range(filename, krmin, krmax):
    data = find_dispersion_relation_range(filename, krmin, krmax)

    plot(data['roots'][0,:].real, data['roots'][0,:].imag,)
    plot(data['roots'][1,:].real, data['roots'][1,:].imag,)
    plot(data['roots'][2,:].real, data['roots'][2,:].imag,)
    plot(data['roots'][3,:].real, data['roots'][3,:].imag,)

def plot_dispersion_relation_range(filename, krmin, krmax):
    data = find_dispersion_relation_range(filename, krmin, krmax)
    
    subplot(2,1,1)

    plot(data['krs'], data['roots'][0,:].real)
    plot(data['krs'], data['roots'][1,:].real)
    plot(data['krs'], data['roots'][2,:].real)
    plot(data['krs'], data['roots'][3,:].real)
    xlabel("kr [1/cm]")
    ylabel("Growth rate [1/sec]")
    axhline(y=0, color='black')

    subplot(2,1,2)

    plot(data['krs'], data['roots'][0,:].imag,)
    plot(data['krs'], data['roots'][1,:].imag,)
    plot(data['krs'], data['roots'][2,:].imag,)
    plot(data['krs'], data['roots'][3,:].imag,)
    xlabel("kr [1/cm]")
    ylabel("Real frequency [1/sec]")


def find_kr(filename, mode, component='vt', localized=1):
    #Use this with caution on the magnetic field components! Boundary
    #conditions might give them DC offsets, which would effect the averaging
    #below. B_r seems especially afflicted.
    
    ncfile = netcdf.netcdf_file(filename, 'r')

    
    r = ncfile.variables['r'][:]
    mag = ncfile.variables[component][mode,:,0]
    phase = ncfile.variables[component][mode,:,1]

    ncfile.close()
    del ncfile
    
    phase = unwrap(phase)

    if (localized == 1):
        #We will use a area-weighted, magnitude-weighted average of the change
        #in phase per unit radius squared, since k**2 is typically what's most
        #important for the damping
        #\int{dr 2\pi r A |\partial \phi / \partial r|^2} / \int{dr 2 \pi rA}
        #where A is the magnitude at a point i
        #Yielding eventually
        #\sum_{i}A_i r_i (\phi_{i+1}-\phi{i})^2/(r[i+1]-r[i]) /
        #\sum_{i}(r_{i+1}-r_{i}) r_i A_{i}

        num = 0.0
        denom = 0.0

        for i in range(0,r.size-1):
            num = num + mag[i]*r[i]*(phase[i+1]-phase[i])**2/(r[i+1]-r[i])
            denom = denom + mag[i]*r[i]*(r[i+1]-r[i])

        kr = sqrt(num/denom)
    
        #Average wavelength = 2\pi/kr
        l = 2*pi/kr

        #print "kr = " + str(kr) + ", \lambda = " + str(l)
    else:
        #If this isn't a localized mode, we should be able to pick out
        #k_r relatively easily by just looking at the peak in the power
        #spectrum. We exclude a couple of points near k = 0 to avoid
        #issues with DC offset or meaningless slow variation.
        k, fourier = get_mode_fft(filename, mode)
        power = fourier*fourier.conjugate()/(2.0*pi)
        midpoint = size(power)/2
        kr = k[power[midpoint+1:-1].argmax() + midpoint + 1]
    
    return kr

    
    
