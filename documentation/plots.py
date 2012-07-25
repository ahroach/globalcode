import sys
sys.path.append("/home/MRI/globalcode/python")
import readnetcdf
import matplotlib.pyplot as pyplot
import matplotlib.ticker
import numpy
import scipy.io.netcdf as netcdf
import localdispersionrelation

benchmarkdir = "/home/MRI/globalcode/documentation/benchmarks/"

def plot_growth_rate_convergence(filename="growth_rate_convergence.eps"):
    datapath = benchmarkdir + "goodman_mri_m0_arpack_resscan/"
    grs = []
    grs_full = []
    n = []
    n_full = []

    ncfile = netcdf.netcdf_file(datapath + "m0_100.nc", 'r')
    grs.append(ncfile.variables['lambda'][0][0])
    n.append(100)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_200.nc", 'r')
    grs.append(ncfile.variables['lambda'][0][0])
    n.append(200)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_400.nc", 'r')
    grs.append(ncfile.variables['lambda'][0][0])
    n.append(400)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_500.nc", 'r')
    grs.append(ncfile.variables['lambda'][0][0])
    n.append(500)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_600.nc", 'r')
    grs.append(ncfile.variables['lambda'][0][0])
    n.append(600)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_800.nc", 'r')
    grs.append(ncfile.variables['lambda'][0][0])
    n.append(800)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_1000.nc", 'r')
    grs.append(ncfile.variables['lambda'][0][0])
    n.append(1000)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_2000.nc", 'r')
    grs.append(ncfile.variables['lambda'][0][0])
    n.append(2000)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_4000.nc", 'r')
    grs.append(ncfile.variables['lambda'][0][0])
    n.append(4000)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_8000.nc", 'r')
    grs.append(ncfile.variables['lambda'][0][0])
    n.append(8000)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_16000.nc", 'r')
    grs.append(ncfile.variables['lambda'][0][0])
    n.append(16000)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_32000.nc", 'r')
    grs.append(ncfile.variables['lambda'][0][0])
    n.append(32000)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_100_full.nc", 'r')
    grs_full.append(ncfile.variables['lambda'][0][0])
    n_full.append(100)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_200_full.nc", 'r')
    grs_full.append(ncfile.variables['lambda'][0][0])
    n_full.append(200)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_400_full.nc", 'r')
    grs_full.append(ncfile.variables['lambda'][0][0])
    n_full.append(400)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_500_full.nc", 'r')
    grs_full.append(ncfile.variables['lambda'][0][0])
    n_full.append(500)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_600_full.nc", 'r')
    grs_full.append(ncfile.variables['lambda'][0][0])
    n_full.append(600)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_800_full.nc", 'r')
    grs_full.append(ncfile.variables['lambda'][0][0])
    n_full.append(800)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_1000_full.nc", 'r')
    grs_full.append(ncfile.variables['lambda'][0][0])
    n_full.append(1000)
    ncfile.close()

    ncfile = netcdf.netcdf_file(datapath + "m0_2000_full.nc", 'r')
    grs_full.append(ncfile.variables['lambda'][0][0])
    n_full.append(2000)
    ncfile.close()

    relerr = abs((numpy.array(grs[:-1])-grs[-1])/grs[-1])
    relerr_full = abs((numpy.array(grs_full)-grs[-1])/grs[-1]) 

    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    ax.loglog()
    ax.plot(n_full, relerr_full, "rs-", ms=9, label="Full solution")
    ax.plot(n[:-1], relerr, 'b*-', label="ARPACK")
    ax.legend(loc='upper right')
    ax.set_xlabel("N (number of grid cells)")
    ax.set_ylabel(r"$|\gamma_N - \gamma_{32000}|/\gamma_{32000}$")
    fig.savefig(filename, bbox_inches='tight', pad_inches=0.01)

def plot_timings(filename="timings.eps"):
    numcells = [50, 100, 200, 400, 800,
                1600, 3200, 6400, 9600]
    t_full = [0.560, 4.080, 34.50, 295.331, 3519.7, 
              37061, 281511, numpy.nan, numpy.nan]
    m_full = [102, 110, 161, 350, 1102,
              4106, 16114, numpy.nan, numpy.nan]
    t_10 = [0.117, 3.410, 3.327, 2.733, 5.357,
            8.967, 60, 156, 264]
    m_10 = [numpy.nan, 99, 100, 101, 106,
            115, 132, 166, 200]
    t_100 = [numpy.nan, 3.957, 4.567, 4.784, 8.886,
             15.621, 68.77, 155.2, 265]
    m_100 = [numpy.nan, 99, 107, 114, 125,
             150, 183, 318, 428]
    t_1000 = [numpy.nan, numpy.nan, numpy.nan, numpy.nan, 473.7,
              646.4, 987.8, 1628, 2426]
    m_1000 = [numpy.nan, numpy.nan, numpy.nan, numpy.nan, 590,
              674, 1071, 1856, 2641]

    fig = pyplot.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax1.plot(numcells, t_full, "b*-", label="Full solution")
    ax1.plot(numcells, t_10, "rd-", label="ARPACK, 10 eigenvalues")
    ax1.plot(numcells, t_100, "gs-", label="ARPACK, 100 eigenvalues")
    ax1.plot(numcells, t_1000, "kp-", label="ARPACK, 1000 eigenvalues")
    ax1.set_xlim(30,10000)
    ax1.loglog()
    ax1.set_ylabel("Run time [s]")

    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(numcells, m_full, "b*-", label="Full solution")
    ax2.plot(numcells, m_10, "rd-", label="ARPACK, 10 values")
    ax2.plot(numcells, m_100, "gs-", label="ARPACK, 100 values")
    ax2.plot(numcells, m_1000, "kp-", label="ARPACK, 1000 values")
    ax2.set_xlim(30,10000)
    ax2.loglog()
    ax2.set_ylabel("Memory used [MB]")
    ax2.set_xlabel("Number of grid cells")
    ax2.legend(loc='upper left')
    
    ax1.axhline(60, color='k')
    ax1.text(35, 60*1.15, "Trip to water cooler")

    ax1.axhline(1800, color='k')
    ax1.text(35, 1800*1.10, "Lunch break")

    ax1.axhline(28800, color='k')
    ax1.text(35, 28800*1.15, "Overnight")

    fig.savefig(filename, bbox_inches='tight', pad_inches=0.01)


def plot_goodman_mri_m0_eigenvalues(filename="goodman_mri_m0_eigenvalues.eps"):
    datapath = benchmarkdir + "goodman_mri/"
    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)

    readnetcdf.plot_eigenvalues(datapath + "m0.nc")

    ax.grid(b='off')
    ax.axvline(0, color='k')
    ax.axhline(0, color='k')

    fig.savefig(filename, bbox_inches='tight', pad_inches=0.01)


def plot_goodman_mri_m1_eigenvalues(filename="goodman_mri_m1_eigenvalues.eps"):
    datapath = benchmarkdir + "goodman_mri/"
    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)

    readnetcdf.plot_eigenvalues(datapath + "m1.nc")

    ax.grid(b='off')
    ax.axvline(0, color='k')
    ax.axhline(0, color='k')

    fig.savefig(filename, bbox_inches='tight', pad_inches=0.01)


def plot_goodman_mri_insulating_mode(filename="goodman_mri_insulating_mode.eps"):
    datapath = benchmarkdir + "goodman_mri/"

    ncfile = netcdf.netcdf_file(datapath + "m0.nc", 'r')

    r = ncfile.variables['r'][:]
    
    vr = numpy.real(ncfile.variables['vr'][0, :, 0] *
                    numpy.exp(1j*ncfile.variables['vr'][0, :, 1]))

    vt = numpy.real(ncfile.variables['vt'][0, :, 0] *
                    numpy.exp(1j*ncfile.variables['vt'][0, :, 1]))

    vz = numpy.real(ncfile.variables['vz'][0, :, 0] *
                    numpy.exp(1j*ncfile.variables['vz'][0, :, 1]))

    br = numpy.real(ncfile.variables['br'][0, :, 0] *
                    numpy.exp(1j*ncfile.variables['br'][0, :, 1]))

    bt = numpy.real(ncfile.variables['bt'][0, :, 0] *
                    numpy.exp(1j*ncfile.variables['bt'][0, :, 1]))

    bz = numpy.real(ncfile.variables['bz'][0, :, 0] *
                    numpy.exp(1j*ncfile.variables['bz'][0, :, 1]))

    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    
    ax.plot(r, vr/3.0, label=r"$\phi_r\times1/3$")
    ax.plot(r, vt, label=r"$\phi_\theta$")
    ax.plot(r, vz*0.07, label=r"$\phi_z\times0.07$")
    ax.plot(r, br, label=r"$\beta_r$")
    ax.plot(r, bt*5, label=r"$\beta_\theta\times5$")
    ax.plot(r, bz, label=r"$\beta_z$")
    ax.set_xlim(ncfile.r1, ncfile.r2)
    ax.set_xlabel(r"$r$ [cm]")
    ax.legend(loc='lower right', ncol=2)

    ncfile.close()
    fig.savefig(filename, bbox_inches='tight', pad_inches=0.01)


def plot_goodman_mri_m0_stable_modes(filename="goodman_mri_m0_stable_modes.eps"):
    fig = pyplot.figure(figsize=(8,10))
    fig.subplots_adjust(left=0.15)
    datapath = benchmarkdir + "goodman_mri/"
    readnetcdf.plot_all_components(datapath + "m0.nc", 20, showpoints=0)
    readnetcdf.plot_all_components(datapath + "m0.nc", 50, showpoints=0)
    readnetcdf.plot_all_components(datapath + "m0.nc", 100, showpoints=0)

    majorformatter = matplotlib.ticker.FormatStrFormatter("%0.1e")
    for ax in fig.axes:
        ax.yaxis.set_major_formatter(majorformatter)
        ax.set_title("")

    fig.savefig(filename, bbox_inches='tight', pad_inches=0.01)

def plot_goodman_mri_m1_stable_modes(filename="goodman_mri_m1_stable_modes.eps"):
    fig = pyplot.figure(figsize=(8,10))
    fig.subplots_adjust(left=0.15)
    datapath = benchmarkdir + "goodman_mri/"
    readnetcdf.plot_all_components(datapath + "m1.nc", 20, showpoints=0)
    readnetcdf.plot_all_components(datapath + "m1.nc", 50, showpoints=0)
    readnetcdf.plot_all_components(datapath + "m1.nc", 100, showpoints=0)

    majorformatter = matplotlib.ticker.FormatStrFormatter("%0.1e")
    for ax in fig.axes:
        ax.yaxis.set_major_formatter(majorformatter)
        ax.set_title("")

    fig.savefig(filename, bbox_inches='tight', pad_inches=0.01)

def plot_narrowgap_alfven_pm1(filename="narrowgap_alfven_pm1.eps"):
    kmin = 100
    kmax = 5000
    #Calculate ks from the global code
    datafile = benchmarkdir + "narrowgap_alfven/m1.nc"
    ncfile = netcdf.netcdf_file(datafile, 'r')
    modes = range(ncfile.variables['lambda'][:,0].size)
    k = numpy.zeros(ncfile.variables['lambda'][:,0].size)

    for mode in modes:
        mode_attrs = readnetcdf.get_mode_attributes_mequalzero(datafile, mode)
        k[mode] = mode_attrs['k']

    #Find ks from the local dispersion relation
    localk = numpy.linspace(kmin, kmax, 200)
    omegas = numpy.zeros([localk.size, 5], dtype=complex)
    findr = localdispersionrelation.find_roots
    omegabar = numpy.sqrt(ncfile.variables['omega'][0] *
                          ncfile.variables['omega'][-1])
    zetabar = (2*(ncfile.r2**2*ncfile.variables['omega'][-1] -
                  ncfile.r1**2*ncfile.variables['omega'][0]) /
               ((ncfile.r2**2 - ncfile.r1**2)*omegabar))
    for i in range(0, localk.size):
        omegas[i,:] = findr(B=ncfile.B0,
                            kr=localk[i],
                            kz=ncfile.kz,
                            m=ncfile.m,
                            nu=ncfile.nu,
                            eta=ncfile.eta,
                            rho=ncfile.rho,
                            r1=ncfile.r1,
                            r2=ncfile.r2,
                            Omega=omegabar,
                            zeta=2)

    omegas = localdispersionrelation.sort_roots(omegas)

    fig=pyplot.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2, sharex=ax1)
    
    ax1.plot(localk, omegas[:,0].imag, 'k-')
    ax1.plot(localk, omegas[:,1].imag, 'k-')
    ax2.plot(localk, omegas[:,0].real, 'k-')
    ax2.plot(localk, omegas[:,1].real, 'k-')
    ax2.plot(localk, omegas[:,2].real, 'k-')
    ax2.plot(localk, omegas[:,3].real, 'k-')
    ax2.plot(localk, omegas[:,4].real, 'k-')
    ax1.plot(k, ncfile.variables['lambda'][:,0], 'r.')
    ax2.plot(k, ncfile.variables['lambda'][:,1], 'r.')

    ax1.set_xlim(kmin, kmax)
    ax2.set_xlim(kmin, kmax)
    
    ax1.set_ylabel(r"Re{$\gamma$}")
    ax2.set_ylabel(r"Im{$\gamma$}")
    ax2.set_xlabel(r"$k_r$ [1/cm]")

    ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    ax2.ticklabel_format(style='sci', scilimits=(0,0), axis='y')

    ncfile.close()
         
    fig.savefig(filename, bbox_inches='tight', pad_inches=0.01)


def plot_narrowgap_alfven_modes(filename="narrowgap_alfven_modes.eps"):
    fig = pyplot.figure(figsize=(8,10))
    fig.subplots_adjust(left=0.15)
    datapath = benchmarkdir + "narrowgap_alfven/"
    readnetcdf.plot_all_components(datapath + "m1.nc", 0, showpoints=0)
    readnetcdf.plot_all_components(datapath + "m1.nc", 22, showpoints=0)
    readnetcdf.plot_all_components(datapath + "m1.nc", 50, showpoints=0)
    readnetcdf.plot_all_components(datapath + "m1.nc", 100, showpoints=0)

    majorformatter = matplotlib.ticker.FormatStrFormatter("%0.1e")
    for ax in fig.axes:
        ax.yaxis.set_major_formatter(majorformatter)
        ax.set_title("")

    fig.savefig(filename, bbox_inches='tight', pad_inches=0.01)


def plot_narrowgap_alfven_pmpoint1(filename="narrowgap_alfven_pmpoint1.eps"):
    kmin = 100
    kmax = 2000
    #Calculate ks from the global code
    datafile = benchmarkdir + "narrowgap_alfven_pm0.1/m1.nc"
    ncfile = netcdf.netcdf_file(datafile, 'r')
    modes = range(ncfile.variables['lambda'][:,0].size)
    k = numpy.zeros(ncfile.variables['lambda'][:,0].size)

    for mode in modes:
        mode_attrs = readnetcdf.get_mode_attributes_mequalzero(datafile, mode)
        k[mode] = mode_attrs['k']

    #Find ks from the local dispersion relation
    localk = numpy.linspace(kmin, kmax, 200)
    omegas = numpy.zeros([localk.size, 5], dtype=complex)
    findr = localdispersionrelation.find_roots
    omegabar = numpy.sqrt(ncfile.variables['omega'][0] *
                          ncfile.variables['omega'][-1])
    zetabar = (2*(ncfile.r2**2*ncfile.variables['omega'][-1] -
                  ncfile.r1**2*ncfile.variables['omega'][0]) /
               ((ncfile.r2**2 - ncfile.r1**2)*omegabar))
    for i in range(0, localk.size):
        omegas[i,:] = findr(B=ncfile.B0,
                            kr=localk[i],
                            kz=ncfile.kz,
                            m=ncfile.m,
                            nu=ncfile.nu,
                            eta=ncfile.eta,
                            rho=ncfile.rho,
                            r1=ncfile.r1,
                            r2=ncfile.r2,
                            Omega=omegabar,
                            zeta=2)

    omegas = localdispersionrelation.sort_roots(omegas)

    fig=pyplot.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2, sharex=ax1)
    
    ax1.plot(localk, omegas[:,0].imag, 'k-')
    ax1.plot(localk, omegas[:,1].imag, 'k-')
    ax1.plot(localk, omegas[:,2].imag, 'k-')
    ax1.plot(localk, omegas[:,3].imag, 'k-')
    ax2.plot(localk, omegas[:,0].real, 'k-')
    ax2.plot(localk, omegas[:,1].real, 'k-')
    ax2.plot(localk, omegas[:,2].real, 'k-')
    ax2.plot(localk, omegas[:,3].real, 'k-')
    ax1.plot(k, ncfile.variables['lambda'][:,0], 'r.')
    ax2.plot(k, ncfile.variables['lambda'][:,1], 'r.')

    ax1.set_xlim(kmin, kmax)
    ax1.set_ylim(-150000, 0)
    ax2.set_xlim(kmin, kmax)

    ax1.set_ylabel(r"Re{$\gamma$}")
    ax2.set_ylabel(r"Im{$\gamma$}")    
    ax2.set_xlabel(r"$k_r$ [1/cm]")

    ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    ax2.ticklabel_format(style='sci', scilimits=(0,0), axis='y')

    ncfile.close()
         
    fig.savefig(filename, bbox_inches='tight', pad_inches=0.01)


def plot_narrowgap_inertial(filename="narrowgap_inertial.eps"):
    kmin = 100
    kmax = 5000
    #Calculate ks from the global code
    datafile = benchmarkdir + "narrowgap_inertial_highk/m0.nc"
    ncfile = netcdf.netcdf_file(datafile, 'r')
    modes = range(ncfile.variables['lambda'][:,0].size)
    k = numpy.zeros(ncfile.variables['lambda'][:,0].size)

    for mode in modes:
        mode_attrs = readnetcdf.get_mode_attributes_mequalzero(datafile, mode)
        k[mode] = mode_attrs['k']

    #Find ks from the local dispersion relation
    localk = numpy.linspace(kmin, kmax, 200)
    omegas = numpy.zeros([localk.size, 5], dtype=complex)
    findr = localdispersionrelation.find_roots
    omegabar = numpy.sqrt(ncfile.variables['omega'][0] *
                          ncfile.variables['omega'][-1])
    zetabar = (2*(ncfile.r2**2*ncfile.variables['omega'][-1] -
                  ncfile.r1**2*ncfile.variables['omega'][0]) /
               ((ncfile.r2**2 - ncfile.r1**2)*omegabar))
    print ncfile.kz
    for i in range(0, localk.size):
        omegas[i,:] = findr(B=ncfile.B0,
                            kr=localk[i],
                            kz=ncfile.kz,
                            m=ncfile.m,
                            nu=ncfile.nu,
                            eta=ncfile.eta,
                            rho=ncfile.rho,
                            r1=ncfile.r1,
                            r2=ncfile.r2,
                            Omega=omegabar,
                            zeta=2)

    omegas = localdispersionrelation.sort_roots(omegas)

    fig=pyplot.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2, sharex=ax1)
    
    ax1.plot(localk, omegas[:,0].imag, 'k-')
    ax1.plot(localk, omegas[:,1].imag, 'k-')
    ax1.plot(localk, omegas[:,2].imag, 'k-')
    ax1.plot(localk, omegas[:,3].imag, 'k-')
    ax1.plot(localk, omegas[:,4].imag, 'k-')
    ax2.plot(localk, omegas[:,0].real, 'k-')
    ax2.plot(localk, omegas[:,1].real, 'k-')
    ax2.plot(localk, omegas[:,2].real, 'k-')
    ax2.plot(localk, omegas[:,3].real, 'k-')
    ax2.plot(localk, omegas[:,4].real, 'k-')
    ax1.plot(k, ncfile.variables['lambda'][:,0], 'r.')
    ax2.plot(k, ncfile.variables['lambda'][:,1], 'r.')

    ax1.set_xlim(kmin, kmax)
    ax2.set_xlim(kmin, kmax)

    ax1.set_ylabel(r"Re{$\gamma$}")
    ax2.set_ylabel(r"Im{$\gamma$}")    
    ax2.set_xlabel(r"$k_r$ [1/cm]")

    ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    ax2.ticklabel_format(style='sci', scilimits=(0,0), axis='y')

    ncfile.close()
         
    fig.savefig(filename, bbox_inches='tight', pad_inches=0.01)


def plot_narrowgap_inertial_modes(filename="narrowgap_inertial_modes.eps"):
    fig = pyplot.figure(figsize=(8,10))
    fig.subplots_adjust(left=0.15)
    datapath = benchmarkdir + "narrowgap_inertial_highk/"
    readnetcdf.plot_all_components(datapath + "m1.nc", 11, showpoints=0)
    readnetcdf.plot_all_components(datapath + "m1.nc", 12, showpoints=0)
    readnetcdf.plot_all_components(datapath + "m1.nc", 99, showpoints=0)
    readnetcdf.plot_all_components(datapath + "m1.nc", 100, showpoints=0)


    majorformatter = matplotlib.ticker.FormatStrFormatter("%0.1e")
    for ax in fig.axes:
        ax.yaxis.set_major_formatter(majorformatter)
        ax.set_title("")

    fig.savefig(filename, bbox_inches='tight', pad_inches=0.01)
