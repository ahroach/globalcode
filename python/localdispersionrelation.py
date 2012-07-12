import math
import cmath
from numpy import *
import scipy
import matplotlib.pyplot as pyplot

#For comparison to Ji and Goodman: k=0.44, zeta=0.07, Omega=186.13, m=0,
#theta=0.785, eta=2000, rho=6.0
#So kr = ksin\theta = 0.311
# and kz = kcos\theta = 0.311

def find_roots_simple(B=0.0, kr=0.0, kz=0.0, Omega=0.0, zeta=0.0, m=1.0,
                      nu=2.98e-3, eta=2.57e3, rho=6.36, r1=7.06, r2=20.3):
    """Solve the local dispersion relation."""
    pi = math.pi

    ktheta = m/(0.5*(r2 + r1))
    k = sqrt(kr**2 + kz**2 + ktheta**2)
    va = B*1.0/math.sqrt(4.0*pi*rho)
    wa = kz*va #Alfven frequency \omega_A

    #gn = \gamma_n = \gamma + \eta k^2 = i\omega + \eta k**2
    gn = scipy.poly1d([1, eta*k**2])
    gv = scipy.poly1d([1, nu*k**2])
    
    dr = (gn*gv + wa**2)**2 + 2*Omega**2*zeta*(kz**2/k**2)*gn**2 \
         + 2*Omega**2*(zeta-2)*(kz**2/k**2)*wa**2
    
    dr_roots = roots(dr)
    
    #Convert to \omega
    #Since \gamma t = -i \omega t, \omega = i \gamma
    
    dr_roots_return = sort(1j*dr_roots)
    
    return dr_roots_return



def find_roots(B=0.0, kr=0.0, kz=0.0, Omega=0.0, zeta=0.0, m=0,
               nu=2.98e-3, eta=2.57e3, rho=6.36, r1=7.06, r2=20.3):

    #Parameters are in cgs units (Gauss, cm, 1/cm, cm^2/sec, gm/cm^3, etc)
    #We will be solving for the Doppler shifted frequency
    #\bar{\omega} = \omega - m\Omega
    
    pi = math.pi
    
    ktheta = m/(0.5*(r2+r1))
    k = sqrt(kr**2 + kz**2 + ktheta**2)
    va = B*1.0/math.sqrt(4.0*pi*rho)
    wa = kz*va #Alfven frequency \omega_A
    wr = (zeta-2)*Omega*kr*ktheta/k**2 #Rossby wave frequency \omega_R
    kfrac = kz/k
    
    #wn = \omega_\eta = \bar{\omega} - i\eta k^2
    #wv = \omega_\nu = \bar{\omega} - i\nu k^2
    
    wn = scipy.poly1d([1, -1j*eta*k**2])
    wv = scipy.poly1d([1, -1j*nu*k**2])

    #Note that for all of these terms are start off with the operations
    #on the poly1d objects, because otherwise I think I've seen instances
    #where the operator overloading fails, and multiplies, adds, and
    #exponentiations don't treat the objects correctly.

    #Term 1: \omega_{\eta}(\omega_A^2 - \omega_{\eta}\omega_{\nu})^2
    
    term1 = wn*(-wn*wv + wa**2)**2
    
    #Term 2: -2\zeta\Omega^2(k_z^2/k^2)\omega_{\eta}^3
    
    term2 = -wn**3*2.0*zeta*Omega**2*kfrac**2
    
    #Term 3: 2(\zeta-2)\Omega^2\omega_A^2(k_z^2/k^2)\omega_{\eta}
    
    term3 = wn*2.0*(zeta-2)*Omega**2*wa**2*kfrac**2
    
    #Term 4: i*\omega_R*(\omega_A^2 - \omega_{\nu}\omega_{\eta})*
    #(\omega_A^2 - \omega_{\eta}^2)
    #where \omega_R = (\zeta-2)\Omega k_r k_{\theta} / k^2
    
    term4 = (-wv*wn + wa**2)*(-wn**2 + wa**2)*1j*wr
    
    dr = term1 + term2 + term3 + term4

    #Minus sign here, because this whole thing was derived with
    #a bizarre assumed dependence as exp(\omega t - k\cdot x), when
    #a reasonable right-going wave is exp(k\cdot x - \omega t).
    
    dr_roots = -roots(dr)

    return sort(dr_roots)

def sort_roots(results):
    #Sorts the results to give smooth curves when the results are plotted
    LARGE=1.0e12 #Very large constant

    #Make an array to store the sorted results
    sorted_results = zeros(results.shape, dtype=complex)
    #We assume that the B or k is varying in the 0th-dimension
    solutions = results.shape[0]
    #And the roots are stored in the 1st dimension
    roots = results.shape[1]

    #Start off by just copying the first row of roots.
    sorted_results[0,:] = results[0,:]

    #Now for each of the remaining rows...
    for i in range(1, solutions):
        #Copy in the roots from results
        tmpvector = results[i,:]
        #And for all of the columns except the last in the sorted results
        for k in range(0,roots-1):
            distance=LARGE
            goodelem=5000

            #Look through the unclaimed results and see which is closest
            #to the root in the same column but in the previous row
            #of the sorted results
            
            for j in range(0,tmpvector.size):
                if (abs(tmpvector[j]-sorted_results[i-1,k]) < distance):
                    distance = abs(tmpvector[j]-sorted_results[i-1,k])
                    goodelem = j

            #Once we've found it, copy that over to the sorted_results array
            sorted_results[i,k] = tmpvector[goodelem]

            #And remove that element from the list of possibilities,
            #so noone else claims it
            tmpvector = delete(tmpvector, goodelem)

        #And once we get to the last column, just copy in the last
        #remaining root
        sorted_results[i,roots-1] = tmpvector[0]

    #Done!
    return sorted_results


def plot_vs_kz(B=0.0, kr=0.0, kz=logspace(-3,2,1000), Omega=0.0, zeta=0.0,
               m=0, nu=2.98e-3, eta=2.57e3, rho=6.36, r1=7.06, r2=20.3,
               logscalex=1):

    omegas = zeros([kz.size, 5], dtype=complex)
    for i in range(0, kz.size):
        omegas[i,:] = find_roots(B=B, kr=kr, kz=kz[i], Omega=Omega,
                                 zeta=zeta, m=m, nu=nu, eta=eta, rho=rho,
                                 r1=r1, r2=r2)
    
    omegas = sort_roots(omegas)

    fig = pyplot.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)

    for i in range(0,5):
        ax1.plot(kz, omegas[:,i].real)
        ax2.plot(kz, omegas[:,i].imag)
    
    ax1.set_xlabel(r"Re[$\omega$]")
    ax2.set_ylabel(r"Im[$\omega$]")
    ax2.set_xlabel(r"$k_z$ [rad/cm]")
    ax2.autoscale(axis='y', tight=True)
    ax2.set_ylim(bottom=-1)
    if(logscalex):
        ax1.set_xscale('log')
        ax2.set_xscale('log')
