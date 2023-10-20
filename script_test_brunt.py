from starter1 import *
import yt
import volavg
import fourier_tools_py3.fourier_filter as Filter

sim = '4_1'
frame = 35

if 'ds' not in dir():
    print('load and cg')
    ds = yt.load("/data/cb1/Projects/P49_EE_BB/%s/DD%04d/data%04d"%(sim, frame, frame))
    print('get cg')
    cg = ds.covering_grid(0, [0.0]*3, [512]*3)
    dds = cg.dds
    rho_full = cg["density"].v
    rho = volavg.volavg(rho_full,rank=3,refine_by=2)

if 'fft' not in dir():
    rho2=rho.sum(axis=0)

    fft = np.fft.fftn( rho )
    power=fft*np.conjugate(fft)
    ff = Filter.FourierFilter(power)
    power_1d3 = np.array([power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
    power_1d3 /= rho.size
    Nzones = np.array([ff.get_shell(bin).sum() for bin in range(ff.nx)])
    kspace=ff.get_shell_k()

    fft2 = np.fft.fftn( rho2 )
    power2=fft2*np.conjugate(fft2)
    ff2 = Filter.FourierFilter(power2)
    power_1d2 = np.array([power2[ff2.get_shell(bin)].sum() for bin in range(ff2.nx)])
    power_1d2 /= rho2.size
    Nzones2 = np.array([ff2.get_shell(bin).sum() for bin in range(ff2.nx)])
    kspace2a=ff2.get_shell_k()

if 0:

    fig,ax=plt.subplots(1,1)
    if np.abs(power_1d2.imag).sum() > 1e-16:
        print("Imaginary Power")
    if np.abs(power_1d3.imag).sum() > 1e-16:
        print("Imaginary Power")
    ax.plot(kspace2a, power_1d2.real,c='r', label='P2d')

    ax.plot(kspace, power_1d3.real,c='g', label='P3d')
    ax.plot(kspace2a, kspace2a*power_1d2.real,c='b', label = 'k P2d')
    ax.legend(loc=0)

    ax.set(yscale='log',xscale='log')
    fig.savefig('%s/PigPen/test1'%(os.environ['HOME']))


if 1:
    #try brunt.  Works to a factor of 2, not bad.
    sigma_rho = (rho**2).sum().real
    #sigma_fft = (power_1d3).sum().real
    sigma_fft = (power).sum().real/power.size
    print("rho %0.2e fft %0.2e ratio %0.2e"%(sigma_rho,sigma_fft,(sigma_rho-sigma_fft)/sigma_rho))

    sigma_col = (rho2**2).sum()
    Rinv = ( kspace2a*power_1d2).sum()/(power_1d2).sum()
    sigma_Brunt = sigma_col*Rinv
    print( sigma_Brunt/sigma_rho)

if 1:
    #Why the factor of 2?
    s2=( kspace2a*power_1d2).sum()
    s3=( power_1d3).sum()
    print(s2,s3,s2/s3)

