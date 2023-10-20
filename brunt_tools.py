
from starter1 import *
import yt
import volavg
import fourier_tools_py3.fourier_filter as Filter

def get_cubes(sim,frame):
    ds = yt.load("/data/cb1/Projects/P49_EE_BB/%s/DD%04d/data%04d"%(sim, frame, frame))
    print('get cg')
    cg = ds.covering_grid(0, [0.0]*3, [512]*3)
    dds = cg.dds
    rho_full = cg["density"].v
    rho = volavg.volavg(rho_full,rank=3,refine_by=2)
    return rho_full, rho

def plot_fft(ftool, outname=None,ax=None):
    if ax is None:
        fig,ax=plt.subplots(1,1)
    if ftool.done2:
        ax.plot(ftool.k2d, ftool.power_1d2.real,c='r', label='P2d')
    if ftool.done3:
        ax.plot(ftool.k3d, ftool.power_1d3.real,c='g', label='P3d')

    if ftool.done2 and ftool.done3:
        ax.plot(ftool.k2d,ftool.k2d*ftool.power_1d2.real,c='b', label = 'k P2d')
    ax.legend(loc=0)

    ax.set(yscale='log',xscale='log')
    if ax is not None:
        fig.savefig(outname)

class fft_tool():
    def __init__(self,rho):
        self.rho=rho
        self.rho2=None
        self.done2=False
        self.done3=False

    def do3(self):
        print('3d fourier transform')
        self.fft3 = np.fft.fftn( self.rho )
        self.power=self.fft3*np.conjugate(self.fft3)
        #self.power/=self.power.size
        ff = Filter.FourierFilter(self.power)
        self.power_1d3 = np.array([self.power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
        self.power_1d3 /= self.rho.size
        self.Nzones = np.array([ff.get_shell(bin).sum() for bin in range(ff.nx)])
        self.k3d=ff.get_shell_k()
        self.done3=True

    def brunt_sigma(self):
        sigma_3d = (self.rho**2).sum()
        sigma_col = (self.rho2**2).sum()
        Rinv = (self.k2d*self.power_1d2).sum()/self.power_1d2.sum()
        sigma_brunt = sigma_col*Rinv
        return sigma_brunt, sigma_3d

    def do2(self,projax=0):
        print('2d fourier transform')
        if self.rho2 is None:
            print('MAKE NEW PROJECTION')
            self.rho2=self.rho.sum(axis=projax)
            Nz = self.rho.shape[projax]
            self.rho2/=Nz
        self.fft2 = np.fft.fftn( self.rho2 )
        self.power2=self.fft2*np.conjugate(self.fft2)
        #self.power2/=self.power2.size
        ff2 = Filter.FourierFilter(self.power2)
        self.power_1d2 = np.array([self.power2[ff2.get_shell(bin)].sum() for bin in range(ff2.nx)])
        self.power_1d2 /= self.rho2.size
        self.Nzones2 = np.array([ff2.get_shell(bin).sum() for bin in range(ff2.nx)])
        self.k2d=ff2.get_shell_k()
        self.done2=True

    def apodize1(self,projax=0):
        self.rho2=self.rho.sum(axis=projax)
        shape=np.array(self.rho2.shape)
        baseshape=(1.*shape).astype('int')
        base = np.zeros(baseshape)
        start = baseshape//2-shape//2

        base[start[0]:(start[0]+shape[0]), start[1]:(start[1]+shape[1])] = self.rho2

        if 0:
            #things that don't quite work
            self.rho2 = base
            from scipy.signal import general_gaussian
            from scipy.signal import convolve2d
            window = np.outer(general_gaussian(baseshape[0],6,3),general_gaussian(baseshape[0],6,3))
            #np.roll(window,baseshape[0]//2,axis=0)
            #np.roll(window,baseshape[0]//2,axis=1)
        if 0:
            #actually work ok
            window = np.zeros(baseshape)
            window[0:3,0:3]=1
        if 1:
            #works pretty well.
            x = np.arange(baseshape[0])
            sigma_conv=2
            g = np.exp(-x**2/(2*sigma_conv**2))**6
            window = np.outer(g,g)

        window/=window.sum()
        self.window=window
        

        #self.rho2 = scipy.convolve(base, window)
        #self.rho2 = convolve2d(base, window)
        #convolve the window function with the base
        a = np.fft.fftn(base)
        b = np.fft.fftn(window)
        c = a*b
        self.rho2 = np.fft.ifftn(c)
        q=np.abs(self.rho2.imag).sum()/self.rho2.imag.size
        if q>1e-13:
            print("!!!!!!!!!!!!!!!!!Imaginary",q)
        self.rho2=self.rho2.real

