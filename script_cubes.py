from starter1 import *
import yt
import volavg
import fourier_tools_py3.fourier_filter as Filter
import brunt_tools as bt
reload(bt)
plt.close('all')
plotdir='%s/PigPen/'%(os.environ['HOME'])



#The simulation and frame from the turbulence suite
sim = '4_1'
frame = 35

#get cubes; rho_full is straight off disk.  rho is downsampled by 2
if 'ds' not in dir():
    print('load and cg')
    rho_full, rho = bt.get_cubes(sim,frame)

if 1:
    #A bunch of cubes.
    #Test the brunt method.
    Nx,Ny,Nz = rho_full.shape
    sigma_3d=[]
    sigma_br=[]

    #Remove this to get actually random.
    np.random.seed(90210)
    Rsize=128
    for ncube in range(20):
        print("try",ncube)
        start = (np.random.random(3)*(Nx-128)).astype('int')
        stop = start + Rsize

        region = rho_full[start[0]:stop[0], start[1]:stop[1], start[2]:stop[2]]

        #make the tool
        ftr = bt.fft_tool(region)

        #sand down the corners
        ftr.apodize1()

        #do the fourier transform in 2d
        ftr.do2()
        #brunt method
        sig_br, sig_3d = ftr.brunt_sigma()
        sigma_3d.append(sig_3d)
        sigma_br.append(sig_br.real)

    fig,ax=plt.subplots(1,1)
    ax.scatter(sigma_3d, sigma_br)
    ax.set(xlabel='3d',ylabel='brunt')
    minmin=min(sigma_3d+sigma_br)
    maxmax=max(sigma_3d+sigma_br)
    ax.plot([minmin,maxmax],[minmin,maxmax])
    ax.plot([minmin,maxmax],[0.5*minmin,0.5*maxmax])
    ax.set(xscale='log',yscale='log', xlim=[minmin,maxmax],ylim=[minmin,maxmax])
    fig.savefig('%s/cubes'%plotdir)

