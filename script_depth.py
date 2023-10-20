from starter1 import *
import yt
import volavg
import fourier_tools_py3.fourier_filter as Filter
import brunt_tools as bt
reload(bt)
plt.close('all')
plotdir='%s/PigPen/'%(os.environ['HOME'])




#sim = '4_1'
#frame = 35
sim='1_1'
frame=31

if 'rho' not in dir():
    print('load and cg')
    rho_full, rho = bt.get_cubes(sim,frame)

if 1:
    #A bunch of cubes.
    #Test the brunt method.
    Nx,Ny,Nz = rho_full.shape
    sigma_3d=[]
    sigma_br=[]
    np.random.seed(90210)
    Rsize=512
    depth_list=[0,1,2,3,4,5,6]#,7,8]
    fig,ax=plt.subplots(3,3,figsize=(8,4))
    axlist=ax.flatten()

    for nextra in depth_list:
        #start = (np.random.random(3)*(Nx-128)).astype('int')
        #start = np.array([64]*3)
        start = np.array([Nx//2 - Rsize//2]*3)

        start[0] =0
        stop = start + Rsize
        stop[0]=Rsize/16
        if start[0]<0:
            print("BORK")
            start[0]=0
        stop[0] +=  nextra*(Rsize/6)
        print(stop[0])
        if stop[0] >= Nx:
            print("BORK2")
            stop[0] = Nx-1

        region = rho_full[start[0]:stop[0], start[1]:stop[1], start[2]:stop[2]]
        axlist[nextra].imshow(region.sum(axis=1))
        ftr = bt.fft_tool(region)
        ftr.apodize1()
        ftr.do2()
        sig_br, sig_3d = ftr.brunt_sigma()
        print("%0.2e %0.2e"%(sig_3d, sig_br.real))
        sigma_3d.append(sig_3d)
        sigma_br.append(sig_br.real)
    fig.savefig('%s/dumb'%plotdir)
    sigma_3d=np.array(sigma_3d)
    sigma_br=np.array(sigma_br)
    fig,ax=plt.subplots(1,1)
    ax.scatter(depth_list, sigma_3d/sigma_br)
    #ax.set(xscale='log',yscale='log')
    fig.savefig('%s/depth'%plotdir)

