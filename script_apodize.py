from starter1 import *
import yt
import volavg
import fourier_tools_py3.fourier_filter as Filter
import brunt_tools as bt
plt.close('all')
plotdir='%s/PigPen/'%(os.environ['HOME'])




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

#if 'ftool' not in dir():
#    ftool = fft_tool(rho)
#    ftool.do3()
#    ftool.do2(projax=0)


if 1:
    Nx,Ny,Nz = rho_full.shape
    #start = (np.random.random(3)*(Nx-128)).astype('int')
    start = np.array([64]*3)
    stop = start + 128
    print(start,stop)


    region = rho_full[start[0]:stop[0], start[1]:stop[1], start[2]:stop[2]]
    ftr = bt.fft_tool(region)
    ftr.apodize1()
    ftr.do2()

    ft1 = bt.fft_tool(region)
    ft1.do2()

    fig,ax=plt.subplots(3,2)
    cmap_name='jet'
    cmap = copy.copy(mpl.cm.get_cmap(cmap_name))
    cmap.set_under('w')
    norm = mpl.colors.Normalize(vmin=ft1.rho2.min(),vmax=ft1.rho2.max())
    imargs={'origin':'lower','interpolation':'nearest', 'cmap':cmap, 'norm':norm}
    ax[0][0].imshow( ft1.rho2, **imargs)
    ax[0][1].plot( ft1.k2d, ft1.power_1d2.real)
    ax[0][1].set(xscale='log',yscale='log')
    ax[1][0].imshow( ftr.rho2, **imargs)
    ax[0][1].plot( ftr.k2d, ftr.power_1d2.real)
    ax[1][1].set(xscale='log',yscale='log')

    norm = mpl.colors.Normalize(vmin=1e-8,vmax=1)
    ax[2][0].imshow( ftr.window,cmap=cmap,norm=norm )
    fig.savefig('%s/apod2'%plotdir)


