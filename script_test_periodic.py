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
#sim = ms_ma ms = 123456 ma = half,1,2
frame=31

#get cubes; rho_full is straight off disk.  rho is downsampled by 2
if 'rho' not in dir():
    print('load and cg')
    rho_full, rho = bt.get_cubes(sim,frame)

if 'ftool' not in dir() or True:
    a=128
    b=a+128
    c=35
    d=c+128
    e=19
    f=e+128
    array=rho[a:b,c:d,e:f]
    array=rho
    ftool = bt.fft_tool(array)
    #ftool.do3()
    ftool.do2_periodic(projax=0)
    ftool.do2(projax=0)
    ftool2=bt.fft_tool(array)
    ftool2.apodize1(projax=0)
    ftool2.do2(projax=0)
if 1:

    #bt.plot_fft(ftool,"%s/ffts_%s_%04d"%(plotdir,sim,frame))
    #bt.plot_fft2(ftool,"%s/fft_periodic_%s_%04d"%(plotdir,sim,frame), apod=ftool2)
    bt.plot_fft3(ftool,outname="%s/fft_periodic_%s_%04d"%(plotdir,sim,frame), apod=ftool2)


if 0:
    #try brunt.  Works to a factor of 2, not bad.
    sigma_rho = (ftool.rho**2).sum().real/ftool.rho.size
    sigma_fft = (ftool.power_1d3).sum().real
    #sigma_fft = (ftool.power).sum().real/ftool.power.size**2
    print("rho %0.2e fft %0.2e ratio %0.2e"%(sigma_rho,sigma_fft,(sigma_fft)/sigma_rho))

    sigma_col = ((ftool.rho2)**2).sum().real/ftool.rho2.size
    sigma_2d_fft = ftool.power_1d2.sum().real
    print("col %0.2e fft %0.2e"%(sigma_col, sigma_2d_fft))
    Rinv = ( ftool.k2d*ftool.power_1d2).sum()/(ftool.power_1d2).sum()*256
    Rinv_actual = ftool.power_1d3.sum()/ftool.power_1d2.sum()
    sigma_Brunt = sigma_col*Rinv
    sigma_Brunt_actual = sigma_col*Rinv_actual
    print( sigma_Brunt/sigma_rho)
    print( sigma_Brunt_actual/sigma_rho)
