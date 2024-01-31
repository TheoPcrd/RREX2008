"""
TP 18/08/2023
Create time series for each month including :
-Tracer percentile evolution
-Mld evolution 
-Tracer concetration per layer

"""


import matplotlib
matplotlib.use('Agg')

import sys 
sys.path.append('/home2/datahome/tpicard/python/Python_Modules_RREX2008/')

######################
#import useful modules
######################

from Modules import *
from Modules_gula import *
from matplotlib.offsetbox import AnchoredText
import cartopy.crs as ccrs
import cartopy
from datetime import date, timedelta, datetime
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import scipy.stats as stats 
import netCDF4 as nc4

dir_save = "/home2/datawork/tpicard/M2_intership/RREX2008_FT_3h/time_evolution/"
month = sys.argv[1]
binz = np.arange(-600,6,3)

t = 0 #Global time

space_time = 1

nc = nc4.Dataset(dir_save+'tracer_percentile_{0}.nc'.format(month),'w') # CHANGE 
nc.createDimension('time',None) # unlimited dimension
nc.createDimension('binz',binz.shape[0]-1) #Number of bins on the vertical axe
nc.createVariable('tpas','f',('time','binz'))
nc.createVariable('Q','f',('time','binz'))
nc.createVariable('time','f','time')
nc.createVariable('mld_d_median','f','time')
nc.createVariable('mld_d_90p','f','time')
nc.createVariable('mld_d_10p','f','time')
nc.createVariable('tracer_99p','f','time')
nc.createVariable('tracer_50p','f','time')
nc.createVariable('tracer_85p','f','time')
nc.createVariable('tracer_90p','f','time')
nc.createVariable('tracer_95p','f','time')
nc.close()

for time in range(0,240,space_time):

    simu_name = 'rrex2008_avg_3h_{0}'.format(month)
    simul = load(simul = simu_name +' [0,810,0,1010,[1,200,1]] ' + format(time), light=False, output=False)
    tpas = var('tpas01',simul).data
    depth = 0
    mld_d= var('hbls_rho',simul).data
    z_r,z_w = tools.get_depths(simul)
    dz_r = np.diff(z_w, axis=-1)
    pnpm = 1./np.tile(simul.pm*simul.pn,(200,1,1))
    pnpm = np.transpose(pnpm,(1,2,0))

    ##############################################################
    # Define velume matrix
    ########################################################

    vol = np.multiply(pnpm,dz_r)
    Q = np.multiply(tpas,vol)

    ##########################################
    # Statistics part
    ##########################################

    Q_tot = np.sum(Q)
    bin_sum, bin_edges, binnumber = stats.binned_statistic(z_r[:,:,50:].ravel(), Q[:,:,50:].ravel(), statistic='sum', bins=binz)
    np.percentile(bin_sum*100/Q_tot,10)
    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width/2

    bin_sum_int = np.zeros(binz.shape[0])

    for i in range(len(bin_sum)):
        bin_sum_int[i] = np.sum(bin_sum[0:i]*100/Q_tot)

    percentile_99_index = np.abs(bin_sum_int - 1).argmin()
    percentile_90_index = np.abs(bin_sum_int - 10).argmin()
    percentile_50_index = np.abs(bin_sum_int - 50).argmin()
    percentile_95_index = np.abs(bin_sum_int - 5).argmin() 

    tracer_99p = bin_centers[percentile_99_index]
    tracer_50p = bin_centers[percentile_50_index]
    tracer_90p = bin_centers[percentile_90_index]
    tracer_95p = bin_centers[percentile_95_index]

    ##########################################
    # Statistics part
    ##########################################

    z_r,z_w = tools.get_depths(simul)
    bin_sum_tpas, bin_edges, binnumber = stats.binned_statistic(z_r[:,:,50:].ravel(), tpas[:,:,50:].ravel(), statistic='mean', bins=binz)
    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width/2

    #MLD_D INTERVAL

    mld_d_list = mld_d.flatten()
    mld_d_median = np.median(mld_d_list)
    mld_d_10p = np.percentile(mld_d_list,10)
    mld_d_90p = np.percentile(mld_d_list,90)

    ##########################################
    # Copy variables in nc file 
    ##########################################

    nc = nc4.Dataset(dir_save+'tracer_percentile_{0}.nc'.format(month),'r+') # CHANGE
    nc.variables['time'][t] = time
    nc.variables['mld_d_median'][t] = mld_d_median
    nc.variables['mld_d_10p'][t] = mld_d_10p
    nc.variables['mld_d_90p'][t] = mld_d_90p
    nc.variables['tracer_99p'][t] = tracer_99p
    nc.variables['tracer_50p'][t] = tracer_50p
    nc.variables['tracer_90p'][t] = tracer_90p
    nc.variables['tracer_95p'][t] = tracer_95p
    nc.variables['tpas'][t,:] = bin_sum_tpas
    nc.variables['Q'][t,:] = bin_sum
    t = t+1
    nc.close()

