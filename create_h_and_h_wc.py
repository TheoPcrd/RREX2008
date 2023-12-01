#CREATION OF JPDF DATABASE

import matplotlib
matplotlib.use('Agg')

import sys 
sys.path.append('/home2/datahome/tpicard/python/Python_Modules_RREX2008_FT/')
#str_month = sys.argv[0]

######################
# import useful modules
# #####################

from Modules import *
from Modules_gula import *
from matplotlib.offsetbox import AnchoredText
import cartopy.crs as ccrs
import cartopy
from datetime import date, timedelta, datetime 
import netCDF4
import scipy.stats as st
import matplotlib as mpl
import matplotlib.colors as colors
import netCDF4 as nc4
from Tools_jpdf import *

# tstart = 0
# tend = 232 = 29 days

def create_H_month(month,t_start,t_end,gap_list,relative_mld=True): 
    
    binx = np.arange(-6,6.05,0.05)
    biny = np.arange(0,6.025,0.025)
    
    t = 0

    H_jpdf = np.zeros((len(binx)-1,len(biny)-1,len(gap_list),t_end-t_start))
    H_wc = np.zeros((len(binx)-1,len(biny)-1,len(gap_list),t_end-t_start))
    H_w = np.zeros((len(binx)-1,len(biny)-1,len(gap_list),t_end-t_start))
        
    simu_name = 'rrex2008_avg_3h_{0}'.format(month)
    #test tstart
    simul = load(simul = simu_name+' [100,700,100,900,[1,200,1]] ' + format(t_start), light=False, output=True)
    u = var('u',simul,depths=[0]).data    
    #test tend
    simul = load(simul = simu_name+' [100,700,100,900,[1,200,1]] ' + format(t_end), light=False, output=True)
    u = var('u',simul,depths=[0]).data
    
    for time in range(t_start,t_end):
            
            if relative_mld: # Si on prend les niveau par rapport a la mld
                
                print('gap relative to mld')
                mld_d = var('hbls_rho',simul).data
                mld_d = np.repeat(mld_d[:, :, np.newaxis], 200, axis=2)
            else :
                print('gap relative to z')
                
            simul = load(simul = simu_name+' [100,700,100,900,[1,200,1]] ' + format(time), light=False, output=False)
            #mld calcul

            tpas = var('tpas01',simul).data
            u = var('u',simul,depths=[0]).data
            v = var('v',simul,depths=[0]).data
            w = var('w',simul).data
            #w = tools.nanbnd(var('w',simul).data) # SI FICHIER HIS
            #w = np.nan_to_num(w)
            strain =  tools.get_strain(u,v,simul.pm,simul.pn) / simul.f
            vrt =  tools.psi2rho(tools.get_vrt(u,v,simul.pm,simul.pn) / tools.rho2psi(simul.f))
            z_r,z_w = tools.get_depths(simul)

            g = 0
            for gap in gap_list:
                if relative_mld:
                    #print("gap < 1")
                    index_mld= np.argmin(np.abs(z_r+mld_d-gap*mld_d),axis=2) # Index of mld for each point
                else:
                    #print("gap > 1")
                    index_mld= np.argmin(np.abs(z_r+gap),axis=2) # Index of mld for each point
                    
                w_mld = np.zeros(z_r.shape[0:2])
                tpas_mld = np.zeros(z_r.shape[0:2])

                for i in range(w_mld.shape[0]):
                     for j in range(w_mld.shape[1]):

                        w_mld[i,j] = w[i,j,index_mld[i,j]]
                        tpas_mld[i,j] = tpas[i,j,index_mld[i,j]]
                            
                wC = np.multiply(tpas_mld,w_mld)
                
                #JPDF
                H_jpdf_t, xedges, yedges, _ = st.binned_statistic_2d(vrt.ravel(),strain.ravel(),vrt.ravel(),\
                                'count', bins=[binx, biny])
                H_jpdf_t = np.rot90(H_jpdf_t); H_jpdf_t = np.flipud(H_jpdf_t)
                H_jpdf[:,:,g,t] = H_jpdf_t

                #JPDF conditionned wc
                H_wc_d_t, xedges, yedges, _ = st.binned_statistic_2d(vrt.ravel(),strain.ravel(),wC.ravel(),\
                'sum', bins=[binx, biny])
                H_wc_d_t = np.rot90(H_wc_d_t); H_wc_d_t = np.flipud(H_wc_d_t)
                H_wc[:,:,g,t] = H_wc_d_t
                
                #JPDF conditionned w
                H_w_d_t, xedges, yedges, _ = st.binned_statistic_2d(vrt.ravel(),strain.ravel(),w_mld.ravel(),\
                'mean', bins=[binx, biny])
                H_w_d_t = np.rot90(H_w_d_t); H_w_d_t = np.flipud(H_w_d_t)
                H_w[:,:,g,t] = H_w_d_t
                g=g+1
                
                
            t=t+1

    return(H_jpdf,H_wc,H_w)

def create_H_nc(month,t_start,t_end,gap_list,relative_mld = True):
    
    binx = np.arange(-6,6.05,0.05)
    biny = np.arange(0,6.025,0.025)
    
    (H_jpdf,H_wc,H_w) = create_H_month(month,t_start,t_end,gap_list,relative_mld)
    
    nb_t = t_end-t_start
    if relative_mld:
        nc = nc4.Dataset('/home2/datawork/tpicard/M2_intership/RREX2008_FT_3h/JPDF_filtred_domain/JPDF_{0}_mld.nc'.format(month),'w')
    else:
        nc = nc4.Dataset('/home2/datawork/tpicard/M2_intership/RREX2008_FT_3h/JPDF_filtred_domain/JPDF_{0}_z.nc'.format(month),'w')
        
    nc.createDimension('binx', len(binx))
    nc.createDimension('biny', len(biny))
    nc.createDimension('Hx', len(binx)-1)
    nc.createDimension('Hy', len(biny)-1)
    nc.createDimension('z', len(gap_list))
    nc.createDimension('time',nb_t)
    nc.createVariable('H_jpdf', 'f4', ('Hx', 'Hy','z','time'))
    nc.createVariable('H_wc', 'f4', ('Hx', 'Hy','z','time'))
    nc.createVariable('H_w', 'f4', ('Hx', 'Hy','z','time'))
    nc.createVariable('gap_list', 'f4', ('z'))
    nc.variables['H_jpdf'][:,:,:,:] = H_jpdf
    nc.variables['H_wc'][:,:,:,:] = H_wc
    nc.variables['H_w'][:,:,:,:] = H_w
    nc.variables['gap_list'][:] = gap_list
    nc.close()


'''
def create_diag_nc(month,t_start,t_end):

    simul = load(simul = 'rrex2008_TS_{0} [0,802,0,1002,[1,200,1]] '.format(month) + format(time1), light=False, output=True)
    
    tpas01_rate = var('tpas01_rate_mld',simul).data
    
    tpas01_Vadv = np.zeros((len(binx)-1,len(biny)-1,len(gap_list),t_end-t_start))
    tpas01_vmix = np.zeros((len(binx)-1,len(biny)-1,len(gap_list),t_end-t_start))
    tpas01_rate = np.zeros((len(binx)-1,len(biny)-1,len(gap_list),t_end-t_start))
    tpas01_hmix = np.zeros((len(binx)-1,len(biny)-1,len(gap_list),t_end-t_start))
    tpas01_Hadv = np.zeros((len(binx)-1,len(biny)-1,len(gap_list),t_end-t_start))
    
    for time in range(t_start,t_end):
        
        simul = load(simul = 'rrex2008_TS_{0} [0,802,0,1002,[1,200,1]] '.format(month) + format(time1), light=False, output=True)

        tpas01_Vadv = var('tpas01_Vadv_mld',simul).data
        tpas01_vmix = var('tpas01_vmix_mld',simul).data
        tpas01_rate = var('tpas01_rate_mld',simul).data
        tpas01_hmix = var('tpas01_hmix_mld',simul).data
        tpas01_Hadv = var('tpas01_Hadv_mld',simul).data
    
    (H_jpdf,H_wc) = create_H_month(month,t_start,t_end,gap_list)
    nb_t = t_end-t_start
    nc = nc4.Dataset('/home2/datawork/tpicard/M2_intership/RREX2008_FT_6h/JPDF_{0}.nc'.format(month),'w')
    nc.createDimension('binx', len(binx))
    nc.createDimension('biny', len(biny))
    nc.createDimension('Hx', len(binx)-1)
    nc.createDimension('Hy', len(biny)-1)
    nc.createDimension('z', len(gap_list))
    nc.createDimension('time',nb_t)
    nc.createVariable('H_jpdf', 'f4', ('Hx', 'Hy','z','time'))
    nc.createVariable('H_wc', 'f4', ('Hx', 'Hy','z','time'))
    nc.variables['H_jpdf'][:,:,:,:] = H_jpdf
    nc.variables['H_wc'][:,:,:,:] = H_wc
    nc.close()
'''

