# 22,Yatir Forest,DarkSeaGreen1,#C1FFC1
# 23,Redwood Forest,Maroon,#800000

import numpy as np
import os
import socket
import xarray as xr
import holoviews as hv
from holoviews import opts
import geoviews as gv
import geoviews.feature as gf
import pandas as pd

import panel as pn

def get_vdim(ds, varname):
    """find the name of the vertical dimension for an xarray variable
    """
    vertical_dim = [d for d in ds[varname].dims if 'bottom_top' in d]
    if any(vertical_dim):
        if len(vertical_dim) > 1:
            raise(ValueError('{} has more than one '
                             'vertical dimension'.format(this_var)))
        else:
            vertical_dim = vertical_dim[0]
    return(vertical_dim)


def get_min_max(ds, varname, hour, zlev):

    vdim = get_vdim(ds, varname)
    idx_run = xr.DataArray(['yatir wet', 'yatir dry'], dims=['WRFrun'])
    idx_dict = {'WRFrun': idx_run,
                'hour': hour}
    if vdim != []:
        idx_dict[vdim] = zlev

    vmin = ds[varname].sel(idx_dict).min().values.tolist()
    vmax = ds[varname].sel(idx_dict).max().values.tolist()
    # vmin = ds[varname].min()
    # vmax = ds[varname].max()
    # print('min, max:', vmin, vmax)
    return((vmin, vmax))

def merge_yatir_fluxes_landuse(fname_ctl='ctl_run_d03_diag_latest.nc',
                               fname_yatir='yatir_run_d03_diag_latest.nc'):
    """merge WRF fluxes and landuse into single xarray dataset
    """
    if 'MacBook' in socket.gethostname():
        cscratch_path = os.path.join(os.path.join('/', 'Users', 'tim',
                                                  'work', 'Data',
                                                  'SummenWRF', 'yatir'))
    elif 'cori' in socket.gethostname():
        cscratch_path = os.path.join('/', 'global', 'cscratch1', 'sd',
                                     'twhilton', 'yatir_output_collected')
    ctlday = WRF_daily_daylight_avg(os.path.join(cscratch_path, fname_ctl))
    ytrday = WRF_daily_daylight_avg(os.path.join(cscratch_path, fname_yatir))
    # landuse_data = yatir_landuse_to_xarray()

    # ytrday = ytrday.assign(
    #     {'LU_INDEX':
    #      landuse_data['d03'].sel(WRFrun='ytr')['LU_INDEX']})
    ctlday, ytrday = (this_dataset.assign(
        {# 'LU_INDEX':
         # landuse_data['d03'].sel(WRFrun=this_key)['LU_INDEX'],
         # 'LANDUSEF':
         # landuse_data['d03'].sel(WRFrun=this_key)['LANDUSEF'],
         'height_agl_stag':
         this_dataset['zstag'] - this_dataset['ter']
        })
                      for (this_dataset, this_key) in
                      zip((ctlday, ytrday), ('ctl', 'ytr')))
    ds_diff =  (ctlday - ytrday)  #.assign_coords({'WRFrun': 'control - Yatir'})

    return(ctlday, ytrday, ds_diff)

def set_attributes_for_plotting(ds):
    """set attrs of xarray dataset for plotting
    """
    ds['XLONG'].attrs['long_name'] = 'Longitude'
    ds['XLONG'].attrs['units'] = 'deg E'
    ds['XLAT'].attrs['long_name'] = 'Latitude'
    ds['XLAT'].attrs['units'] = 'deg N'
    ds['height_agl_stag'].attrs['long_name'] = 'height above ground level'
    ds['XLAT'].attrs['units'] = 'm'


    try:
        for k, v in ds.data_vars.items():
            ds[k].attrs['long_name'] = v.attrs['description']
    except KeyError:
        ds[k].attrs['long_name'] = k
        print(('variable {} has no attribute \'description\','
               ' setting long_name to {}'.format(k, k)))

    return(ds)

def yatir_WRF_to_xarray(fname):
    """read Yatir forest WRF output to xarray

    returns an xarray suitable for plotting with
    [GeoViews](http://geoviews.org/user_guide/)
    """
    ds = xr.open_dataset(fname)
    for dim in ['XLAT', 'XLONG']:
        ds[dim] = ds[dim].isel(Time=0)
    return(ds)


def WRF_daily_daylight_avg(fname):
    """read Yatir forest WRF output to xarray containing daily means

    returns an xarray suitable for plotting with
    [GeoViews](http://geoviews.org/user_guide/)
    """
    ds = yatir_WRF_to_xarray(fname)
    # is_daytime = ds['SWDOWN'] > 0.1
    # ds_day_mean = ds.where(is_daytime, drop=True).groupby('XTIME.hour').mean(keep_attrs=True)
    for this_var in ds.data_vars:
        vertical_dim = [d for d in ds[this_var].dims if 'bottom_top' in d]
        if any(vertical_dim):
            if len(vertical_dim) > 1:
                raise(ValueError('{} has more than one '
                                 'vertical dimension'.format(this_var)))
            else:
                vertical_dim = vertical_dim[0]
            # keep only the surface value
            # bottom_level_only = ds[this_var].sel({vertical_dim: 0})
            # ds[this_var] = bottom_level_only
    ds_day_mean = ds.groupby('XTIME.hour').mean(keep_attrs=True)
    return(ds_day_mean)

def yatir_landuse_to_xarray():
    """parse land use data for Yatir, Control runs for domains d02 and d03

    RETURNS:
      dict keyed by ['d02', 'd03']; values are xarray.DataSet objects
      containing land use data concatenated on new dimension WRFrun
    """
    if 'MacBook' in socket.gethostname():
        dir_path = os.path.join('/', 'Users', 'tim', 'work',
                                'Data', 'SummenWRF', 'yatir')
    elif 'cori' in socket.gethostname():
        dir_path = os.path.join('/', 'global', 'cscratch1', 'sd',
                                'twhilton', 'yatir_land_use')
    ctable = get_IGBP_modMODIS_21Category_PFTs_table()
    land_cat_names = list(ctable['long_name'])
    land_cat_names = [x if x != 'BareGroundTundra' else 'Yatir'
                      for x in land_cat_names]
    dict_runs = {}
    for WRFdomain in ['d02', 'd03']:
        dict_runs[WRFdomain] = xr.concat((xr.open_dataset(
            os.path.join(dir_path, 'land_data_{run}_{dom}.nc').format(
                run=WRFrun, dom=WRFdomain)).squeeze()
                                          for WRFrun in ['ctl', 'ytr']),
                                         dim='WRFrun')
        # remove the hyphen from z dimension name.  '-' also being an
        # operator messes up assign_coords()
        # dict_runs[WRFdomain] = dict_runs[WRFdomain].rename(
        #     {'z-dimension0021': 'zdimension0021'})

        # assign integral coordinate values to spatial coordinate
        # variables
        new_coords = {this_var: range(dict_runs[WRFdomain][this_var].size) for
                      this_var in ['south_north', 'west_east']}
        new_coords = {**new_coords,
                      'WRFrun': ['ctl', 'ytr'],
                      'z-dimension0021': land_cat_names,
                      'lat': (('west_east', 'south_north'),
                              dict_runs[WRFdomain]['CLAT'][0, ...].values),
                      'lon': (('west_east', 'south_north'),
                              dict_runs[WRFdomain]['CLONG'][0, ...].values)}
        dict_runs[WRFdomain] = dict_runs[WRFdomain].assign_coords(new_coords)
        dict_runs[WRFdomain] = dict_runs[WRFdomain].rename(
            {'z-dimension0021': 'PFT'})
    return(dict_runs)


if __name__ == '__main__':
    test_get_landuse_to_xarray = False
    test_get_data_file = False
    test_get_xarray = False
    test_get_xarray_daily = False
    test_merge = False
    test_postprocessed = True

    if test_get_landuse_to_xarray:
        dict_runs = yatir_landuse_to_xarray()

    if test_get_data_file:
        roughness_data_dir = os.path.join('/', 'Users',
                                          'tim', 'work',
                                          'Data', 'SummenWRF',
                                          'yatir', 'roughness_experiments')
        fnames = {'ctl': 'energyfluxes_control_d03.nc',
                  'yatir': 'energyfluxes_yatir_run_d03.nc',
                  'vegparm_z0_10': 'energyfluxes_vegparm_z0_10_d03.nc',
                  'landuse_z0_500': 'energyfluxes_landuse_z0_500_d03.nc'}
        data_paths = {k: os.path.join(roughness_data_dir, v)
                      for k, v in fnames.items()}
        dsgv, dsxr = get_data_file(data_paths)

    if test_get_xarray:
        # LHctl = yatir_to_xarray(os.path.join('/', 'Users', 'tim', 'work',
        #                                      'Data', 'SummenWRF', 'yatir',
        #                                      'LH_d03_yatirZ50.nc'),
        #                         varname='LH',
        #                         groupname='ctl',
        #                         timerange=10)
        # LHytr = yatir_to_xarray(os.path.join('/', 'Users', 'tim', 'work',
        #                                      'Data', 'SummenWRF', 'yatir',
        #                                      'LH_d03_yatirZ50.nc'),
        #                         varname='LH',
        #                         groupname='yatirZ050',
        #                         timerange=10)

        ctl = yatir_WRF_to_xarray('/Users/tim/work/Data/SummenWRF/yatir/fluxes_ctl_run_d03.nc')
        ytr = yatir_WRF_to_xarray('/Users/tim/work/Data/SummenWRF/yatir/fluxes_yatir_run_d03.nc')

    if test_get_xarray_daily:
        ctl = yatir_WRF_to_xarray('/Users/tim/work/Data/SummenWRF/yatir/fluxes_ctl_run_d03.nc')
        ytr = yatir_WRF_to_xarray('/Users/tim/work/Data/SummenWRF/yatir/fluxes_yatir_run_d03.nc')
        ctlday = WRF_daily_daylight_avg('/Users/tim/work/Data/SummenWRF/yatir/fluxes_ctl_run_d03.nc')
        ytrday = WRF_daily_daylight_avg('/Users/tim/work/Data/SummenWRF/yatir/fluxes_yatir_run_d03.nc')
        dims = define_dims(ctl)

    if test_merge:
        ctlday, ytrday, ctl_minus_ytr = merge_yatir_fluxes_landuse()
        (ctlday_TP,
         ytrday_TP,
         ctl_minus_ytr_TP) = merge_yatir_fluxes_landuse(
             fname_ctl='ctl_run_d03_diag_TP.nc',
             fname_yatir='yatir_run_d03_diag_TP.nc')
        ctlall = xr.merge((ctlday, ctlday_TP))

    if test_postprocessed:
        (ctlday_TP,
         ytrday_TP,
         ctl_minus_ytr_TP) = merge_yatir_fluxes_landuse(
             fname_ctl='/global/cscratch1/sd/twhilton/yatir_output_collected/wetsoil/ctl_run_d03_diag_TP_VWCx2.nc',
             fname_yatir='/global/cscratch1/sd/twhilton/yatir_output_collected/wetsoil/yatir_run_d03_diag_TP_VWCx2.nc')
             # fname_ctl='ctl_d03_postprocessed.nc',
             # fname_yatir='ytr_d03_postprocessed.nc')
