# 22,Yatir Forest,DarkSeaGreen1,#C1FFC1
# 23,Redwood Forest,Maroon,#800000

import os
import socket
import xarray as xr


def get_vdim(ds, varname):
    """find the name of the vertical dimension for an xarray variable
    """
    vertical_dim = [d for d in ds[varname].dims if 'bottom_top' in d]
    if any(vertical_dim):
        if len(vertical_dim) > 1:
            raise(ValueError('{} has more than one '
                             'vertical dimension'.format(varname)))
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
    cscratch_path = os.getcwd()
    ctlday = WRF_daily_daylight_avg(os.path.join(cscratch_path, fname_ctl))
    ytrday = WRF_daily_daylight_avg(os.path.join(cscratch_path, fname_yatir))
    ctlday, ytrday = (this_dataset.assign(
        {'height_agl_stag': this_dataset['zstag'] - this_dataset['ter']}
    )
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
    for this_var in ds.data_vars:
        vertical_dim = [d for d in ds[this_var].dims if 'bottom_top' in d]
        if any(vertical_dim):
            if len(vertical_dim) > 1:
                raise(ValueError('{} has more than one '
                                 'vertical dimension'.format(this_var)))
            else:
                vertical_dim = vertical_dim[0]
    ds_day_mean = ds.groupby('XTIME.hour').mean(keep_attrs=True)
    return(ds_day_mean)


def yatir_landuse_to_xarray():
    """parse land use data for Yatir, Control runs for domains d02 and d03

    RETURNS:
      dict keyed by ['d02', 'd03']; values are xarray.DataSet objects
      containing land use data concatenated on new dimension WRFrun
    """
    dict_runs = {}
    for WRFdomain in ['d02', 'd03']:
        dict_runs[WRFdomain] = xr.concat((xr.open_dataset(
            os.path.join(os.getcwd(), 'land_data_{run}_{dom}.nc').format(
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
                      this_var in ['south_north', 'west_east',
                                   'z-dimension0021']}
        new_coords = {**new_coords,
                      'WRFrun': ['ctl', 'ytr'],
                      'lat': (('west_east', 'south_north'),
                              dict_runs[WRFdomain]['CLAT'][0, ...].values),
                      'lon': (('west_east', 'south_north'),
                              dict_runs[WRFdomain]['CLONG'][0, ...].values)}
        dict_runs[WRFdomain] = dict_runs[WRFdomain].assign_coords(new_coords)
        dict_runs[WRFdomain] = dict_runs[WRFdomain].rename(
            {'z-dimension0021': 'PFT'})
    return(dict_runs)
