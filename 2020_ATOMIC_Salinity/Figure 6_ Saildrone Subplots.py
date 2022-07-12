"""
Created Mar 29 2021, Author: Kashawn Hall
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmocean.cm as cmo
from glob import glob

'''HYCOM Data Pull'''
data_dir1 = './data/HYCOM/'
files = glob(data_dir1+'*nc4')
hycom = xr.open_mfdataset(data_dir1+'*.nc4', concat_dim='time').isel(depth=0)

#this next section removes duplicate timesteps
_, index = np.unique(hycom['time'], return_index=True)
hycom2 = hycom.isel(time=index)

#change saildrone coordinates to match with hycom (0-359:-180-179)
hycom_lon = hycom2.assign_coords(longitude=(((hycom2.lon + 180) % 360) - 180))
hycom2 = hycom_lon.swap_dims({'lon': 'longitude'})

#remove nans from hycom data
filled = hycom2.chunk({'time': -1}).interpolate_na(dim="time", method="nearest", fill_value='extrapolate')
filled2 = filled.interpolate_na(dim="lat", method="nearest", fill_value='extrapolate')
filled3 = filled2.interpolate_na(dim="lon", method="nearest", fill_value='extrapolate')
filled3_clip = filled3.sel(time=slice('2020-02-14T00:00:00.000000000', '2020-02-21T00:00:00.000000000'))
HYCOM_avg = filled3_clip.mean(dim='time')
HYCOM_8avg = filled3.sel(time=('2020-02-17T12:00:00.000000000'))
# print(HYCOM_8avg)


'''Saildrone Data Pull'''
ds_1026 = xr.open_dataset('./data/saildrone-gen_5-atomic_eurec4a_2020-sd1026-20200117T000000-20200302T235959-1_minutes-v1.1595708344687.nc')  # Importing SD 1026
ds_1026 = ds_1026.isel(trajectory=0).swap_dims({'obs': 'time'})  # Switching dimensions from "obs" to "time"
ds_1026_dot = ds_1026.sel(time=slice('2020-02-18T00:00:00'))  # Slicing SD1026 time to showcase the SD in the middle of the fresh tongue
ds_1026_clip = ds_1026.sel(time=slice('2020-02-16T00:00:00', '2020-02-20T00:00:00.000000000'))  # Slicing SD1026 time to showcase the SD traversing the fresh tongue


'''SMAP Data Pull'''
SMAP_JPL = xr.open_dataset('./data/JPL_SMAP_8DAYS_V5.0.nc')  # Importing JPL 8 day average
SMAP_RSS = xr.open_dataset('./data/RSS_SMAP_8day_running_v04.0.nc')  # Importing RSS 8 day average
# print(SMAP_JPL)

'''Extents'''
dx, dy = 3.05, 3.05
x1, x2 = ds_1026.longitude.min().data - dx, ds_1026.longitude.max().data + dx  # Setting X extents based on the max SD extents
y1, y2 = ds_1026.latitude.min().data - dy, ds_1026.latitude.max().data + dy  # Setting Y extents based on the max SD extents

'''Clipping the HYCOM data to the extents to limit the color bar range'''
HYCOM_clipped = HYCOM_8avg.where((filled3_clip.lat >= y1) & (filled3_clip.lat <= y2) & (filled3_clip.longitude >= x1) & (filled3_clip.longitude <= x2), drop=True)
# print(HYCOM_clipped)
sss_HYCOM = HYCOM_clipped.salinity.values  # Setting the HYCOM Salinity variable
lat_HYCOM = HYCOM_clipped.lat.data  # Setting the HYCOM Latitude variable
lon_HYCOM = HYCOM_clipped.longitude.data  # Setting the HYCOM Longitude variable

'''Clipping the SMAP data to the extents to limit the color bar range'''
SMAP_JPL_clipped = SMAP_JPL.where((SMAP_JPL.latitude >= y1-0.125) & (SMAP_JPL.latitude <= y2+0.125) & (SMAP_JPL.longitude >= x1-0.125) & (SMAP_JPL.longitude <= x2+0.250), drop=True)
# print(SMAP_JPL_clipped)
sss_JPL = SMAP_JPL_clipped.smap_sss.values  # Setting the JPL Salinity variable
anc_JPL = SMAP_JPL_clipped.anc_sss.values  # Setting the JPL HYCOM Salinity variable
lat_JPL = SMAP_JPL_clipped.latitude.data  # Setting the JPL Latitude variable
lon_JPL = SMAP_JPL_clipped.longitude.data  # Setting the JPL Longitude variable

SMAP_RSS_clipped = SMAP_RSS.where((SMAP_RSS.lat >= y1-0.125) & (SMAP_RSS.lat <= y2+0.125) & (SMAP_RSS.lon >= 298.625-1.125) & (SMAP_RSS.lon <= 313.375+1.250), drop=True)
# print(SMAP_RSS_clipped)
sss_RSS = SMAP_RSS_clipped.sss_smap.values  # Setting the RSS Salinity variable
sss_RSS_40km = SMAP_RSS_clipped.sss_smap_40km.values  # Setting the RSS 40km Salinity variable
lat_RSS = SMAP_RSS_clipped.lat.data   # Setting the RSS Latitude variable
lon_RSS = SMAP_RSS_clipped.lon.data  # Setting the RSS Longitude variable
# print(len(lon_RSS), len(lat_RSS))

'''Plotting'''
cmap = cmo.haline  # Colormap choice
norm = mpl.colors.Normalize(vmin=31.85, vmax=37, clip=False)  # Normalizing colorbar between all subplots
feature = cf.NaturalEarthFeature(name='land', category='physical', scale='10m', edgecolor='#000000', facecolor='#FFFFFF')  # Import land features
v = np.linspace(30.15, 38.28, 9, endpoint=True)
v = [round(num, 2) for num in v]

# Plotting RSS SMAP Salinity
fig = plt.figure()
plt.title('02-17-2020 8-day Averages', weight='bold')
ax1 = plt.subplot(2, 2, 1, projection=ccrs.PlateCarree())
ax1.set_extent([-58.6, -54.6, 9, 11.5])
fd1 = ax1.contourf(lon_RSS, lat_RSS, sss_RSS, 60, cmap=cmap, norm=norm)
ax1.scatter(ds_1026_clip.longitude, ds_1026_clip.latitude, c=ds_1026_clip.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(), cmap=cmap, label='Saildrone 1026', norm=norm)
dot = plt.plot(ds_1026_dot.longitude.data[-1], ds_1026_dot.latitude.data[-1],  markersize=5, marker='o', color='white')
ax1.set_xticks([-58.6, -56.6, -54.6], crs=ccrs.PlateCarree())
ax1.set_yticks([9, 10.25, 11.5], crs=ccrs.PlateCarree())
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)
ax1.text(-55.5, 11.3, '(a) RSS70', fontsize=10, weight='bold', color='black')
ax1.coastlines(resolution='10m')

# Plotting RSS 40km SMAP Salinity
ax2 = plt.subplot(2, 2, 2, projection=ccrs.PlateCarree())
ax2.set_extent([-58.6, -54.6, 9, 11.5])
fd2 = ax2.contourf(lon_RSS, lat_RSS, sss_RSS_40km, 60, cmap=cmap, norm=norm)
ax2.scatter(ds_1026_clip.longitude, ds_1026_clip.latitude, c=ds_1026_clip.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(), cmap=cmap, label='Saildrone 1026', norm=norm)
dot = plt.plot(ds_1026_dot.longitude.data[-1], ds_1026_dot.latitude.data[-1],  markersize=5, marker='o', color='white')
ax2.set_xticks([-58.6, -56.6, -54.6], crs=ccrs.PlateCarree())
ax2.set_yticks([9, 10.25, 11.5], crs=ccrs.PlateCarree())
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax2.xaxis.set_major_formatter(lon_formatter)
ax2.yaxis.set_major_formatter(lat_formatter)
ax2.text(-55.5, 11.3, '(b) RSS40', fontsize=10, weight='bold')
ax2.coastlines(resolution='10m')

# Plotting JPL SMAP Salinity
ax3 = plt.subplot(2, 2, 3, projection=ccrs.PlateCarree())
ax3.set_extent([-58.6, -54.6, 9, 11.5])
fd3 = ax3.contourf(lon_JPL, lat_JPL, sss_JPL, 60, cmap=cmap, norm=norm)
ax3.scatter(ds_1026_clip.longitude, ds_1026_clip.latitude, c=ds_1026_clip.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(), cmap=cmap, label='Saildrone 1026', norm=norm)
dot = plt.plot(ds_1026_dot.longitude.data[-1], ds_1026_dot.latitude.data[-1],  markersize=5, marker='o', color='white')
ax3.set_xticks([-58.6, -56.6, -54.6], crs=ccrs.PlateCarree())
ax3.set_yticks([9, 10.25, 11.5], crs=ccrs.PlateCarree())
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax3.xaxis.set_major_formatter(lon_formatter)
ax3.yaxis.set_major_formatter(lat_formatter)
ax3.text(-55.25, 11.3, '(b) JPL', fontsize=10, weight='bold')
ax3.coastlines(resolution='10m')

#Plotting HYCOM Salinity
ax4 = plt.subplot(2, 2, 4, projection=ccrs.PlateCarree())
ax4.add_feature(feature)
ax4.set_extent([-58.6, -54.6, 9, 11.5])
ax4.pcolormesh(lon_HYCOM, lat_HYCOM, sss_HYCOM, cmap=cmap, norm=norm)
ax4.scatter(ds_1026_clip.longitude, ds_1026_clip.latitude, c=ds_1026_clip.SAL_SBE37_MEAN, s=.15, transform=ccrs.PlateCarree(), cmap=cmap, label='Saildrone 1026', norm=norm)
dot = plt.plot(ds_1026_dot.longitude.data[-1], ds_1026_dot.latitude.data[-1],  markersize=5, marker='o', color='white')
ax4.set_xticks([-58.6, -56.6, -54.6], crs=ccrs.PlateCarree())
ax4.set_yticks([9, 10.25, 11.5], crs=ccrs.PlateCarree())
plt.xticks(fontweight='semibold')
plt.yticks(fontweight='semibold')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax4.xaxis.set_major_formatter(lon_formatter)
ax4.yaxis.set_major_formatter(lat_formatter)
ax4.text(-55.55, 11.3, '(d) HYCOM', fontsize=10, weight='bold')
ax4.coastlines(resolution='10m')


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.8, 0.11, 0.02, 0.77])  # Set Location of colorbar
cb = fig.colorbar(fd3, cax=cbar_ax, ticks=v, orientation="vertical")  # Create colorbar
cb.ax.set_ylabel('Salinity (PSU)', fontsize=10, fontweight='bold')  # Set Label of colorbar
cb.ax.set_yticklabels(labels=cb.ax.get_yticklabels(), weight='bold')
figManager = plt.get_current_fig_manager()
figManager.window.state('zoomed')
# plt.savefig('./data/figures/Satellite_subplots_boxes_zoomed.png', bbox_inches='tight')
plt.show()
