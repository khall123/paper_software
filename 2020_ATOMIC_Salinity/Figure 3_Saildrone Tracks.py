"""
Created Jul 02 2022, Author: Kashawn Hall
"""

import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.patheffects as pe

ds_1026 = xr.open_dataset('./data/saildrone-gen_5-atomic_eurec4a_2020-sd1026-20200117T000000-20200302T235959-1_minutes-v1.1589306725934.nc')
ds_1060 = xr.open_dataset('./data/saildrone-gen_5-atomic_eurec4a_2020-sd1060-20200117T000000-20200302T235959-1_minutes-v1.1589306886594.nc')
ds_1061 = xr.open_dataset('./data/saildrone-gen_5-atomic_eurec4a_2020-sd1061-20200117T000000-20200302T235959-1_minutes-v1.1589307121602.nc')
ds_1026 = ds_1026.isel(trajectory=0).swap_dims({'obs': 'time'})
ds_1060 = ds_1060.isel(trajectory=0).swap_dims({'obs': 'time'})
ds_1061 = ds_1061.isel(trajectory=0).swap_dims({'obs': 'time'})
ds_1026_start = ds_1026.sel(time=slice('2020-01-17T00:00:00'))
ds_1026_end = ds_1026.sel(time=slice('2020-03-02T23:59:00'))
ds_1060_start = ds_1060.sel(time=slice('2020-01-17T00:00:00'))
ds_1060_end = ds_1060.sel(time=slice('2020-03-02T23:59:00'))
ds_1061_start = ds_1061.sel(time=slice('2020-01-17T00:00:00'))
ds_1061_end = ds_1061.sel(time=slice('2020-03-02T23:59:00'))

gebco = xr.open_dataset("./data/gebco_2022_n14.3013_s5.9228_w-60.4674_e-47.5498.nc")
# print(gebco)
gebco_lon = gebco.lon.values
gebco_lat = gebco.lat.values
gebco_ele = gebco.elevation.values

dx, dy = 1.0, 1.5
x1, x2 = ds_1060.longitude.min().data - dx, ds_1060.longitude.max().data + dx
y1, y2 = ds_1060.latitude.min().data - dy, ds_1060.latitude.max().data + dy

'''Plotting'''
cmap = 'gray' #cmo.haline
norm = mpl.colors.Normalize(vmin=-5000, vmax=0, clip=True)
feature = cf.NaturalEarthFeature(name='land', category='physical', scale='10m', edgecolor='#000000', facecolor='#FFFFFF')

# Plotting SD Tracks
fig = plt.figure()
ax1 = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())

mesh = ax1.pcolormesh(gebco_lon, gebco_lat, gebco_ele, cmap=cmap, norm=norm)
ax1.plot(ds_1026.longitude, ds_1026.latitude, c='red', label='Saildrone 1026')
ax1.plot(ds_1060.longitude, ds_1060.latitude, c='gold', label='Saildrone 1060')
ax1.plot(ds_1061.longitude, ds_1061.latitude, c='blue', label='Saildrone 1061')
ax1.plot(ds_1026_start.longitude.data[-1], ds_1026_start.latitude.data[-1],  markersize=7, marker='o', color='red', label='Saildrone 1026 Start')
ax1.plot(ds_1060_start.longitude.data[-1], ds_1060_start.latitude.data[-1],  markersize=7, marker='o', color='gold', label='Saildrone 1060 Start')
ax1.plot(ds_1061_start.longitude.data[-1], ds_1061_start.latitude.data[-1],  markersize=7, marker='o', color='blue', label='Saildrone 1061 Start')
ax1.plot(ds_1026_end.longitude.data[-1], ds_1026_end.latitude.data[-1],  markersize=7, marker='D', color='red', label='Saildrone 1026 End')
ax1.plot(ds_1060_end.longitude.data[-1], ds_1060_end.latitude.data[-1],  markersize=7, marker='D', color='gold', label='Saildrone 1060 End')
ax1.plot(ds_1061_end.longitude.data[-1], ds_1061_end.latitude.data[-1],  markersize=7, marker='D', color='blue', label='Saildrone 1061 End')
ax1.set_extent([x1, x2, y1, y2])
ax1.legend(loc='upper right', ncol=3)
ax1.text(-59.35, 13.0939, 'Barbados', fontsize=12, weight='bold', color='white', path_effects=[pe.withStroke(linewidth=1, foreground="black")])
ax1.text(-49.5, 6.5, 'South\nAmerica', zorder=10, fontsize=12, weight='bold', color='white', path_effects=[pe.withStroke(linewidth=1, foreground="black")])
ax1.coastlines(resolution='10m')
ax1.add_feature(cf.LAND, zorder=100, edgecolor='k', facecolor='black')
ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
ax1.yaxis.set_major_locator(plt.NullLocator())
ax1.xaxis.set_major_formatter(plt.NullFormatter())

fig.subplots_adjust(bottom=0.2)
cbar_ax = fig.add_axes([0.2415, 0.14, 0.542, 0.02])
cb = fig.colorbar(mesh, cax=cbar_ax, orientation="horizontal")
cb.ax.set_xlabel('Seafloor Depth', fontsize=10, fontweight='semibold')
figManager = plt.get_current_fig_manager()
figManager.window.state('zoomed')
# plt.savefig('./data/Saildrone_Tracks_Gebco.png', bbox_inches='tight')
plt.show()
