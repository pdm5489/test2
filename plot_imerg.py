import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import netCDF4 as nc
import numpy as np
import cartopy


# open the files and merge them into a single dataset
ds = xr.open_dataset('E:\\sattelite\\imerg-merging - version 1\\imerg_2000_2021.nc')
def wrap_lat(data, mode=""):
    if mode == 90:
        if data.lat.values[0] < 0:
            data = data.reindex(lat=list(reversed(data.lat)))
            print("latitude axis switched from -90 to 90 ----> 90 to -90")
        else:
            data = data
    if mode == -90:
        if data.lat.values[0] > 0:
            data = data.reindex(lat=list(reversed(data.lat)))
            print("latitude axis switched from 90 to -90 ----> -90 to 90")
        else:
            data = data

    return data

ds=wrap_lat(ds,mode=-90)
comp = dict(zlib=True, complevel=6)
encoding = {var: comp for var in ds.data_vars}
ds.to_netcdf("E:\\sattelite\\imerg-merging - version 1\\comp_transpos_imerg_2000_2021.nc", encoding=encoding, unlimited_dims='time')

p_year = ds.groupby("time.year").sum("time")
trend = p_year['precipitationCal'].polyfit(dim="year", deg=1)

print('salam')

fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=(18, 7))
#ax.set_extent([-10, 60, 30, 65])
# ax.add_feature(cartopy.feature.OCEAN, zorder=4)
ax.add_feature(cartopy.feature.BORDERS)
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS, zorder=1,  linewidth=0.3)
g = plt.pcolormesh(trend.lon, trend.lat, trend['polyfit_coefficients'][0, :, :],
                   cmap='bwr_r', shading='auto', transform=ccrs.PlateCarree())
cbar = plt.colorbar(g, extend='both')
cbar.set_label('mm/year', fontsize=18)
g.set_clim(-20, 20)
plt.title('yearly trend 2000-2021', fontsize=15, style='italic', horizontalalignment='center')
#plt.suptitle('Land-Ocean fraction - 0.5-degree',  fontsize=19, horizontalalignment='center',  x=0.5) #y and x needed as you have adjusted the subplot size already.
# plt.savefig(path_figs + "/groundwater.png", dpi=1600, bbox_inches='tight')
plt.show()