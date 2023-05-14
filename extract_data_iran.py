import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs

# Open the dataset
ds = xr.open_dataset('E:\\imerg-new\\imerg_2000_2021.nc')

# Define the bounding box coordinates
lat_min, lat_max = 40, 65
lon_min, lon_max = 20, 40

# Clip the dataset to the bounding box
ds_clip = ds.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))

pre_yearly = ds_clip.groupby("time.year").sum("time")
trend = pre_yearly['precipitationCal'].polyfit(dim="year", deg=1)

fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=(18, 7))
#ax.set_extent([44, 63, 24, 41])
# ax.add_feature(cartopy.feature.OCEAN, zorder=4)
ax.add_feature(cartopy.feature.BORDERS)
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS, zorder=1,  linewidth=0.3)
g = plt.pcolormesh(trend.lat, trend.lon, trend['polyfit_coefficients'][0, :, :],
                   cmap='bwr_r', shading='auto', transform=ccrs.PlateCarree())
cbar = plt.colorbar(g)
cbar.set_label('mm/year', fontsize=18)
g.set_clim(-50, 50)
plt.title('yearly trend 2000-2021', fontsize=15, style='italic', horizontalalignment='center')
#plt.suptitle('Land-Ocean fraction - 0.5-degree',  fontsize=19, horizontalalignment='center',  x=0.5) #y and x needed as you have adjusted the subplot size already.
# plt.savefig(path_figs + "/groundwater.png", dpi=1600, bbox_inches='tight')
plt.show()




def calculate_trend(da, var):
    """
    calculate linear trend

    Parameters
    ----------
    da:   dataset
    var:   variable name

    Returns:   The yearly trend
    -------

    """
    import xarray as xr
    WSC_Yearly = da.groupby("time.year").mean("time")
    # WSC_Yearly = da
    p = WSC_Yearly[var].polyfit(dim="year", deg=1)
    fit = xr.polyval(WSC_Yearly["year"], p.polyfit_coefficients)
    return p



