import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs



with xr.open_dataset('C:\\Users\\ASUS\\Downloads\\G3P_v1.11\\IRAN_DATA_SWE_1km.nc') as ds:
    print("data ready")

print(1)
#convert monthly to yearly mean
SWE_Yearly = ds.groupby("time.year").mean("time")
#calculate linear yearly trend
trend = SWE_Yearly['swe'].polyfit(dim="year", deg=1)

fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=(18, 7))
#ax.set_extent([-10, 60, 30, 65])
# ax.add_feature(cartopy.feature.OCEAN, zorder=4)
ax.add_feature(cartopy.feature.BORDERS)
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS, zorder=1,  linewidth=0.3)
g = plt.pcolormesh(trend.lon, trend.lat, trend['polyfit_coefficients'][0, :, :],
                   cmap='bwr_r', shading='auto', transform=ccrs.PlateCarree())
cbar = plt.colorbar(g)
cbar.set_label('mm/year', fontsize=18)
g.set_clim(-1, 1)
plt.title('yearly trend 2002-2020', fontsize=15, style='italic', horizontalalignment='center')
#plt.suptitle('Land-Ocean fraction - 0.5-degree',  fontsize=19, horizontalalignment='center',  x=0.5) #y and x needed as you have adjusted the subplot size already.
# plt.savefig(path_figs + "/groundwater.png", dpi=1600, bbox_inches='tight')
plt.show()
