conda install -c conda-forge xarray dask netCDF4 bottleneck
---------------------------
downscale
import os

# Define the input and output file names
input_file = 'C:\\Users\\ASUS\\Downloads\\G3P_v1.11\\downscale\\swe_2020_Ano_0.5_degree_Gauss_200km_land_masked_GravIS.nc'
output_file = 'C:\\Users\\ASUS\\Downloads\\G3P_v1.11\\downscale\\swe_2020_Ano_0.1_degree_Gauss_200km_land_masked_GravIS.nc'

# Define the command to downscale the data
command = f'ncks -d lat,,,2 -d lon,,,2 {input_file} {output_file}'

# Run the command
os.system(command)
------------------------------
merging
import glob
ds = xr.merge([xr.open_dataset(f) for f in glob.glob('C:\\Users\\ASUS\\Downloads\\G3P_v1.11\\SWE\\*.nc')])
-----------------------------
plot trend
import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs


# open the files and merge them into a single dataset
ds = xr.open_mfdataset("C:\\Users\\ASUS\\Downloads\\G3P_v1.11\\SWE\\*.nc")

# select the variable
x = ds["swe"]

# calculate the mean over time
x_mean = x.mean(dim="time")


fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.Robinson()), figsize=(18, 7))
#ax.set_extent([44, 63, 24, 41])     #iran cordinate
# ax.add_feature(cartopy.feature.OCEAN, zorder=4)
ax.add_feature(cartopy.feature.BORDERS)
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS, zorder=1,  linewidth=0.3)

plt.title('swe', fontsize=15, style='italic', horizontalalignment='center')

# plot the trend
x_mean.plot()


# show the plot
plt.show()
-------------------------------------

extract mask iran cordinate

# Load the NetCDF file
data = xr.open_dataset('C:\\Users\\ASUS\\Downloads\\G3P_v1.11\\downscale\\swe_2020_Ano_0.5_degree_Gauss_200km_land_masked_GravIS.nc')

# Extract data at specific coordinates
ds_iran = data.sel(lat=slice(25.0, 40.0), lon=slice(44.0, 63.0))

--------------------------------------
out to nc file
ds = xr.open_dataset('C:\\Users\\ASUS\\Downloads\\G3P_v1.11\\downscale\\swe_2020_Ano_0.5_degree_Gauss_200km_land_masked_GravIS.nc')
ds_out = ds.interp(lat=720, lon=1440)
ds_out.to_netcdf('C:\\Users\\ASUS\\Downloads\\G3P_v1.11\\downscale\\swe_2020_downscale.nc')
print(ds_out)
--------------------------------
nimbus-7 code plot
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import netCDF4 as nc
import numpy as np

nc_file = nc.Dataset(
    'C:\\Users\\ASUS\\Downloads\\nimbus-7\\GlobSnow_v3.0_monthly_SWE\\201805_northern_hemisphere_monthly_swe_0.25grid.nc')
x = nc_file.variables['x'][:].data
y = nc_file.variables['y'][:].data
z = nc_file.variables['swe'][:].data

fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.Robinson()), figsize=(18, 7))
ax.set_extent([44, 63, 24, 41])     #iran cordinate
# ax.add_feature(cartopy.feature.OCEAN, zorder=4)
ax.coastlines()
ax.add_feature(cfeature.LAND,
               edgecolor='black',facecolor="none",linewidth=1)
#ax.add_feature(cartopy.feature.BORDERS, zorder=1,  linewidth=0.3)
g = plt.pcolormesh(y, x, z,
                   shading='auto')
cbar = plt.colorbar(g)
plt.title('swe', fontsize=15, style='italic', horizontalalignment='center')
#plt.suptitle('Land-Ocean fraction - 0.5-degree',  fontsize=19, horizontalalignment='center',  x=0.5) #y and x needed as you have adjusted the subplot size already.
# plt.savefig(path_figs + "/groundwater.png", dpi=1600, bbox_inches='tight')
plt.show()
---------------------------------------------------
nasa
from subprocess import Popen
from getpass import getpass
import platform
import os
import shutil

urs = 'urs.earthdata.nasa.gov'    # Earthdata URL to call for authentication
prompts = ['Enter NASA Earthdata Login Username \n(or create an account at urs.earthdata.nasa.gov): ',
           'Enter NASA Earthdata Login Password: ']

homeDir = os.path.expanduser("~") + os.sep

with open(homeDir + '.netrc', 'w') as file:
    file.write('machine {} login {} password {}'.format(urs, getpass(prompt=prompts[0]), getpass(prompt=prompts[1])))
    file.close()
with open(homeDir + '.urs_cookies', 'w') as file:
    file.write('')
    file.close()
with open(homeDir + '.dodsrc', 'w') as file:
    file.write('HTTP.COOKIEJAR={}.urs_cookies\n'.format(homeDir))
    file.write('HTTP.NETRC={}.netrc'.format(homeDir))
    file.close()

print('Saved .netrc, .urs_cookies, and .dodsrc to:', homeDir)

# Set appropriate permissions for Linux/macOS
if platform.system() != "Windows":
    Popen('chmod og-rw ~/.netrc', shell=True)
else:
    # Copy dodsrc to working directory in Windows
    shutil.copy2(homeDir + '.dodsrc', os.getcwd())
    print('Copied .dodsrc to:', os.getcwd())
-------------------------------
extract iran
with xr.open_dataset('E:\\imerg-new\\imerg_2000_2021.nc') as ds:
    lon = ds.lon
    lat = ds["lat"]
lat_bnds, lon_bnds = [80, 10], [10, 80]
iran = ds.sel(lat=slice(*lat_bnds), lon=slice(*lon_bnds))
pre_yearly = iran.groupby("time.year").sum("time")
trend = pre_yearly['precipitationCal'].polyfit(dim="year", deg=1)
--------------------------
extract data from netcdf file
import netCDF4
import pandas as pd
import openpyxl

# Open the NetCDF file
nc = netCDF4.Dataset('E:\\sattelite\\era5-land-snow and temp\\compress_era5_land_snow_temp_2000_2022.nc')

# Extract variable of interest
var = nc.variables['latitude']

# Convert data to a Pandas DataFrame
df = pd.DataFrame(var[:])

# Save data as an Excel file
df.to_excel('E:\\sattelite\\era5-land-snow and temp\\latitude_compress_era5_land_snow_temp_2000_2022.xlsx', index=False)

# Save data as a text file
#df.to_csv('output.txt', index=False, sep='\t')
------------------
transpos and compress
ds = xr.open_dataset('E:\\sattelite\\imerg-merging - version 1\\imerg_2000_2021.nc')
ds = ds.transpose('time', 'bnds', 'lat', 'lon')
comp = dict(zlib=True, complevel=6)
encoding = {var: comp for var in ds.data_vars}
ds.to_netcdf("E:\\sattelite\\imerg-merging - version 1\\comp_transpos_imerg_2000_2021.nc", encoding=encoding, unlimited_dims='time')
-----------
fix coordinate
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
