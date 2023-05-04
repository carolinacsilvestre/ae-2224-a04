import matplotlib.colors
import numpy
import numpy as np
from RF_plots_season import path_chooser
import glob                                     # Dynamic file name loading
import numpy as np                              # Array processing
from netCDF4 import Dataset                     # NETCDF file handling
import tqdm
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, SpectralClustering
import os
import numpy as np
from mpl_toolkits.basemap import Basemap        # For map plotting
import matplotlib.cm as cm


def read_files(f_string, verbose=True):
    """
    It reads in the net fluxes for each of the 28 emission points (3 months) and returns a list of 28
    net fluxes

    :param f_string: the file path to the folder containing the files
    :param verbose: If True, prints out the file names that are being read in, defaults to True
    (optional)
    """
    # Read in file names based on f_string variable
    filenames_all = sorted(glob.glob(f_string))  # Get all file names in f_string
    if verbose:
        pbar = tqdm.tqdm(total=len(filenames_all))

    # Variables declaration
    time = []  # Attila
    plat = []  # Attila
    plon = []  # Attila
    ppress = []  # Attila in Pa

    for file in filenames_all:
        if verbose:
            pbar.update(1)
        if 'attila' in file:
            data = Dataset(file, 'r')
            # print('\n')
            # print('File loaded: ')
            # print(file)

            # Latitudes and longitudes
            lats = data.variables['lat'][:]
            lons_0to36 = data.variables['lon'][:]  # Varies from 0 to 360 deg
            lons_18to18 = data.variables['lon'][:]  # Will be converted later to -180 to 180 deg

            # Time
            temp = data.variables['time'][:]
            time.append(temp)

            # Air parcel longitudinal position
            temp = data.variables['PLON'][:]
            plon.append(temp)

            # Air parcel latitudinal position
            temp = data.variables['PLAT'][:]
            plat.append(temp)

            # Air parcel pressure altitude
            temp = data.variables['PPRESS'][:]
            ppress.append(temp)

    # Concatenation of variables, lists become multi-dimensional numpy arrays
    time = np.concatenate(time, axis=0)
    plon = np.concatenate(plon, axis=0)
    plat = np.concatenate(plat, axis=0)
    ppress = np.concatenate(ppress, axis=0) / 100  # Convert to hPa

    # Convert longitude range from [0, 360] to [-180, 180]
    plon = (plon + 180) % 360 - 180  # Air parcel longitudinal coordinates
    lons_18to18 = (lons_18to18 + 180) % 360 - 180  # EMAC grid longitudes
    if verbose:
        pbar.close()
    return ppress,time

season = "winter"
altitude = 300
# Below are the quick examples
f_string = path_chooser(season, altitude)

# construct median for every EP
altitude_list, time_list = read_files(f_string)
median_list = []
for EP in range(1,29):
    median = np.median(altitude_list[:,(EP-1)*50:EP*50],axis=1)
    median_list.append(median)
median_list = np.array(median_list)

def clustering(inputdata,n_clusters=6):
    kmeans = KMeans(n_clusters=n_clusters)
    prediction = kmeans.fit_predict(inputdata)
    return prediction

results = clustering(median_list,3)
#make subplot
fig, axs = plt.subplots(1, 2, figsize=(10, 5))  # Create a figure with two subplots
fig.suptitle('Altitude Median vs Time and KMeans Clustering Results')  # Add a title to the figure

import matplotlib.colors as mcolor
color_map = {i: cm.tab20(i) for i in range(20)}

color = []
for i in range(28):
    color.append(mcolor.to_hex(color_map.get(results[i])))
color_grid = numpy.reshape(color,(7, 4),order="F")

cluster_list = np.array(results).reshape(7, 4,order="F")
lat = np.linspace(85, 25, 7)  # define x as an array with 4 elements
lon = np.linspace(-115, -55, 4)  # define y as an array with 7 elements

# Define map projection and settings
# For more info: https://matplotlib.org/basemap/users/cyl.html
mp = Basemap(projection='cyl',  # equidistant cylindrical projection
             llcrnrlon=-135,
             llcrnrlat=10,
             urcrnrlon=-35,
             urcrnrlat=90,
             resolution='i', ax=axs[0])  # h=high, f=full, i=intermediate, c=crude

# Shift the fluxes from [0,360] to [-180,180]
# net_flx_EP_shft, lons_shft = shiftgrid(180.,global_net_flx[1], lons_0to36,start=False)

# Format the lat and lon arrays for map graphing,
# makes lat array a lat x lon array and same for lon array
# lon, lat = np.meshgrid(lons_shft, lats)
x, y = mp(lon, lat)

# Choose the settings for the coastlines, countries, meridians...
mp.drawcoastlines(linewidth=0.2)
mp.drawcountries(linewidth=0.2)

# Draw lon lines every 20º
meridians = mp.drawmeridians(np.arange(-180, 200, 20), labels=[False, False, False, True], linewidth=0.2, fontsize=10)
# Draw lat lines every 20º
mp.drawparallels(np.arange(-90, 110, 20), labels=[True, False, False, True], linewidth=0.2, fontsize=10)

# Plot the flux on the map
for n in range(4):
    for m in range(7):
        mp.scatter(x[n], y[m],s=100,color=color_grid[m,n])


for i in range(28):
    # Plot the KMeans clustering results in the second subplot
    axs[1].plot(time_list, median_list[i], color=color[i],label="EP"+str(i+1))
axs[1].set_title('KMeans Clustering Results')  # Add a title to the second subplot

os.chdir('C:/Users/twanv/OneDrive - Delft University of Technology/Bsc - 2/Q3/Project/plots')

plt.savefig("KMeans"+season+str(altitude)+".png", dpi= 500)
plt.show()  # Display the figure"""

""""
# Assume median_list is a NumPy array of shape (n, m)
n, m = median_list.shape
print(n,m)
# Convert median_list to a list of streamlines, where each streamline is a list of 3D coordinates
streamlines = [median_list[i, :].reshape(-1, 1) for i in range(n)]
print(len(streamlines[0]),len(streamlines))
streamlines = Streamlines(streamlines)

# Print the resulting streamline


# Define the QuickBundles algorithm with a given threshold
qb = QuickBundles(threshold=0.1,max_nb_clusters=10)

# Cluster the streamlines
clusters = qb.cluster(streamlines)
print(clusters)




print("Nb. clusters:", len(clusters))
print("Cluster sizes:", map(len, clusters))
print("Small clusters:", clusters < 2)
print("Streamlines indices of the first cluster:\n", clusters[0].indices)
print("Centroid of the last cluster:\n", clusters[-1].centroid)

# Plot the altitude median vs time in the first subplot
for median in median_list:
    axs[0].plot(time_list, median)
axs[0].set_title('Altitude Median vs Time')  # Add a title to the first subplot

import matplotlib.cm as cm

color_map = {i: cm.tab20(i) for i in range(20)}

for i in range(28):
    color = color_map.get(results[i], cm.tab20(0))
# Plot the KMeans clustering results in the second subplot
    axs[1].plot(time_list, median_list[i], color=color,label="EP"+str(i+1))
axs[1].set_title('QuickBundels Clustering Results')  # Add a title to the second subplot
axs[1].legend()

os.chdir('C:/Users/twanv/OneDrive - Delft University of Technology/Bsc - 2/Q3/Project/plots')

plt.savefig("KMeans"+season+str(altitude)+".png", dpi= 500)
plt.show()  # Display the figure"""



