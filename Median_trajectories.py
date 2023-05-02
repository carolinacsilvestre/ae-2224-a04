import numpy as np
from RF_plots_season import path_chooser
import glob                                     # Dynamic file name loading
import numpy as np                              # Array processing
from netCDF4 import Dataset                     # NETCDF file handling
import tqdm
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

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


# Below are the quick examples
f_string = path_chooser("summer", 250)

altitude_list, time_list = read_files(f_string)
median_list = []
for EP in range(1,29):
    median = np.median(altitude_list[:,(EP-1)*50:EP*50],axis=1)
    #plt.plot(time_list, median)
    median_list.append(median)
#plt.show()
median_list = np.array(median_list)
print(median_list.shape)



color_map = {0: 'blue', 1: 'red', 2: 'green', 3: 'black', 4: 'yellow'}
"""
for i in range(28):
    color = color_map.get(results[i], 'blue')
    plt.plot(time_list,median_list[i], color=color)
plt.show()"""

kmeans = KMeans(n_clusters=5)
kmeans.fit(median_list)
results = kmeans.predict(median_list)
print(results)

fig, axs = plt.subplots(1, 2, figsize=(10, 5))  # Create a figure with two subplots
fig.suptitle('Altitude Median vs Time and KMeans Clustering Results')  # Add a title to the figure

# Plot the altitude median vs time in the first subplot
for median in median_list:
    axs[0].plot(time_list, median)
axs[0].set_title('Altitude Median vs Time')  # Add a title to the first subplot

# Plot the KMeans clustering results in the second subplot
for i in range(28):
    color = color_map.get(results[i], 'blue')
    axs[1].plot(time_list, median_list[i], color=color)
axs[1].set_title('KMeans Clustering Results')  # Add a title to the second subplot

plt.show()  # Display the figure





