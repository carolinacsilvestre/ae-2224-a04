'''Taking 1 emission point. For each of the 50 air parcels: define time window (first 45d), finding rate of descent,
   average mixing ratio, plot average mixing ratio against rate of descent. Do statistical analysis'''

'''Sample script for reading in data and visualizing the
Lagrangian air parcel trajectories along with the net radiative fluxes
arising from a short-term increase in ozone.'''

# Import libraries

# =============================================================================
# LOADING DATA FROM NETCDF (.NC) FILES
# =============================================================================

# USER INPUT - File path
# f_string = 'C:/Users/Carolina Silvestre\Desktop\dataproject*' #Insert file path to input data, do not forget wildcard

import scipy.stats
from matplotlib.animation import FuncAnimation
import glob  # Dynamic file name loading
import numpy as np  # Array processing
from netCDF4 import Dataset  # NETCDF file handling
import matplotlib.pyplot as plt  # Plotting
from mpl_toolkits.basemap import Basemap  # For map plotting
from mpl_toolkits.basemap import shiftgrid  # For shifting longitudes
import matplotlib.colors  # To create new colorbar
import math

load_data_from = str('Winter300')

f_string = 'C:/Users/alexm/AE2224/DATA_ANALYSIS/' + load_data_from + '/*' 

if len(sorted(glob.glob(f_string))) == 0:
    f_string = 'C:/Users/Carolina Silvestre/Desktop/dataproject/' + load_data_from + '/*'
if len(sorted(glob.glob(f_string))) == 0:
    f_string = 'D:/Python safe/all test dat/*'
if len(sorted(glob.glob(f_string))) == 0:
    f_string = 'E:/all data/' + load_data_from + '/*'
if len(sorted(glob.glob(f_string))) == 0:
    f_string = 'D:/Python safe/Rest test data/'+ load_data_from +'/*'

# USER INPUT - Switches to determine which data types should be loaded
attila_switch = True
o3tracer_switch = True
rad_fluxes_switch = False

# Read in file names based on f_string variable
filenames_all = sorted(glob.glob(f_string))  # Get all file names in f_string
print('\n')  # Move to next line, improve readability
print('Files in input folder: ')
print('\n')  # Check that all files present in folder are being printed
print(filenames_all)

# Variables declaration
time = []  # Attila
plat = []  # Attila
plon = []  # Attila
ppress = []  # Attila in Pa
airO3_001 = []  # O3lg in mol/mol for emission point 1
rad_flx_SW = []  # Holds 1 specific call for SW radiative fluxes
rad_flx_LW = []  # Holds 1 specific call for LW radiative fluxes
rad_flx_SW_02 = []  # SW flux for call 02
rad_flx_LW_02 = []  # LW flux for call 02
global_net_flx = []  # Holds all net fluxes for the 28 EPs (3 months)

# Load relevant variables from NETCDF files
# Positions of air parcels
if attila_switch == True:
    for file in filenames_all:
        if 'attila' in file:
            data = Dataset(file, 'r')
            print('\n')
            print('File loaded: ')
            print(file)

            #Latitudes and longitudes
            lats = data.variables['lat'][:]
            lons_0to36 = data.variables['lon'][:]  # Varies from 0 to 360 deg
            # Will be converted later to -180 to 180 deg
            lons_18to18 = data.variables['lon'][:]

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
    ppress = np.concatenate(ppress, axis=0)/100  # Convert to hPa
    #
    # Convert longitude range from [0, 360] to [-180, 180]
    plon = (plon + 180) % 360 - 180  # Air parcel longitudinal coordinates
    lons_18to18 = (lons_18to18 + 180) % 360 - 180  # EMAC grid longitudes

# O3 data along air parcel trajectories
if o3tracer_switch == True:
    for file in filenames_all:
        if 'O3lg' in file:
            data = Dataset(file, 'r')
            print('\n')
            print('File loaded: ')
            print(file)

            # O3 along air parcels for each emission point
            # Only the first 1400 trajectories
            temp = data.variables['airO3_001'][:, :1400]
            airO3_001.append(temp)

            # NOTE: do not forget there are 28 emission points in total!

    # Concatenation of variables, lists become multi-dimensional numpy arrays
    airO3_001 = np.concatenate(airO3_001, axis=0)

# Radiative fluxes corresponding to O3 increase, specifically for 250hPa
if rad_fluxes_switch == True:
    for file in filenames_all:

        # Get call 2 separately, which will be subtracted from calls 3-30
        if file.find('fluxes_tp') != -1 and file.find('EP02') != -1:
            # fluxes_tp_NAmerica_July2014_EP02
            data = Dataset(file, 'r')
            # print('\n')
            # print('File loaded: ')
            # print(file)

            # Radiative fluxes
            rad_flx_SW_02 = data.variables['flxs_tp'][:]
            rad_flx_LW_02 = data.variables['flxt_tp'][:]

        # Start at 3 since first two calls are not emission points
        for ep in range(3, 31):
            # If there is no match, output is -1.
            if file.find('fluxes_tp') != -1 and file.find('EP'+str(ep).zfill(2)) != -1:
                data = Dataset(file, 'r')
                # print('\n')
                # print('File loaded: ')
                # print(file)

                # Radiative fluxes
                rad_flx_SW = data.variables['flxs_tp'][:]
                rad_flx_LW = data.variables['flxt_tp'][:]

                # Calculate net flux and append
                global_net_flx.append((rad_flx_LW+rad_flx_SW) -
                                      (rad_flx_LW_02+rad_flx_SW_02))

    # Delete unnecessary variables
    del rad_flx_LW, rad_flx_SW, rad_flx_LW_02, rad_flx_SW_02

# Fluxes from the VISO submodel, for 200 and 300 hPa
# if rad_fluxes_switch == True:
#     #Start at 3 since first two calls are not emission points
#     for ep in range(3,31):
#         for file in filenames_all:
#             if 'viso' in file:
#                 data = Dataset(file,'r')
#                 print('\n')
#                 print('File loaded: ')
#                 print(file)
#
#                 #Store call 2 containing only background O3 fluxes
#                 temp = data.variables['flxs_tp_02'][:]
#                 rad_flx_SW_02.append(temp)
#
#                 temp = data.variables['flxt_tp_02'][:]
#                 rad_flx_LW_02.append(temp)
#
#                 #Other emission point
#                 temp = data.variables['flxs_tp_'+str(ep).zfill(2)][:]
#                 rad_flx_SW.append(temp)
#
#                 temp = data.variables['flxt_tp_'+str(ep).zfill(2)][:]
#                 rad_flx_LW.append(temp)
#
#         #Concatenate the SW and LW flux data for all 3 months
#         rad_flx_SW = np.concatenate(rad_flx_SW, axis=0)
#         rad_flx_LW = np.concatenate(rad_flx_LW, axis=0)
#         rad_flx_SW_02 = np.concatenate(rad_flx_SW_02, axis=0)
#         rad_flx_LW_02 = np.concatenate(rad_flx_LW_02, axis=0)
#
#         #Calculate and append net fluxes for each of the 28 emission points
#         global_net_flx.append((rad_flx_LW+rad_flx_SW)-(rad_flx_LW_02+rad_flx_SW_02))
#
#         #Clear the radiative flux list holding the 3-month period data for a given call
#         rad_flx_SW = []
#         rad_flx_LW = []
#         rad_flx_SW_02 = []
#         rad_flx_LW_02 = []
#
#         print('Net flux for EP'+str(ep-2)+' loaded...')
#
#     del rad_flx_LW, rad_flx_SW, rad_flx_LW_02, rad_flx_SW_02
#
# #Delete unnecessary variables
# if attila_switch or o3tracer_switch or rad_fluxes_switch and not '250' in f_string:
#     del temp, data


def scatterplot(x,y, nrow, ncol):
    ax, fig = plt.subplots(nrows = nrow, ncols = ncol)

    if nrow == 1 and ncol == 1:

        plt.scatter(x,y)

    elif nrow > 1 or ncol > 1:
        r = 0
        c = 0
        while c <= ncol:
            ax[r,c].scatter(x,y)
            c = c + 1
            if r == nrow:
                break
        if c == ncol:
            r = r + 1
            c = 0

    plt.show()

    return 0

# =============================================================================
# PLOT TYPE 1 - VERTICAL EVOLUTION OF LAGRANGIAN AIR PARCELS
# =============================================================================

# Requires ATTILA air parcel trajectory location data

activate_plot1 = False  # Activation of vertical position plot##


if attila_switch == True and activate_plot1 == True:

    # Set up axis object for plotting the map
    fig, ax = plt.subplots()

    # Adjust dimensions of map plot
    fig.set_figheight(8)
    fig.set_figwidth(14)

    parcel = 29  # Parcel ID, 0 means first.

    # Scatter plot settings
    plt.scatter(time, ppress[:, parcel], s=30, marker='o', color='green')
    # print('1', ppress[:, parcel])
    # Defining range of axes
    plt.xticks(np.arange(0, 110, 10), np.arange(0, 110, 10), fontsize=20)
    plt.yticks(np.arange(0, 1200, 200), np.arange(0, 1200, 200), fontsize=20)

    # Inversion of vertical axis
    plt.gca().invert_yaxis()  # The higher the pressure, the closer to the surface

    # Labeling of axes
    plt.xlabel(
        'Time elapsed since emission \n [Days]', fontsize=22, weight='bold')
    plt.ylabel(
        'Air parcel pressure altitude \n [hPa]', fontsize=22, weight='bold')

    # Add plot title
    plt.title('Air Parcel ID='+str(parcel), fontsize=26, weight='bold')

    # To save the plot, give it a name, format and dpi represents the resolution.
    # Typically, for publications, dpi=300.
    plt.tight_layout()  # Ensure all parts of the plot will show after saving
    plt.savefig("air_parcel_ID"+str(parcel)+".png", format="png", dpi=300)

    plt.show()
    plt.close()


# =============================================================================
# PLOT TYPE 2 - HORIZONTAL EVOLUTION OF LAGRANGIAN AIR PARCELS (ON MAP)
# =============================================================================

activate_plot2 = False  # Activation of trajectory plot ##

# Requires ATTILA air parcel trajectory location data
if attila_switch == True and activate_plot2 == True:

    parcel2 = list(range(100, 140))  # Parcel ID, 0 means first.
    # print('parcel2', parcel2)
    # Set up axis object for plotting the map
    fig, ax = plt.subplots()  # Subplots are useful for drawing multiple plots together

    # Adjust dimensions of map plot
    fig.set_figheight(8)
    fig.set_figwidth(14)

    # Define map projection and settings
    # For more info: https://matplotlib.org/basemap/users/cyl.html
    mp = Basemap(projection='cyl',  # equidistant cylindrical projection
                 llcrnrlon=-180,
                 llcrnrlat=-90,
                 urcrnrlon=180,
                 urcrnrlat=90,
                 resolution='i', ax=ax)  # h=high, f=full, i=intermediate, c=crude

    # Format the lat and lon arrays for map graphing,
    # makes lat array a lat x lon array and same for lon array
    lon, lat = np.meshgrid(lons_18to18, lats)
    x, y = mp(lon, lat)

    # Choose the settings for the coastlines, countries, meridians...
    mp.drawcoastlines(linewidth=0.2)
    mp.drawcountries(linewidth=0.2)

    meridians = mp.drawmeridians(np.arange(-180, 200, 20),
                                 labels=[False, False, False, True],
                                 linewidth=0.2, fontsize=10)  # Draw lon lines every 20º

    mp.drawparallels(np.arange(-90, 110, 20),
                     labels=[True, False, False, True],
                     linewidth=0.2, fontsize=10)  # Draw lat lines every 20º

    mp.fillcontinents(color='lightgray')

    color = list(matplotlib.colors._colors_full_map.values())

    # Plot a Lagrangian air parcel with parcel ID given by "parcel2"

    # ax.scatter(plon[:,parcel2], parcel2, plat[:parcel2], parcel2, s=20, marker='o', color=color[])

    for i in range(0, len(parcel2)):
        # for j in range(0, len(parcel2[i])):
        ax.scatter(plon[:, parcel2[i]], plat[:, parcel2[i]], s=20, marker='o', color=color[i],
                   zorder=2)

    # ani = FuncAnimation(fig, animate, interval=100, frames=30, repeat = True)
    # ani.save('animation.gif')

    # Plot start and end points with an "S" and "F" respectively.
    ax.scatter(plon[0, parcel2], plat[0, parcel2], s=140, marker='$S$', color='red',
               zorder=2)

    ax.scatter(plon[-1, parcel2], plat[-1, parcel2], s=140, marker='$F$', color='red',
               zorder=2)

    # Save and close the map plot
    plt.savefig("air_parcel_ID"+str(parcel2)+"_map.png", format="png", dpi=300)
    plt.show()
    plt.close()


# =============================================================================
# PLOT TYPE 3 - VERTICAL EVOLUTION OF LAGRANGIAN AIR PARCELS W/ COLORBAR
# =============================================================================

# Requires ATTILA air parcel trajectory locatio and O3 data


############# INPUT SECTION ############

time_window = 40  # Time window (days) for calculating the RoD ##
# emission_point = 25                                                          ## The number of emission point (choose from 1 to 28) ##

#########################################

# Time window splitted into an array, for future plotting and finding elements ##
time_window_arr = np.arange(0, time_window, 0.25)
# Total number of time increments until time window (the increament 0.25 is given in the data, cannot change) ##
number_of_t = int(time_window / 0.25)


#########Empty Array Definitions##########
# An array representing the rate of descent of the 50 parcels in location N(input), it should contain 50 elements##
RoD_arr = np.array([])
# An array representing the average mixing ratio of the 50 parcels ##
MR_arr = np.array([])
RoD_average_arr = np.array([])
MR_average_arr = np.array([])
End_point_arr = np.array([])



###########################################################################################################################################
########################################################### average trajectory method #####################################################



for emission_point in range(1, 29):

    # A loop covering all 50 parcels in one emission location ##
    for i in range(((emission_point-1) *50), (emission_point)*50):

        ppress_temp = ppress[:, i]
        ppress_temp1 = ppress_temp[0: number_of_t]
        # Find where the minimum altitude A.K.A. maximum pressure (that's why the max in the function) ##
        min = int(np.where(ppress_temp1 == np.max(ppress_temp1))[0])

        if min == 0:

            min = int(np.where(ppress_temp1 == np.min(ppress_temp1))[0])

        # Find out the time corresponding to the minimum altitude ##
        time_at_minimum = time_window_arr[min]
        RoD = (- ppress_temp1[0] + ppress_temp1[min]) / time_at_minimum
        # Append the R.O.D. of each single parcel into an array ##
        RoD_arr = np.append(RoD_arr, RoD)
        # Mixing ratio of a single parcel, expressed as an array W.R.T. time window ##
        mr_one_parcel = airO3_001[:, i][0:number_of_t]
        # Average mixing ratio of a single parcel throughout the time window ##
        average_mr_one_parcel = np.average(mr_one_parcel)
        # Append the mixing ratio ##
        MR_arr = np.append(MR_arr, average_mr_one_parcel)

    

    ##############OLD STUFF FOR AVG################
    MR_average = np.average(MR_arr)
    RoD_average = np.average(RoD_arr)

    MR_average = np.average(MR_arr)
    RoD_average = np.average(RoD_arr)

    RoD_average_arr = np.append(RoD_average_arr, RoD_average)
    MR_average_arr = np.append(MR_average_arr, MR_average)



#############################################################################################################################################
#################################################### Median trajectory method ###############################################################



plot_median_trajectory = False 

median_trajectory_matrix = np.array([])
RoD_arr_for_median = np.array([])
MR_arr_for_median = np.array([])
MR_average_arr_median = np.array([])

fig, ax = plt.subplots(ncols = 2, nrows = 1)
colors = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"]
cmap = matplotlib.colors.ListedColormap(colors)
ax[0].invert_yaxis()
plt.suptitle('Median Trajectory Method')


for emission_point in range(1, 29):
    median_trajectory = np.array([])
    MR_arr_for_median = np.array([])
    
    median_trajectory_matrix = ppress[:, ((emission_point-1) * 50):(emission_point * 50)]
    MR_matrix_for_median = airO3_001[:, ((emission_point-1) * 50):(emission_point * 50)]

    # print(np.shape(median_trajectory_matrix))
    for k in range(0, number_of_t):
        median_vector = median_trajectory_matrix[k,:]
        median_MR_vector = MR_matrix_for_median[k,:]
        median_MR = np.median(median_MR_vector)
        median_value = np.median(median_vector)
        median_trajectory = np.append(median_trajectory, median_value)
        MR_arr_for_median = np.append(MR_arr_for_median, median_MR)

    
    MR_average_median = np.average(MR_arr_for_median)

    MR_average_arr_median = np.append(MR_average_arr_median, MR_average_median)
    
    
    # print('dadasfoa', np.shape(MR_average_arr_median))

    min = int(np.where(median_trajectory == np.max(median_trajectory))[0])
    if min == 0:
        min = int(np.where(median_trajectory == np.min(median_trajectory))[0])
    time_at_minimum = time_window_arr[min]
    RoD = (- median_trajectory[0] + median_trajectory[min]) / time_at_minimum
    RoD_arr_for_median = np.append(RoD_arr_for_median, RoD)

    # print(np.shape(time_window_arr), np.shape(MR_average_arr_median))
    ax[0].scatter(time_window_arr, median_trajectory)
    ax[0].set_title(load_data_from)
    ax[0].set_xlabel('Days')
    ax[0].set_ylabel('Altitude')
# plt.colorbar()
MR_average_arr_median = np.delete(MR_average_arr_median, 0)
RoD_arr_for_median = np.delete(RoD_arr_for_median, 0)
ax[1].set_title('RoD vs. Mr')
ax[1].set_xlabel('RoD')
ax[1].set_ylabel('Average MR')
ax[1].scatter(RoD_arr_for_median, (MR_average_arr_median * 10E9), s = 2)
fig.subplots_adjust(wspace = 0.3)

if plot_median_trajectory == True:
    ccp, pp = scipy.stats.pearsonr(RoD_arr_for_median, MR_average_arr_median * 10E15)
    # print("Median, Pearson correlation coefficient + p-value: ", str(ccp), ", ", str(pp))
    ccs, ps = scipy.stats.spearmanr(RoD_arr_for_median, MR_average_arr_median * 10E15)
    # print("Median, Spearman correlation coefficient + p-value: ", str(ccs), ", ", str(ps))
    cck, pk = scipy.stats.kendalltau(RoD_arr_for_median, MR_average_arr_median * 10E15)
    # print("Median, Kendall correlation coefficient + p-value: ", str(cck), ", ", str(pk))
    plt.show()
else: 
    plt.close()


# print('2', np.shape(MR_average_arr_median))
# print('1', np.shape(RoD_arr_for_median))



############################################################################################################################################
################################################ Closest trajectory method #################################################################


plot_closest_trajectory = True


index_array = np.array([])
RoD_arr_new_method = np.array([])

for emission_point in range(1, 29):

    
    trajectory_matrix = np.copy(ppress[:,((emission_point-1) * 50):(emission_point * 50)])
    trajectory_matrix_needed = np.transpose(trajectory_matrix[0: number_of_t])
    # print('shape of trajectory matrix needed', np.shape(trajectory_matrix_needed))
    traj_matrix_1emission_point_needed = np.transpose(trajectory_matrix[0:number_of_t])
    # print('2',np.min(ppress[:,29]))
    median_trajectory = np.array([])

    for k in range(0, number_of_t):
        median_vector = trajectory_matrix[k, :]
        median_value = np.median(median_vector)
        median_trajectory = np.append(median_trajectory, median_value)
        traj_matrix_1emission_point_needed[:, k] = np.abs(traj_matrix_1emission_point_needed[:, k] - median_trajectory[k])

    distance_array = np.array([])

    for i in range(0,50):
        distance_array = np.append(distance_array, np.sum(traj_matrix_1emission_point_needed[i, :]))
    index = int(np.where(distance_array == np.min(distance_array))[0])
    index_array = np.append(index_array, int((index + 50 * (emission_point-1))))
    # print(index_array)
    # for i in range(0, len(index_array)):

# print(index_array)

MR_avg_new_method_arr = np.array([])
for i in index_array:
    MR_new_method = np.copy(airO3_001[:, 0:1400])
    MR_new_method_needed = np.transpose(MR_new_method[0:number_of_t])
    MR_avg_new_method = np.average(MR_new_method_needed[int(i),:])
    MR_avg_new_method_arr = np.append(MR_avg_new_method_arr, MR_avg_new_method)


    new_trajectory_matrix = np.transpose(np.copy(ppress[:, 0:1400])[0: number_of_t])
    new_trajectory = new_trajectory_matrix[int(i), :]
    # new_trajectory_matrix_needed = new_trajectory_matrix[0: number_of_t]
    # print('gdsfas', len(new_trajectory))

    min = int(np.where(new_trajectory == np.max(new_trajectory))[0])
    if min == 0:
        min = int(np.where(new_trajectory == np.min(new_trajectory))[0])

    time_at_minimum = time_window_arr[min]
    RoD = (- new_trajectory[0] + new_trajectory[min]) / time_at_minimum
    RoD_arr_new_method = np.append(RoD_arr_new_method, RoD)

print('Rate of descent of 28 trajectories', RoD_arr_new_method)
MR_avg_new_method_arr = np.delete(MR_avg_new_method_arr, 0)
RoD_arr_new_method = np.delete(RoD_arr_new_method, 0)

# print('1', np.shape(MR_avg_new_method_arr))
# print('2', np.shape(RoD_arr_new_method))
# print(len(index_array))
# print(type(index_array))
# print('1', np.min(ppress[:, 29]))

fig, ax = plt.subplots(nrows = 1, ncols = 1)
#plt.suptitle('Closest Trajectory Method', fontweight = 'bold')
ax.invert_yaxis()

colors = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"]
# plot_colors = 
cmap = matplotlib.colors.ListedColormap(colors)
bounds = [0, 15, 30, 45, 60, 75]
# ax[1].scatter(RoD_arr_new_method, MR_avg_new_method_arr * 10E9, s = 2)
# ax[1].set_title('RoD vs. MR', fontweight = 'bold')
# ax[1].set_xlabel('RoD', fontweight = 'bold')
# ax[1].set_ylabel('Average MR', fontweight = 'bold')
norm = matplotlib.colors.Normalize(vmin=0, vmax=0.2)
fig.subplots_adjust(wspace = 0.4)
for i in index_array:
    traj = ppress[:, int(i)]
    traj_needed = traj[0 : number_of_t]
    # ax, fig = plt.subplots()
    # fig.invert_yaxis()
    sc = ax.scatter(time_window_arr, traj_needed, s = 3, cmap = 'viridis', norm = norm)
    # ax.set_title(load_data_from, fontweight = 'bold')
    ax.set_xlabel('Time Elapsed Since Emission [Days]', fontweight = 'bold', fontsize = 12)
    ax.set_ylabel('Air Paarcel Pressure Altitude [hPa]', fontweight = 'bold', fontsize = 12)
# fig.colorbar(sc, ticks=bounds, extend='both')
plt.savefig('28 Closest trajectory from ' + load_data_from, dpi = 512)


if plot_closest_trajectory == True:
    ccp, pp = scipy.stats.pearsonr(RoD_arr_new_method, MR_avg_new_method_arr * 10E9)
    print("New method, Pearson correlation coefficient + p-value: ", str(ccp), ", ", str(pp))
    ccs, ps = scipy.stats.spearmanr(RoD_arr_new_method, MR_avg_new_method_arr * 10E9)
    print("New method, Spearman correlation coefficient + p-value: ", str(ccs), ", ", str(ps))
    cck, pk = scipy.stats.kendalltau(RoD_arr_new_method, MR_avg_new_method_arr * 10E9)
    print("New method, Kendall correlation coefficient + p-value: ", str(cck), ", ", str(pk))
    plt.show()
else:
    plt.close()


fig, ax = plt.subplots()
plt.scatter(RoD_arr_new_method, MR_avg_new_method_arr * 10E9, s = 10)
# ax.set_title('RoD vs. MR[hPa/D]', fontweight = 'bold')
ax.set_xlabel('RoD [hPa/D]', fontweight = 'bold', fontsize = 12)
dot_multiply = '\u00B7'
superscript_minus_one = '\u207B\u00B9'
subscript_three = '\u2083'
ax.set_ylabel('Mean O' + subscript_three + ' Volume Mixing Ratio [nmol' + dot_multiply + 'mol' + superscript_minus_one + ']', fontweight = 'bold', fontsize = 11)
plt.savefig('RoD vs MR ' + load_data_from, dpi = 512)
plt.show()





#################################################################################################################################################
#################################################################################################################################################

    # print('there', np.shape(distance_array))
    # print('there', np.shape(traj_matrix_1emission_point_needed))
    # median_trajectory = np.array([])]
    # distance_50_traj = np.array([])
    # for i in range(((emission_point-1) * 50),(emission_point * 50)):
        
    #     individual_trajectory = ppress[:, i]
    #     individual_trajectory_needed = individual_trajectory[0 : number_of_t]
    #     # print('2', len(individual_trajectory_needed))
    #     # for m in range(0, 50):

    #     for m in range(0,number_of_t):
    #         distance_arr_one_traj = np.array([])
    #         distance = np.abs(individual_trajectory_needed[m] - median_trajectory[m])
    #         distance_arr_one_traj = np.append(distance_arr_one_traj, distance)
        
    #         distance_50_traj = np.append(distance_50_traj, np.sum(distance_arr_one_traj))
    # print('shitshow', len(distance_50_traj))
    # closest_traj = np.where(distance_50_traj == np.min(distance_50_traj))
        # print(closest_traj)
            # total_distance = np.sum(distance_arr)
            # total_distance_arr = np.append(total_distance)
    
    
    # print('show time', len(distance_50_traj))
    # print('verify', emission_point)

# print(len(median_trajectory))
# print('this is it',median_trajectory)

# ccp, pp = scipy.stats.pearsonr(RoD_average_arr, MR_average_arr)
# print("Average, Pearson correlation coefficient + p-value: ", str(ccp), ", ", str(pp))
# ccs, ps = scipy.stats.spearmanr(RoD_average_arr, MR_average_arr)
# print("Average, Spearman correlation coefficient + p-value: ", str(ccs), ", ", str(ps))
# cck, pk = scipy.stats.kendalltau(RoD_average_arr, MR_average_arr)
# print("Average, Kendall correlation coefficient + p-value: ", str(cck), ", ", str(pk))

# plot_parcel = False
# if plot_parcel == True:
#     fig, ax = plt.subplots()                                                        ## Plot MR W.R.T. RoD ##
#     fig.set_figheight(8)
#     fig.set_figwidth(15)
#     ax.grid(True)
#     # ax.set_aspect('equal')
#     # ax.set_xlim([0, 30])
#     # ax.set_ylim([0, 6E-8])
#     # ax.set_xticks(np.arange(0,30,2.5))
#     # ax.set_yticks(np.arange(0,6E-8,3E-9))
#     ax.scatter(RoD_arr, MR_arr)
#     ax.set_xlabel('Rate of descent of 50 parcels')
#     ax.set_ylabel('Mixing ratio of 50 parcels')
#     ax.set_title('MR with respect to RoD, emission location' + str(emission_point))
#     plt.show()
#     plt.close()

# plot_emission_point = False 
# if plot_emission_point == True: 
#     fig, ax = plt.subplots()  # Plot MR W.R.T. RoD ##
#     fig.set_figheight(8)
#     fig.set_figwidth(15)
#     ax.grid(True)
#     ax.scatter(RoD_average_arr, MR_average_arr * 10E9)
#     ax.set_xlabel('Emission point rate of descent [hPa/day]')
#     ax.set_ylabel('Emission point ozone mixing ratio [nmol/mol]')
#     ax.set_title("July 2014, 250hPa, 40 day window")
#     plt.savefig('Mr vs RoD, July 2014 40D')
#     plt.show()
#     plt.close()

activate_plot3 = False  # activation of vertical location plot with colorbar##

if attila_switch == True and o3tracer_switch == True and activate_plot3 == True:

    # Set up axis object for plotting the map
    fig, ax = plt.subplots(nrows= 7, ncols= 4)

    # Adjust dimensions of map plot
    fig.set_figheight(10)
    fig.set_figwidth(20)

    # parcel3 = 0  # Parcel ID, 0 means first.

    # Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
    colors = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"]
    cmap = matplotlib.colors.ListedColormap(colors)
    bounds = [0, 15, 30, 45, 60, 75] 

    cmap.set_under("w")
    cmap.set_over("crimson")

    norm = matplotlib.colors.Normalize(vmin=0, vmax=75)
    m = 0
    n = 0

    # Scatter plot command

    plot_vertical_trajectory = False
    if plot_vertical_trajectory == True: 
        for emission_points in range(1,29): 
            for parcel3 in range((emission_points - 1) * 50,emission_points * 50): 
                array = airO3_001[:, parcel3]
                # print('here', array[1])
                median = np.median(array)
                # average = np.average(array)
            #     min_filtered = np.min(filtered_array)
                power = np.abs(math.log10(median))
                # print(power)
                sc = ax[m,n].scatter(time, ppress[:, parcel3], s=2, marker='o',
                        c=airO3_001[:, parcel3] * (10 ** 9),
                        cmap=cmap, norm=norm, linewidth=1)
            ax[m,n].set_title(str(emission_points))
            ax[m,n].invert_yaxis()    
            cb = fig.colorbar(sc, ticks=bounds, extend='both')
            cb.ax.tick_params(labelsize=10)
            # cb.set_label(
            # label="O$_3$ Mixing Ratio [nmol·mol$^{-1}$]", size=10, weight='bold')
            n = n + 1
            if n == 4:
                m = m + 1
                n = 0
            if m == 7:
                break
            
    # print('lennnnnnn', len(airO3_001[:,parcel3]))

    # Pressure altitude increases towards the surface, reading convention
    
    # Set spacing and sizing of axes tickmarks
    # ax.set_xticks(np.arange(0, 110, 10))
    # ax.set_yticks(np.arange(0, 1200, 200))

    # ax.xaxis.set_tick_params(labelsize=20) 
    # ax.yaxis.set_tick_params(labelsize=20) 

    # Labeling of axes and plot title
    # ax.set_xlabel(
    #     "Time elapsed since emission \n [Days]", fontsize=22, weight='bold')
    # ax.set_ylabel(
    #     "Air parcel pressure altitude \n [hPa]", fontsize=22, weight='bold')
    # ax.set_title("Air parcel with colorbar - vertical",
    #              fontsize=24, weight='bold')

    # Define colorbar features
    
    # Adjust colorbar tickmark size
    # cb.ax.tick_params(labelsize=18)

    # Label the colorbar

    # Save and close the map plot
        plt.tight_layout()  # Ensure all parts of the plot will show after saving
        # plt.savefig("air_parcel_ID"+str(parcel3)+"_vertical_colorbar.png",format="png", dpi=300)
        plt.show()
        plt.close()


####################################################################################################################################################
####################################################################################################################################################

original_plot = False 

if attila_switch == True and o3tracer_switch == True and original_plot == True:

    #Set up axis object for plotting the map
    fig, ax = plt.subplots()
    
    #Adjust dimensions of map plot
    fig.set_figheight(8)
    fig.set_figwidth(14)
    
    # parcel3 = 1378 #Parcel ID, 0 means first.
    
    #Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
    colors = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"]
    cmap= matplotlib.colors.ListedColormap(colors)
    bounds = [0, 15, 30, 45, 60, 75]
    
    cmap.set_under("w")
    cmap.set_over("crimson")
    
    norm= matplotlib.colors.Normalize(vmin=0,vmax=75)
    
    #Scatter plot command
    for parcel3 in range(1350, 1400):
        sc = ax.scatter(time, ppress[:,parcel3], s=30, marker='o', 
                    c=airO3_001[:,parcel3]*1E09,
                    cmap=cmap,norm=norm ,linewidth=1)
    
    #Pressure altitude increases towards the surface, reading convention
    ax.invert_yaxis()
    
    #Set spacing and sizing of axes tickmarks
    ax.set_xticks(np.arange(0,110,10))
    ax.set_yticks(np.arange(0,1200,200))
    
    ax.xaxis.set_tick_params(labelsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    
    #Labeling of axes and plot title
    ax.set_xlabel("Time elapsed since emission \n [Days]" , fontsize=22, weight='bold')
    ax.set_ylabel("Air parcel pressure altitude \n [hPa]", fontsize=22, weight='bold')
    # ax.set_title("Emission point 1 parcels with colorbar, t=90 days", fontsize=24, weight='bold')
    
    #Define colorbar features
    cb = fig.colorbar(sc, ticks=bounds, extend='both')
    
    #Adjust colorbar tickmark size
    cb.ax.tick_params(labelsize=18)
    
    #Label the colorbar
    cb.set_label(label="O$_3$ Mixing Ratio [nmol·mol$^{-1}$]",size=18,weight='bold')
    
    #Save and close the map plot
    plt.tight_layout() #Ensure all parts of the plot will show after saving
    plt.savefig("air_parcel_ID_0-50A_vertical_colorbar.png",
                format="png",dpi=300)
    plt.show()
    plt.close()

    # plt.show()
    # plt.close()

# =================================================================================
# PLOT TYPE 4 - HORIZONTAL EVOLUTION OF LAGRANGIAN AIR PARCELS (ON MAP) W/ COLORBAR
# =================================================================================

# Requires ATTILA air parcel trajectory locatio and O3 data

activate_plot4 = False  # Activation mixing ratio plot ##

if attila_switch == True and o3tracer_switch == True and activate_plot4 == True:

     # parcel4 Parcel ID, 0 means first.
    parcel4 = 24
    # Set up axis object for plotting the map
    fig, ax = plt.subplots()  # Subplots are useful for drawing multiple plots together

    # Adjust dimensions of map plot
    fig.set_figheight(8)
    fig.set_figwidth(14)

    # Define map projection and settings
    # For more info: https://matplotlib.org/basemap/users/cyl.html
    mp = Basemap(projection='cyl',  # equidistant cylindrical projection
                 llcrnrlon=-180,
                 llcrnrlat=-90,
                 urcrnrlon=180,
                 urcrnrlat=90,
                 resolution='i', ax=ax)  # h=high, f=full, i=intermediate, c=crude

    # Format the lat and lon arrays for map graphing,
    # makes lat array a lat x lon array and same for lon array
    lon, lat = np.meshgrid(lons_18to18, lats)
    x, y = mp(lon, lat)

    # Choose the settings for the coastlines, countries, meridians...
    mp.drawcoastlines(linewidth=0.2)
    mp.drawcountries(linewidth=0.2)

    meridians = mp.drawmeridians(np.arange(-180, 200, 20),
                                 labels=[False, False, False, True],
                                 linewidth=0.2, fontsize=10)  # Draw lon lines every 20º

    mp.drawparallels(np.arange(-90, 110, 20),
                     labels=[True, False, False, True],
                     linewidth=0.2, fontsize=10)  # Draw lat lines every 20º

    mp.fillcontinents(color='lightgray')

    # Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
    colors = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"] 
    cmap = matplotlib.colors.ListedColormap(colors)
    bounds = [0, 15, 30, 45, 60, 75]

    cmap.set_under("w")
    cmap.set_over("crimson")

    norm = matplotlib.colors.Normalize(vmin=0, vmax=75)


    #################### Array definitions ########################################
    print('plon',np.shape(plon))
    print('plat', np.shape(plat))
    plon_median = np.array([])
    plat_median = np.array([])
    airO3_001_median_arr = np.array([])
    ##################### Determining the median trajectories #####################
    # fig, ax = plt.subplots(ncols = 2, nrows = 1)
    # ax[0].invert_yaxis()

    for emission_point in range(1, 29):
        median_trajectory = np.array([])
        MR_arr_for_median_hor = np.array([])
        
        plat_median_matrix = plat[:, ((emission_point-1) * 50):(emission_point * 50)]
        plon_median_matrix = plon[:, ((emission_point-1) * 50):(emission_point * 50)]
        airO3_001_median_arr = airO3_001[:, ((emission_point-1) * 50):(emission_point * 50)]

        # print(np.shape(median_trajectory_matrix))
        for k in range(0, number_of_t):
            plon_median_vector = plon_median_matrix[k,:]
            plat_median_vector = plat_median_matrix[k,:]

            airO3_001_median_vector = airO3_001_median_arr[k,:]
            airO3_001_median = np.median(airO3_001_median_vector)
            plon_median_value = np.median(plon_median_vector)
            plat_median_value = np.median(plat_median_vector)
            plon_median = np.append(plon_median, plon_median_value)
            plat_median = np.append(plat_median, plat_median_value)
            MR_arr_for_median_hor = np.append(MR_arr_for_median_hor, airO3_001_median)

        
        MR_average_median = np.average(MR_arr_for_median)

        MR_average_arr_median = np.append(MR_average_arr_median, MR_average_median)


    ##################### Plotting the median trajectories ########################
    for emission_point in range(0,29):  
        # Plot a Lagrangian air parcel with parcel ID given by "parcel4"
        # sc = ax.scatter(plon[:, parcel4], plat[:, parcel4], s=20, marker='o',
        #                 c=airO3_001[:, parcel4]*1E09, cmap=cmap, norm=norm, zorder=2)

        # # Plot starting point with an "S", "+4" is added to avoid overlay of letter on point
        # ax.scatter(plon[0, parcel4]+4, plat[0, parcel4], s=140, marker='$S$', color='black',
        #         zorder=2)
        
        sc = ax.scatter(plon_median[:, emission_point], plat_median[:, emission_point], s=20, marker='o',
                        c=airO3_001_median[:, emission_point]*1E09, cmap=cmap, norm=norm, zorder=2)

        # Plot starting point with an "S", "+4" is added to avoid overlay of letter on point
        ax.scatter(plon_median[0, emission_point]+4, plat_median[0, emission_point], s=140, marker='$S$', color='black',
                zorder=2)

    # Define colorbar features
    cb = fig.colorbar(sc, ticks=bounds, extend='both',
                      orientation='vertical', fraction=0.04, pad=0.03)

    # Adjust colorbar tickmark size
    cb.ax.tick_params(labelsize=14)

    # Label the colorbar
    cb.set_label(
        label="O$_3$ Mixing Ratio [nmol·mol$^{-1}$]", size=14, weight='bold')

    # Save and close the map plot
    plt.savefig("air_parcel_ID"+str(parcel4) +
                "_map_colorbar.png", format="png", dpi=300)
    plt.show()
    plt.close()
    ########################## Statistical analysis of the median trajectories and MR #############################


# =============================================================================
# PLOT TYPE 5 - Net radiative fluxes from short-term ozone increase (Single EP)
# =============================================================================

# Requires lat,lon from ATTILA data files and fluxes

activate_plot5 = False

if attila_switch == True and rad_fluxes_switch == True and activate_plot5 == True:

    # Set up axis object for plotting the map
    fig, ax = plt.subplots()  # Subplots are useful for drawing multiple plots together

    # Adjust dimensions of map plot
    fig.set_figheight(8)
    fig.set_figwidth(14)

    # Define map projection and settings
    # For more info: https://matplotlib.org/basemap/users/cyl.html
    mp = Basemap(projection='cyl',  # equidistant cylindrical projection
                 llcrnrlon=-180,
                 llcrnrlat=-90,
                 urcrnrlon=180,
                 urcrnrlat=90,
                 resolution='i', ax=ax)  # h=high, f=full, i=intermediate, c=crude

    # Shift the fluxes from [0,360] to [-180,180]
    net_flx_EP_shft, lons_shft = shiftgrid(180., global_net_flx[0],
                                           lons_0to36, start=False)

    # Format the lat and lon arrays for map graphing,
    # makes lat array a lat x lon array and same for lon array
    lon, lat = np.meshgrid(lons_shft, lats)
    x, y = mp(lon, lat)

    # Choose the settings for the coastlines, countries, meridians...
    mp.drawcoastlines(linewidth=0.2)
    mp.drawcountries(linewidth=0.2)

    meridians = mp.drawmeridians(np.arange(-180, 200, 20),
                                 labels=[False, False, False, True],
                                 linewidth=0.2, fontsize=10)  # Draw lon lines every 20º

    mp.drawparallels(np.arange(-90, 110, 20),
                     labels=[True, False, False, True],
                     linewidth=0.2, fontsize=10)  # Draw lat lines every 20º

    # Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
    colors = ["#ffffff", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"]
    cmap = matplotlib.colors.ListedColormap(colors)
    bounds = [0, 0.5, 1, 1.5, 2, 2.5]

    cmap.set_under("w")
    cmap.set_over("red")

    norm = matplotlib.colors.Normalize(vmin=0, vmax=2.5)

    # Time-average for the first emission point
    time_avg_flux = np.mean(net_flx_EP_shft, axis=0)

    # Plot the flux on the map
    sc2 = mp.pcolor(x, y, time_avg_flux*1000,
                    cmap=cmap, norm=norm, shading='auto')

    # Define colorbar features
    cb = fig.colorbar(sc2, ticks=bounds, extend='both',
                      orientation='horizontal', fraction=0.052,
                      pad=0.065)

    # Adjust colorbar tickmark size
    cb.ax.tick_params(labelsize=14)

    # Label the colorbar
    cb.set_label(
        label="Radiative Forcing from Short-term O$_3$ [mW·m$^{-2}$]", size=14, weight='bold')

    # Save and close the map plot
    plt.savefig("rad_fluxes_map_example.png", format="png", dpi=300)
    plt.show()
    plt.close()
