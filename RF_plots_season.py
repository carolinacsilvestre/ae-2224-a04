# Import libraries
import glob                                     # Dynamic file name loading
import numpy as np                              # Array processing
from netCDF4 import Dataset                     # NETCDF file handling
import matplotlib.pyplot as plt                 # Plotting
from mpl_toolkits.basemap import Basemap        # For map plotting
from mpl_toolkits.basemap import shiftgrid      # For map plotting
from matplotlib.lines import Line2D             # For map plotting
import matplotlib.colors                        # For map plotting
import os
import time
import pandas as pd                             #For exporting data to csv

start = time.time()
# =============================================================================
# LOADING DATA FROM NETCDF (.NC) FILES
# =============================================================================
# USER INPUT - File path
# Insert file path to folder containing all input data, do not forget wildcard


def path_chooser(season, altitude):
    """
    Returns a file path as a string based on the input parameters.

    Parameters:
    season (str): The season, either "summer" or "winter"
    altitude (int): The altitude, either 200, 250, or 300

    Returns:
    str: A file path representing a specific file or files in the "Data" directory

    """
    # Check the season and altitude parameters to determine the file path
    if season == "summer":
        if altitude == 200:
            f_string = 'C:/Users/twanv/OneDrive - Delft University of Technology/Bsc - 2/Q3/Project/Data/200_July/*'
        if altitude == 250:
            f_string = 'C:/Users/twanv/OneDrive - Delft University of Technology/Bsc - 2/Q3/Project/Data/250_July/*'
        if altitude == 300: 
            f_string = 'C:/Users/twanv/OneDrive - Delft University of Technology/Bsc - 2/Q3/Project/Data/300_July/*' 
    else: 
        if altitude == 200:
            f_string = 'C:/Users/twanv/OneDrive - Delft University of Technology/Bsc - 2/Q3/Project/Data/200_Jan/*'
        if altitude == 250:
            f_string = 'C:/Users/twanv/OneDrive - Delft University of Technology/Bsc - 2/Q3/Project/Data/250_Jan/*'
        if altitude == 300:
            f_string = 'C:/Users/twanv/OneDrive - Delft University of Technology/Bsc - 2/Q3/Project/Data/300_Jan/*'
    
    # Return the file path as a string
    return f_string


def read_files(f_string, emission_point=None, RF=True, Atilla=False, O3=False, verbose=True):
    """
    It reads in the net fluxes for each of the 28 emission points (3 months) and returns a list of 28
    net fluxes
    
    :param f_string: the file path to the folder containing the files
    :param verbose: If True, prints out the file names that are being read in, defaults to True
    (optional)
    :return: A list of 28 net fluxes, each corresponding to a different emission point.
    """
    # Read in file names based on f_string variable
    filenames_all = sorted(glob.glob(f_string)) # Get all file names in f_string
    if verbose:
        print('Files in input folder: ')
        print(filenames_all)

    if RF:
        # Variables declaration
        rad_flx_SW = []         # Holds 1 specific call for SW radiative fluxes
        rad_flx_LW = []         # Holds 1 specific call for LW radiative fluxes
        rad_flx_SW_02 = []      # SW flux for call 02
        rad_flx_LW_02 = []      # LW flux for call 02
        global_net_flx = []     # Holds all net fluxes for the 28 EPs (3 months)

        # Radiative fluxes corresponding to O3 increase, specifically for 250hPa
        if '250' in f_string:
            for file in filenames_all:
                # Get call 2 separately, which will be subtracted from calls 3-30
                if file.find('fluxes_tp') != -1 and file.find('EP02') != -1:
                    # fluxes_tp_NAmerica_July2014_EP02
                    data = Dataset(file, 'r')
                    if verbose:
                        print('\n')
                        print('File loaded: ')
                        print(file)

                    # Radiative fluxes
                    rad_flx_SW_02 = data.variables['flxs_tp'][:]
                    rad_flx_LW_02 = data.variables['flxt_tp'][:]
                    
                # Start at 3 since first two calls are not emission points
                for ep in range(3, 31):
                    # If there is no match, output is -1.
                    if file.find('fluxes_tp') != -1 and file.find('EP'+str(ep).zfill(2)) != -1:
                        data = Dataset(file, 'r')
                        if verbose:
                            print('\n')
                            print('File loaded: ')
                            print(file)

                        # Radiative fluxes
                        rad_flx_SW = data.variables['flxs_tp'][:]
                        rad_flx_LW = data.variables['flxt_tp'][:]
                        # Calculate net flux and append
                        global_net_flx.append((rad_flx_LW+rad_flx_SW) - (rad_flx_LW_02+rad_flx_SW_02))
                        
                        # Delete unnecessary variables
                        del rad_flx_LW, rad_flx_SW

        # Fluxes from the VISO submodel, for 200 and 300 hPa
        if not '250' in f_string:
            # Start at 3 since first two calls are not emission points
            for ep in range(3,31):
                for file in filenames_all:
                    if 'viso' in file:
                        data = Dataset(file,'r')
                        if verbose:
                            print('\n')
                            print('File loaded: ')
                            print(file)

                        
                        
                        #Store call 2 containing only background O3 fluxes
                        temp = data.variables['flxs_tp_02'][:]
                        rad_flx_SW_02.append(temp)
                        temp = data.variables['flxt_tp_02'][:]
                        rad_flx_LW_02.append(temp)
                        
                        #Other emission point
                        temp = data.variables['flxs_tp_'+str(ep).zfill(2)][:]
                        rad_flx_SW.append(temp)
                        
                        temp = data.variables['flxt_tp_'+str(ep).zfill(2)][:]
                        rad_flx_LW.append(temp)        
                #Concatenate the SW and LW flux data for all 3 months
                rad_flx_SW = np.concatenate(rad_flx_SW, axis=0)
                rad_flx_LW = np.concatenate(rad_flx_LW, axis=0) 
                rad_flx_SW_02 = np.concatenate(rad_flx_SW_02, axis=0)
                rad_flx_LW_02 = np.concatenate(rad_flx_LW_02, axis=0)
                
                #Calculate and append net fluxes for each of the 28 emission points
                global_net_flx.append((rad_flx_LW+rad_flx_SW)-(rad_flx_LW_02+rad_flx_SW_02))
                
                #Clear the radiative flux list holding the 3-month period data for a given call
                rad_flx_SW = []
                rad_flx_LW = []
                rad_flx_SW_02 = []
                rad_flx_LW_02 = []
                
                if verbose:
                    print('Net flux for EP'+str(ep-2)+' loaded...')
                
            del rad_flx_LW, rad_flx_SW, rad_flx_LW_02, rad_flx_SW_02

        #Delete unnecessary variables
        if not '250' in f_string:
            del temp, data
        return global_net_flx
    
    if Atilla:
        #Variables declaration
        time = [] #Attila
        plat = [] #Attila
        plon = [] #Attila
        ppress = [] #Attila in Pa

        for file in filenames_all:
            if 'attila' in file:
                data = Dataset(file,'r')
                #print('\n')
                #print('File loaded: ')
                #print(file)
                
                #Latitudes and longitudes
                lats = data.variables['lat'][:]
                lons_0to36 = data.variables['lon'][:] #Varies from 0 to 360 deg
                lons_18to18 = data.variables['lon'][:] #Will be converted later to -180 to 180 deg
                
                #Time
                temp = data.variables['time'][:]
                time.append(temp)
                
                #Air parcel longitudinal position
                temp = data.variables['PLON'][:]
                plon.append(temp)
                
                #Air parcel latitudinal position
                temp = data.variables['PLAT'][:]
                plat.append(temp)
                
                #Air parcel pressure altitude
                temp = data.variables['PPRESS'][:]
                ppress.append(temp)
    
        #Concatenation of variables, lists become multi-dimensional numpy arrays
        time = np.concatenate(time, axis=0)
        plon = np.concatenate(plon, axis=0)
        plat = np.concatenate(plat, axis=0)
        ppress = np.concatenate(ppress, axis=0)/100 #Convert to hPa
        
        #Convert longitude range from [0, 360] to [-180, 180]
        plon = (plon + 180) % 360 -180 #Air parcel longitudinal coordinates
        lons_18to18 = (lons_18to18 + 180) % 360 -180 #EMAC grid longitudes
        return lats ,lons_0to36,plon,plat
    if O3:
        airO3_001 = []
        key = "airO3_" + str(emission_point).zfill(3)
        #O3 data along air parcel trajectories
        for file in filenames_all:
            if 'O3lg' in file:
                data = Dataset(file,'r')
                #print('\n')
                #print('File loaded: ')
                #print(file)
                
                #O3 along air parcels for each emission point
                temp = data.variables[key][:,:1400] #Only the first 1400 trajectories
                airO3_001.append(temp)
                
                #NOTE: do not forget there are 28 emission points in total!
        
        #Concatenation of variables, lists become multi-dimensional numpy arrays
        airO3_001 = np.concatenate(airO3_001, axis=0)
        return airO3_001


def calculate_values(global_net_flx):
    """
    It takes a 3D array (global_net_flx) and returns a 3D array (flux_list) that is the sum of the mean
    of each 2D array in the 3D array
    
    :param global_net_flx: a list of 28 arrays, each array is a 2D array of shape (time, position)
    :return: lat, lon, flux_list
    """
    flux_list = []
    for n in range(28):
        flux = global_net_flx[n]
        flux_time_avg = np.mean(flux,axis=0)
        flux_list.append(np.mean(flux_time_avg)*10**(-3)*A_section)
    flux_list = np.array(flux_list).reshape(7, 4,order="F") 
    lat = np.linspace(85, 25, 7)  # define x as an array with 4 elements
    lon = np.linspace(-115, -55, 4)  # define y as an array with 7 elements
    max = np.max(flux_list)
    return lat, lon, flux_list,max


def plot_RF_all_points(lat,lon,flux_list,colors="Reds",show=False, save=True, dpi=300,max=10,season="winter",altitude=300):
    """
    This function takes in a list of latitudes, a list of longitudes, a list of fluxes, and a color
    scheme, and plots the fluxes on a map
    
    :param lat: latitude array
    :param lon: longitude array
    :param flux_list: a list of the fluxes at each grid point
    :param colors: the color scheme for the map, defaults to Reds (optional)
    :param show: if True, the plot will be shown. If False, the plot will be saved to a file, defaults
    to True (optional)
    """

    # Set up axis object for plotting the map
    fig, ax = plt.subplots() # Subplots are useful for drawing multiple plots together
        
    # Adjust dimensions of map plot
    fig.set_figheight(8)
    fig.set_figwidth(14)
        
    # Define map projection and settings
    # For more info: https://matplotlib.org/basemap/users/cyl.html
    mp = Basemap(projection = 'cyl', #equidistant cylindrical projection
                            llcrnrlon = -135,
                            llcrnrlat = 10,
                            urcrnrlon = -35,
                            urcrnrlat = 90,
                            resolution = 'i', ax=ax) # h=high, f=full, i=intermediate, c=crude
        
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
    meridians = mp.drawmeridians(np.arange(-180, 200, 20), labels=[False,False,False,True], linewidth=0.2, fontsize=10)
    # Draw lat lines every 20º
    mp.drawparallels(np.arange(-90, 110, 20), labels=[True, False, False, True], linewidth=0.2, fontsize=10)
        
        
    norm= matplotlib.colors.Normalize(vmin=0,vmax=max)
    # Plot the flux on the map
    sc2 = mp.pcolor(x, y, flux_list, cmap='Reds',shading='auto',norm=norm)
    # Define colorbar features
    cb = fig.colorbar(sc2, extend='both', orientation='horizontal', fraction=0.052, pad=0.065)
        
    # Adjust colorbar tickmark size
    cb.ax.tick_params(labelsize=14)
        
    # Label the colorbar
    cb.set_label(label="Average Radiative Forcing (3 months) in "+season+ " at " +str(altitude)+"hPa in [W]",size=14,weight='bold',labelpad=15)
        
    # Save and close the map plot
    if save:
        plt.savefig("RFmap"+season+str(altitude)+".png", dpi=dpi)
    if show:
        plt.show()
    plt.close()


def plot_RF_per_point(global_net_flx,lons_0to36,lats,EP):
    # Set up axis object for plotting the map
    fig, ax = plt.subplots() #Subplots are useful for drawing multiple plots together
    
    #Adjust dimensions of map plot
    fig.set_figheight(8)
    fig.set_figwidth(14)
    
    #Define map projection and settings
    #For more info: https://matplotlib.org/basemap/users/cyl.html
    mp = Basemap(projection = 'cyl', #equidistant cylindrical projection
                        llcrnrlon = -180,
                        llcrnrlat = -90,
                        urcrnrlon = 180,
                        urcrnrlat = 90,
                        resolution = 'i', ax=ax) #h=high, f=full, i=intermediate, c=crude
    
    #Shift the fluxes from [0,360] to [-180,180]
    net_flx_EP_shft, lons_shft = shiftgrid(180.,global_net_flx[EP], 
                                        lons_0to36,start=False)
    
    #Format the lat and lon arrays for map graphing, 
    #makes lat array a lat x lon array and same for lon array
    lon, lat = np.meshgrid(lons_shft, lats)
    x, y = mp(lon, lat)
    
    #Choose the settings for the coastlines, countries, meridians...
    mp.drawcoastlines(linewidth=0.2)
    mp.drawcountries(linewidth=0.2)

    # Draw lon lines every 20º
    meridians = mp.drawmeridians(np.arange(-180,200,20), labels=[False,False,False,True], linewidth=0.2, fontsize=10)
    # Draw lat lines every 20º
    mp.drawparallels(np.arange(-90,110,20), labels=[True, False, False, True], linewidth=0.2, fontsize=10)
    
    # Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
    colors = ["#ffffff", "#fec44f", "#d95f0e", "#e34a33", "#b30000"]
    bounds = [0, 0.5, 1, 1.5, 2, 2.5]
    cmap= matplotlib.colors.LinearSegmentedColormap.from_list(bounds,colors)
    
    cmap.set_under("w")
    cmap.set_over("red")
    
    norm= matplotlib.colors.Normalize(vmin=0,vmax=2.5)
    
    # Time-average for the first emission point
    time_avg_flux = np.mean(net_flx_EP_shft, axis=0)
    
    # Plot the flux on the map
    sc2 = mp.pcolor(x, y, time_avg_flux*1000,
                    cmap=cmap, norm=norm, shading='auto')
    
    # Define colorbar features
    cb = fig.colorbar(sc2, #ticks=bounds, extend='both', 
                    orientation='horizontal',fraction=0.052, 
                    pad=0.065)
    
    # Adjust colorbar tickmark size
    cb.ax.tick_params(labelsize=14)
    
    # Label the colorbar
    cb.set_label(label="Radiative Forcing from Short-term O$_3$ [mW·m$^{-2}$]",size=14,weight='bold')
    
    # Save and close the map plot
    plt.savefig("rad_fluxes_"+str(altitude)+season +"ep" +str(EP+1) + ".png",format="png",dpi=300)
    if show:
        plt.show()
    plt.close()


def plot_big(lst:list,show:bool=True,save:bool=False):
    # First loop over all the files to establish a maximum and save all the values for later plotting
    max_list = []
    flux_list_all = []
    for season, altitude in lst:   
        path = path_chooser(season,altitude)
        flux = read_files(path, verbose=verbose)
        lat, lon, flux_list, maximum = calculate_values(flux)
        max_list.append(maximum)
        flux_list_all.append(flux_list)

    # plot the values saved earlier
    for i in range(len(lst)):
        season, altitude = lst[i]
        plot_RF_all_points(lat,lon,flux_list_all[i],max=max(max_list),season=season,altitude=altitude,show=show,save=save)
    """"
    # Used for saving the data to csv
    lst = ['summer200', 'summer250', 'summer300','winter200', 'winter250', 'winter300']
    df = pd.DataFrame(flux_list_all[0].reshape(28,order='F'), columns=[lst[0]])
    for i in range(1, len(flux_list_all)):
        df[lst[i]] = flux_list_all[i].reshape(28,order='F')
    print(df)
    df.to_csv('data.csv')"""


def plot_small(EP,season,altitude,show=True):
    path = path_chooser(season,altitude)
    flux = read_files(path, verbose=verbose,RF=True,Atilla=False)
    lats ,lons, plon,plat = read_files(path, verbose=verbose,RF=False,Atilla=True)
    plot_RF_per_point(flux,lons,lats,EP)


def plot_overlay(EP,global_net_flx,lons_0to36,lats,plon,plat,airO3_001,cut=30,plot_below_cut=False,c1='black',c2='green'):
    # Requires lat,lon from ATTILA data files and fluxes
    # Set up axis object for plotting the map
    emission_point = EP
    # Set up axis object for plotting the map
    fig, ax = plt.subplots() # Subplots are useful for drawing multiple plots together
    
    # Adjust dimensions of map plot
    fig.set_figheight(8)
    fig.set_figwidth(14)
    
    # Define map projection and settings
    # For more info: https://matplotlib.org/basemap/users/cyl.html
    mp = Basemap(projection = 'cyl', #equidistant cylindrical projection
                         llcrnrlon = -180,
                         llcrnrlat = -90,
                         urcrnrlon = 180,
                         urcrnrlat = 90,
                         resolution = 'i', ax=ax) #h=high, f=full, i=intermediate, c=crude
    
    #Shift the fluxes from [0,360] to [-180,180]
    net_flx_EP_shft, lons_shft = shiftgrid(180.,global_net_flx[emission_point-1], lons_0to36,start=False)
    
    #Format the lat and lon arrays for map graphing, 
    #makes lat array a lat x lon array and same for lon array
    lon, lat = np.meshgrid(lons_shft, lats)
    x, y = mp(lon, lat)
    
    #Choose the settings for the coastlines, countries, meridians...
    mp.drawcoastlines(linewidth=0.2)
    mp.drawcountries(linewidth=0.2)

    # Draw lon lines every 20º
    meridians = mp.drawmeridians(np.arange(-180,200,20), labels=[False,False,False,True], linewidth=0.2, fontsize=10)
    # Draw lat lines every 20º
    mp.drawparallels(np.arange(-90,110,20), labels=[True,False,False,True], linewidth=0.2, fontsize=10)
    
    # Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
    colors = ["#ffffff", "#fec44f", "#d95f0e", "#e34a33", "#b30000"]
    bounds = [0, 0.5, 1, 1.5, 2, 2.5]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(bounds,colors)
    cmap.set_under("w")
    cmap.set_over("red")
    
    norm = matplotlib.colors.Normalize(vmin=0,vmax=2.5)
    
    # Time-average for the first emission point
    time_avg_flux = np.mean(net_flx_EP_shft, axis=0)
    
    # Plot the flux on the map
    for i in range(50):
        # Plot a Lagrangian air parcels for emission point 'EP'
        ax.scatter(plon[:,(emission_point-1)*50+i][airO3_001[:,(emission_point-1)*50+i]*1E09>cut], plat[:,(emission_point-1)*50+i][airO3_001[:,(emission_point-1)*50+i]*1E09>cut], s=20, marker='o', color=c1,
               zorder=2,alpha = 0.05)
        if plot_below_cut:
            ax.scatter(plon[:,(emission_point-1)*50+i][airO3_001[:,(emission_point-1)*50+i]*1E09<30], plat[:,(emission_point-1)*50+i][airO3_001[:,(emission_point-1)*50+i]*1E09<30], s=20, marker='o', color=c2,zorder=2,alpha = 0.025)

    # Add the legend, and set the handler map to use the change_alpha function

    if plot_below_cut:
            legend_elements = [ Line2D([0], [0], markerfacecolor=c1,marker='o', color='w', label="mixing ratio above "+str(cut)+"%",markersize=15),
                        Line2D([0], [0], markerfacecolor=c2,marker='o', color='w', label="mixing ratio below "+str(cut)+"%", markersize=15),]
    else:
        legend_elements = [ Line2D([0], [0], markerfacecolor=c1,marker='o', color='w', label="mixing ratio above "+str(cut)+"%",markersize=15),]

    # Create the RF on the map
    plt.legend(handles=legend_elements)

    sc2 = mp.pcolor(x, y, time_avg_flux*1000,
                    cmap=cmap, norm=norm, shading='auto')
    # Define colorbar features
    cb = fig.colorbar(sc2, ticks=bounds, extend='both', 
                      orientation='horizontal',fraction=0.052, 
                      pad=0.065)
    # Adjust colorbar tickmark size
    cb.ax.tick_params(labelsize=14)
    # Label the colorbar
    cb.set_label(label="Radiative Forcing from Short-term O$_3$ [mW·m$^{-2}$]",size=14,weight='bold')
    
    # Save and close the map plot
    plt.savefig("rad_fluxes_map_overlay" + str(emission_point) + ".png",format="png",dpi=300)
    plt.show()
    plt.close()

# Changing the current working directory to the plots' folder.
os.chdir('C:/Users/twanv/OneDrive - Delft University of Technology/Bsc - 2/Q3/Project/plots')

""""
A_earth = 5.1e8 #km^2
frac = 2.8**2/(360*180)
A_section = A_earth*frac*10**6 #m^2
"""
A_section = 1

# Parameters for saving and showing the plots
verbose = False
show = False
save = True
lst = []

# Looping through the seasons and altitudes and plotting the fluxes for each of them.
for season in ["summer","winter"]:
    for altitude in [200,250,300]:
        lst.append([season,altitude])

plot_big(lst, show, save)

# plot overlay of position on map with radiative forcing
# select point and dataset:
# season = "summer"
# altitude = 250
# emission_point = 12

# import data
# f_string = path_chooser(season,altitude)
# global_net_flux = read_files(f_string,RF=True,Atilla=False,O3=False,verbose=False)
# lats , lons_0to36,plon,plat = read_files(f_string,RF=False,Atilla=True,O3=False,verbose=False)
# airO3 = read_files(f_string,emission_point=emission_point,RF=False,Atilla=False,O3=True,verbose=False)

# plot data
# plot_overlay(emission_point,global_net_flux,lons_0to36,lats,plon,plat,airO3,cut=15,plot_below_cut=True)

print(time.time() - start)
