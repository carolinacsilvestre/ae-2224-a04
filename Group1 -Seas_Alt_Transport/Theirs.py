'''Sample script for reading in data and visualizing the
Lagrangian air parcel trajectories along with the net radiative fluxes
arising from a short-term increase in ozone.'''  

#Import libraries
import glob #Dynamic file name loading
import numpy as np #Array processing
from netCDF4 import Dataset #NETCDF file handling
import matplotlib.pyplot as plt #Plotting
from mpl_toolkits.basemap import Basemap #For map plotting
from mpl_toolkits.basemap import shiftgrid #For shifting longitudes
import matplotlib.colors #To create new colorbar
import matplotlib.cm as cm

# =============================================================================
# LOADING DATA FROM NETCDF (.NC) FILES
# =============================================================================

#USER INPUT - File path
f_string = r"C:\Users\moheb\Desktop\Q3_Proj (Group Git)\*" #Insert file path to folder containing all input data, do not forget wildcard

#USER INPUT - Switches to determine which data types should be loaded
attila_switch = False
o3tracer_switch = False
rad_fluxes_switch = True
parcel_id = 249
#Read in file names based on f_string variable
filenames_all = sorted(glob.glob(f_string)) #Get all file names in f_string
print('Files in input folder: ')
print(filenames_all)


#Variables declaration
time = [] #Attila
plat = [] #Attila
plon = [] #Attila
ppress = [] #Attila in Pa
airO3_001 = [] #O3lg in mol/mol for emission point 1
rad_flx_SW = [] #Holds 1 specific call for SW radiative fluxes
rad_flx_LW = [] #Holds 1 specific call for LW radiative fluxes
rad_flx_SW_02 = [] #SW flux for call 02
rad_flx_LW_02 = [] #LW flux for call 02
global_net_flx = [] #Holds all net fluxes for the 28 EPs (3 months)

#Load relevant variables from NETCDF files
#Positions of air parcels
if attila_switch == True:
    for file in filenames_all:
        if 'attila.nc' in file:
            data = Dataset(file,'r')
            print('\n')
            print('File loaded: ')
            print(file)
            
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


#O3 data along air parcel trajectories
if o3tracer_switch == True:
    for file in filenames_all:
        if 'O3lg' in file:
            data = Dataset(file,'r')
            print('\n')
            print('File loaded: ')
            print(file)
            
            #O3 along air parcels for each emission point
            temp = data.variables['airO3_001'][:,:1400] #Only the first 1400 trajectories
            airO3_001.append(temp)
            
            #NOTE: do not forget there are 28 emission points in total!
    
    #Concatenation of variables, lists become multi-dimensional numpy arrays
    airO3_001 = np.concatenate(airO3_001, axis=0)
    
#Radiative fluxes corresponding to O3 increase, specifically for 250hPa
if rad_fluxes_switch == True:# and '250' in f_string:
    for file in filenames_all:
        
        #Get call 2 separately, which will be subtracted from calls 3-30
        if file.find('fluxes_tp') != -1 and file.find('EP02') != -1:
            #fluxes_tp_NAmerica_July2014_EP02
            data = Dataset(file,'r')
            print('\n')
            print('File loaded: ')
            print(file)
            
            #Radiative fluxes
            rad_flx_SW_02 = data.variables['flxs_tp'][:]
            rad_flx_LW_02 = data.variables['flxt_tp'][:]
            
        #Start at 3 since first two calls are not emission points
        for ep in range(3,31):
            #If there is no match, output is -1.
            if file.find('fluxes_tp') != -1 and file.find('EP'+str(ep).zfill(2)) != -1:
                data = Dataset(file,'r')
                print('\n')
                print('File loaded: ')
                print(file)
                
                #Radiative fluxes
                rad_flx_SW = data.variables['flxs_tp'][:]
                rad_flx_LW = data.variables['flxt_tp'][:]
                #Calculate net flux and append
                global_net_flx.append((rad_flx_LW+rad_flx_SW)- #\
                                      (rad_flx_LW_02+rad_flx_SW_02))
                
    #Delete unnecessary variables
    del rad_flx_LW, rad_flx_SW, rad_flx_LW_02, rad_flx_SW_02

#Fluxes from the VISO submodel, for 200 and 300 hPa
# if rad_fluxes_switch == True and not '250' in f_string:
#     #Start at 3 since first two calls are not emission points
#     for ep in range(3,31):
#         for file in filenames_all:
#             if 'viso' in file:
#                 data = Dataset(file,'r')
#                 print('\n')
#                 print('File loaded: ')
#                 print(file)
                
#                 #Store call 2 containing only background O3 fluxes
#                 temp = data.variables['flxs_tp_02'][:]
#                 rad_flx_SW_02.append(temp)
#                 temp = data.variables['flxt_tp_02'][:]
#                 rad_flx_LW_02.append(temp)
                
#                 #Other emission point
#                 temp = data.variables['flxs_tp_'+str(ep).zfill(2)][:]
#                 rad_flx_SW.append(temp)
                
#                 temp = data.variables['flxt_tp_'+str(ep).zfill(2)][:]
#                 rad_flx_LW.append(temp)        
#         #Concatenate the SW and LW flux data for all 3 months
#         rad_flx_SW = np.concatenate(rad_flx_SW, axis=0)
#         rad_flx_LW = np.concatenate(rad_flx_LW, axis=0) 
#         rad_flx_SW_02 = np.concatenate(rad_flx_SW_02, axis=0)
#         rad_flx_LW_02 = np.concatenate(rad_flx_LW_02, axis=0)
        
#         #Calculate and append net fluxes for each of the 28 emission points
#         global_net_flx.append((rad_flx_LW+rad_flx_SW)-(rad_flx_LW_02+rad_flx_SW_02))
        
#         #Clear the radiative flux list holding the 3-month period data for a given call
#         rad_flx_SW = []
#         rad_flx_LW = []
#         rad_flx_SW_02 = []
#         rad_flx_LW_02 = []
        
#         print('Net flux for EP'+str(ep-2)+' loaded...')
        
#     del rad_flx_LW, rad_flx_SW, rad_flx_LW_02, rad_flx_SW_02

# #Delete unnecessary variables
# if attila_switch or o3tracer_switch or rad_fluxes_switch and not '250' in f_string:
#     del temp, data
"""
# =============================================================================
# PLOT TYPE 1 - VERTICAL EVOLUTION OF LAGRANGIAN AIR PARCELS
# =============================================================================

#Requires ATTILA air parcel trajectory location data
if attila_switch == True:

    #Set up axis object for plotting the map
    fig, ax = plt.subplots()
    
    #Adjust dimensions of map plot
    fig.set_figheight(8)
    fig.set_figwidth(14)
    
    parcel = parcel_id #Parcel ID, 0 means first.
    
    #Scatter plot settings
    plt.scatter(time, ppress[:,parcel], s=30, marker='o', color='green')
    
    #Defining range of axes
    plt.xticks(np.arange(0,110,10), np.arange(0,110,10), fontsize=20)
    plt.yticks(np.arange(0,1200,200), np.arange(0,1200,200), fontsize=20)
    
    #Inversion of vertical axis
    plt.gca().invert_yaxis() #The higher the pressure, the closer to the surface
    
    #Labeling of axes
    plt.xlabel('Time elapsed since emission \n [Days]', fontsize=22, weight='bold')
    plt.ylabel('Air parcel pressure altitude \n [hPa]', fontsize=22, weight='bold')
    
    #Add plot title
    plt.title('Air Parcel ID='+str(parcel), fontsize=26, weight='bold')
    
    #To save the plot, give it a name, format and dpi represents the resolution.
    #Typically, for publications, dpi=300.
    plt.tight_layout() #Ensure all parts of the plot will show after saving
    plt.savefig("air_parcel_ID"+str(parcel)+".png",format="png",dpi=300)
    
    plt.show()
    plt.close()



# =============================================================================
# PLOT TYPE 3 - VERTICAL EVOLUTION OF LAGRANGIAN AIR PARCELS W/ COLORBAR
# =============================================================================

#Requires ATTILA air parcel trajectory locatio and O3 data
if attila_switch == True and o3tracer_switch == True:

    #Set up axis object for plotting the map
    fig, ax = plt.subplots()
    
    #Adjust dimensions of map plot
    fig.set_figheight(8)
    fig.set_figwidth(14)
    
    parcel3 = parcel_id #Parcel ID, 0 means first
    
    #Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
    colors = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"]
    cmap= matplotlib.colors.ListedColormap(colors)
    bounds = [0, 15, 30, 45, 60, 75]
    
    cmap.set_under("w")
    cmap.set_over("crimson")
    
    norm= matplotlib.colors.Normalize(vmin=0,vmax=75)
    
    #Scatter plot command
    sc = ax.scatter(time, ppress[:,parcel3], s=30, marker='o', c=airO3_001[:,parcel3]*1E09, cmap=cmap,norm=norm ,linewidth=1)
    
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
    ax.set_title("Air parcel with colorbar - vertical", fontsize=24, weight='bold')
    
    #Define colorbar features
    cb = fig.colorbar(sc, ticks=bounds, extend='both')
    
    #Adjust colorbar tickmark size
    cb.ax.tick_params(labelsize=18)
    
    #Label the colorbar
    cb.set_label(label="O$_3$ Mixing Ratio [nmol·mol$^{-1}$]",size=18,weight='bold')
    
    #Save and close the map plot
    plt.tight_layout() #Ensure all parts of the plot will show after saving
    plt.savefig("air_parcel_ID"+str(parcel3)+"_vertical_colorbar.png",
                format="png",dpi=300)
    plt.show()
    plt.close()

# =============================================================================
# PLOT TYPE 2 - HORIZONTAL EVOLUTION OF LAGRANGIAN AIR PARCELS (ON MAP)
# =============================================================================

#Requires ATTILA air parcel trajectory location data
if attila_switch == True:

    parcel2 = parcel_id #Parcel ID, 0 means first.
    
    #Set up axis object for plotting the map
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
    
    #Format the lat and lon arrays for map graphing, 
    #makes lat array a lat x lon array and same for lon array
    lon, lat = np.meshgrid(lons_18to18, lats)
    x, y = mp(lon, lat)
    
    #Choose the settings for the coastlines, countries, meridians...
    mp.drawcoastlines(linewidth=0.2)
    mp.drawcountries(linewidth=0.2)
    
    meridians = mp.drawmeridians(np.arange(-180,200,20), 
                         labels=[False,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lon lines every 20º
    
    mp.drawparallels(np.arange(-90,110,20), 
                         labels=[True,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lat lines every 20º
    
    mp.fillcontinents(color='lightgray')
    
    #Plot a Lagrangian air parcel with parcel ID given by "parcel2"
    ax.scatter(plon[:,parcel2], plat[:,parcel2], s=20, marker='o', color='green',
               zorder=2)
    
    #Plot start and end points with an "S" and "F" respectively.
    ax.scatter(plon[0,parcel2], plat[0,parcel2], s=140, marker='$S$', color='red', #s is the marker size
               zorder=2)
    
    ax.scatter(plon[-1,parcel2], plat[-1,parcel2], s=140, marker='$F$', color='red',
               zorder=2)
    
    #Save and close the map plot
    plt.savefig("air_parcel_ID"+str(parcel2)+"_map.png",format="png",dpi=300)
    plt.show()
    plt.close()

# =============================================================================
# PLOT TYPE 6 - HORIZONTAL EVOLUTION OF LAGRANGIAN ALL AIR PARCELS OVERLAY (ON MAP)
# =============================================================================

#Requires ATTILA air parcel trajectory location data
if attila_switch == True:
    
    Number_IDs = 28
    #Set up axis object for plotting the map
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
    
    #Format the lat and lon arrays for map graphing, 
    #makes lat array a lat x lon array and same for lon array
    lon, lat = np.meshgrid(lons_18to18, lats)
    x, y = mp(lon, lat)
    
    #Choose the settings for the coastlines, countries, meridians...
    mp.drawcoastlines(linewidth=0.2)
    mp.drawcountries(linewidth=0.2)
    
    meridians = mp.drawmeridians(np.arange(-180,200,20), 
                         labels=[False,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lon lines every 20º
    
    mp.drawparallels(np.arange(-90,110,20), 
                         labels=[True,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lat lines every 20º
    
    mp.fillcontinents(color='lightgray')
    
    colors = cm.rainbow(np.linspace(0,1,Number_IDs))
    for ID in range(Number_IDs):
        parcel2 = ID #Parcel ID, 0 means first.
        #Plot a Lagrangian air parcel with parcel ID given by "parcel2"
        ax.scatter(plon[:,parcel2], plat[:,parcel2], s=20, marker='o', color=colors[ID-1],
               zorder=2,alpha = 0.3)
    
        #Plot start and end points with an "S" and "F" respectively.
        #ax.scatter(plon[0,parcel2], plat[0,parcel2], s=140, marker='$S$', color='red', #s is the marker size
               #zorder=2)
    
        #ax.scatter(plon[-1,parcel2], plat[-1,parcel2], s=140, marker='$F$', color='red',
               #zorder=2)
    
    #Save and close the map plot
    plt.savefig("air_parcel_ID"+str(parcel2)+"_map.png",format="png",dpi=300)
    plt.show()
    plt.close()


# =================================================================================
# PLOT TYPE 4 - HORIZONTAL EVOLUTION OF LAGRANGIAN AIR PARCELS (ON MAP) W/ COLORBAR
# =================================================================================

#Requires ATTILA air parcel trajectory locatio and O3 data
if attila_switch == True and o3tracer_switch == True:

    parcel4 = parcel_id #Parcel ID, 0 means first.
    
    #Set up axis object for plotting the map
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
    
    #Format the lat and lon arrays for map graphing, 
    #makes lat array a lat x lon array and same for lon array
    lon, lat = np.meshgrid(lons_18to18, lats)
    x, y = mp(lon, lat)
    
    #Choose the settings for the coastlines, countries, meridians...
    mp.drawcoastlines(linewidth=0.2)
    mp.drawcountries(linewidth=0.2)
    
    meridians = mp.drawmeridians(np.arange(-180,200,20), 
                         labels=[False,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lon lines every 20º
    
    mp.drawparallels(np.arange(-90,110,20), 
                         labels=[True,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lat lines every 20º
    
    mp.fillcontinents(color='lightgray')
    
    #Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
    colors = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"]
    cmap= matplotlib.colors.ListedColormap(colors)
    bounds = [0, 15, 30, 45, 60, 75]
    
    cmap.set_under("w")
    cmap.set_over("crimson")
    
    norm= matplotlib.colors.Normalize(vmin=0,vmax=75)
    
    #Plot a Lagrangian air parcel with parcel ID given by "parcel4"
    sc = ax.scatter(plon[:,parcel4], plat[:,parcel4], s=20, marker='o',
               c=airO3_001[:,parcel4]*1E09,cmap=cmap,norm=norm,zorder=2)
    
    #Plot starting point with an "S", "+4" is added to avoid overlay of letter on point
    ax.scatter(plon[0,parcel4]+4, plat[0,parcel4], s=140, marker='$S$', color='black',
               zorder=2)
    
    #Define colorbar features
    cb = fig.colorbar(sc, ticks=bounds, extend='both', 
                      orientation='vertical',fraction=0.04, pad=0.03)
    
    #Adjust colorbar tickmark size
    cb.ax.tick_params(labelsize=14)
    
    #Label the colorbar
    cb.set_label(label="O$_3$ Mixing Ratio [nmol·mol$^{-1}$]",size=14,weight='bold')
    
    #Save and close the map plot
    plt.savefig("air_parcel_ID"+str(parcel4)+"_map_colorbar.png",format="png",dpi=300)
    plt.show()
    plt.close()
 """   
# =============================================================================
# PLOT TYPE 5 - Net radiative fluxes from short-term ozone increase (Single EP)
# =============================================================================

#Requires lat,lon from ATTILA data files and fluxes
if attila_switch == True and rad_fluxes_switch == True:
    
    #Set up axis object for plotting the map
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
    net_flx_EP_shft, lons_shft = shiftgrid(180.,global_net_flx[5], 
                                           lons_0to36,start=False)
    
    #Format the lat and lon arrays for map graphing, 
    #makes lat array a lat x lon array and same for lon array
    lon, lat = np.meshgrid(lons_shft, lats)
    x, y = mp(lon, lat)
    
    #Choose the settings for the coastlines, countries, meridians...
    mp.drawcoastlines(linewidth=0.2)
    mp.drawcountries(linewidth=0.2)
    
    meridians = mp.drawmeridians(np.arange(-180,200,20), 
                         labels=[False,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lon lines every 20º
    
    mp.drawparallels(np.arange(-90,110,20), 
                         labels=[True,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lat lines every 20º
    
    #Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
    colors = ["#ffffff", "#fec44f", "#d95f0e", "#e34a33", "#b30000"]
    cmap= matplotlib.colors.ListedColormap(colors)
    bounds = [0, 0.5, 1, 1.5, 2, 2.5]
    
    cmap.set_under("w")
    cmap.set_over("red")
    
    norm= matplotlib.colors.Normalize(vmin=0,vmax=2.5)
    
    #Time-average for the first emission point
    time_avg_flux = np.mean(net_flx_EP_shft, axis=0)
    
    #Plot the flux on the map
    sc2 = mp.pcolor(x, y, time_avg_flux*1000,
                    cmap=cmap, norm=norm, shading='auto')
    
    #Define colorbar features
    cb = fig.colorbar(sc2, ticks=bounds, extend='both', 
                      orientation='horizontal',fraction=0.052, 
                      pad=0.065)
    
    #Adjust colorbar tickmark size
    cb.ax.tick_params(labelsize=14)
    
    #Label the colorbar
    cb.set_label(label="Radiative Forcing from Short-term O$_3$ [mW·m$^{-2}$]",size=14,weight='bold')
    
    #Save and close the map plot
    plt.savefig("rad_fluxes_map_example.png",format="png",dpi=300)
    plt.show()
    plt.close()
"""
# =============================================================================
# PLOT TYPE 7 - Net radiative fluxes from short-term ozone increase with airparcel trajectory (Single EP)
# =============================================================================



#Requires lat,lon from ATTILA data files and fluxes
if attila_switch == True and rad_fluxes_switch == True:
    emission_point = 1
    #Set up axis object for plotting the map
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
    net_flx_EP_shft, lons_shft = shiftgrid(180.,global_net_flx[emission_point-1], 
                                           lons_0to36,start=False)
    
    #Format the lat and lon arrays for map graphing, 
    #makes lat array a lat x lon array and same for lon array
    lon, lat = np.meshgrid(lons_shft, lats)
    x, y = mp(lon, lat)
    
    #Choose the settings for the coastlines, countries, meridians...
    mp.drawcoastlines(linewidth=0.2)
    mp.drawcountries(linewidth=0.2)
    
    meridians = mp.drawmeridians(np.arange(-180,200,20), 
                         labels=[False,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lon lines every 20º
    
    mp.drawparallels(np.arange(-90,110,20), 
                         labels=[True,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lat lines every 20º
    
    #Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
    colors = ["#ffffff", "#fec44f", "#d95f0e", "#e34a33", "#b30000"]
    cmap= matplotlib.colors.ListedColormap(colors)
    bounds = [0, 0.5, 1, 1.5, 2, 2.5]
    
    cmap.set_under("w")
    cmap.set_over("red")
    
    norm= matplotlib.colors.Normalize(vmin=0,vmax=2.5)
    
    #Time-average for the first emission point
    time_avg_flux = np.mean(net_flx_EP_shft, axis=0)
    
    #Plot the flux on the map

    
    #for ID in range(50):
        #Plot a Lagrangian air parcel with parcel ID given by "parcel2"
        #ax.scatter(plon[:,(emission_point-1)*50+ID], plat[:,(emission_point-1)*50+ID], s=20, marker='o', color='black',
               #zorder=2,alpha = 0.05)
    
    sc2 = mp.pcolor(x, y, time_avg_flux*1000,
                    cmap=cmap, norm=norm, shading='auto')
    #Define colorbar features
    cb = fig.colorbar(sc2, ticks=bounds, extend='both', 
                      orientation='horizontal',fraction=0.052, 
                      pad=0.065)
    
    #Adjust colorbar tickmark size
    cb.ax.tick_params(labelsize=14)
    
    #Label the colorbar
    cb.set_label(label="Radiative Forcing from Short-term O$_3$ [mW·m$^{-2}$]",size=14,weight='bold')
    
    #Save and close the map plot
    plt.savefig("rad_fluxes_map_" + str(emission_point) + ".png",format="png",dpi=300)
    plt.show()
    plt.close()

# =================================================================================
# PLOT TYPE 8 - HORIZONTAL EVOLUTION OF LAGRANGIAN AIR PARCELS (ON MAP) W/ COLORBAR
# =================================================================================

#Requires ATTILA air parcel trajectory locatio and O3 data
if attila_switch == True and o3tracer_switch == True:
    emission_point = 1
    #Set up axis object for plotting the map
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
    
    #Format the lat and lon arrays for map graphing, 
    #makes lat array a lat x lon array and same for lon array
    lon, lat = np.meshgrid(lons_18to18, lats)
    x, y = mp(lon, lat)
    
    #Choose the settings for the coastlines, countries, meridians...
    mp.drawcoastlines(linewidth=0.2)
    mp.drawcountries(linewidth=0.2)
    
    meridians = mp.drawmeridians(np.arange(-180,200,20), 
                         labels=[False,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lon lines every 20º
    
    mp.drawparallels(np.arange(-90,110,20), 
                         labels=[True,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lat lines every 20º
    
    mp.fillcontinents(color='lightgray')
    
    #Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
    colors = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"]
    cmap= matplotlib.colors.ListedColormap(colors)
    bounds = [0, 15, 30, 45, 60, 75]
    
    cmap.set_under("w")
    cmap.set_over("crimson")
    
    norm= matplotlib.colors.Normalize(vmin=0,vmax=75)
    for i in range(50):
        #Plot a Lagrangian air parcel with parcel ID given by "parcel4"
        sc = ax.scatter(plon[:,(emission_point-1)*50+i][airO3_001[:,(emission_point-1)*50+i]*1E09>=30], plat[:,(emission_point-1)*50+i][airO3_001[:,(emission_point-1)*50+i]*1E09>30], s=20, marker='o',
               c=airO3_001[:,(emission_point-1)*50+i][airO3_001[:,(emission_point-1)*50+i]*1E09>30]*1E09,cmap=cmap,norm=norm,zorder=2,alpha=0.15)
    #Define colorbar features
    cb = fig.colorbar(sc, ticks=bounds, extend='both', 
                      orientation='vertical',fraction=0.04, pad=0.03)
    
    #Adjust colorbar tickmark size
    cb.ax.tick_params(labelsize=14)
    
    #Label the colorbar
    cb.set_label(label="O$_3$ Mixing Ratio [nmol·mol$^{-1}$]",size=14,weight='bold')
    
    #Save and close the map plot
    plt.savefig("emission point"+str((emission_point))+"_map_colorbar.png",format="png",dpi=300)
    plt.show()
    plt.close()


from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# Define the map projection and size
m = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-90, urcrnrlon=180, urcrnrlat=90, resolution='l')
fig = plt.figure(figsize=(10, 6))

x = plon[:,1]
y = plat[:,1]
print(x,y)
print(len(plon[:,1]))

# Define the animation function
def animate(i):
    plt.cla()  # clear the previous frame
    
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary()
    m.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180,180,60),labels=[0,0,0,1])
    for n in range(i):
        if n == i+1:
            m.scatter(x[n], y[n], latlon=True, s=10, color='r')
        else:
            m.scatter(x[n], y[n], latlon=True, s=10, color='r',alpha=0.5)
    print(i)
    plt.title('Animated Map')
    return

# Create the animation
anim = FuncAnimation(fig, animate, frames=len(plon[:,1]), interval=200)

# Save the animation as a GIF
anim.save('earthquake_map2.gif', writer='pillow', fps=5)
"""


#plot of RF result
flux_list = []
for n in range(28):
    flux = global_net_flx[n]
    flux_time_avg = np.mean(flux,axis=0)
    flux_list.append(np.sum(flux_time_avg))
flux_list = np.array(flux_list).reshape(7, 4,order="F") 
print(flux_list)
lat = np.linspace(85,25,7)  # define x as an array with 4 elements
lon = np.linspace(-115,-55,4)  # define y as an array with 7 elements
#plt.pcolormesh(lon,lat,flux_list, cmap='RdYlGn_r')
#plt.colorbar()
#plt.show()



#Set up axis object for plotting the map
fig, ax = plt.subplots() #Subplots are useful for drawing multiple plots together
    
#Adjust dimensions of map plot
fig.set_figheight(8)
fig.set_figwidth(14)
    
#Define map projection and settings
#For more info: https://matplotlib.org/basemap/users/cyl.html
mp = Basemap(projection = 'cyl', #equidistant cylindrical projection
                         llcrnrlon = -135,
                         llcrnrlat = 10,
                         urcrnrlon = -35,
                         urcrnrlat = 90,
                         resolution = 'i', ax=ax) #h=high, f=full, i=intermediate, c=crude
    
#Shift the fluxes from [0,360] to [-180,180]
#net_flx_EP_shft, lons_shft = shiftgrid(180.,global_net_flx[1], lons_0to36,start=False)
    
#Format the lat and lon arrays for map graphing, 
#makes lat array a lat x lon array and same for lon array
#lon, lat = np.meshgrid(lons_shft, lats)
x, y = mp(lon, lat)
    
#Choose the settings for the coastlines, countries, meridians...
mp.drawcoastlines(linewidth=0.2)
mp.drawcountries(linewidth=0.2)
    
meridians = mp.drawmeridians(np.arange(-180,200,20), 
                         labels=[False,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lon lines every 20º
    
mp.drawparallels(np.arange(-90,110,20), 
                         labels=[True,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lat lines every 20º
    
#Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
colors = ["#ffffff", "#fec44f", "#d95f0e", "#e34a33", "#b30000"]
cmap= matplotlib.colors.ListedColormap(colors)
    
cmap.set_under("w")
cmap.set_over("red")
    
   
#Plot the flux on the map
sc2 = mp.pcolor(x, y, flux_list,
                    cmap='RdYlGn_r',shading='auto')
    
#Define colorbar features
cb = fig.colorbar(sc2, extend='both', 
                     orientation='horizontal',fraction=0.052, 
                      pad=0.065)
    
#Adjust colorbar tickmark size
cb.ax.tick_params(labelsize=14)
    
#Label the colorbar
cb.set_label(label="Radiative Forcing due to emissions from every emission point",size=14,weight='bold')
    
    #Save and close the map plot
plt.show()
plt.close()