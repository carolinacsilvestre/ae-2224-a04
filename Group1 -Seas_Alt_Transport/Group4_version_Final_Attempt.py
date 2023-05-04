#Import libraries
import glob #Dynamic file name loading
import numpy as np #Array processing
from netCDF4 import Dataset #NETCDF file handling
import matplotlib.pyplot as plt #Plotting
from mpl_toolkits.basemap import Basemap #For map plotting
from mpl_toolkits.basemap import shiftgrid #For shifting longitudes
import matplotlib.colors #To create new colorbar

# =============================================================================
# LOADING DATA FROM NETCDF (.NC) FILES
# =============================================================================

#USER INPUT - File path
#foldernamelist = ["C:/Users/31683/Desktop/project data/Summer200/*","C:/Users/31683/Desktop/project data/Summer250/*","C:/Users/31683/Desktop/project data/Summer300/*","C:/Users/31683/Desktop/project data/Winter200/*","C:/Users/31683/Desktop/project data/Winter250/*","C:/Users/31683/Desktop/project data/Winter300/*"]
foldernamelist = ["C:/Users/joren/Documents/project data/Winter/*","C:/Users/joren/Documents/project data/Summer/*"]
#Joren foldernamelist = ["C:/Users/joren/Documents/project data/Winter","C:/Users/joren/Documents/project data/Summer"]

#USER INPUT - Switches to determine which data types should be loaded
attila_switch = True
o3tracer_switch = False
rad_fluxes_switch = False

graphs28 = False
graphsperlat = True
MedianTrajectories = False
MedianTrajectoriesAll = False

step = 20
columns = 360/step
rows = 180/step

def TrendMap(EmissionPoint):

    parcel2A = np.arange(EmissionPoint*50,(EmissionPoint+1)*50-1) #[0]#
    TrendMapPlot = np.zeros((int(rows),int(columns)))
    
    for parcel2Ai in parcel2A:
        for point in range(len(plon[:,parcel2Ai])):
            for Nlong in range(int(columns)):
                if -180+(Nlong*step) <= plon[point,parcel2Ai] and plon[point,parcel2Ai] <= -180+((Nlong+1)*step):
                    for Nlat in range(int(rows)):
                        if 90-(Nlat*step) >= plat[point,parcel2Ai] and plat[point,parcel2Ai] >= 90-((Nlat+1)*step):
                            TrendMapPlot[Nlat][Nlong]+=1
    print(np.sum(TrendMapPlot))
    TrendMapPlot=np.flip(TrendMapPlot,0)
    return TrendMapPlot


for f_string in foldernamelist:

    print(f_string)
    #Read in file names based on f_string variable
    filenames_all = sorted(glob.glob(f_string)) #Get all file names in f_string
    print('\n') #Move to next line, improve readability
    print('Files in input folder: ')
    print('\n') #Check that all files present in folder are being printed
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
            if 'attila' in file:
                data = Dataset(file,'r')
                print('\n')
                print('File loaded: ')
                print(file)
                
                #Latitudes and longitudes
                lats = data.variables['lat'][:]
                lons_0to36 = data.variables['lon'][:] #Varies from 0 to 360 deg
                lons_18to18 = data.variables['lon'][:] #Will be converted later to -180 to 180 deg
                
                #changed
                #Time
                temp = data.variables['time'][:]
                time.append(temp)
            # print(temp)
                
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
    if rad_fluxes_switch == True: #and '250' in f_string:
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
                    global_net_flx.append((rad_flx_LW+rad_flx_SW)- \
                                        (rad_flx_LW_02+rad_flx_SW_02))
                    
        #Delete unnecessary variables
        del rad_flx_LW, rad_flx_SW, rad_flx_LW_02, rad_flx_SW_02

        #Delete unnecessary variables
    if attila_switch or o3tracer_switch or rad_fluxes_switch and not '250' in f_string:
        del temp, data


    #####################################################################################################################
    # MAKE  28 GRAPHS PLOT #
    
    #################################################################################################################
    if graphs28:
        fig, axs = plt.subplots(nrows=7, ncols=4, figsize=(32, 32))
        fig.tight_layout()
        axs = axs.transpose()
        for i, ax in enumerate(axs.flat):
            flux_list=TrendMap(i)
            lat = np.linspace(-90,90,int(rows)+1)  # define x as an array with 4 elementss
            lon = np.linspace(-180,180,int(columns)+1)
            #Set up axis object for plotting the map
            #fig, ax = plt.subplots(7, ncols=4, figsize=(16, 16)) #Subplots are useful for drawing multiple plots together=nrows= 7, ncols=4, figsize=(16, 16)
            
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
            #net_flx_EP_shft, lons_shft = shiftgrid(180.,global_net_flx[1], lons_0to36,start=False)
                
            #Format the lat and lon arrays for map graphing, 
            #makes lat array a lat x lon array and same for lon array
            #lon, lat = np.meshgrid(lons_shft, lats)
            x, y = mp(lon, lat)
                
            #Choose the settings for the coastlines, countries, meridians...
            mp.drawcoastlines(linewidth=0.2)
            mp.drawcountries(linewidth=0.2)
                

            #space between numbers
        ## meridians = mp.drawmeridians(np.arange(-180,200,45), 
            ##                        labels=[False,False,False,True], 
            ##                       linewidth=0.2, fontsize=10) #Draw lon lines every 20º
                
            ##mp.drawparallels(np.arange(-90,110,20), 
            ##                        labels=[True,False,False,True], 
            ##                        linewidth=0.2, fontsize=10) #Draw lat lines every 20º
                
            #Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
            colors = ["#ffffff", "#fec44f", "#d95f0e", "#e34a33", "#b30000"]
            cmap= matplotlib.colors.ListedColormap(colors)
                
            cmap.set_under("w")
            cmap.set_over("red")
                

            #Plot the flux on the map
            sc2 = mp.pcolor(x, y, flux_list, cmap='hot_r',shading='auto')
                
            ##Define colorbar features
            #cb = fig.colorbar(sc2, extend='both', 
            #                   orientation='horizontal',fraction=0.052, 
            #                  pad=0.065)
                
            ##Adjust colorbar tickmark size
            #cb.ax.tick_params(labelsize=14)
                
            ##Label the colorbar
            #cb.set_label(label="Heat map of all airparcels from the first emission point",size=14,weight='bold')
                
                #Save and close the map plot




            ax.set_ylabel(f'E. P. {i}', loc='top')
            #ax.set_xlabel('XLabel', loc='left')
            #cbar = fig.colorbar(sc)
            #cbar.set_label("ZLabel", loc='top')





        #Define colorbar features
        ##cb = fig.colorbar(sc2, extend='both', 
                            #  orientation='horizontal',fraction=0.052, 
                            #  pad=0.065)
                
        #Adjust colorbar tickmark size
        ##cb.ax.tick_params(labelsize=14)
                
        #Label the colorbar
        ##cb.set_label(label="Heat map of all airparcels from the first emission point",size=14,weight='bold',loc='right')

        plt.show()
        plt.close()

#####################################################################################################################
    # MAKE  GRAPHS PER LAT  #
    
#################################################################################################################
    if graphsperlat:
        fig, axs = plt.subplots(nrows=7, ncols=1, figsize=(32, 32))
        fig.tight_layout()
        for i, ax in enumerate(axs.flat):
            flux_list=TrendMap(i)
            flux_list+=TrendMap(i+7)
            flux_list+=TrendMap(i+14)
            flux_list+=TrendMap(i+21)

            lat = np.linspace(-90,90,int(rows)+1)  # define x as an array with 4 elementss
            lon = np.linspace(-180,180,int(columns)+1)
            #Set up axis object for plotting the map
            #fig, ax = plt.subplots(7, ncols=4, figsize=(16, 16)) #Subplots are useful for drawing multiple plots together=nrows= 7, ncols=4, figsize=(16, 16)
            
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
            #net_flx_EP_shft, lons_shft = shiftgrid(180.,global_net_flx[1], lons_0to36,start=False)
                
            #Format the lat and lon arrays for map graphing, 
            #makes lat array a lat x lon array and same for lon array
            #lon, lat = np.meshgrid(lons_shft, lats)
            x, y = mp(lon, lat)
                
            #Choose the settings for the coastlines, countries, meridians...
            mp.drawcoastlines(linewidth=0.2)
            mp.drawcountries(linewidth=0.2)
                

            #space between numbers
        ## meridians = mp.drawmeridians(np.arange(-180,200,45), 
            ##                        labels=[False,False,False,True], 
            ##                       linewidth=0.2, fontsize=10) #Draw lon lines every 20º
                
            ##mp.drawparallels(np.arange(-90,110,20), 
            ##                        labels=[True,False,False,True], 
            ##                        linewidth=0.2, fontsize=10) #Draw lat lines every 20º
                
            #Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
            colors = ["#ffffff", "#fec44f", "#d95f0e", "#e34a33", "#b30000"]
            cmap= matplotlib.colors.ListedColormap(colors)
                
            cmap.set_under("w")
            cmap.set_over("red")
                

            #Plot the flux on the map
            sc2 = mp.pcolor(x, y, flux_list, cmap='hot_r',shading='auto')
                
            ##Define colorbar features
            #cb = fig.colorbar(sc2, extend='both', 
            #                   orientation='horizontal',fraction=0.052, 
            #                  pad=0.065)
                
            ##Adjust colorbar tickmark size
            #cb.ax.tick_params(labelsize=14)
                
            ##Label the colorbar
            #cb.set_label(label="Heat map of all airparcels from the first emission point",size=14,weight='bold')
                
                #Save and close the map plot




            ax.set_ylabel(f'E. P. {i}', loc='top')
            #ax.set_xlabel('XLabel', loc='left')
            #cbar = fig.colorbar(sc)
            #cbar.set_label("ZLabel", loc='top')





        #Define colorbar features
        ##cb = fig.colorbar(sc2, extend='both', 
                            #  orientation='horizontal',fraction=0.052, 
                            #  pad=0.065)
                
        #Adjust colorbar tickmark size
        ##cb.ax.tick_params(labelsize=14)
                
        #Label the colorbar
        #graphtitle = 
        cb.set_label(label="Heat map of all airparcels from the first emission point",size=14,weight='bold',loc='right')

        plt.show()
        plt.close()

#####################################################################################################################
    # MAKE  plot/calc median  #
    
#################################################################################################################
    if MedianTrajectories == True:
            lonmedian = np.median(plon[:,:50],axis = 1)
            latmedian = np.median(plat[:,:50],axis = 1)
            lon180to180median = np.median(lons_18to18[:50])

            parcel2 = 25 #Parcel ID, 0 means first.
        
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
            lon, lat = np.meshgrid(lon180to180median, latmedian)
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
            
            ax.scatter(plon[:,:50], plat[:,:50], s=20, marker='o', color='green',
                    zorder=2)
            ax.scatter(lonmedian[:], latmedian[:], s=20, marker='o', color='red',
                    zorder=2)
            
            #Plot start and end points with an "S" and "F" respectively.
            ax.scatter(lonmedian[0], latmedian[0], s=140, marker='$S$', color='orange',
                    zorder=2)
            
            ax.scatter(lonmedian[-1], latmedian[-1], s=140, marker='$F$', color='orange',
                    zorder=2)
            
            #Save and close the map plot
            plt.savefig("air_parcel_ID"+str(parcel2)+"_map.png",format="png",dpi=300)
            plt.show()
            plt.close()


#####################################################################################################################
    # MAKE  plot/calc median ALL  #
    
#################################################################################################################

    if MedianTrajectoriesAll == True:
            
            parcel2 = 25 #Parcel ID, 0 means first.
            
            #Set up axis object for plotting the map
            fig, ax = plt.subplots() #Subplots are useful for drawing multiple plots together
            
            #Adjust dimensions of map plot
            fig.set_figheight(8)
            fig.set_figwidth(14)
            
            for k in range(0,28):
                lonmedian = np.median(plon[:,:50*k],axis = 1)
                latmedian = np.median(plat[:,:50*k],axis = 1)
                lon180to180median = np.median(lons_18to18[:50*k])

                
                
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
                lon, lat = np.meshgrid(lon180to180median, latmedian)
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
                

                ax.scatter(lonmedian[:], latmedian[:], s=20, marker='o', color='red',
                        zorder=2)
                
                #Plot start and end points with an "S" and "F" respectively.
                
                ax.scatter(lonmedian[0], latmedian[0], s=140, marker='$S$', color='orange',
                        zorder=2)
                
                ax.scatter(lonmedian[-1], latmedian[-1], s=140, marker='$F$', color='orange',
                        zorder=2)
                
            #Save and close the map plot
            #plt.savefig("air_parcel_ID"+str(parcel2)+"_map.png",format="png",dpi=300)
            plt.show()
            plt.close()