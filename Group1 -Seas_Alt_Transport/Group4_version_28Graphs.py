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
f_string =r"C:\Users\moheb\Desktop\DATA_PROJ_Q3\Summer\250hpa\*"#C:\Users\moheb\Desktop\Q3_Proj (Group Git)\*" #'P:/AE2224I_GroupA4/250hPa/NAmerica/201407/*' #Insert file path to input data, do not forget wildcard

#USER INPUT - Switches to determine which data types should be loaded
attila_switch = True
o3tracer_switch = False
rad_fluxes_switch = False

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

            print("HIIIIIIIII")
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
'''''  
#Fluxes from the VISO submodel, for 200 and 300 hPa
if rad_fluxes_switch == True and not '250' in f_string:
    #Start at 3 since first two calls are not emission points
    for ep in range(3,31):
        for file in filenames_all:
            if 'viso' in file:
                data = Dataset(file,'r')
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
        
        print('Net flux for EP'+str(ep-2)+' loaded...')
        
    del rad_flx_LW, rad_flx_SW, rad_flx_LW_02, rad_flx_SW_02
'''''
#Delete unnecessary variables
if attila_switch or o3tracer_switch or rad_fluxes_switch and not '250' in f_string:
    del temp, data



###########

#parcel10=0

#print(plon[0,parcel10],plat[0,parcel10])
#why  40 not correct????
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

#hello
print(TrendMap(0))





###########
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
    
    parcel = 25 #Parcel ID, 0 means first.
    

    setOFcolor=np.array(range(len(ppress[:,parcel])))
    print(setOFcolor)
    #Scatter plot settings

    #ADDED
    plt.scatter(time, ppress[:,parcel], s=30, marker='o', c=ppress[:,parcel]) #color='green')
   
  
   #BEFORE ADDED plt.scatter(time, ppress[:,parcel], s=30, marker='o', color='green')
    #END ADDED

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
# PLOT TYPE 2 - HORIZONTAL EVOLUTION OF LAGRANGIAN AIR PARCELS (ON MAP)
# =============================================================================

#Requires ATTILA air parcel trajectory location data
if attila_switch == True:


    #for i in range(50):

   # parcel2 = 25 #Parcel ID, 0 means first.
    

    #ADDED LATER
     #Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
   ## colors = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"]
   ## cmap= matplotlib.colors.ListedColormap(colors)
   ## bounds = [0, 15, 30, 45, 60, 75]
    #END ADDED LATER

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
  #before ADEED  ax.scatter(plon[:,parcel2], plat[:,parcel2], s=20, marker='o', color='green',
       #        zorder=2) #added c=setOFcolor
    
    for i in range(1):
        parcel2 = i
        ax.scatter(plon[:,parcel2], plat[:,parcel2], s=20, marker='o', c=ppress[:,parcel],
                zorder=2)

        #print(plon[:,parcel2].shape)

        ax.scatter(plon[-1,parcel2], plat[-1,parcel2], s=140, marker='$F$', color='red',
                zorder=2)
        
        #Plot start and end points with an "S" and "F" respectively.
        ax.scatter(plon[0,parcel2], plat[0,parcel2], s=140, marker='$S$', color='red',
                zorder=2)

#original
  ##  parcel2=50
 ##   ax.scatter(plon[:,parcel2], plat[:,parcel2], s=20, marker='o', c=ppress[:,parcel],
  ##              zorder=2)
  ##  ax.scatter(plon[-1,parcel2], plat[-1,parcel2], s=140, marker='$F$', color='red',
  ##             zorder=2)


  ##  #Plot start and end points with an "S" and "F" respectively.
  ##  ax.scatter(plon[0,parcel2], plat[0,parcel2], s=140, marker='$S$', color='red',
  ##             zorder=2)


#added later
     #Plot a Lagrangian air parcel with parcel ID given by "parcel2"
    ##ax.scatter(plon[:,parcel22], plat[:,parcel2], s=20, marker='o', color='green',
     ##          zorder=2)
    
    #Plot start and end points with an "S" and "F" respectively.
   ## ax.scatter(plon[0,parcel22], plat[0,parcel2], s=140, marker='$S$', color='red',
    ##           zorder=2)
    
   ## ax.scatter(plon[-1,parcel22], plat[-1,parcel2], s=140, marker='$F$', color='red',
     ##          zorder=2)
#end added later

    #Save and close the map plot
    plt.savefig("air_parcel_ID"+str(parcel2)+"_map.png",format="png",dpi=300)
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
    
    parcel3 = 25 #Parcel ID, 0 means first.
    
    #Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
    colors = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"]
    cmap= matplotlib.colors.ListedColormap(colors)
    bounds = [0, 15, 30, 45, 60, 75]
    
    cmap.set_under("w")
    cmap.set_over("crimson")
    
    norm= matplotlib.colors.Normalize(vmin=0,vmax=75)
    
    #Scatter plot command
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

# =================================================================================
# PLOT TYPE 4 - HORIZONTAL EVOLUTION OF LAGRANGIAN AIR PARCELS (ON MAP) W/ COLORBAR
# =================================================================================

#Requires ATTILA air parcel trajectory locatio and O3 data
if attila_switch == True and o3tracer_switch == True:

    parcel4 = 25 #Parcel ID, 0 means first.
    
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
    net_flx_EP_shft, lons_shft = shiftgrid(180.,global_net_flx[0], 
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
    colors = ["#ffffff", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"]
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


'''
##################
import numpy as np
import matplotlib.pyplot as plt
  
data = TrendMap(0)
#Adjust dimensions of map plot

# we can use differenrt cmaps (initial one =>'autumn' )
m=plt.imshow( data , cmap = 'hot' , interpolation = 'nearest',aspect="auto" )




#3
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
#3

  
plt.title( "2-D Heat Map" )
plt.show()
'''


############################################################28 by 28 fig############################################
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
        


    vmin=0
    vmax=700
    #Plot the flux on the map
    sc2 = mp.pcolor(x, y, flux_list, cmap='hot_r',shading='auto',vmin=vmin, vmax=vmax)#,)
        
    ##Define colorbar features
    #cb = fig.colorbar(sc2, extend='both', 
     #                   orientation='horizontal',fraction=0.052, 
      #                  pad=0.065)
        
    ##Adjust colorbar tickmark size
    #cb.ax.tick_params(labelsize=14)
        
    ##Label the colorbar
    #cb.set_label(label="Heat map of all airparcels from the first emission point",size=14,weight='bold')
        
        #Save and close the map plot




    ax.set_ylabel(f'E. P. {i+1}', loc='top')
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
############################################################end############################################



'''
#################################Single Heat Map (airparcel trajectory can be added)####################################### 
fig, ax = plt.subplots()

flux_list=TrendMap(0)
print(f"this is wdsfdpso{flux_list}")
lat = np.linspace(-90,90,int(rows)+1)  # define x as an array with 4 elementss
lon = np.linspace(-180,180,int(columns)+1)
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
    

x, y = mp(lon, lat)
    
#Choose the settings for the coastlines, countries, meridians...
mp.drawcoastlines(linewidth=0.2)
mp.drawcountries(linewidth=0.2)
    
#Set up custom colorbar, colors may be chosen with the help from colorbrewer2.org
colors = ["#ffffff", "#fec44f", "#d95f0e", "#e34a33", "#b30000"]
cmap= matplotlib.colors.ListedColormap(colors)
    




cmap.set_under("w")
cmap.set_over("red")
    
vmin=0
vmax=300
#Plot the flux on the map
sc2 = mp.pcolor(x, y, flux_list, cmap='hot_r',shading='auto',vmin=vmin, vmax=vmax)
    
meridians = mp.drawmeridians(np.arange(-180,200,20), 
                         labels=[False,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lon lines every 20º
    
mp.drawparallels(np.arange(-90,110,20), 
                         labels=[True,False,False,True], 
                         linewidth=0.2, fontsize=10) #Draw lat lines every 20º


######draws the parcels themselves on top of them w start and end point and color bar#######

for i in range(50):
        parcel2 = i
        ax.scatter(plon[:,parcel2], plat[:,parcel2], s=20, marker='o', c=ppress[:,parcel],
                zorder=2)

        #print(plon[:,parcel2].shape)

        ax.scatter(plon[-1,parcel2], plat[-1,parcel2], s=140, marker='$F$', color='red',
                zorder=2)
        
        #Plot start and end points with an "S" and "F" respectively.
        ax.scatter(plon[0,parcel2], plat[0,parcel2], s=140, marker='$S$', color='red',
                zorder=2)
        

cb = fig.colorbar(sc2, extend='both', 
                        orientation='horizontal',fraction=0.052, 
                        pad=0.065)
        
#Adjust colorbar tickmark size
cb.ax.tick_params(labelsize=14)
        
#Label the colorbar
cb.set_label(label="Heat map of airparcel 0 from the first emission point (E.P. 0)",size=14,weight='bold')

######END draws the parcels themselves on top of them w start and end point and color bar#######

plt.show()
plt.close()

################################# END ####################################### 
'''


print("DONE!")
