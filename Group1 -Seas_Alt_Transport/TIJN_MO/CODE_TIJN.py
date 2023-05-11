 #Import libraries
import glob #Dynamic file name loading
import numpy as np #Array processing
from netCDF4 import Dataset #NETCDF file handling
import matplotlib.pyplot as plt #Plotting
from mpl_toolkits.basemap import Basemap #For map plotting
from mpl_toolkits.basemap import shiftgrid #For shifting longitudes
import matplotlib.colors #To create new colorbar
import matplotlib.ticker as mticker


# =============================================================================
# LOADING DATA FROM NETCDF (.NC) FILES
# =============================================================================

#USER INPUT - File path
f_string =r"C:\Users\moheb\Desktop\Q3_Proj (Group Git)\*"#""C:\Users\twanv\OneDrive - Delft University of Technology\Bsc - 2\Q3\Project\Data\250_july\*" #'P:/AE2224I_GroupA4/250hPa/NAmerica/201407/*' #Insert file path to input data, do not forget wildcard

#USER INPUT - Switches to determine which data types should be loaded
attila_switch = True
o3tracer_switch = True
rad_fluxes_switch = False
load_03_all_eps = True

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
airO3_008 = [] #O3lg in mol/mol for emission point 7
airO3_015 = [] #O3lg in mol/mol for emission point 15
airO3_022 = [] #O3lg in mol/mol for emission point 22
rad_flx_SW = [] #Holds 1 specific call for SW radiative fluxes
rad_flx_LW = [] #Holds 1 specific call for LW radiative fluxes
rad_flx_SW_02 = [] #SW flux for call 02
rad_flx_LW_02 = [] #LW flux for call 02
global_net_flx = [] #Holds all net fluxes for the 28 EPs (3 months)
labels = [f"{lat+20} <-> {lat}" for lat in range(70,-110,-20)]

#Load relevant variables from NETCDF files
#Positions of air parcels
if attila_switch == True:
    for file in filenames_all:
        if 'attila.nc' in file:
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
if o3tracer_switch == True and load_03_all_eps:
    for file in filenames_all:
        if 'O3lg' in file:
            data = Dataset(file,'r')
            print('\n')
            print('File loaded: ')
            print(file)
            
            #O3 along air parcels for each emission point
            temp = data.variables['airO3_001'][:,:1400] #Only the first 1400 trajectories
            airO3_001.append(temp)
            temp = data.variables['airO3_008'][:,:1400] #Only the first 1400 trajectories
            airO3_008.append(temp)
            temp = data.variables['airO3_015'][:,:1400] #Only the first 1400 trajectories
            airO3_015.append(temp)
            temp = data.variables['airO3_022'][:,:1400] #Only the first 1400 trajectories
            airO3_022.append(temp)
            
            #NOTE: do not forget there are 28 emission points in total!
    
    #Concatenation of variables, lists become multi-dimensional numpy arrays
    airO3_001 = np.concatenate(airO3_001, axis=0)
    airO3_008 = np.concatenate(airO3_008, axis=0)
    airO3_015 = np.concatenate(airO3_015, axis=0)
    airO3_022 = np.concatenate(airO3_022, axis=0)
    
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



###########

#parcel10=0

#print(plon[0,parcel10],plat[0,parcel10])
#why  40 not correct????
step = 20
columns = 360/step
rows = 180/step

def TrendMap(EmissionPoint):

    parcel2A = np.arange(EmissionPoint*50,(EmissionPoint+1)*50) #[0]#
    TrendMapPlot = np.zeros((int(rows),int(columns)))
    
    for parcel2Ai in parcel2A:
        for point in range(len(plon[:,parcel2Ai])):
            for Nlong in range(int(columns)):
                if -180+(Nlong*step) <= plon[point,parcel2Ai] and plon[point,parcel2Ai] <= -180+((Nlong+1)*step):
                    for Nlat in range(int(rows)):
                        if 90-(Nlat*step) >= plat[point,parcel2Ai] and plat[point,parcel2Ai] >= 90-((Nlat+1)*step):
                            if EmissionPoint == 0 and airO3_001[point,parcel2Ai] > 0.2e-9:
                                TrendMapPlot[Nlat][Nlong]+=1
                            if EmissionPoint == 7 and airO3_008[point,parcel2Ai] > 0.2e-9:
                                TrendMapPlot[Nlat][Nlong]+=1
                            if EmissionPoint == 14 and airO3_015[point,parcel2Ai] > 0.2e-9:
                                TrendMapPlot[Nlat][Nlong]+=1
                            if EmissionPoint == 21 and airO3_022[point,parcel2Ai] > 0.2e-9:
                                TrendMapPlot[Nlat][Nlong]+=1
    #print(np.sum(TrendMapPlot))
    #TrendMapPlot=np.flip(TrendMapPlot,0)
    return TrendMapPlot

print("ghello")
#print(TrendMap(0))
#print(TrendMap(1))
#Total_trendmap=np.zeros(TrendMap(0).shape)

"""
#print(TrendMap(0)+TrendMap(1))
for i in range(28):
    Total_trendmap=Total_trendmap+TrendMap(i)
#print(Total_trendmap)

#Main_trendmap=TrendMap(0)
#print(Main_trendmap)

New_main_trandmap=np.sum(Total_trendmap,axis=1)
#print(New_main_trandmap)
plt.barh(np.flip(labels), np.flip(New_main_trandmap))
plt.title("All emission points, all altitudes [summer]")
plt.show()
"""

#####first ALT (0,7,14, 21)
customrange=np.array([0,7,14,21])
Total_trendmap = np.zeros((9,18))
for i in (customrange):
    Total_trendmap=Total_trendmap+TrendMap(i)
print(Total_trendmap)
Total_trendmap = 100*Total_trendmap/np.sum(Total_trendmap)


New_main_trandmap=np.sum(Total_trendmap,axis=1)
#print(New_main_trandmap)
print(f"st row of emission points, all altitudes [summer]{New_main_trandmap}")
plt.barh(np.flip(labels), np.flip(New_main_trandmap))
plt.title("1st row of emission points, all altitudes [summer]")
plt.xticks(ticks=[n for n in range(0,60,10)], labels=[f'{n}%' for n in range(0,60,10)])
plt.savefig(f"With O3 filter for points {customrange}")
plt.show()
###########
#####first ALT (1,8,15, 22)

"""
customrange=np.array([1,8,15, 22])
for i in (customrange):
    Total_trendmap=Total_trendmap+TrendMap(i)
#print(Total_trendmap)


New_main_trandmap=np.sum(Total_trendmap,axis=1)
#print(New_main_trandmap)

labels = ["Lat 70 - 50", "Lat 70-50", "Lat 50-30", "Lat 30-10", "Lat 10- -10", "Lat -10 - -30","Lat -30- -50","Lat -50- -70","Lat 70-90"]
plt.barh(np.flip(labels), np.flip(New_main_trandmap))
plt.title("2nd row of emission points, all altitudes [summer]")
plt.show()



################
customrange=np.array([2,9,16, 23])
for i in (customrange):
    Total_trendmap=Total_trendmap+TrendMap(i)

New_main_trandmap=np.sum(Total_trendmap,axis=1)

labels = ["Lat 70 - 50", "Lat 70-50", "Lat 50-30", "Lat 30-10", "Lat 10- -10", "Lat -10 - -30","Lat -30- -50","Lat -50- -70","Lat 70-90"]
plt.barh(np.flip(labels), np.flip(New_main_trandmap))
plt.title("3rd row of emission points, all altitudes [summer]")
plt.show()
##################
customrange=np.array([3,10,17, 24])
for i in (customrange):
    Total_trendmap=Total_trendmap+TrendMap(i)
#print(Total_trendmap)


New_main_trandmap=np.sum(Total_trendmap,axis=1)
#print(New_main_trandmap)

labels = ["Lat 90-70", "Lat 70-50", "Lat 50-30", "Lat 30-10", "Lat 10- -10", "Lat -10 - -30","Lat -30- -50","Lat -50- -70","Lat 70-90"]
plt.barh(np.flip(labels), np.flip(New_main_trandmap))
plt.title("4th row of emission points, all altitudes [summer]")
plt.show()
#####################
customrange=np.array([4,11,18, 25])
for i in (customrange):
    Total_trendmap=Total_trendmap+TrendMap(i)
#print(Total_trendmap)


New_main_trandmap=np.sum(Total_trendmap,axis=1)
#print(New_main_trandmap)

labels = ["Lat 90-70", "Lat 70-50", "Lat 50-30", "Lat 30-10", "Lat 10- -10", "Lat -10 - -30","Lat -30- -50","Lat -50- -70","Lat 70-90"]
plt.barh(np.flip(labels), np.flip(New_main_trandmap))
plt.title("5th row of emission points, all altitudes [summer]")
plt.show()
###################
customrange=np.array([5,12,19, 26])
for i in (customrange):
    Total_trendmap=Total_trendmap+TrendMap(i)
#print(Total_trendmap)


New_main_trandmap=np.sum(Total_trendmap,axis=1)
#print(New_main_trandmap)

labels = ["Lat 90-70", "Lat 70-50", "Lat 50-30", "Lat 30-10", "Lat 10- -10", "Lat -10 - -30","Lat -30- -50","Lat -50- -70","Lat 70-90"]
plt.barh(np.flip(labels), np.flip(New_main_trandmap))
plt.title("6th row of emission points, all altitudes [summer]")
plt.show()
###################
customrange=np.array([6,13,20, 27])
for i in (customrange):
    Total_trendmap=Total_trendmap+TrendMap(i)
#print(Total_trendmap)


New_main_trandmap=np.sum(Total_trendmap,axis=1)
#print(New_main_trandmap)

labels = ["Lat 90, 70", "Lat 70, 50", "Lat 50, 30", "Lat 30, 10", "Lat 10, -10", "Lat -10, -30","Lat -30, -50","Lat -50, -70","Lat -70, -90"]
plt.barh(np.flip(labels), np.flip(New_main_trandmap))
plt.title("7th row of emission points, all altitudes [summer]")
plt.show()
"""


print("DONE!")