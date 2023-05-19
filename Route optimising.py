import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import numpy as np

data = pd.read_excel("RF.xlsx" ,engine='openpyxl',sheet_name='Data')
RF = np.array(data[data.columns[5]]).reshape(7, 4,order="F")
routes = pd.read_excel("C:/Users/twanv/OneDrive - Delft University of Technology/Bsc - 2/Q3/Project/Data/routes.xlsx",sheet_name=None,header=0,engine='openpyxl')

def route_points(name_route:str or list,path="C:/Users/twanv/OneDrive - Delft University of Technology/Bsc - 2/Q3/Project/Data/routes.xlsx"):
    def deg_to_dec(deg_array:np.array):
        dec_array = [int(deg_array[i][0]) + int(deg_array[i][1]) / 60 + int(deg_array[i][2]) / 3600 for i in range(len(deg_array))]
        return dec_array
    route_data = routes[name_route]
    indices = route_data.columns
    x_coor, y_coor = np.array(route_data[indices[2]]) , np.array(route_data[indices[1]])
    x_coor, y_coor = [x_coor[i][2:].replace("'", '').replace('°','.').split('.') for i in range(len(x_coor))] , [y_coor[i][2:].replace("'", '').replace('°','.').split('.') for i in range(len(y_coor))]
    x_coor_dec , y_coor_dec = deg_to_dec(x_coor) , deg_to_dec(y_coor)
    return np.array(x_coor_dec) , np.array(y_coor_dec) , route_data[indices[0]]

def interpolate_data(RF_array:np.array):
    y = np.linspace(-85,-25, 7)  # define x as an array with 4 elements
    x = np.linspace(-115, -55, 4)  # define y as an array with 7 elements

    #create interpolant
    data = RF.T
    interp = RegularGridInterpolator((x, y), data, bounds_error=False, fill_value=None)
    return interp

def interpolated_plot(interp,refine:int=25):
    xx = np.linspace(-125, -45, 4*refine)
    yy = np.linspace(15,90, 7*refine)
    X, Y = np.meshgrid(xx, yy, indexing='ij')

    # Set up axis object for plotting the map
    fig, ax = plt.subplots()  # Subplots are useful for drawing multiple plots together

    # Adjust dimensions of map plot
    fig.set_figheight(8)
    fig.set_figwidth(14)

    # Define map projection and settings
    # For more info: https://matplotlib.org/basemap/users/cyl.html
    mp = Basemap(projection='cyl',  # equidistant cylindrical projection
                 llcrnrlon=-135,
                 llcrnrlat=10,
                 urcrnrlon=-35,
                 urcrnrlat=90,
                 resolution='i', ax=ax)  # h=high, f=full, i=intermediate, c=crude
    x, y = mp(xx, yy)

    # Choose the settings for the coastlines, countries, meridians...
    mp.drawcoastlines(linewidth=0.2)
    mp.drawcountries(linewidth=0.2)

    # Draw lon lines every 20º
    meridians = mp.drawmeridians(np.linspace(-115, -55, 4), labels=[False, False, False, True], linewidth=0.2,
                                 fontsize=10)
    # Draw lat lines every 20º
    mp.drawparallels(np.linspace(25,85, 7), labels=[True, False, False, True], linewidth=0.2, fontsize=10)

    # Plot the flux on the map
    sc2 = mp.pcolor(x, y, interp((X, -Y)).T, shading='auto',cmap='Reds')
    mp.colorbar(sc2)

interp = interpolate_data(RF)
interpolated_plot(interp)

names = ["JFK-LAX" ]#, "JFK-ORD" ,"LAX-ORD","ATL-MCO","LAS-LAX"]
for name in names:
    x_route, y_route, waypoints = route_points(name)
    plt.plot(-x_route,y_route,label=name)
plt.legend()
plt.show()





