import matplotlib.pyplot as plt
import scipy.interpolate
from mpl_toolkits.basemap import Basemap
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import numpy as np

dataset = 3

data = pd.read_excel("RF.xlsx" ,engine='openpyxl',sheet_name='Data')
RF = np.array(data[data.columns[dataset]]).reshape(7, 4,order="F")
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

def plot_rf(interp=None,refine:int=5):
    xx = np.linspace(-115, -55, 4)
    yy = np.linspace(85,25, 7)
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
                 llcrnrlat=25,
                 urcrnrlon=-70,
                 urcrnrlat=55,
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
    sc2 = mp.pcolor(x, y, RF , shading='auto',cmap='binary') #interp((X, -Y)).T
    mp.colorbar(sc2,location='bottom')

plot_rf()

def total_rf_per_route(route, RF,n_points):
    # grid
    x_grid = np.linspace(-115, -55, 4)
    y_grid = np.linspace(85, 25, 7)

    x_discrete = np.linspace(x_route[0], x_route[-1], n_points)
    y_discrete = route(x_discrete)

    RF_impact = 0

    for i in range(n_points):
        for n in range(4):
            lon = x_grid[n]
            for m in range(7):
                lat = y_grid[m]
                if ((lon - 10 < -x_discrete[i] < lon + 10) and (lat - 5 < y_discrete[i] < lat + 5)):
                    RF_impact += RF[m,n] / n_points

    #plt.scatter(-x_discrete, f_inter(x_discrete))
    return RF_impact

impact_list = []
names = ["JFK-LAX", "JFK-ORD" ,"LAX-ORD","ATL-MCO","LAS-LAX"]
for name in names:
    x_route, y_route, waypoints = route_points(name)
    f_inter = scipy.interpolate.interp1d(x_route,y_route,'linear')
    impact_list.append( total_rf_per_route(f_inter,RF,100))
    plt.plot(-x_route, y_route, label=name)
plt.title(data.columns[dataset])
plt.legend()
plt.savefig(f"route_visualization{data.columns[dataset]}.png")
plt.show()
plt.close()

plt.bar(names,impact_list)
plt.savefig(f"Barplot{data.columns[dataset]}.png")
plt.show()






