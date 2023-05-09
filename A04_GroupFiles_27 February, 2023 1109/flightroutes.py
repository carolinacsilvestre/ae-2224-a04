import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as Basemap
import numpy as np 
from mpl_toolkits.basemap import shiftgrid #For shifting longitudes
import matplotlib.colors #To create new colorbar


df = pd.read_excel('C:/Users/Carolina Silvestre/Desktop/dataproject/routes.xlsx', sheet_name='JFK-LAX')
