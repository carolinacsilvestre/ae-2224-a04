import numpy as np
import matplotlib.pyplot as plt
  
data = np.random.random(( 12 , 12 ))
plt.imshow( data , cmap = 'autumn' , interpolation = 'nearest', aspect=10 )
  
plt.title( "2-D Heat Map" )
plt.show()