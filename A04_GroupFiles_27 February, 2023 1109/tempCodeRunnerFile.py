colors = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"]
    cmap= matplotlib.colors.ListedColormap(colors)
    bounds = [0, 15, 30, 45, 60, 75]
    
    cmap.set_under("w")
    cmap.set_over("crimson")
    
    norm= matplotlib.colors.Normalize(vmin=0,vmax=75)