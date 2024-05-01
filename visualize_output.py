import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
colors = [
    '#FF0000',  # Red
    '#00FF00',  # Green
    '#0000FF',  # Blue
    '#F033FF',  # Magenta
    '#33FFF6',  # Cyan
    '#FFC300',  # Yellow
    '#DAF7A6',  # Light Green
    '#581845',  # Purple
    '#C70039',  # Crimson
    '#FF5733',  # Red Orange
    '#900C3F',  # Dark Red
    '#FFC0CB',  # Pink
    '#9C27B0',  # Deep Purple
    '#3F51B5',  # Indigo
    '#03A9F4'   # Light Blue
]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
scatters = []
groups = {}

min_x = float('inf')
max_x = float('-inf')
min_y = float('inf')
max_y = float('-inf')
min_z = float('inf')
max_z = float('-inf')
with open('output.txt', 'r') as file:
    for line in file:
        parts = line.strip().split("(")
        
        
        group = parts[0]
       
        
        
        if group not in groups:
            groups[group] = {'xs': [], 'ys': [], 'zs': []}
        
        for i in range (1, len(parts)):
            parts_valid = parts[i].split(")")[0]
            coords = parts_valid.split(", ")
            x = float(coords[0].strip())
            y = float(coords[1].strip())
            z = float(coords[2].strip())
            groups[group]['xs'].append(x)
            groups[group]['ys'].append(y)
            groups[group]['zs'].append(z)
            min_x = min(min_x, x)
            max_x = max(max_x, x)
            min_y = min(min_y, y)
            max_y = max(max_y, y)
            min_z = min(min_z, z)
            max_z = max(max_z, z)



def update(frame):
    ax.clear()
    ax.set_title(f"Visualization of {list(groups.keys())[frame]}")
    coords = groups[list(groups.keys())[frame]]
    # Plot each coordinate with the color corresponding to its order across all groups
    # mean_x = sum(coords['xs']) / len(coords['xs'])
    # mean_y = sum(coords['ys']) / len(coords['ys'])
    # mean_z = sum(coords['zs']) / len(coords['zs'])
    for i in range(len(coords['xs'])):
        ax.scatter(coords['xs'][i], coords['ys'][i], coords['zs'][i], color=colors[i % len(colors)], label=f'Coordinate {i+1}')

    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    
    # Setting dynamic limits based on the data
    # ax.set_xlim(min(coords['xs']) - mean_x, max(coords['xs']) + mean_x)
    # ax.set_ylim(min(coords['ys']) - mean_y, max(coords['ys']) + mean_y)
    # ax.set_zlim(min(coords['zs']) - mean_z, max(coords['zs']) + mean_z)

    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.set_zlim(min_z, max_z)

    
    # Add a legend if there are multiple coordinates
    if len(coords['xs']) > 1:
        ax.legend()

    return scatters

ani = FuncAnimation(fig, update, frames=len(groups), repeat=False, interval=1000)

plt.show()