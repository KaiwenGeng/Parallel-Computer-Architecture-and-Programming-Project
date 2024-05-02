import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation 
import re
def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]
plt.ioff()
# colors = [
#     '#FF0000',  # Red
#     '#00FF00',  # Green
#     '#0000FF',  # Blue
#     '#F033FF',  # Magenta
#     '#33FFF6',  # Cyan
#     '#FFC300',  # Yellow
#     '#DAF7A6',  # Light Green
#     '#581845',  # Purple
#     '#C70039',  # Crimson
#     '#FF5733',  # Red Orange
#     '#900C3F',  # Dark Red
#     '#FFC0CB',  # Pink
#     '#9C27B0',  # Deep Purple
#     '#3F51B5',  # Indigo
#     '#03A9F4'   # Light Blue
# ]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
scatters = []
groups = {}
color_maps = []
min_x = float('inf')
max_x = float('-inf')
min_y = float('inf')
max_y = float('-inf')
min_z = float('inf')
max_z = float('-inf')
with open('position_result.txt', 'r') as file:
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
sorted_group = list(groups.keys())
sorted_group.sort(key=natural_sort_key)
norm = plt.Normalize(-0.05, 0.05)
cmap = plt.get_cmap('coolwarm')
with open('weight_charge_result.txt', 'r') as file:
    for line in file:
        parts = line.strip().split("(")
        for i in range (1, len(parts)):
            parts_valid = parts[i].split(")")[0]
            weight_charge = parts_valid.split(", ")
            weight = float(weight_charge[0].strip())
            charge = float(weight_charge[1].strip())
            color_maps.append((weight, charge))


def update(frame):
    ax.clear()
    ax.set_title(f"Visualization of {sorted_group[frame]}")
    coords = groups[sorted_group[frame]]
    for i in range(len(coords['xs'])):
        w,c = color_maps[i]
        color = cmap(norm(c))
    
        # Scatter plot for each coordinate with size scaled by weight and color by charge
        ax.scatter(coords['xs'][i], coords['ys'][i], coords['zs'][i], color=color, 
                s=w*10,  # Scale size by weight; adjust the factor as needed
                label=f'Coordinate {i+1}')
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.set_zlim(min_z, max_z)
    
    # Add a legend if there are multiple coordinates
    # if len(coords['xs']) > 1:
    #     ax.legend()
    return scatters
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.set_zlabel('Z Coordinate')
ani = FuncAnimation(fig, update, frames=len(groups), repeat=False, interval=1000)
writergif = animation.PillowWriter(fps=12)
ani.save('filename.gif',writer=writergif)
plt.show()