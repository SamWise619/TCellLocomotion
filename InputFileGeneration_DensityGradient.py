# Sami Alawadhi, David M. Rutkowski (2025) 
# Vavylonis Group
# Department of Physics, Lehigh University
# If you use any part of the scripts in this package, please cite:
# TBD

import math
import random

box_lx = 434.78
box_ly = 434.78

eqBondLength = 0.2730435

nucleus_num_points = 100
membrane_num_points = 200
nucleus_radius = (0.5*eqBondLength)/(math.sin(math.pi/nucleus_num_points)) #radius of nuc
membrane_radius = (0.5*eqBondLength)/(math.sin(math.pi/membrane_num_points)) #radius of CM
nucleus_membrane_num_points=nucleus_num_points + membrane_num_points # total number of beads for the cell

print(nucleus_radius)
print("diameter = ", 2*membrane_radius)
print("Area_Cell = ", math.pi*membrane_radius*membrane_radius)

#change amount of obstacles total.
num_excludedVolPoints = 3200 # total number of obstacle beads
initial_point_pos = [0.0, 0.0]

angle_division = 2.0*math.pi/float(nucleus_num_points)
point_list = []

distance_between_points = 0.0

for i in range(0, nucleus_num_points):
    del_x = math.sin(i*angle_division) * nucleus_radius
    del_y = nucleus_radius - math.cos(i*angle_division) * nucleus_radius
    
    current_point_pos = [initial_point_pos[0]+del_x, initial_point_pos[1]-del_y]
    
    if i == 1:
        curr_dist = math.sqrt(del_x**2 + del_y**2)
        distance_between_points = curr_dist

    # Shifting nucleus to avoid symmetry with obstacles.
    point_list.append([0, current_point_pos[0] - 0.45, current_point_pos[1] + 9.8])

print("dist between points: " + str(distance_between_points))

angle_division = 2.0*math.pi / float(membrane_num_points) #place the particles around the cortex perimeter

for i in range(0, membrane_num_points):
    del_x = math.sin(i*angle_division) * membrane_radius
    del_y = membrane_radius - math.cos(i*angle_division) * membrane_radius
    
    current_point_pos = [initial_point_pos[0]+del_x, (1+initial_point_pos[1])-del_y]
    
    if i == 1:
        curr_dist = math.sqrt(del_x**2 + del_y**2)
        distance_between_points = curr_dist

    # Shifting cell to avoid symmetry with obstacles.
    point_list.append([0, current_point_pos[0] - 0.45, current_point_pos[1] + 13.2])

print("dist between points: " + str(distance_between_points))

# periodic obstacles
region_x_ev_points = box_lx # max x
min_x = 0

region_y_ev_points = box_ly - nucleus_radius*2.0 # max y
min_y = 0

x_0 = -box_lx/2
y_0 = -17.39
pts_in_row = 28 
f = 0.5
d_Excl = 0.2608696
d_x = f * (box_lx/pts_in_row - d_Excl) + d_Excl
d_y = math.sqrt(d_x**2 - (d_x**2)/4)
current_row = 0
for i in range(0, num_excludedVolPoints):
    #Grid of increasing obstacle density:
    point_list.append([1, x_0, y_0])
    x_0 = x_0 + d_x
    if x_0 > box_lx / 2:
        if current_row % 2 == 0:
            x_0 = -box_lx/2 + random.random()*d_x
            y_0 = y_0 - d_y
        else:
            x_0 = -box_lx/2 + random.random()*d_x
            y_0 = y_0 - d_y
        current_row = current_row + 1
        f = f - 0.008 # narrowing horizontal spacing in each row - corresponds to f_n = 1 - n*epsilon in paper (Methods)

        #d_x corresponds to dy_n in the paper (Methods) because the grid was rotated in the simulations to show the cell migrating through the density grid in a horizontal manner.
        d_x = f * (box_lx/pts_in_row - d_Excl) + d_Excl

        #d_y corresponds to dx_n in the paper (Methods) because the grid was rotated in the simulations to show the cell migrating through the density grid in a horizontal manner.
        d_y = (math.sqrt(3)/2)*d_x

# input file generation to be inputted into input.txt (positions, bonds, angles):
output_file_name = 'Positions.xyz'

with open(output_file_name, 'w') as fp:
    fp.write('{0}\n'.format(len(point_list)))
    fp.write('t=0.0\n')
    
    for i in range(0, len(point_list)):
        fp.write('{0} {1} {2} {3} {4}\n'.format(point_list[i][0], i, point_list[i][1], point_list[i][2], 0.0))
     
     
outputbnd_file_name = 'Bonds.bnd'

with open(outputbnd_file_name, 'w') as fp:
    fp.write('{0}\n'.format(nucleus_membrane_num_points))
    fp.write('t=0.0\n')
    
    for i in range(0, nucleus_num_points-1):
        fp.write('{0} {1} {2}\n'.format(0, i, i+1))
    #connect last and first beads
    fp.write('{0} {1} {2}\n'.format(0, 0, nucleus_num_points-1))

    for i in range(nucleus_num_points, nucleus_membrane_num_points-1):
        fp.write('{0} {1} {2}\n'.format(1, i, i+1))
    #connect last and first beads
    fp.write('{0} {1} {2}\n'.format(1, nucleus_num_points, nucleus_membrane_num_points-1))
        
        
outputbnd_file_name = 'Angles.ang'

with open(outputbnd_file_name, 'w') as fp:
    fp.write('{0}\n'.format(nucleus_membrane_num_points))
    fp.write('t=0.0\n')
    
    for i in range(0, nucleus_num_points-2):
        fp.write('{0} {1} {2} {3}\n'.format(0, i, i+1, i+2))

    fp.write('{0} {1} {2} {3}\n'.format(0, nucleus_num_points-2, nucleus_num_points-1, 0))
    fp.write('{0} {1} {2} {3}\n'.format(0, nucleus_num_points-1, 0, 1))
    
    for i in range(nucleus_num_points, nucleus_membrane_num_points-2):
        fp.write('{0} {1} {2} {3}\n'.format(0, i, i+1, i+2))

    fp.write('{0} {1} {2} {3}\n'.format(0, nucleus_membrane_num_points-2, nucleus_membrane_num_points-1, nucleus_num_points))
    fp.write('{0} {1} {2} {3}\n'.format(0,  nucleus_membrane_num_points-1, nucleus_num_points, nucleus_num_points+1))
    