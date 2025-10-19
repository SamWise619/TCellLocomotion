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

angle_division = 2.0*math.pi / float(membrane_num_points)

for i in range(0, membrane_num_points):
    del_x = math.sin(i*angle_division) * membrane_radius
    del_y = membrane_radius - math.cos(i*angle_division) * membrane_radius
    
    current_point_pos = [initial_point_pos[0]+del_x, (1+initial_point_pos[1])-del_y]
    
    if i == 1:
        curr_dist = math.sqrt(del_x**2 + del_y**2)
        distance_between_points = curr_dist

    # shifted cell to avoid symmetry with obstacles.
    point_list.append([0, current_point_pos[0] - 0.45, current_point_pos[1] + 13.2])

print("dist between points: " + str(distance_between_points))

current_row = 0
COM = [0.0, 0.0] # center of cell.
for i in range(100, point_list.__len__()):
    COM[0] += point_list[i][1] # value in 2D array
    COM[1] += point_list[i][2]
COM_x = COM[0]/membrane_num_points
COM_y = COM[1]/membrane_num_points
#Exclude obstacles from inside the cell:
# Obstacles are placed at random locations according to a uniform distribution around the simulated cell placed in the middle of the simulation domain:
while(point_list.__len__() < num_excludedVolPoints+nucleus_membrane_num_points):
    rand_x = box_lx*random.random() + COM_x - 0.5*box_lx
    rand_y = box_ly*random.random() + COM_y - 0.5*box_ly

    distance_x = abs(rand_x - COM_x)
    distance_y = abs(rand_y - COM_y)
    distance_mag = math.sqrt(distance_x**2 + distance_y**2)

    if distance_mag >= 1.2*membrane_radius:
        point_list.append([1, rand_x, rand_y])

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
