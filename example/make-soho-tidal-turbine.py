# Creates a turbine blade
import sys
import numpy
import glob
# Import other self written modules:
import TGAM as tgam

# Blade cross section data files:
filenames = glob.glob('soho_aerof*.txt')
sorted_filenames = sorted(filenames)
# Blade geometry options:
# Total length of blade:
factor = 1.0
blade_length = factor*0.340 #mm length, and a radius of 400mm
# Array of coordinates where we want cross-sections:
extrusion_length = [-factor*(0 + i*0.02) for i in range(18)]
extrusion_length = numpy.around(extrusion_length, decimals=2)
R = 0.400
c_over_R = [0.2425, 0.2105, 0.185, 0.1588, 0.141, 0.1263, 0.115, 0.1053, 0.0978, 0.0908, 0.0845, 0.079, 0.0733, 0.0685, 0.065, 0.0625, 0.06, 0.0578]
extrusion_scale = [factor*cR * R for cR in c_over_R]
extrusion_scale = numpy.around(extrusion_scale, decimals=4)
pitch_angles = [27.3, 23.1, 20.0, 16.3, 13.8, 12.6, 11.0, 9.5, 8.2, 7.4, 6.8, 5.8, 5.0, 4.5, 4.0, 3.6, 3.3, 3.1]

# Turbine options:
num_blades = 3
turbine_type = 'windturbine'
spin_direction = 'clockwise'
# Hub settings:
blade_hub_connection_circle_radius=factor*0.015
blade_hub_connection_distance=factor*0.02
hub_radius=factor*0.050
ellipse_tip_distance=factor*0.06
hub_length = factor*(0.577 + (0.577*0.25))
hub_height = factor*(2-1.06)
base_radius = factor*0.04
base_length = factor*1.06
base_direction = 'y' # accepts: '+y', '-y', 1, -1
blade_base_distance = factor*0.577
# Meshing options:
meshingscheme = 'tetmesh'
sizingfunction = 3
meshingoption = 'auto'
# Turbine array options:
num_turbines = 1
turbine_placement = [[0,0,0]]



# Initialize class to generate the wanted geometry/mesh:
bladegenerator = tgam.TGAM(
                            sorted_filenames, dimension=3, rotation=pitch_angles, 
                            extrusion_scale=extrusion_scale, extrusion_length=extrusion_length,
                            meshingscheme=meshingscheme, sizingfunction=sizingfunction, 
                            num_blades=num_blades, turbine_type=turbine_type, spin_direction=spin_direction,
                            blade_hub_connection_circle_radius=blade_hub_connection_circle_radius,
                            blade_hub_connection_distance=blade_hub_connection_distance, hub_radius=hub_radius, 
                            ellipse_tip_distance=ellipse_tip_distance, hub_length=hub_length, hub_height=hub_height, 
                            base_radius=base_radius, base_length=base_length, base_direction=base_direction,
                            blade_base_distance=blade_base_distance, meshingoption=meshingoption,
                            num_turbines=num_turbines, turbine_placement=turbine_placement
                          )

# Generate the cubit journal file to generate the entire turbine array geometry/mesh files:
bladegenerator.generate_geometry()
