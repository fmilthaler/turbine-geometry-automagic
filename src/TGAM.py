#File: TGAM.py
#Copyright (C) 2013 Frank Milthaler.
#
#This file is part of turbine-geometry-automagic.
#
#turbine-geometry-automagic is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#turbine-geometry-automagic is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with turbine-geometry-automagic. If not, see <http://www.gnu.org/licenses/>.

# This python class automagically writes a "journal" script which contains all the
# Cubit commands necessary to generate the geometry and mesh of a (wind/tidal) turbine.
# The generated geometry/mesh is simplified such that it does not replace the hard work
# associated with a full CAD model with all its bearing etc, but is useful when doing
# numerical modelling. 
# The user can set various parameters of the turbine, e.g. height of the hub, 
# twist angles of the blade, length of the blade, and of course the aerofoil 
# cross sections of a blade. For a full list of parameters, check the documentation,
# or examples.
import sys
import os
from math import sqrt, pi, sin, cos, asin


class TGAM:

  def __init__(self, textfilename, dimension=3, rotation=[0.0], blade_placement=None, extrusion_scale=[1.0], extrusion_length=[0.0], blade_crosssection_type=None, blade_curve_type=None, num_blades=None, blade_hub_connection_circle_radius=None, blade_hub_connection_distance=None, hub_radius=None, ellipse_tip_distance=None, hub_length=None, hub_height=None, base_length=None, base_direction='-y', base_radius=None, blade_base_distance=None, meshingscheme='tetmesh', sizingfunction=5, meshingoption='auto', cad_package='cubit', turbine_type=None, num_turbines=1, turbine_placement=None, spin_direction='counterclockwise'):
      # Blade geometry settings:
      if (not (isinstance(textfilename, list) or isinstance(textfilename, str))):
          print "ERROR: the input argument 'textfilename' with the coordinate of the aerofoil section(s) was not given as a string or list."
          raise SystemExit()
      self.textfilename = textfilename
      self.dimension = dimension # default is 3, which will construct an entire turbine, 2 will only draw the blade cross section
      self.global_scale = extrusion_scale[0]
      self.rotation = rotation
      self.extrusion_scale = extrusion_scale
      self.extrusion_length = extrusion_length
      if (dimension == 3 and num_blades is None): self.num_blades = 3
      elif (dimension == 2 and num_blades is None): self.num_blades = 1
      else: self.num_blades = num_blades
      self.cad_package = cad_package

      # Check the cross section curve type for the blade:
      if (blade_crosssection_type is None):
          self.blade_crosssection_type = 'linear'
      else:
          if (not (blade_crosssection_type == 'linaer' or blade_crosssection_type == 'spline')):
              print "ERROR: The 'blade_crosssection_type' must be 'linear' or 'spline'. Exiting..."
              raise SystemExit()
          else:
              if (blade_crosssection_type == 'spline'):
                  print "Warning: The spline curve is currently not fully supported for connecting the points of a cross section."
                  print "For now, use 'linear' which is the default option."
                  raise SystemExit()
              self.blade_crosssection_type = blade_crosssection_type
      # Check the curve connection type for the blade:
      if (blade_curve_type is None):
          self.blade_curve_type = 'linear'
      else:
          if (not (blade_curve_type == 'linear' or blade_curve_type == 'spline')):
              print "ERROR: The 'blade_curve_type' must be 'linear' or 'spline'. Exiting..."
              raise SystemExit()
          else:
              self.blade_curve_type = blade_curve_type
      # For 2D geometries, blade cross sections:
      if (blade_placement is None):
          self.blade_placement = [[0,0,0]]
      else:
          self.blade_placement = blade_placement
      # Connecting volume blade - hub:
      if (blade_hub_connection_circle_radius is None): self.blade_hub_connection_circle_radius = self.global_scale*0.25
      else: self.blade_hub_connection_circle_radius = blade_hub_connection_circle_radius
      if (blade_hub_connection_distance is None): self.blade_hub_connection_distance = abs(extrusion_length[-1]*0.1)# 10% of blade's length
      else: self.blade_hub_connection_distance = blade_hub_connection_distance
      # For the ellipsoid volume, the tip of the hub:
      if (hub_radius is None): self.hub_radius = self.global_scale*0.45
      else: self.hub_radius = hub_radius
      if (ellipse_tip_distance is None): self.ellipse_tip_distance = self.global_scale*0.75
      else: self.ellipse_tip_distance = ellipse_tip_distance
      if (hub_length is None): self.hub_length = self.global_scale*3
      else: self.hub_length = hub_length
      if (hub_height is None): self.hub_height = abs(extrusion_length[-1] * 3.5) # default value, 3times as high as each blade is long!
      else: self.hub_height = abs(hub_height)
      if (base_length is None): self.base_length = self.hub_height
      else: self.base_length = base_length
      self.base_direction = base_direction
      if (base_radius is None): self.base_radius = self.hub_radius
      elif (base_radius <= self.hub_radius): self.base_radius = base_radius
      else:
          print "=================================================================================================="
          print "ERROR: The radius of the base structure must be less than or equal to the hub radius"
          raise SystemExit()
      # If base_direction is -y or -1, then the base_length must be equal to hub_height
      if (self.base_direction == '-y' or self.base_direction == -1 or self.base_direction == '-1'):
          if (not self.hub_height == self.base_length):
              print "ERROR: The base structure is extruded from the rotor level downwards, but the parameters of 'hub_height' and _base_length' are unequal."
              print "These parameters must be equal, otherwise the base structure is above the ground level (above 0, 0 being the ground)"
              raise SystemExit()
      # Meshing settings:
      self.meshingscheme = meshingscheme
      self.sizingfunction = sizingfunction
      self.meshingoption = meshingoption
      # Geometry settings: Wind or tidal turbine
      if (not (turbine_type=='windturbine' or turbine_type=='tidalturbine')):
          self.turbine_type = None
      else:
          self.turbine_type = turbine_type
      self.num_turbines = num_turbines
      if (self.dimension == 3):
          if (turbine_placement is None):
              self.turbine_placement = [0, 0, 0]
          else:
              self.turbine_placement = turbine_placement
      else:
          self.turbine_placement = turbine_placement

      # Terminate program if given arguments are of unequal length:
      if (self.dimension == 3 and (len(rotation) != len(extrusion_scale) or len(rotation) != len(extrusion_length) or len(extrusion_scale) != len(extrusion_length))):
          print "=================================================================================================="
          print "ERROR, given arguments are not of the same length!\nAborting program..."
          print "len(rotation) = ", len(rotation)
          print len(extrusion_length)
          print len(extrusion_scale)
          raise SystemExit()
      # Or, if using different aerofoils per blade, also check if the list of aerofoils is of equal lengths to the other lists:
      if (isinstance(textfilename,list)):
          if (len(textfilename) != len(rotation)):
              print "=================================================================================================="
              print "ERROR, given arguments are not of the same length!\nAborting program..."
              raise SystemExit()

      # Give error if extrusion length is in the wrong z-direction:
      if (self.dimension == 3 and spin_direction == 'clockwise'):
          if (True in (t > 0 for t in extrusion_length)):
              print "=================================================================================================="
              print "Error: The blade must be extruded in the -z direction for clockwise spinning turbine. Provide "
              print "negative real numbers for the input argument 'extrusion_length', or switch to counterclockwise "
              print "rotation."
              raise SystemExit()
      elif (self.dimension == 3 and spin_direction == 'counterclockwise'):
          if (True in (t < 0 for t in extrusion_length)):
              print "=================================================================================================="
              print "Error: The blade must be extruded in the -z direction for clockwise spinning turbine. Provide "
              print "negative real numbers for the input argument 'extrusion_length', or switch to counterclockwise "
              print "rotation."
              raise SystemExit()
      # If rotation angle/twist is unacceptable, give the user an error:
      if (True in (t<-90 or t>90 for t in rotation)):
          print "=================================================================================================="
          print "Error: What are you trying to do? The twist angles must be real numbers between -90 and 90."
          raise SystemExit()
      # Give error message, if an unacceptable combination of base radius, hub length and blade_hub_connection_circle_radius was given:
      if (self.dimension == 3 and self.hub_length - self.blade_hub_connection_circle_radius*2*1.5 - self.base_radius*2 <= self.base_radius):
          print "=================================================================================================="
          print "Error: The base of the turbine cannot be placed, please increase the 'hub_length' or decrease either the 'hub_radius' or 'blade_hub_connection_circle_radius'."
          raise SystemExit()
      # If however the distance between the blade and base structure is specified, check again:
      if (self.dimension == 3 and not blade_base_distance is None):
          if (self.hub_length - self.blade_hub_connection_circle_radius*1.5 - self.base_radius - blade_base_distance <= 0.0 or blade_base_distance - self.blade_hub_connection_circle_radius*1.5 - self.base_radius*1.5 <= 0):
              print "=================================================================================================="
              print "Error: The given distance between blades and the base structure is not large enough based on the given hub_length or blade_hub_connection_circle_radius and base_radius."
              raise SystemExit()
      self.blade_base_distance = blade_base_distance
      # Error if hub_height is too short:
      if (self.dimension == 3 and self.hub_height <= extrusion_length[-1]):
          print "=================================================================================================="
          print "Error, your blades will hit the ground, set your hub_height sufficiently high."
          print "in order to NOT smash your blades into pieces!"
      # Error for wrong meshing option:
      if (not (self.meshingoption == 'auto' or self.meshingoption == 'adaptive' or self.meshingoption == 'size')):
          print "=================================================================================================="
          print "Error, wrong meshing option found. Choose one of the following:"
          print "  auto"
          print "  adaptive"
          print "  size"
          raise SystemExit()
      # Exit if both, blade_placement and turbine_placement are given:
      if (not (blade_placement is None) and not (turbine_placement is None)):
          print "=================================================================================================="
          print "Error: Both 'blade_placement' and 'turbine_placement' are given."
          print "For 2D cross sections/2D geometries, specify 'blade_placement'"
          print "and for 3D geometries (turbines) specify 'turbine_placement'"
          raise SystemExit()
      # Exit if 3D and gmsh was selected:
      if (cad_package != 'cubit' and dimension == 3):
          print "=================================================================================================="
          print "Error, currently only cubit is supported to automatically generate a 3D geometry/mesh"
          print "Choose: cad_package='cubit'"
          raise SystemExit()
      if (self.dimension == 2 and ((self.num_blades != len(self.blade_placement)) or self.num_blades != len(self.rotation) )):
          print "=================================================================================================="
          print "Error, given number of blades is unequal to the length of the"
          print "corresponding 'blade_placement' or 'rotation' vector"
          #raise SystemExit()
      if (self.dimension == 3 and not (spin_direction == 'clockwise' or spin_direction == 'counterclockwise')):
          print "=================================================================================================="
          print "Error, the input argument 'spin_direction' must be one of the following:"
          print "  * 'clockwise' or"
          print "  * 'counterclockwise'"
          print "without the single quotes."
          raise SystemExit()
      else:
          self.spin_direction = spin_direction

      # Create class variable that will hold the cubit/gmsh commands:
      self.cadlines = []
      # A variable for the blade design name:
      self.geometryname = ''
      # And a variable for the number of points in the given textfile:
      self.num_2d_coords = 0
      self.num_cs = len(extrusion_length)-1
      # And the coordinates itself:
      #self.coordinates = []
      # Depending on which cad software was specified, set the comment character:
      if (cad_package.lower() == 'cubit'):
          self.commentchar = '#'
          self.cubitcmdopen = 'cubit.cmd("'
          self.cubitcmdclose = '")'
      elif (cad_package.lower() == 'gmsh'):
          self.commentchar = '//'
      else:
          print "Error: No valid cad_package given, only 'cubit' or 'gmsh' are currently allowed!"
          raise SystemExit()



  def list_id_increase(self, l, num_entries):
  # Function to increase a list by given number of entries:
  # Eases the trouble we are facing when keeping track of vertex/curve/surface/volume lists below:
    """ This routine takes in a list 'l' and adds 'num_entries' to it, whereas they are increased by 1
         Input: 
           l: list of integers
           num_entries: Integer, that many entries will be appended to l, each entry is by 1 greater
             than the previous entry.
         Output:
           l: list of integers
    """
    for i in range(num_entries):
        l.append(l[-1]+1)
    return l


  def query_entity_id_cubit(self, var_name, type, offset=0):
      """ This routine writes command to query vertex/curve/surface/volume IDs 
          and stores them in given variable names. These can be used in the online
          operation to manipulate the geometry, e.g. move a certain surface.
          Input:
           var_name: String of the variable name that should be used in the online operation
           type: String of either 'vertex', 'curve', 'surface', 'volume
           offset: optional integer argument. By default the last ID number is stored in the 
             variable 'var_name'.
      """
      if (not (type in ['vertex', 'curve', 'surface', 'volume'])):
          print "FATAL ERROR: The function 'query_entity_id_cubit' was called with an invalid "
          print "  value for the input argument 'type'."
          print "  The value of 'type' is: "+type
          print "Please send this errormessage with your setup files to the developer."
          raise SystemExit()
      if (offset>=0): sign='+'
      else: sign=''
      # aprepro command:
      self.cadlines.append(self.cubitcmdopen+'${'+var_name+'=Id(\''+type+'\') '+sign+str(offset)+' }'+self.cubitcmdclose)
      self.cadlines.append(var_name+' = int(cubit.get_aprepro_value_as_string(\''+var_name+'\'))')


  ########################################################################
  # Extract 2D coordinates and information of the initial cross section: #
  ########################################################################
  def get_initial_cross_section_coordinates(self, filename, global_scale=None):
      """ This routines reads in the 2D coordinates of the cross section
          of a blade. It is assumed the first line of the given file contains
          the name of the blade design. The coordinates are returned in form
          of a list.
          Input:
           filename: String of the filename with the coordinates of the blade's
             cross section
           global_scale: Real number, that is applied as a scaling factor to the
             given 2D geometry in the file 'filename'.
          Output:
           coordinates: 2 dimensional list of coordinates of blade's cross section
      """
      if (global_scale is None):
          global_scale = self.global_scale

      # Start processing the file:
      infile = open(filename, 'r')
      # Assume that the first line always contains the name of the geometry:
      self.geometryname = infile.readline().strip()
      # Start with header of cad file:
      self.cadlines.append(self.commentchar+' Geometry: '+self.geometryname+'\n')
      # Now remove special characters from geometryname, such that we can use it for a filename:
      self.geometryname = self.geometryname.replace('=','').replace('.','-').replace(' ','_')
      if (self.cad_package == 'gmsh'):
          self.cadlines.append(self.commentchar+' Element edge length:')
          self.cadlines.append('el='+str(self.sizingfunction)+';\n')
      self.cadlines.append(self.commentchar+' Starting geometry, points:')
      # Extract coordinates from Text file:
      coordinates = self.get_cross_section_coordinates(filename, global_scale=global_scale)
      # Return coordinates:
      return coordinates


  ################################################
  # Extract 2D coordinates of the cross section: #
  ################################################
  def get_cross_section_coordinates(self, filename, global_scale=None):
      """ This routines reads in the 2D coordinates of the cross section
          of a blade. It is assumed the first line of the given file contains
          the name of the blade design. The coordinates are returned in form
          of a list.
          Input:
           filename: String of the filename with the coordinates of the blade's
             cross section
           global_scale: Real number, that is applied as a scaling factor to the
             given 2D geometry in the file 'filename'.
          Output:
           coordinates: 2 dimensional list of coordinates of blade's cross section
      """
      if (global_scale is None):
          global_scale = self.global_scale

      # Start processing the file:
      infile = open(filename, 'r')
      # Check if the first line contains data or a comment/text:
      firstline = infile.readline().strip()
      try:
          if (len(firstline.split()) == 2):
              float(firstline.split()[0])
              float(firstline.split()[1])
              # if all that was successful, close the file, and open it again:
              infile.close()
              infile = open(filename, 'r')
          else:
              raise Exception
      except:
          None # nothing to do here

      # Start processing the content of the file:
      # Extract coordinates from Text file:
      # Variables for point ids:
      pointid = 0
      coordinates = []
      # Process text file and get coordinates:
      for line in infile:
          # Getting rid of newline character at the end:
          line = line.strip()
          if (line == '' or line.isspace()):
              continue
          # increase pointid:
          pointid += 1
          coord = []
          for linecoord in line.split():
              # Get coordinate bit by bit, and apply global scaling factor to it:
              coord.append(float(linecoord)*global_scale)
          # At this point, all coordinates from text have been assembled,
          # now add as many zero coordinates for the missing dimensions:
          for i in range(3-len(line.split())):
              coord.append(0.0)
          # Now three coordinates for this line have been assembled, add them
          # to the global list of coordinates:
          coordinates.append(coord)
      # Done processing text input file:
      infile.close()
      # Deleting last entry in coordinates only, iff it has the same coordinates as the first entry:
      if (coordinates[0] == coordinates[-1]): del coordinates[-1];
      # Give the user a warning, if too many coordinates were found:
      if (len(coordinates) >= 99):
          print "-----------------------------------------------------------------------------"
          print "| WARNING: More than 98 coordinates were found for the blade cross-section. |"
          print "| This most likely will cause problems for the Bladegenerator.              |"
          print "| If the turbine could not be generated, try to reduce the number of        |"
          print "|  coordinates for the cross-section to be less than 98 points.             |"
          print "-----------------------------------------------------------------------------"
      return coordinates


  def write_2d_vertices(self, coordinates, vertex_list, placement, rotation):
    """ A higher-level interface to write 2D coordinates of a blade cross section:
        Input:
          coordinates: list of coordinates for a cross section of a blade.
          vertex_list: list of vertices generated in cubit
          placement: 1D-list/vector by which the points should
            be translated
        Output:
          vertex_list: list of vertices generated in cubit
    """
    if (self.cad_package == 'cubit'):
        vertex_list = self.write_2d_vertices_cubit(coordinates, vertex_list)
    else:
        for i in range(self.num_blades):
            vertex_list = self.write_2d_vertices_gmsh(coordinates, vertex_list, placement[i], rotation[i])
    return vertex_list


  def write_2d_vertices_cubit(self, coordinates, vertex_list):
    """ This routine writes the cubit commands to set the coordinates
        of the 2D cross section of the blade.
        Input:
          coordinates: list of coordinates for a cross section of a blade.
          vertex_list: list of vertices generated in cubit
        Output:
          vertex_list: list of vertices generated in cubit
    """
    i = 0
    for coord in coordinates:
        # Form new cadline:
        cadline = 'create vertex '
        for x in coord:
            cadline = cadline+str(x)+' '
        # Add assembled string for defining a point to cadlines:
        self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
        # Increase counter and add it (as vertex id) to vertex_list:
        i += 1
        vertex_list.append(i)
    return vertex_list


  def write_2d_vertices_gmsh(self, coordinates, vertex_list, placement, rotation):
    """ This routine writes the cubit commands to set the coordinates
        of the 2D cross section of the blade.
        Input:
          coordinates: list of coordinates for a cross section of a blade.
          vertex_list: list of vertices generated in cubit
          placement: 1D-list/vector by which the points should
            be translated
          rotation: Angle by which the blade should be rotated
        Output:
          vertex_list: list of vertices generated in cubit
    """
    # Variable for point ids:
    if (len(vertex_list) == 0): pointid = 0
    else: pointid = vertex_list[-1]
    if (rotation%360.0 != 0.0):
        # Get center point of blade:
        M = self.get_blade_centre_point(coordinates)
        # It is already scaled, but we still have to move it:
        M[0] = M[0] + placement[0]
        M[1] = M[1] + placement[1]
        M[2] = M[2] + placement[2]
    for coord in coordinates:
        # increase pointid:
        pointid += 1
        cadline = 'Point ('+str(pointid)+') = {'
        # Here we extract the x/y/z coordinates, apply the wanted scaling
        # and placement:
        geo_coord = []
        for i in range(len(coord)):
            scaled_moved_coord = float(coord[i]) + placement[i]
            geo_coord.append(scaled_moved_coord)
        if (rotation%360.0 != 0.0):
            # Get rotated coordinates:
            geo_coord = self.rotate_blade_point(rotation, M, geo_coord)
        for i in range(len(geo_coord)):
            cadline = cadline+str(round(geo_coord[i], 5))+', ' # rounding to 5 decimals
        # At this point, all coordinates have been assembled,
        # now add as many zero coordinates for the missing dimensions:
        for i in range(3-len(coord)):
            cadline = cadline+str(0.0)+', '
        # At this point, 3 coordinates for this point have been assembled, now
        # append the desired edge length to it:
        cadline = cadline+'el};'
        self.cadlines.append(cadline)
        # Append pointid to vertex_list:
        vertex_list.append(pointid)
    return vertex_list


  def get_blade_centre_point(self, coordinates):
    """
        This subroutines computes the centre point of given coordinates
    """
    M = []
    maxx = 0; minx = 0;
    maxy = 0; miny = 0;
    maxz = 0; minz = 0;
    for coord in coordinates:
        if (coord[0] > maxx):
            maxx = coord[0]
        elif (coord[0] < minx):
            minx = coord[0]
        if (coord[1] > maxy):
            maxy = coord[1]
        elif (coord[1] < miny):
            miny = coord[1]
        if (coord[2] > maxz):
            maxz = coord[2]
        elif (coord[2] < minz):
            minz = coord[2]
    mx = (maxx-minx)/2.0 + minx
    my = (maxy-miny)/2.0 + miny
    mz = (maxz-minz)/2.0 + minz
    M.append(mx); M.append(my); M.append(mz)
    return M


  def rotate_blade_point(self, dphi, M, coord):
    """ This subroutine rotates a given point around the rotationl 
        axis in M by a given angle.
         Input:
           dphi: Floating point number of the angle [in Degrees]
           M: vector of the rotational axis around which the coordinates
             are rotated
           coord: Vector of a point, which is rotated around M
         Output:
           coord: The new coordinates of the given point
    """
    oldcoord = []
    for i in range(len(coord)):
        oldcoord.append(coord[i])
    oldx = oldcoord[0]; oldy = oldcoord[1];
    mx = M[0]; my = M[1]
    # Converting angle 'dphi' from degree into radian:
    dphi = dphi * pi/180.0
    # First: compute current relative position of the given point:
    px = oldx-mx; py = oldy-my
    # The radius (length):
    r = sqrt(px*px + py*py)
    # The current angle in the polar coordinate system:
    if (r == 0.0): # New positions are old positions
        nx = oldx;
        ny = oldy;
    else:
        if (px>=0):
            phi = asin(py/r)
        else:
            phi = pi - asin(py/r)
        # New coordinates:
        nx = r * cos( phi+dphi ) + mx;
        ny = r * sin( phi+dphi ) + my;
    if (len(oldcoord) == 2):
        coord = [nx, ny]
    elif (len(oldcoord) == 3):
        coord = [nx, ny, oldcoord[2]]
    return coord


  ############################################################
  # Creating intermediate cross sections (points) of blade:  #
  ############################################################
  def create_scaled_cross_sections_cubit(self, coordinates, extrusion_scale, extrusion_length, vertex_list):
    """ This routines assembles the cubit commands for creating the scaled cross sections at
        different z-coordinates (depths).
         Input:
           coordinates: list of coordinates of the original cross section.
           extrusion_scale: scaling factors.
           extrusion_length: z-coordinate of the cross sections.
           vertex_list: list of vertices generated in cubit
         Output:
           vertex_list: list of vertices generated in cubit
    """
    # To make a copy of the surface, create new vertices:
    self.cadlines.append('\n'+self.commentchar+' Generate another aerofoil cross section, by either using the same aerofoil design or different one (if specified):')
    self.cadlines.append(self.commentchar+" Here we create the vertices again, because the number of the vertices is not definitely defined when copying the surface by using Cubit's function:")
    for coord in coordinates:
        # Form new cadline:
        cadline = 'create vertex '
        for i in range(len(coord)):
            if (i == 2):
                # For z-coordinate, use given location:
                cadline = cadline+str(extrusion_length)
            else:
                cadline = cadline+str(coord[i]*extrusion_scale)+' '
        # Add assembled string for defining a point to cadlines:
        self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
        # Add new vertex id to vertex_list:
        vertex_list.append(vertex_list[-1] + 1)

#    # Create one more, that is scaled down to the length equal to the diameter of the circle that then sits on the hub:
#    self.cadlines.append('\n'+self.commentchar+' Make copy of surface, by creating the same points:')
#    self.cadlines.append(self.commentchar+" This surface will be between the blade itself and the circle that will sit on the hub of the turbine:")
#    for coord in coordinates:
#          # Form new cadline:
#          cadline = 'create vertex '
#          for i in range(len(coord)):
#              if (i == 2):
#                  # For z-coordinate, use given location:
#                  cadline = cadline+str(self.blade_hub_connection_distance*0.8-self.hub_radius)
#              else:
#                  cadline = cadline+str(coord[i]*(2*self.blade_hub_connection_circle_radius/self.global_scale))+' '
#          # Add assembled string for defining a point to cadlines:
#          self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
#          # Add new vertex id to vertex_list:
#          vertex_list.append(vertex_list[-1] + 1)
    return vertex_list


  ###############################################################
  # Creating lines/curves and surfaces for blade cross section: #
  ###############################################################
  def create_lines_surface_cross_sections(self, num_points, num_cs, vertex_list, curve_list, surface_list):
    """ A higher-level interface to generate lines (linear/spline) to connect the created points:
        Input:
           num_points: Integer, number of points of original cross section
           num_cs: Integer, number of additional cross sections per blade
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
           surface_list: list of surface generated in cubit
         Output:
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
           surface_list: list of surface generated in cubit
    """
    if (self.cad_package == 'cubit'):
        if (self.blade_crosssection_type == 'linear'):
            (vertex_list, curve_list) = self.create_linear_lines_cross_sections_cubit(num_points, num_cs, vertex_list, curve_list)
            surface_list = self.create_surface_cross_sections_cubit(num_points, num_cs, surface_list)
        #elif (self.blade_crosssection_type == 'spline'):
        #    (vertex_list, curve_list) = self.create_spline_lines_cross_sections_cubit(num_points, num_cs, vertex_list, curve_list)
        #    # Here we used a spline curve to connect points of the blade cross sections, thus the 'num_points' is ALWAYS
        #    # 2, but in the linear case, it is the number of points/coordinates from the data/text file.
        #    surface_list = self.create_surface_cross_sections_cubit(2, num_cs, surface_list)
    else:
        for i in range(self.num_blades):
            curve_list = self.create_linear_lines_cross_sections_gmsh(i+1, num_points, vertex_list, curve_list)
            surface_list = self.create_surface_cross_sections_gmsh(i+1, num_points, surface_list)
    return (vertex_list, curve_list, surface_list)


  def create_linear_lines_cross_sections_cubit(self, num_points, num_cs, vertex_list, curve_list):
    """ This routine connects all points of each cross section 
        with linear lines and assembles the corresponding commands for cubit
         Input:
           num_points: Integer, number of points of original cross section
           num_cs: Integer, number of additional cross sections per blade
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
         Output:
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
    """
    # Now assemble Cubit journal entries to define lines, assume that 
    # the given coordinates in the textfile were immediate neightbours
    # from one line to the next
    self.cadlines.append('\n'+self.commentchar+' Define Lines/Curves:')
    for i in range(1, num_points+1):
        # Linear lines for the cross section:
        if (i < num_points):
            self.cadlines.append(self.cubitcmdopen+'create curve vertex '+str(i)+' '+str(i+1)+self.cubitcmdclose)
        else:
            # For the last point, close the loop:
            self.cadlines.append(self.cubitcmdopen+'create curve vertex '+str(i)+' 1'+self.cubitcmdclose)
        # Add id to curve_list:
        curve_list.append(i)
    curve_list.append(1); curve_list.append(2)

    # Do the same for the extruded surfaces and the surface between the blade and the circle sitting on the hub:
    for j in range(1, num_cs+1):
      self.cadlines.append('\n'+self.commentchar+' Define Lines/Curves for the extruded surface j='+str(j)+':')
      for i in range(1, num_points+1):
          if (i < num_points):
              self.cadlines.append(self.cubitcmdopen+'create curve vertex '+str(num_points*j+i)+' '+str(num_points*j+i+1)+self.cubitcmdclose)
          else:
              # For the last point, close the loop:
              self.cadlines.append(self.cubitcmdopen+'create curve vertex '+str(num_points*j+i)+' '+str(num_points*j+1)+self.cubitcmdclose)
          # Add id to curve_list:
          curve_list.append(num_points*j+i)
    
    # Stupid cubit creates new vertices on top of old vertices when creating lines....:
    for i in range(len(vertex_list)):
        vertex_list[i] = vertex_list[-1] + i
    # Plus one:
    vertex_list.append(vertex_list[-1] + 1)
    return vertex_list, curve_list


  def create_spline_lines_cross_sections_cubit(self, num_points, num_cs, vertex_list, curve_list):
    """ This routine connects all points of each cross section 
        with a spline curve, plus one linear line to close the curve.
         Input:
           num_points: Integer, number of points of original cross section
           num_cs: Integer, number of additional cross sections per blade
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
         Output:
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
    """
    # Create a spline curve instead of having linear lines between vertices:
    self.cadlines.append('\n'+self.commentchar+' Define a Spline Curve:')
    cadline = 'create curve spline vertex '
    for i in range(1, num_points+1):
        cadline = cadline+' '+str(i)
    self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
    # Now make one more line, to connect last vertex with first vertex:
    self.cadlines.append(self.cubitcmdopen+'create curve vertex '+str(1)+' '+str(num_points)+self.cubitcmdclose)
    # Adding id to curve_list:
    curve_list.append(1); curve_list.append(2)
    # Thus increase vertex_list accordingly:
    vertex_list = self.list_id_increase(vertex_list, 2)

    # Do the same for the additional blade cross sections:
    for j in range(1, num_cs+1):
      self.cadlines.append('\n'+self.commentchar+' Define Spline curve for the extruded surface j='+str(j)+':')
      cadline = 'create curve spline vertex '
      for i in range(1, num_points+1):
          cadline = cadline+' '+str(num_points*j+i)
      self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
      # Now make one more line, to connect last vertex with first vertex:
      self.cadlines.append(self.cubitcmdopen+'create curve vertex '+str(num_points*j+1)+' '+str(num_points*j+num_points)+self.cubitcmdclose)
      # Adding id to curve_list:
      curve_list = self.list_id_increase(curve_list, 2)
      # Thus increase vertex_list accordingly:
      vertex_list = self.list_id_increase(vertex_list, 2)
    return vertex_list, curve_list


  def create_linear_lines_cross_sections_gmsh(self, blade, num_points, vertex_list, curve_list):
    """ Drawing lines to connect points created for GMSH
         Input:
           blade: Integer of the i-th blade that we draw lines for
           num_points: Integer, number of points of original cross section
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
         Output:
           curve_list: list of curves generated in cubit
    """
    if (len(curve_list) == 0): curveid = 0
    else: curveid = curve_list[-1]
    # First pointid for this blade:
    firstpointid = (blade-1) * num_points + 1
    # Curve id offset:
    curveidoffset = curveid
    # Now assemble gmsh entries to define lines, assume that 
    # the given coordinates in the textfile were immediate neightbours
    # from one line to the next
    self.cadlines.append('\n'+self.commentchar+' Lines:')
    for i in range(1, num_points+1):
        if (i < num_points):
            self.cadlines.append('Line ('+str(i+curveidoffset)+') = {'+str(i+curveidoffset)+', '+str(i+curveidoffset+1)+'};')
        else:
            # For the last point, close the loop:
            self.cadlines.append('Line ('+str(i+curveidoffset)+') = {'+str(i+curveidoffset)+', '+str(firstpointid)+'};')
        # Add line id to 'curve_list':
        curve_list.append(i+curveidoffset)
    return curve_list


  def create_surface_cross_sections_cubit(self, num_points, num_cs, surface_list):
    """ This routines assembles cubit commands for creating surfaces for
        all generated cross sections
         Input:
           num_points: Integer, number of points of original cross section
           num_cs: Integer, number of additional cross sections per blade
           surface_list: list of surface generated in cubit
         Output:
           surface_list: list of surface generated in cubit
    """
    self.cadlines.append('\n'+self.commentchar+' Creating the corresponding surface:')
    # For linear lines, do:
    cadline = 'create surface curve'
    for i in range(1, num_points+1):
        cadline = cadline+' '+str(i)
    self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
    # Adding surface ids:
    surface_list.append(1)
    # Same for the extruded surface and the surface between the blade and the circle sitting on the hub:
    for j in range(1, num_cs+1):
      self.cadlines.append('\n'+self.commentchar+' Creating the surface of the extruded surface j='+str(j)+':')
      cadline = 'create surface curve'
      for i in range(1, num_points+1):
          cadline = cadline+' '+str(num_points*j+i)
      self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
      # Adding surface ids:
      surface_list.append(j+1)
    return surface_list


  def create_surface_cross_sections_gmsh(self, blade, num_points, surface_list):
    """ This routines assembles surfaces in GMSH
         Input:
           blade: Integer of the i-th blade that we create the surface for
           num_points: Integer, number of points of original cross section
           surface_list: list of surface generated in cubit
         Output:
           surface_list: list of surface generated in cubit
    """
    # Curve id offset:
    curveidoffset = (blade-1) * num_points
    self.cadlines.append('\n'+self.commentchar+' Surface of blade number '+str(blade)+':')
    cadline = 'Line Loop ('+str(blade)+') = {'
    for i in range(1, num_points+1):
        if (i < num_points):
            cadline = cadline+str(i+curveidoffset)+', '
        else:
            cadline = cadline+str(i+curveidoffset)+'};'
    self.cadlines.append(cadline)
    self.cadlines.append('Plane Surface('+str(blade)+') = {'+str(blade)+'};')
    surface_list.append(blade)
    return surface_list


  def set_physical_ids(self, blade, num_points, curve_list, surface_list):
    """ This subroutine sets the physical ids of the blade cross section
         Input:
           blade: Integer of the i-th blade that we create the surface for
           num_points: Integer, number of points of original cross section
           curve_list: list of curves generated in cubit
           surface_list: list of surface generated in cubit
    """
    self.cadlines.append('\n'+self.commentchar+' Physical IDs:')
    curveidoffset = (blade-1) * num_points
    cadline = 'Physical Line('+str(blade)+') = {'
    for i in range(1, num_points+1):
        if (i < num_points):
            cadline = cadline+str(i+curveidoffset)+', '
        else:
            cadline = cadline+str(i+curveidoffset)+'};'
    self.cadlines.append(cadline)
    self.cadlines.append('Physical Surface('+str(blade)+') = {'+str(surface_list[blade-1])+'};')


  def apply_twist_to_cross_sections_cubit(self,rotation=None):
    """ This routine applies the rotation specified by the user to each
        cross section of the blade.
         Input:
          rotation: list of floats with the twist angle defining the twist of each
            cross section of a blade.
    """
    if (rotation is None):
        rotation = self.rotation
    # Applying correct pitch angle for all aerofoil cross sections:
    # Be careful, we'll first only rotate the extruded surfaces by the difference of twist compared to the initial surface (nearest to the hub),
    # then later on, once the blade volume is created, we'll apply the twist of the initial surface to the volume to get the twist angles
    # of the entire blade correct:
    if (self.spin_direction == 'clockwise'): sign = -1;
    elif (self.spin_direction == 'counterclockwise'): sign = 1;
    twist = [sign*rotation[0] - sign*i for i in rotation]
#    if (self.spin_direction == 'clockwise'): sign = 1;
#    elif (self.spin_direction == 'counterclockwise'): sign = -1;
    for i in range(1, len(rotation)):
        self.cadlines.append('\n'+self.commentchar+' Now rotating surface i='+str(i+1)+':')
        self.cadlines.append(self.cubitcmdopen+'rotate surface '+str(i+1)+' angle '+str(-sign*twist[i])+' about z'+self.cubitcmdclose)# include_merged


  ###############################################################
  # With the extruded surface, its vertices, and curves         #
  # create connecting curves and surface between both surfaces: #
  ###############################################################
  def create_blade_surface_cubit(self, curve_type, num_points, num_cs, vertex_list, curve_list, surface_list):
    """ This routine creates spline curves and based on those a surface to define
        the blade's surface.
         Input:
           curve_type: String, either 'spline' or 'linear'. Type of line connecting
             the intermediate surfaces of the blade.
           num_points: Integer, number of points of original cross section
           num_cs: Integer, number of additional cross sections per blade
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
           surface_list: list of surface generated in cubit
         Output:
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
           surface_list: list of surface generated in cubit
    """
    # Create curves to connect the points of different cross sections:
    if (curve_type == 'linear'):
        # As we are going to use cubit's lofted volume function for this, we'll return straight away:
        return (vertex_list, curve_list, surface_list)
    elif (not (curve_type == 'spline')):
        print "ERROR: Unrecognized type of curve given. Will exit..."
        raise SystemExit()
    else: # This is, if curve_type == 'spline'
        self.cadlines.append('\n'+self.commentchar+' Creating Curves to connect all end and intermediate surfaces:')
        for i in range(1, num_points+1):
            # Here we have to create the spline from the vertex on the last created surface (that is the one between blade and circle sitting on the hub), then the vertex on the surface of the original geometry of the blade, and then the other extruded surfaces...
            cadline = 'create curve spline vertex '+str(i)
            for j in range(1, num_cs+1):
                cadline = cadline + ' '+str(num_points*j+i)
            # Now we have reached the last surface, so add the assembled command to cadlines:
            self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
            # Stupid cubit creates new vertices on top of old vertices when creating these lines....:
            vertex_list = self.list_id_increase(vertex_list, 2)
            curve_list = self.list_id_increase(curve_list, 1)

        # Now create surfaces of the above generated curves:
        self.cadlines.append('\n'+self.commentchar+' Creating the surfaces for the connecting curves of front and back surfaces:')
        for i in range(1, num_points+1):
            if (i < num_points):
                # curve on first cross section, curve on last cross section, first connecting curve and second connecting curve
                cadline = 'create surface curve '+str(i)+' '+str(num_points*num_cs+i)+' '+str(num_points*(num_cs+1)+i)+' '+str(num_points*(num_cs+1)+i+1)
            else:
                cadline = 'create surface curve '+str(i)+' '+str(num_points*num_cs+i)+' '+str(num_points*(num_cs+1)+i)+' '+str(num_points*(num_cs+1)+1)
            # Adding ids to surface_list:
            surface_list = self.list_id_increase(surface_list, 1)
            # Incredibly stupid cubit..., for each surface we create here, it creates 3 new curves
            curve_list = self.list_id_increase(curve_list, 3)
            # And even better, for each curve, three new vertices are created.... Could it get any better? I'm sure there is a reason for that, but what a pain for automatisation...
            vertex_list = self.list_id_increase(vertex_list, 3)
            self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)

        # Now remove the last entry in vertex_list, as the first created surface only created two new vertices, instead of three..... No words...
        delv = 1
        for i in range(delv):
            del vertex_list[-1]
    return (vertex_list, curve_list, surface_list)


  ########################################################
  # Create connecting surface to hub of turbine:         #
  ########################################################
  def create_connection_to_hub_cubit(self, vertex_list, curve_list, surface_list, volume_list, num_points, num_cs, extrusion_length=None):
    """ This routine creates a surface to connect the blade with the hub
         Input:
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
           surface_list: list of surface generated in cubit
           volume_list: List of integers with the volume ids
           num_points: Integer, number of points of original cross section
           num_cs: Integer, number of additional cross sections per blade
           extrusion_length: list of z-coordinates of the cross sections.
         Output:
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
           surface_list: list of surface generated in cubit
           volume_list: List of integers with the volume ids
    """
    if (extrusion_length is None):
        extrusion_length = self.extrusion_length
    self.cadlines.append('\n'+self.commentchar+' Creating a circle below the blade, that connects the blade with the hub:')
    # First, create some variables for the circle that we'll create:
    blade_hub_connection_circle_radius = self.blade_hub_connection_circle_radius
    blade_hub_connection_distance = self.blade_hub_connection_distance
    circle_surface_id = surface_list[-1] + 1
    var_circle_surface_id_name = 'circle_surface_id'
    self.query_entity_id_cubit(var_circle_surface_id_name, 'surface', offset=1)

    # Updating id lists:
    vertex_list.append(vertex_list[-1] + 1)
    curve_list.append(curve_list[-1] + 1)
    surface_list.append(circle_surface_id)
    # Now draw the cylinder and move it into position:
    self.cadlines.append(self.cubitcmdopen+'create surface circle radius '+str(blade_hub_connection_circle_radius)+' zplane'+self.cubitcmdclose)
    if (self.extrusion_length[-1] > 0):
        self.cadlines.append(self.cubitcmdopen+'move surface "+str('+var_circle_surface_id_name+')+" x '+str(blade_hub_connection_circle_radius)+' z '+str(-blade_hub_connection_distance)+' include_merged'+self.cubitcmdclose)
    elif (self.extrusion_length[-1] < 0):
        self.cadlines.append(self.cubitcmdopen+'move surface "+str('+var_circle_surface_id_name+')+" x '+str(blade_hub_connection_circle_radius)+' z '+str(blade_hub_connection_distance)+' include_merged'+self.cubitcmdclose)
    # Now create lofted volume from cylinder and original blade cross section:
    self.cadlines.append(self.cubitcmdopen+'create volume loft surface 1 "+str('+var_circle_surface_id_name+')+" heal'+self.cubitcmdclose)
    if (self.blade_curve_type == 'spline'):
        hub_connection_volume_id = circle_surface_id+1
    elif (self.blade_curve_type == 'linear'):
        hub_connection_volume_id = volume_list[-1]+2 # plus 2, because 1 intermediate volume was created due to the create surface circle command above!
    volume_list.append(hub_connection_volume_id)

    # Finally, let's move the original blade volume and the connection volume to the hub
    # into the negative z direction, such that the bottom of the cylinder is at z=0:
    self.cadlines.append('\n'+self.commentchar+' Moving blade into +/-z direction:')
    cadline = 'move Volume'
    for volid in volume_list:
        cadline = cadline+' '+str(volid)
    if (self.extrusion_length[-1] > 0):
        cadline = cadline + ' z '+str(blade_hub_connection_distance)+' include_merged'
    elif (self.extrusion_length[-1] < 0):
        cadline = cadline + ' z '+str(-blade_hub_connection_distance)+' include_merged'
    self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)

    # Now apply the pitch angle by rotating the blade-hub- and blade-volume around the surface normal of 'circle_surface_id':
    self.cadlines.append('\n'+self.commentchar+' Applying the blades pitch angle:')
    self.cadlines.append(self.cubitcmdopen+'rotate volume '+str(volume_list[-2])+' '+str(volume_list[-1])+' angle '+str(-90.0+self.rotation[0])+' about normal of surface "+str('+var_circle_surface_id_name+')+" include_merged'+self.cubitcmdclose)
    # Since we rotated around the surface normal of the circle at the bottom of the connecting volume of blade and hub, we don't need to correct the blade's position.

    # Also, delete volume of the cylinder only, that cubit stupidly creates, but is not needed:
    if (volume_list[0]+2 == volume_list[1]):
        self.cadlines.append(self.cubitcmdopen+'delete volume '+str(volume_list[0]+1)+self.cubitcmdclose)
    # While we are at it, also delete the volumes for the intermediate cross sections which are not needed anymore:
    cadline = 'delete volume'
    for i in range(1, num_cs+2):
        cadline = cadline+' '+str(i)
    self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)

    # Due to creating a lofted volume, we have to correct our id lists again...:
    # First, vertices, two new surfaces, with each having as many vertices as the orginal cross section, thus:
    # Same for curves and surfaces:
    vertex_list = self.list_id_increase(vertex_list, num_points)
    curve_list = self.list_id_increase(curve_list, num_points)
    surface_list = self.list_id_increase(surface_list, num_points)
    # But also, two more new curves for each connecting surface (which then will form the volume...):
    vertex_list = self.list_id_increase(vertex_list, num_points)
    curve_list = self.list_id_increase(curve_list, 2*num_points)
    # 2 more surfaces were created, for the front/end:
    surface_list = self.list_id_increase(surface_list, 2)
    # Now we should be up to date with the id lists...
    # Return the assembled list of volume ids:
    return (vertex_list, curve_list, surface_list, volume_list)


  ###########################################
  # Create the volume of the blade:         #
  ###########################################
  def create_blade_volume_cubit(self, vertex_list, curve_list, surface_list, volume_list, num_points, num_cs, num_blades=None):
    """ This routine creates a volume of the blade
         Input:
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
           surface_list: list of surface generated in cubit
           volume_list: List of integers with the volume ids
           num_points: Integer, number of points of original cross section
           num_cs: Integer, number of additional cross sections per blade
           num_blades: Integer, number of blades the turbine should have
         Output:
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
           surface_list: list of surface generated in cubit
           volume_list: List of integers with the volume ids
    """
    if (num_blades is None):
        num_blades = self.num_blades
    # Create Volume of generated geometry:
    self.cadlines.append('\n'+self.commentchar+' Creating a volume of the generated surfaces/geometry:')
    if (self.blade_curve_type == 'spline'):
        cadline = 'create volume surface '+str(1)+' '+str(num_cs+1)
        for i in range(num_cs+2, surface_list[-1]+1):
            cadline = cadline+' '+str(i)
    elif (self.blade_curve_type == 'linear'):
        cadline = 'create volume loft surface '
        for i in range(1, num_cs+2):
            cadline = cadline+' '+str(i)
        # The volume loft operation creates a number of new vertices, curves and surfaces...
        # First, vertices, two new surfaces, with each having as many vertices as the orginal cross section, thus:
        # Same for curves and surfaces:
        vertex_list = self.list_id_increase(vertex_list, num_points)
        curve_list = self.list_id_increase(curve_list, num_points)
        surface_list = self.list_id_increase(surface_list, num_points)
        # But also, two more new curves for each connecting surface (which then will form the volume...):
        vertex_list = self.list_id_increase(vertex_list, num_points)
        curve_list = self.list_id_increase(curve_list, 2*num_points)
        # 2 more surfaces were created, for the front/end:
        surface_list = self.list_id_increase(surface_list, 2)

#    #elif (self.blade_curve_type == 'linear'):
#        # we are going to use cubit's create volume loft surface command, which only requires the surface ids of the intermediate cross sections of the blade!
    # Now add healing command to the Cubit command:
    cadline = cadline+' heal sheet'
    self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
    if (self.blade_curve_type == 'spline'):
        volume_list.append(surface_list[-1])
    elif (self.blade_curve_type == 'linear'):
        volume_list.append(num_cs+2)
        # Then merge the volume with previous volumes:
        cadline = 'merge volume '+str(volume_list[-1])+' with volume'
        for sid in range(1, num_cs+2):
            cadline = cadline+' '+str(sid)
        self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)

    # Now that the volume is done, let's merge all entities, which helps to prevent errors during the meshing operation:
    self.cadlines.append('\n'+self.commentchar+' Imprinting/Merging all volumes:')
    self.cadlines.append(self.cubitcmdopen+'imprint volume all'+self.cubitcmdclose)
    self.cadlines.append(self.cubitcmdopen+'merge volume all'+self.cubitcmdclose)

    # Return the list of volume ids:
    return vertex_list, curve_list, surface_list, volume_list


  def create_copies_of_blade_cubit(self, vertex_list, curve_list, surface_list, volume_list, num_points, num_cs, num_blades=None):
    """ This routine creates copies of the generated blade and rotates the
        copied volumes about the x-axis
         Input:
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
           surface_list: list of surfaces generated in cubit
           volume_list: List of integers with the volume ids
           num_points: Integer, number of points of original cross section
           num_cs: Integer, number of additional cross sections per blade
           num_blades: Integer, number of blades the turbine should have
         Output:
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
           surface_list: list of surfaces generated in cubit
           volume_list: List of integers with the volume ids
    """
    if (num_blades is None):
        num_blades = self.num_blades
    # Now create a copies of the volume and rotate it:
    if (num_blades > 1):
      blade_rotangle = 360.0/num_blades
      volid_string = ' '.join([str(volid) for volid in volume_list])
      self.cadlines.append('\n'+self.commentchar+' Creating '+str(num_blades-1)+' copies of the volume, and rotate these:')
      self.cadlines.append(self.cubitcmdopen+'volume '+volid_string+' copy rotate '+str(blade_rotangle)+' about x repeat '+str(num_blades-1)+self.cubitcmdclose)
      if (self.extrusion_length[-1] > 0):
          self.cadlines.append(self.cubitcmdopen+'move volume '+volid_string+' z '+str(round(self.hub_radius*0.85,5))+' include_merged'+self.cubitcmdclose)
      elif (self.extrusion_length[-1] < 0):
          self.cadlines.append(self.cubitcmdopen+'move volume '+volid_string+' z '+str(-round(self.hub_radius*0.85,5))+' include_merged'+self.cubitcmdclose)
      # Now that n copies of the blade were generated, move those away from 0,0,0
      if (self.extrusion_length[-1] > 0): movefactor = 1.0
      elif (self.extrusion_length[-1] < 0): movefactor = -1.0
      for i in range(2,num_blades+1):
        volid_copy_string = str(volume_list[-1]+1)+' '+str(volume_list[-1]+2)
        z_move = movefactor * (self.hub_radius*0.85) * cos((i-1)*blade_rotangle*(pi/180.0))
        y_move = movefactor * (self.hub_radius*0.85) * (-1.0)*sin((i-1)*blade_rotangle*(pi/180.0))
        # Rounding to 5 decimal points:
        z_move = round(z_move,5)
        y_move = round(y_move,5)
        if (abs(z_move) < 1e-15): z_move = 0.0
        if (abs(y_move) < 1e-15): y_move = 0.0
        self.cadlines.append(self.cubitcmdopen+'move volume '+volid_copy_string+' y '+str(y_move)+' z '+str(z_move)+' include_merged'+self.cubitcmdclose)
        # Append volume number to volume_list:
        volume_list.append(volume_list[-1]+1) # for the blade volume
        volume_list.append(volume_list[-1]+1) # for the volume connecting the blade with the hub

    # Again, we have to update the vertex, curve and surface id lists:
    factor = num_points*4
    # don't start questioning the numbers... TIMES 4 because we have vertices for:
    # 1: original cross section, 2: back cross section, 3: original cross section again for connecting blade with hub, 4: circle that intersects with the hub
    vertex_list = self.list_id_increase(vertex_list, factor*(num_blades-1))
    #####################################
    # Curves:                           #
    #####################################
    factor = num_points*6
    # don't start questioning the numbers... TIMES 6 because we have curves for:
    # 1: original cross section, 2: back cross section, 3: original cross section again for connecting blade with hub, 4: circle that intersects with the hub, 5: curves for connecting original cross section with back cross section, 6: curves for connecting original cross section wit circle that sits in the hub...
    curve_list = self.list_id_increase(curve_list, factor*(num_blades-1))
    #####################################
    # Surfaces:                         #
    #####################################
    factor = num_points*2+4
    # don't start questioning the numbers... TIMES 2 because we have as many surfaces as we have vertices on the original cross section, once for the blade, and once for the volume that connects the blade with the hub...
    # don't start questioning the numbers... PLUS 4 because we have surfaces for:
    # 1: original cross section, 2: back cross section, 3: original cross section again for connecting blade with hub, 4: circle that intersects with the hub.
    surface_list = self.list_id_increase(surface_list, factor*(num_blades-1))
    # Return the list of volume ids:
    return vertex_list, curve_list, surface_list, volume_list


  def create_hub_wind_turbine_cubit(self, vertex_list, curve_list, surface_list, volume_list, group_list, num_points, num_cs, num_blades=None):
    """ This routine creates a standard hub for a wind turbine.
         Input:
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
           surface_list: list of surfaces generated in cubit
           volume_list: List of integers with the volume ids
           group_list: List of lists, first element of the inner list is a string of the 
             turbine's name/number, the next element of the inner list, is a list itself with
             the first element being the group's name, and the second element being a list of 
             the corresponding volume ids
           num_points: Integer, number of points of original cross section
           num_cs: Integer, number of additional cross sections per blade
           num_blades: Integer, number of blades the turbine should have
         Output:
           vertex_list: list of vertices generated in cubit
           curve_list: list of curves generated in cubit
           surface_list: list of surfaces generated in cubit
           volume_list: List of integers with the volume ids
           group_list: List of lists, first element of the inner list is a string of the 
             turbine's name/number, the next element of the inner list, is a list itself with
             the first element being the group's name, and the second element being a list of 
             the corresponding volume ids
    """
    if (num_blades is None):
        num_blades = self.num_blades
    hub_radius = self.hub_radius
    ellipse_tip_distance = self.ellipse_tip_distance
    hub_length = self.hub_length
    hub_height = self.hub_height
    base_length = self.base_length
    base_direction = self.base_direction
    base_radius = self.base_radius
    blade_base_distance = self.blade_base_distance
    # First, create the ellipsoid volume, that rotates with the blades:
    self.cadlines.append('\n'+self.commentchar+' Creating the hub of the turbine:')
    self.cadlines.append(self.commentchar+'First defining some variables:')
    var_hub_ellipsoid_surface_id_name = 'hub_ellipsoid_surface_id'
    var_vertex_circle_centre_id_name = 'vertex_circle_centre_id'
    var_vertex_on_circle_ids_name = 'vertex_on_circle_ids'
    var_vertex_tip_ellipsoid_id_name = 'vertex_tip_ellipsoid_id'
    self.query_entity_id_cubit(var_hub_ellipsoid_surface_id_name, 'surface', offset=1)
    self.query_entity_id_cubit(var_vertex_circle_centre_id_name, 'vertex', offset=1)
    self.cadlines.append(var_vertex_on_circle_ids_name+' = ['+var_vertex_circle_centre_id_name+' + i for i in range(1,5)]')
    self.cadlines.append(var_vertex_tip_ellipsoid_id_name+' = '+var_vertex_on_circle_ids_name+'[-1] + 1')
    # obsolete:
    hub_ellipsoid_surface_id = surface_list[-1] + 1
    vertex_circle_centre_id = vertex_list[-1] + 1
    vertex_on_circle_ids = [i + vertex_list[-1]+1 for i in range(1, 5)]
    vertex_tip_ellipsoid_id = vertex_on_circle_ids[-1] + 1

    # Create 6 vertices, 4 on a circle, 1 for the centre of the circle, and 1 for the tip of the ellipsoid:
    self.cadlines.append(self.cubitcmdopen+'create vertex 0 0 0'+self.cubitcmdclose) # centre of circle
    self.cadlines.append(self.cubitcmdopen+'create vertex 0 '+str(hub_radius)+' 0'+self.cubitcmdclose) # on the circle
    self.cadlines.append(self.cubitcmdopen+'create vertex 0 0 '+str(hub_radius)+''+self.cubitcmdclose) # on the circle
    self.cadlines.append(self.cubitcmdopen+'create vertex 0 '+str(-hub_radius)+' 0'+self.cubitcmdclose) # on the circle
    self.cadlines.append(self.cubitcmdopen+'create vertex 0 0 '+str(-hub_radius)+''+self.cubitcmdclose) # on the circle
    self.cadlines.append(self.cubitcmdopen+'create vertex '+str(-ellipse_tip_distance)+' 0 0'+self.cubitcmdclose) # tip of ellipsoid
    # Now create arcs for the circle:
    # Query curve IDs:
    var_last_curve_id_name = 'last_curve_id'
    self.query_entity_id_cubit(var_last_curve_id_name, 'curve')
    var_circle_curve_ids_name = 'circle_curve_ids'
    self.cadlines.append(var_circle_curve_ids_name+' = ['+var_last_curve_id_name+' + i for i in range(1,5)]')

#    for i in range(len(vertex_on_circle_ids)):
    for i in range(4):
#        if (i == len(vertex_on_circle_ids)-1):
        if (i == 3):
            self.cadlines.append(self.cubitcmdopen+'create curve arc center vertex "+str('+var_vertex_circle_centre_id_name+')+" "+str('+var_vertex_on_circle_ids_name+'['+str(i)+'])+" "+str('+var_vertex_on_circle_ids_name+'[0])+" '+self.cubitcmdclose)
        else:
            self.cadlines.append(self.cubitcmdopen+'create curve arc center vertex "+str('+var_vertex_circle_centre_id_name+')+" "+str('+var_vertex_on_circle_ids_name+'['+str(i)+'])+" "+str('+var_vertex_on_circle_ids_name+'['+str(i+1)+'])+" '+self.cubitcmdclose)
    # OBSOLETE:
    # Updating the id lists again:
    # we created 6 new vertices, and cubit created 4 more due to creating curves:
    vertex_list = self.list_id_increase(vertex_list, 10)
    # we created 4 new curves:
    curve_list = self.list_id_increase(curve_list, 4)
    # Store the circle-curve ids in a seperate list:
    circle_curve_ids = curve_list[-4:]


    # Now create elliptical curves for the tip of the hub:
    var_last_curve_id_name = 'last_curve_id'
    self.query_entity_id_cubit(var_last_curve_id_name, 'curve')
    var_elliptical_curve_ids_name = 'elliptical_curve_ids'
    self.cadlines.append(var_elliptical_curve_ids_name+' = ['+var_last_curve_id_name+' + i for i in range(1,5)]')
    for i in range(4): # as we want to create 4 lines
        self.cadlines.append(self.cubitcmdopen+'create curve vertex "+str('+var_vertex_on_circle_ids_name+'['+str(i)+'])+" vertex "+str('+var_vertex_tip_ellipsoid_id_name+')+" vertex "+str('+var_vertex_circle_centre_id_name+')+" ellipse'+self.cubitcmdclose)
    # OBSOLETE:
    # Again, updating id lists:
    vertex_list = self.list_id_increase(vertex_list, 7) # 7 new vertices were created by cubit
    curve_list = self.list_id_increase(curve_list, 4)


    rotating_hub_vol_id = volume_list[-1] + 1
    # Later on, we need to translate several volumes, so we need to assemble a list for cubit that holds all the relevant volume IDs:
    self.cadlines.append(self.commentchar+' We\'ll need to translate a bunch of volumes later on. Here is where we start assembling a list of relevant volume IDs:')
    var_translate_volume_list_name = 'translate_volume_list'
    var_rotating_hub_vol_id_name = 'rotating_hub_vol_id'
    self.query_entity_id_cubit(var_rotating_hub_vol_id_name, 'volume', offset=1)
    self.cadlines.append(var_translate_volume_list_name+' = []')
    self.cadlines.append(var_translate_volume_list_name+'.append('+var_rotating_hub_vol_id_name+')')

    # Now, create surfaces that will form the ellipsoid with the last 4 curves we created:
    for i in range(4):
        if (i == 3):
            self.cadlines.append(self.cubitcmdopen+'create surface curve "+str('+var_circle_curve_ids_name+'['+str(i)+'])+" "+str('+var_elliptical_curve_ids_name+'['+str(i)+'])+" "+str('+var_elliptical_curve_ids_name+'[0])+" '+self.cubitcmdclose) # first id is the arc, the other two, elliptical 
        else:
            self.cadlines.append(self.cubitcmdopen+'create surface curve "+str('+var_circle_curve_ids_name+'['+str(i)+'])+" "+str('+var_elliptical_curve_ids_name+'['+str(i)+'])+" "+str('+var_elliptical_curve_ids_name+'['+str(i+1)+'])+" '+self.cubitcmdclose)
    # OBSOLETE:
    # Again, updating id lists:
    vertex_list = self.list_id_increase(vertex_list, 7) # again, 7 new vertices were created by cubit
    curve_list = self.list_id_increase(curve_list, 4)
    surface_list = self.list_id_increase(surface_list, 4)
    volume_list = self.list_id_increase(volume_list, 4) # for each surface, cubit also created a new volume...

    # And create one more surface, to close the circle:
    self.cadlines.append(self.cubitcmdopen+'create surface curve "+str('+var_circle_curve_ids_name+'[0])+" "+str('+var_circle_curve_ids_name+'[1])+" "+str('+var_circle_curve_ids_name+'[2])+" "+str('+var_circle_curve_ids_name+'[3])+" '+self.cubitcmdclose)
    # OBSOLETE:
    # Update id lists:
    vertex_list = self.list_id_increase(vertex_list, 4)
    curve_list = self.list_id_increase(curve_list, 4)
    surface_list = self.list_id_increase(surface_list, 1)
    volume_list = self.list_id_increase(volume_list, 1)


    # Create volume from the surfaces:
    # Thus we first need the last five surface IDs:
    var_closing_tiphub_circle_surface_id_name = 'closing_tiphub_circle_surface_id'
    self.query_entity_id_cubit(var_closing_tiphub_circle_surface_id_name, 'surface')
    var_tiphub_surface_ids_name = 'tiphub_surface_ids'
    self.cadlines.append(var_tiphub_surface_ids_name+' = ['+var_closing_tiphub_circle_surface_id_name+' -4 +i for i in range(5)]')
    cadline = 'create volume surface'
    for i in range(1, 6):
        cadline = cadline + '  "+str('+var_tiphub_surface_ids_name+'['+str(-i)+'])+" '
    cadline = cadline + ' heal'
    self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)


    # The ellipsoid is done, now extrude the circle of the ellipsoid into x-direction. For performance sake, the distance of the extrusion should be limited to the diameter of the circle of the blade-hub-connection, as the blades, as well as the ellipsoid and the extruded cylinder we are creating below are going to be the moving parts, the rest can be assumed as static, thus we want to split the mesh into moving and non-moving parts.
    hub_circle_surface_id = surface_list[-1]
    length_blade_on_hub_factor = 1.5 # 50% more than the bottom circle diameter coming from the blade
    length_blade_on_hub = self.blade_hub_connection_circle_radius*2*length_blade_on_hub_factor
#    self.cadlines.append(self.cubitcmdopen+'sweep surface '+str(hub_circle_surface_id)+' direction x 1 distance '+str(length_blade_on_hub)+self.cubitcmdclose)
    self.cadlines.append(self.cubitcmdopen+'sweep surface "+str('+var_closing_tiphub_circle_surface_id_name+')+" direction x 1 distance '+str(length_blade_on_hub)+self.cubitcmdclose)
    # OBSOLETE:
    # Update id lists:
    vertex_list = self.list_id_increase(vertex_list, 5)
    curve_list = self.list_id_increase(curve_list, 10)
    surface_list = self.list_id_increase(surface_list, 4) # 4 new surfaces created!

    # Create a copy of the surface
    self.cadlines.append(self.cubitcmdopen+'Surface "+str('+var_closing_tiphub_circle_surface_id_name+')+" copy'+self.cubitcmdclose)
    # OBSOLETE:
    vertex_list = self.list_id_increase(vertex_list, 4)
    curve_list = self.list_id_increase(curve_list, 4)
    surface_list = self.list_id_increase(surface_list, 1)
    volume_list = self.list_id_increase(volume_list, 1)
    # Save this volume id:
    base_vol_id = volume_list[-1]

    # Extrude that last surface (the copied surface) for the remaining length of the hub:
    # Thus, find out the surface ID first:
    var_last_surface_id_name = 'last_surface_id'
    self.query_entity_id_cubit(var_last_surface_id_name, 'surface')
    hub_end_x_coord = hub_length
    hub_length = hub_length - length_blade_on_hub
#    self.cadlines.append(self.cubitcmdopen+'sweep surface '+str(surface_list[-1])+' direction x 1 distance '+str(hub_length)+self.cubitcmdclose)
    self.cadlines.append(self.cubitcmdopen+'sweep surface "+str('+var_last_surface_id_name+')+" direction x 1 distance '+str(hub_length)+self.cubitcmdclose)
    # Get the volume ID and append it to the list of translated volumes:
    self.cadlines.append(self.commentchar+' Append the list of translated volume IDs:')
    self.query_entity_id_cubit('last_volume_id', 'volume')
    self.cadlines.append(var_translate_volume_list_name+'.append(last_volume_id)')
    # OBSOLETE:
    vertex_list = self.list_id_increase(vertex_list, 4)
    curve_list = self.list_id_increase(curve_list, 8)
    surface_list = self.list_id_increase(surface_list, 5)

    # And create another volume at the back of the hub, that has no sharp edges:
    # Create points on the surface of the back of the hub:
#===================================================================================================================
# OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE :
#    self.cadlines.append(self.cubitcmdopen+'create vertex '+str(hub_end_x_coord)+' 0 0'+self.cubitcmdclose)
#    self.cadlines.append(self.cubitcmdopen+'create vertex '+str(hub_end_x_coord)+' '+str(-hub_radius)+' 0'+self.cubitcmdclose)
#    self.cadlines.append(self.cubitcmdopen+'create vertex '+str(hub_end_x_coord)+' 0 '+str(hub_radius)+self.cubitcmdclose)
#    self.cadlines.append(self.cubitcmdopen+'create vertex '+str(hub_end_x_coord)+' '+str(hub_radius)+' 0'+self.cubitcmdclose)
#    self.cadlines.append(self.cubitcmdopen+'create vertex '+str(hub_end_x_coord)+' 0 '+str(-hub_radius)+self.cubitcmdclose)
    vertex_on_back_hub_circle_centre_id = vertex_list[-1] + 1
    vertex_on_back_hub_ids = [i + vertex_list[-1]+1 for i in range(1, 5)]
#    for i in range(len(vertex_on_back_hub_ids)):
#        if (i == len(vertex_on_back_hub_ids)-1):
#            self.cadlines.append(self.cubitcmdopen+'create curve arc center vertex '+str(vertex_on_back_hub_circle_centre_id)+' '+str(vertex_on_back_hub_ids[i])+' '+str(vertex_on_back_hub_ids[0])+self.cubitcmdclose)
#        else:
#            self.cadlines.append(self.cubitcmdopen+'create curve arc center vertex '+str(vertex_on_back_hub_circle_centre_id)+' '+str(vertex_on_back_hub_ids[i])+' '+str(vertex_on_back_hub_ids[i+1])+self.cubitcmdclose)
    # Updating id lists again and again...
    vertex_list = self.list_id_increase(vertex_list, 9)
    curve_list = self.list_id_increase(curve_list, 4)

    # Store the circle-curve ids in a seperate list:
    curves_on_back_hub_ids = curve_list[-4:]
#    self.cadlines.append(self.cubitcmdopen+'create surface curve '+str(curves_on_back_hub_ids[0])+' '+str(curves_on_back_hub_ids[1])+' '+str(curves_on_back_hub_ids[2])+' '+str(curves_on_back_hub_ids[3])+self.cubitcmdclose)
    # Update id lists:
    surface_list = self.list_id_increase(surface_list, 1)
    volume_list = self.list_id_increase(volume_list, 1)
#===================================================================================================================
    
    # Simply copy the extruded surface again:
    self.cadlines.append(self.cubitcmdopen+'Surface "+str('+var_last_surface_id_name+')+" copy'+self.cubitcmdclose)
    # Now create lists to store relevant vertices/curves IDs of that surface:
    var_hub_end_big_circle_surface_id_name = 'hub_end_big_circle_surface_id'
    var_hub_end_big_circle_curve_ids_name = 'hub_end_big_circle_curve_ids'
    var_hub_end_big_circle_vertex_ids_name = 'hub_end_big_circle_vertex_ids'
    self.query_entity_id_cubit(var_hub_end_big_circle_surface_id_name, 'surface')
    self.query_entity_id_cubit('last_curve_id', 'curve')
    self.cadlines.append(var_hub_end_big_circle_curve_ids_name+' = [last_curve_id -3 +i for i in range(4)]')
    self.query_entity_id_cubit('last_vertex_id', 'vertex')
    self.cadlines.append(var_hub_end_big_circle_vertex_ids_name+' = [last_vertex_id -3 +i for i in range(4)]')
    # Now that we have the IDs of the bigger surface, and its curves and vertices, we can create a copy of the 
    # surface and scale it down by 10%:
    self.cadlines.append(self.cubitcmdopen+'Surface "+str('+var_hub_end_big_circle_surface_id_name+')+" copy scale 0.9 about '+str(hub_end_x_coord)+' 0 0'+self.cubitcmdclose)
    # And query their IDs:
    var_hub_end_small_circle_surface_id_name = 'hub_end_small_circle_surface_id'
    var_hub_end_small_circle_curve_ids_name = 'hub_end_small_circle_curve_ids'
    var_hub_end_small_circle_vertex_ids_name = 'hub_end_small_circle_vertex_ids'
    self.query_entity_id_cubit(var_hub_end_small_circle_surface_id_name, 'surface')
    self.query_entity_id_cubit('last_curve_id', 'curve')
    self.cadlines.append(var_hub_end_small_circle_curve_ids_name+' = [last_curve_id -3 +i for i in range(4)]')
    self.query_entity_id_cubit('last_vertex_id', 'vertex')
    self.cadlines.append(var_hub_end_small_circle_vertex_ids_name+' = [last_vertex_id -3 +i for i in range(4)]')

#===================================================================================================================
# OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE :
    back_hub_surface_id = surface_list[-1]
    # Create points that will form the slightly smaller surface behind the hub:
#    self.cadlines.append(self.cubitcmdopen+'create vertex '+str(hub_end_x_coord)+' 0 0'+self.cubitcmdclose)
#    self.cadlines.append(self.cubitcmdopen+'create vertex '+str(hub_end_x_coord)+' '+str(-hub_radius*0.9)+' 0'+self.cubitcmdclose)
#    self.cadlines.append(self.cubitcmdopen+'create vertex '+str(hub_end_x_coord)+' 0 '+str(hub_radius*0.9)+self.cubitcmdclose)
#    self.cadlines.append(self.cubitcmdopen+'create vertex '+str(hub_end_x_coord)+' '+str(hub_radius*0.9)+' 0'+self.cubitcmdclose)
#    self.cadlines.append(self.cubitcmdopen+'create vertex '+str(hub_end_x_coord)+' 0 '+str(-hub_radius*0.9)+self.cubitcmdclose)
    # Now create arcs for the circle:
    vertex_small_back_hub_circle_centre_id = vertex_list[-1] + 1
    vertex_small_back_hub_ids = [i + vertex_list[-1]+1 for i in range(1, 5)]
#    for i in range(len(vertex_small_back_hub_ids)):
#        if (i == len(vertex_small_back_hub_ids)-1):
#            self.cadlines.append(self.cubitcmdopen+'create curve arc center vertex '+str(vertex_small_back_hub_circle_centre_id)+' '+str(vertex_small_back_hub_ids[i])+' '+str(vertex_small_back_hub_ids[0])+self.cubitcmdclose)
#        else:
#            self.cadlines.append(self.cubitcmdopen+'create curve arc center vertex '+str(vertex_small_back_hub_circle_centre_id)+' '+str(vertex_small_back_hub_ids[i])+' '+str(vertex_small_back_hub_ids[i+1])+self.cubitcmdclose)
    # Updating id lists again and again...
    vertex_list = self.list_id_increase(vertex_list, 9)
    curve_list = self.list_id_increase(curve_list, 4)
    # Store the circle-curve ids in a seperate list:
    curves_small_back_hub_ids = curve_list[-4:]

    # And create one more surface, to close the circle:
#    self.cadlines.append(self.cubitcmdopen+'create surface curve '+str(curves_small_back_hub_ids[0])+' '+str(curves_small_back_hub_ids[1])+' '+str(curves_small_back_hub_ids[2])+' '+str(curves_small_back_hub_ids[3])+self.cubitcmdclose)
    # Update id lists:
    surface_list = self.list_id_increase(surface_list, 1)
    volume_list = self.list_id_increase(volume_list, 1)
#===================================================================================================================

    # Move the circle away from the hub:
    self.cadlines.append(self.cubitcmdopen+'move Surface "+str('+var_hub_end_small_circle_surface_id_name+')+" x '+str(0.1*hub_radius)+' include_merged'+self.cubitcmdclose)
    # Now create arcs from circle on back of hub to the slightly smaller circle behind the hub:
    normal_string = ['0 0 1', '0 1 0', '0 -1 0', '0 0 -1']
    for i in range(4):
        self.cadlines.append(self.cubitcmdopen+'create curve arc center vertex "+str('+var_hub_end_big_circle_vertex_ids_name+'['+str(i)+'])+" "+str('+var_hub_end_small_circle_vertex_ids_name+'['+str(i)+'])+" radius '+str(0.1*hub_radius)+' normal '+normal_string[i]+' '+self.cubitcmdclose)
    # Get the curve IDs:
    var_hub_end_big2small_connect_curve_ids_name = 'hub_end_big2small_connect_curve_ids'
    self.query_entity_id_cubit('last_curve_id', 'curve')
    self.cadlines.append(var_hub_end_big2small_connect_curve_ids_name+' = [last_curve_id -3 +i for i in range(4)]')

    # OBSOLETE:
    # Update id lists:
    vertex_list = self.list_id_increase(vertex_list, 8)
    curve_list = self.list_id_increase(curve_list, 4)

    curve_connect_back_hub_ids = curve_list[-4:]
    self.cadlines.append(self.commentchar+' Creating connecting surfaces between big and small circle at the end of the hub:')
    # Create surfaces:
    iind = [0, 0, 2, 1]; jind = [1, 2, 3, 3]
    for i in range(4):
        self.cadlines.append(self.cubitcmdopen+'create surface curve "+str('+var_hub_end_big_circle_curve_ids_name+'['+str(i)+'])+" "+str('+var_hub_end_small_circle_curve_ids_name+'['+str(i)+'])+" "+str('+var_hub_end_big2small_connect_curve_ids_name+'['+str(iind[i])+'])+" "+str('+var_hub_end_big2small_connect_curve_ids_name+'['+str(jind[i])+'])+" '+self.cubitcmdclose)

    surface_list = self.list_id_increase(surface_list, 4)
    # Create volume:
    self.cadlines.append(self.commentchar+' Make a volume of the small area at the end of the hub:')
    # First, get the surface IDs of the surfaces we just created:
    var_hub_end_big2small_connect_surface_ids_name = 'hub_end_big2small_connect_surface_ids'
    self.query_entity_id_cubit('last_surface_id', 'surface')
    self.cadlines.append(var_hub_end_big2small_connect_surface_ids_name+' = [last_surface_id -3 +i for i in range(4)]')
    cadline = 'create volume surface "+str('+var_hub_end_big_circle_surface_id_name+')+" "+str('+var_hub_end_small_circle_surface_id_name+')+" '
    for i in range(4):
        cadline = cadline + '"+str('+var_hub_end_big2small_connect_surface_ids_name+'['+str(i)+'])+" '
    cadline = cadline + ' heal'
    self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
    # Get the volume ID and append it to the list of translated volumes:
    self.cadlines.append(self.commentchar+' Append the list of translated volume IDs:')
    self.query_entity_id_cubit('last_volume_id', 'volume')
    self.cadlines.append(var_translate_volume_list_name+'.append(last_volume_id)')

    # OBSOLETE:
    vertex_list = self.list_id_increase(vertex_list, 11)
    curve_list = self.list_id_increase(curve_list, 12)
    volume_list = self.list_id_increase(volume_list, 4)

    # Now create the base of the wind turbine:
    self.cadlines.append(self.commentchar+' Now draw the base structure of the turbine:')
    self.cadlines.append(self.cubitcmdopen+'create surface circle radius '+str(base_radius)+' yplane'+self.cubitcmdclose)
    # Get the circle's surface ID:
    var_base_circle_surface_id_name = 'base_circle_surface_id'
    self.query_entity_id_cubit(var_base_circle_surface_id_name, 'surface')
    # OBSOLETE:
    # Updating id lists:
    vertex_list = self.list_id_increase(vertex_list, 1)
    curve_list = self.list_id_increase(curve_list, 1)
    surface_list = self.list_id_increase(surface_list, 1)
    volume_list = self.list_id_increase(volume_list, 1)

    # Move the circle in place:
    if (blade_base_distance is None):
      base_x_coordinate = length_blade_on_hub + (hub_length)*0.5
    else:
      base_x_coordinate = blade_base_distance - 0.5*length_blade_on_hub
    self.cadlines.append(self.cubitcmdopen+'move Surface "+str('+var_base_circle_surface_id_name+')+" x '+str(base_x_coordinate)+' include_merged '+self.cubitcmdclose)
    # Extrude circle into +/-y direction:
    if (base_direction == '-y' or base_direction == -1 or base_direction == '-1'):
      y_factor = -1
    elif (base_direction == 'y' or base_direction == '+y' or base_direction == 1 or base_direction == '1' or base_direction == '+1'):
      y_factor = 1
    self.cadlines.append(self.cubitcmdopen+'sweep surface "+str('+var_base_circle_surface_id_name+')+" vector 0 '+str(y_factor)+' 0 distance '+str(base_length)+self.cubitcmdclose)
    # OBSOLETE:
    vertex_list = self.list_id_increase(vertex_list, 1)
    curve_list = self.list_id_increase(curve_list, 2)
    surface_list = self.list_id_increase(surface_list, 2)

    # Due to meshing problems, move the volumes of the complete hub into the negative x direction a bit, 
    # such that the blades intersect with the center of the surface/volume where they are connected to the hub:
    rotating_hub_vol_id_index = volume_list.index(rotating_hub_vol_id)
#    move_volume_string_list = ' '.join([str(i) for i in volume_list[rotating_hub_vol_id_index:-1]])
    move_volume_distance = -0.5*(length_blade_on_hub-self.blade_hub_connection_circle_radius*2)
    self.cadlines.append(self.commentchar+' Move the hub-volumes a bit in the negative x direction due to meshing problems otherwise:')
    self.cadlines.append('move_volume_string_list = \' \'.join([str(i) for i in '+var_translate_volume_list_name+'])')
    self.cadlines.append('cadline = "move Volume "+move_volume_string_list+" x '+str(move_volume_distance)+' include_merged"')
    self.cadlines.append(self.cubitcmdopen+'"+cadline+"'+self.cubitcmdclose)
#    self.cadlines.append(self.cubitcmdopen+'move Volume '+move_volume_string_list+' x '+str(move_volume_distance)+' include_merged '+self.cubitcmdclose)

    # Get rid of free vertex:
    self.cadlines.append('\n'+self.commentchar+' Deleting a free vertex:')
    self.cadlines.append(self.cubitcmdopen+'delete vertex "+str('+var_vertex_circle_centre_id_name+')+" '+self.cubitcmdclose)

    # Before starting to unite volumes, we should make sure all bodies have topology, 
    # cubit's "imprint" command takes care of that:
    self.cadlines.append('\n'+self.commentchar+' Imprinting all bodies:')
    self.cadlines.append(self.cubitcmdopen+'imprint all'+self.cubitcmdclose)

    # Now that all geometries have been created, let's unite the intersecting volumes on the tip of the hub:
    self.cadlines.append('\n'+self.commentchar+' Uniting volumes of the tip of the hub:')
    cadline = 'unite volume'
    # Let's assemble a list of volume ids, of the volumes we want to unite, those are basically
    # all the volumes, except for the blade-volumes!!!
    unite_volids = []; bladevol_indices = []; bladevol_ids = []
    for n in range(0,num_blades):
        bladevol_indices.append(n*2)
        bladevol_ids.append(volume_list[n*2])
    for i in range(len(volume_list[:rotating_hub_vol_id_index+1])):
        unite_volids.append(volume_list[rotating_hub_vol_id_index-i])
    # Now remove the elements of unite_volids which belong to the blade volumes, as we found out above:
    z = 0 # each time we remove an element from 'unite_volids' we have to increase z, such that we remove the correct second element from 'unite_volids'
    for rmid in bladevol_indices:
        del unite_volids[(len(unite_volids)-1)-rmid+z]
        z = z+1
    # Now we can finally go on and assemble the string with the relevant cubit command:
    for uvolid in unite_volids:
        cadline = cadline + ' '+str(uvolid)
    self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
    # Creating a group for the moving parts:
    group_blade_volids = []
    for bid in bladevol_ids:
        group_blade_volids.append(bid)
    # And the volume id of the tip of the hub:
    group_blade_volids.append(volume_list[rotating_hub_vol_id_index])
    volume_string = 'Volume '
    for vid in group_blade_volids:
        volume_string = volume_string+' '+str(vid)
    self.cadlines.append(self.cubitcmdopen+'group \'blades\' add '+str(volume_string)+self.cubitcmdclose)

    # Unite the last few created volumes that belong to the hub. This is fairly easy, 
    # as the hub ids will always have the same relation: the last two volumes plus the one
    # we saved before under 'base_vol_id'
    cadline = 'unite volume '+str(base_vol_id)+' '+str(volume_list[-2])+' '+str(volume_list[-1])
    self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
    # Now create a group for the base:
    self.cadlines.append(self.cubitcmdopen+'group \'base\' add Volume '+str(base_vol_id)+self.cubitcmdclose)

    # Updating group list:
    group_list.append(['turbine-1', ['blades', group_blade_volids], ['base', [base_vol_id]]])

    return vertex_list, curve_list, surface_list, volume_list, group_list


  def create_block_ids_cubit(self, group_list, blockid_offset):
    """ This routines assembles the cubit commands to create block and 
        sideset ids for the generated geometry.
         Input:
           group_list: List of lists, first element of the inner list is a string of the 
             turbine's name/number, the next element of the inner list, is a list itself with
             the first element being the group's name, and the second element being a list of 
             the corresponding volume ids
           blockid_offset: Integer, the offset for new block ids
         Output:
           group_list: List of lists, first element of the inner list is a string of the 
             turbine's name/number, the next element of the inner list, is a list itself with
             the first element being the group's name, and the second element being a list of 
             the corresponding volume ids
    """
    #self.cadlines.append('\n'+self.commentchar+' Block and Sideset IDs:')
    # For each list in 'group_list' we want to create one block id, such that we can then export different meshes based on the block ids:
    block_id = blockid_offset
    for turbine_group in group_list[-1][1:]: # only consider the last entry (turbine) in group_list, and without the corresponding name/number!
        block_id += 1
        volid = turbine_group[-1] # volume id of that group
        volid_string = ' '.join(str(vid) for vid in volid)
        self.cadlines.append(self.cubitcmdopen+'block '+str(block_id)+' volume '+volid_string+self.cubitcmdclose)


  ################################################################
  # Setting up meshing options, and meshing specific volume ids: #
  ################################################################
  def mesh_volume_ids_cubit(self, volume_ids, meshingscheme=None, sizingfunction=None, meshingoption=None):
    """ This routine assembles the cubit commands to create one mesh of all the volume ids in a given list.
         Input:
           volume_ids: List of integers of volume ids that we want to mesh.
           meshingscheme: String, should be 'tetmesh' by default for Fluidity
           sizingfunction: Float, specifies the sizing of the elements
    """
    if (meshingscheme is None):
        meshingscheme = self.meshingscheme
    if (sizingfunction is None):
        sizingfunction = self.sizingfunction
    if (meshingoption is None):
        meshingoption = self.meshingoption
    self.cadlines.append('\n'+self.commentchar+' Setting up meshing options, and meshing the volume(s)')
    # Create string of 'volume firstvolid secondvolid ...'
    volume_string = 'volume'
    for volid in volume_ids:
        volume_string = volume_string+' '+str(volid)
    # First, set meshing scheme, should always be 'tetmesh' for Fluidity:
    cadline = volume_string+' scheme '+str(meshingscheme)
    self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
    # Set sizing function:
    if (meshingoption == 'auto'):
        self.cadlines.append(self.cubitcmdopen+volume_string+' size auto factor '+str(sizingfunction)+self.cubitcmdclose)
    elif (meshingoption == 'adaptive'):
        self.cadlines.append(self.cubitcmdopen+volume_string+' sizing function type skeleton scale '+str(sizingfunction)+' time_accuracy_level 3'+self.cubitcmdclose) # Just as a default, could be changed manually
    elif (meshingoption == 'size'):
        self.cadlines.append(self.cubitcmdopen+volume_string+' size '+str(sizingfunction)+self.cubitcmdclose)
    # Meshing the volumes of volume_ids:
    self.cadlines.append(self.cubitcmdopen+'mesh '+volume_string+self.cubitcmdclose)


  ########################################################
  # Setting up meshing options, and meshing the turbine: #
  ########################################################
  def mesh_turbine_cubit(self, group_list, meshingscheme=None, sizingfunction=None, meshingoption=None):
    """ This routine assembles the cubit commands to mesh the rest of the turbine geometry (everything except for the blades).
         Input:
           group_list: List of lists, first element of the inner list is a string of the 
             turbine's name/number, the next element of the inner list, is a list itself with
             the first element being the group's name, and the second element being a list of 
             the corresponding volume ids
           meshingscheme: String, should be 'tetmesh' by default for Fluidity
           sizingfunction: Float, specifies the sizing of the elements
    """
    if (meshingscheme is None):
        meshingscheme = self.meshingscheme
    if (sizingfunction is None):
        sizingfunction = self.sizingfunction
    if (meshingoption is None):
        meshingoption = self.meshingoption
    # We want to create as many meshes, as we have defined groups, which should only be 2:
    for turbine in group_list:
        for turbine_group in turbine[1:]:
            volid = turbine_group[-1] # volume id of that group
            self.cadlines.append('\n'+self.commentchar+' Setting up meshing options, and meshing the volume')
            # Create string of 'volume firstvolid secondvolid ...'
            # As we set it up in this way, the volume we want to mesh here, is the last volume id in this turbine_group:
            volume_string = 'volume '+str(volid[-1])
            # First, set meshing scheme, should always be 'tetmesh' for Fluidity:
            cadline = volume_string+' scheme '+str(meshingscheme)
            self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)
            # Set sizing function:
            if (meshingoption == 'auto'):
                self.cadlines.append(self.cubitcmdopen+volume_string+' size auto factor '+str(sizingfunction)+self.cubitcmdclose)
            elif (meshingoption == 'adaptive'):
                self.cadlines.append(self.cubitcmdopen+volume_string+' sizing function type skeleton scale '+str(sizingfunction)+' time_accuracy_level 3'+self.cubitcmdclose) # Just as a default, could be changed manually
            elif (meshingoption == 'size'):
                self.cadlines.append(self.cubitcmdopen+volume_string+' size '+str(sizingfunction)+self.cubitcmdclose)
            # Meshing the volumes of group_list:
            self.cadlines.append(self.cubitcmdopen+'mesh '+volume_string+self.cubitcmdclose)

    # Now that all volumes of the group_list have been meshed, move them:
    self.cadlines.append('\n'+self.commentchar+' Moving meshed turbine in positive y direction:')
    cadline = 'move group'
    groupid = 2 # first created group id in cubit is always 2, 1 is reserved for 'picked'.
    for turbine in group_list:
        for turbine_group in turbine[1:]:
            cadline = cadline+' '+str(groupid)
            groupid = groupid + 1
    cadline = cadline+' y '+str(self.hub_height)+' include_merged'
    self.cadlines.append(self.cubitcmdopen+cadline+self.cubitcmdclose)


  def export_mesh_cubit(self, group_list, filename=None):
    """ This routine assembles the cubit commands export the mesh to the hard drive.
         Input:
           group_list: List of lists, first element of the inner list is a string of the 
             turbine's name/number, the next element of the inner list, is a list itself with
             the first element being the group's name, and the second element being a list of 
             the corresponding volume ids
           filename: String of the filename of the exported mesh file, by default it is the same 
             basename as the input textfile.
    """
    if (filename is None):
        filename = self.geometryname
    filename = filename.lower()
    self.cadlines.append('\n'+self.commentchar+' Exporting generated mesh to disk.')
    # We want to create as many meshes, as we have defined groups:
    blockid = 1
    for turbine in group_list:
        turbine_name = turbine[0]
        for turbine_group in turbine[1:]:
            meshname_postfix = turbine_name+'_'+turbine_group[0] # e.g. "turbine-1_blades"
            volid = turbine_group[-1] # volume id of that group
            # Exporting the mesh:
            self.cadlines.append(self.cubitcmdopen+'undo group begin'+self.cubitcmdclose)
            self.cadlines.append(self.cubitcmdopen+'set large exodus file on'+self.cubitcmdclose)
            cwd = os.getcwd()
            self.cadlines.append(self.cubitcmdopen+'export genesis \''+str(cwd)+'/'+filename.split('.')[0]+'-'+meshname_postfix+'.exo\' dimension 3 block '+str(blockid)+' overwrite'+self.cubitcmdclose)
            self.cadlines.append(self.cubitcmdopen+'undo group end'+self.cubitcmdclose)
            blockid = blockid + 1

  def mesh_gmsh(self, outputfilename):
    """ This routine simply runs gmsh on the written output geo-file.
         Input:
           outputfilename: String of the filename of the generated geo-file
    """
    print "GMSH output:"
    os.system('gmsh -'+str(self.dimension)+' -bin '+outputfilename)
    
  ############################################################
  # Write assembled data to Cubit journal or GMSH geo file:  #
  ############################################################
  def write_generated_geometry(self, filename):
    """ This routines writes the assembled cubit commands to a textfile
        which can be loaded and run by cubit.
         Input:
           filename: String of the filename, which is
             extended by the correct file extension, whether it is a 
             Cubit journal or GMSH geo file.
    """
    # Open a file and write assembled data to that file:
    outfile = open(filename, 'w')
    for cadline in self.cadlines:
        outfile.write(cadline+'\n')
    outfile.close()


  def get_outputfilename(self, inputfilename):
    """ Simply constructing a string based on input arguments that will
        form the output filename
        Input:
         inputfilename: String of the input filename with the 2D coordinates
        Output:
         outputfilename: String of the output filename
    """
    outputfilename = inputfilename.split('.')[0]
    if (self.dimension == 3): postfix = '-3D'
    elif(self.dimension == 2): postfix = '-2D'
    else: print 'Error, invalid dimension given. Should never get here.'; raise SystemExit()
    if (self.cad_package == 'cubit'): fileextension = 'py'
    elif (self.cad_package == 'gmsh'): fileextension = 'geo'
    else: print 'Error, invalid cad_package name given. Should never get here.'; raise SystemExit()
    outputfilename = outputfilename+postfix+'.'+fileextension
    return outputfilename


  def print_next2be_geometry_ids(self, vertex_list, curve_list, surface_list, volume_list):
    # This routine is for debugging purposes. It prints out the next ID numbers of
    # vertices, curves, surfaces and volumes.
#        self.cadlines.append('\n'+self.commentchar+' =============================================')
#        self.cadlines.append(self.cubitcmdopen+self.commentchar+' Test:'+self.cubitcmdclose)
#        self.cadlines.append(self.cubitcmdopen+self.commentchar+' Next vertex ID:  '+str(vertex_list[-1]+1)+self.cubitcmdclose)
#        self.cadlines.append(self.cubitcmdopen+self.commentchar+' Next curve ID:   '+str(curve_list[-1]+1)+self.cubitcmdclose)
#        self.cadlines.append(self.cubitcmdopen+self.commentchar+' Next surface ID: '+str(surface_list[-1]+1)+self.cubitcmdclose)
#        self.cadlines.append(self.cubitcmdopen+self.commentchar+' Next volume ID:  '+str(volume_list[-1]+1)+self.cubitcmdclose)
#        self.cadlines.append(self.cubitcmdopen+self.commentchar+' =============================================\n'+self.cubitcmdclose)
        self.cadlines.append('\nprint "============================================="')
        self.cadlines.append('print "Test:"')
        self.cadlines.append('print "Next vertex ID:  '+str(vertex_list[-1]+1)+'"')
        self.cadlines.append('print "Next curve ID:   '+str(curve_list[-1]+1)+'"')
        self.cadlines.append('print "Next surface ID: '+str(surface_list[-1]+1)+'"')
        self.cadlines.append('print "Next volume ID:  '+str(volume_list[-1]+1)+'"')
        self.cadlines.append('print "=============================================\n"')


  def generate_geometry(self):
    # Copy values of variables from class variables into local variables:
    textfilename = self.textfilename
    dimension = self.dimension
    global_scale = self.global_scale
    blade_placement = self.blade_placement
    rotation = self.rotation
    extrusion_scale = self.extrusion_scale
    extrusion_length = self.extrusion_length
    blade_curve_type = self.blade_curve_type
    num_blades = self.num_blades
    meshingscheme = self.meshingscheme
    meshingoption = self.meshingoption
    sizingfunction = self.sizingfunction
    turbine_type = self.turbine_type
    num_cs = self.num_cs
    num_turbines = self.num_turbines
    turbine_placement = self.turbine_placement
    cad_package = self.cad_package

    # Initialise:
    vertex_list = []
    curve_list = []
    surface_list = []
    volume_list = []
    group_list = []

    
    # If we are generating a python script for Cubit, import modules first:
    #if (cad_package == 'cubit'):
    #  self.cadlines.append(self.commentchar+' Loading modules:')
    #  self.cadlines.append('import cubit '+self.commentchar+' make sure the path of your "cubit.py" is added to PYTHONPATH\n')

    # Start generating the cubit/gmsh commands:
    # Read in textfile and get coordinates:
    if (isinstance(textfilename, list)): aerofoilfilename=textfilename[0]
    elif (isinstance(textfilename, str)): aerofoilfilename=textfilename
    coordinates = self.get_initial_cross_section_coordinates(aerofoilfilename, global_scale=global_scale)
    vertex_list = self.write_2d_vertices(coordinates, vertex_list, blade_placement, rotation)

    # Number of coordinates from textfile:
    num_2d_coords = len(coordinates)
    self.num_2d_coords = num_2d_coords

    # Creating intermediate cross section points:
    if (dimension == 3):
        for j in range(1, len(extrusion_scale)):
            if (isinstance(textfilename, list)):
                aerofoilfilename = textfilename[j]
                coordinates = self.get_cross_section_coordinates(aerofoilfilename, global_scale=1.0)
                # Check if the length of coordinates is equal to the length of previously extracted coordinates:
                if (len(coordinates) != num_2d_coords):
                    print "=================================================================================================="
                    print "len(coordinates) = ", len(coordinates)
                    print "num_2d_coords = ", num_2d_coords
                    print "ERROR: The number of points per aerofoil cross section must be equal to the number of points "
                    print "used for other aerofoil cross sections along a blade!"
                    print "Exiting program..."
                    raise SystemExit()
            # Then create the scaled cross section (with either updated coordinates or previously used coordinates):
            vertex_list = self.create_scaled_cross_sections_cubit(coordinates, extrusion_scale[j], extrusion_length[j], vertex_list)


    # Creating lines (linear or spline) to connect the points of each cross section generated and their surfaces:
    (vertex_list, curve_list, surface_list) = self.create_lines_surface_cross_sections(num_2d_coords, num_cs, vertex_list, curve_list, surface_list)

    if (cad_package == 'gmsh' and dimension == 2):
        # Set physical ids:
        for i in range(num_blades):
            self.set_physical_ids(i+1, num_2d_coords, curve_list, surface_list)

    # Apply twist to cross sections:
    if (dimension == 3):
        self.apply_twist_to_cross_sections_cubit(rotation=rotation)

        # Create blade's surface:
        (vertex_list, curve_list, surface_list) = self.create_blade_surface_cubit(blade_curve_type, num_2d_coords, num_cs, vertex_list, curve_list, surface_list)

        # Create volume of blade and get list of volume ids:
        (vertex_list, curve_list, surface_list, volume_list) = self.create_blade_volume_cubit(vertex_list, curve_list, surface_list, volume_list, num_2d_coords, num_cs, num_blades=num_blades)
        # Now mesh this initial blade volume, such that the mesh is copied with the volume geometry later on:
        self.mesh_volume_ids_cubit([volume_list[-1]], meshingscheme=meshingscheme, sizingfunction=sizingfunction, meshingoption=meshingoption)

        # Create connecting surface to the hub:
        (vertex_list, curve_list, surface_list, volume_list) = self.create_connection_to_hub_cubit(vertex_list, curve_list, surface_list, volume_list, num_2d_coords, num_cs)

        # Copy the blade and rotate it:
        (vertex_list, curve_list, surface_list, volume_list) = self.create_copies_of_blade_cubit(vertex_list, curve_list, surface_list, volume_list, num_2d_coords, num_cs, num_blades=num_blades)

        # Create hub of the turbine:
        if (turbine_type == 'windturbine'):
            (vertex_list, curve_list, surface_list, volume_list, group_list) = self.create_hub_wind_turbine_cubit(vertex_list, curve_list, surface_list, volume_list, group_list, num_2d_coords, num_cs, num_blades=num_blades)

    # Get outputfilename:
    if (isinstance(textfilename, list)): aerofoilfilename=textfilename[0]
    elif (isinstance(textfilename, str)): aerofoilfilename=textfilename
    outputfilename = self.get_outputfilename(aerofoilfilename)

    # Mesh the volumes:
    if (cad_package == 'cubit'):
        if (dimension == 3):
            self.mesh_turbine_cubit(group_list)
            # For as many turbines as specified, copy (if not first turbine) and move the turbine 
            # around to create several turbines and put them to their specified location:
            self.cadlines.append('\n'+self.commentchar+' Finishing off the turbine(s):')
            blockid_offset = 0
            # Assemble string with all relevant volume ids of the first turbine:
            volume_copy_id_string = ' '.join(str(volid) for volid in group_list[0][1][1])
            volume_copy_id_string = volume_copy_id_string+' '+' '.join(str(volid) for volid in group_list[0][2][1])
            for i in range(1, num_turbines+1):
                if (i > 1):
                    self.cadlines.append('\n'+self.commentchar+' Copying and moving turbine 1 to create turbine number '+str(i))
                    # Movement vector components:
                    move_x = turbine_placement[i-1][0]; move_y = turbine_placement[i-1][1]; move_z = turbine_placement[i-1][2]
                    self.cadlines.append(self.cubitcmdopen+'volume '+volume_copy_id_string+' copy move x '+str(move_x)+' y '+str(move_y)+' z '+str(move_z)+self.cubitcmdclose)
                    # Updating volume_list:
                    volume_list = self.list_id_increase(volume_list, 2+num_blades)
                    # Updating group:
                    group_list.append(['turbine-'+str(i), ['blades', volume_list[-num_blades-2:-2+1]], ['base', [volume_list[-1]]]])

                # Setting block ids for turbine volumes:
                self.cadlines.append(self.commentchar+' Creating block ids for turbine number '+str(i))
                self.create_block_ids_cubit(group_list, blockid_offset)
                # Increasing blockid_offset for next turbine:
                blockid_offset = blockid_offset+2

            # Group list and block ids were set correctly at this point, so now we can write the mesh files to the hard drive:
            # Export the mesh files for this turbine:
            self.export_mesh_cubit(group_list, filename=self.geometryname)


    # Most important step:
    # Write cubit journal textfile:
    self.write_generated_geometry(filename=outputfilename)

    if (cad_package == 'gmsh'):
        self.mesh_gmsh(outputfilename)

    # Done:
    print "\n============================================================================"
    print " Conversion from textfile to a "+cad_package.upper()+" readable file was successful."
    print " Start "+cad_package.upper()+", and run the file:\n    "+outputfilename
    line = " through "
    if (cad_package == 'cubit'):
        line = line+"Cubit's 'journal editor'"
    elif (cad_package == 'gmsh'):
        line = line+"GMSH"
    print line
    print " to automatically generate the desired geometry and mesh files."
    print "============================================================================"

