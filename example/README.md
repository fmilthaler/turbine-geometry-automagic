Dependencies
==========================
TGAM depends on the following software:
  * Python (> v. 2.7, see http://python.org)
  * Cubit (tested with v. 13.0, see http://cubit.sandia.gov)

Example
==========================

## PYTHONPATH
First of all, make sure `TGAM.py` is in a directory which is part of `PYTHONPATH`.

## Example files
  * All "soho_aerof*.txt" files contain x,y coordinates of different blade cross sections. These were obtained from [JavaFoil](http://www.mh-aerotools.de/airfoils/javafoil.htm). Moreover, they were saved in a way such that their filenames allowed for sorting them.
  * The file "make-soho-tidal-turbine.py" is user generated and contains settings to generate the turbine(s), e.g. filenames of the blade cross sections, radius of the hub, height of the hub, number of blades, pitch angles, spin direction, and even option for the meshing. Theses options are currently required to be set in a python file, such as "make-soho-tidal-turbine.py", and are then passed to TGAM.

## Running the example
To run the example, run "make-soho-tidal-turbine.py" through Python: `python make-soho-tidal-turbine.py`.
If everything works out fine it results in a file called "soho_aerof01_naca_4826_31point-3D.py". This can then be run through Cubit's Python interface to generate the geometry and mesh of the specified turbine.