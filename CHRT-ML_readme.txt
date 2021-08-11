This repository pertains to an article submitted to XXXXX concerning
the extraction of soundings representing shallow-water bathymetry
from lidar point clouds. The article describes how and where all
1620  .las files in the data base can be accessed.

Avaialble in this repository, is the source Python code for the ML
portion of the CHRT-ML algorithm described in the article. (File
name: Soundings_WrkFlw_ML_CHRT_StreamLined.py)

A prototype implementation of CHRT exists as a source code library
held at the Center for Coastal and Ocean Mapping, and is being
licensed to members of the Industrial Consortium associated with
the Center on a non-exclusive, royalty-free basis. Further details
of the Industrial Consortium and licensing agreements are available
from the authors.

Because all readers may not be able to access CHRT, this repository
includes two CHRT output files ("chrt_aughypos" and "chrt_out") for
two exemplar LiDAR data tiles. The associated .las files can be
downloaded as described in the "Data Availability" statement in
the article. (The naming convention of the .las files is the UTM
easting and northing of the northwestern corner of each .las data
tile.) Once the .las files associated with the CHRT outputs have
been downloaded, interested readers can modify statements in the
file Soundings_WrkFlw_ML_CHRT_StreamLined.py that identify the
location of the CHRT and .las files to run the ML portion of CHRT
ML. 