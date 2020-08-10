## Description 

Repository containing code and ParaView statefiles to execute and visualize output from the Rossby Wave Identification algorithm as detailed in the paper 

An Integrated Geometric and Topological Approach for the Identification and Visual Analysis of Rossby Wave Packets.  
Karran Pandey, Joy Merwin Monteiro, and Vijay Natarajan.  
Monthly Weather Review, 2020.  
https://journals.ametsoc.org/doi/abs/10.1175/MWR-D-20-0014.1

While fully funcitonal code for our algorithm is present in the repository, we are currently still in the process of developing a user friendly API package for convenient use of our algorithm. Hence, please contact us if you wish to use the code in its current form. 

(e-mail: karran13 [AT] gmail.com)

---

## Requirements

To compile and execute this code, you will need installed ParaView binaries (tested with 5.6.0) and python 3.6 along with the following python packages: 

1. vtk
2. scikit-learn
3. networkx
4. numpy/scipy

---

## Instructions for Python version

The python scripts for the pipeline can be found in the Python_Only_version folder. The scripts and their individual input formats are described below: 

1. data_script.py - This script converts an input meridional wind field in .nc format to the desired .vtk format scalar field for further processing. To run the script, the following terminal command is to be executed - 

[ParaView Install Directory]/bin/pvpython [Path to .nc file to be converted] [Desired Output Path] [Repo Directory]/Python_Only_version/

This will store a folder with the name of the nc file in the desired output path containing scalar_field.vtk files for each timestep in the .nc file. 

2. cluster_max.py/cluster_min.py - These scripts take as input the meridional wind field in .vtk format and output the maxima and minima along with their cluster IDs in a vtu format. The vtu can be directly visualized in ParaView. To specify the meridional wind field path, one must modify the following line in the python script - 

file_reader.SetFileName('[path to scalar_field.vtk]')
 
3. graph_computation.py - This script takes as input the clustered maxima and minima in vtu format along with the meridional wind field in vtk format, and outputs the final RWP graph in the vtu format. The output vtu can again be visualized in ParaView. Further it is stored such that subsequent points in the output list constitute a line (for e.g the first and second point in the output list form a line, the third and fourth form a line and so on). This ordered storage allows for easily saving to csv format from within ParaView, or using a csvFileWriter in vtk without losing access to the topology of the graph. To specify input clustered maxima/minima and scalar field, the following paths in the script must be modified: 

file_reader.SetFileName('clustered_max.vtu') - maxima

file_reader.SetFileName('clustered_min.vtu') - minima

file_reader.SetFileName('forecast_bust_0.vtk') - scalar field

The code further contains access to modifiable thresholds - the scalar pruning threshold and the edge threshold as described in our paper. The following lines must be modified to change the thresholds:

scalar_pruned_assoc_graph = scalarThresh(assoc_graph,"Min Scalar Intensity",scalar_thresh,100) - for scalar threshold change the value of scalar_thresh in this line

edge_wt_pruned_graph = scalarThresh(edge_wted_graph,"Updated Dist",edge_thresh,1000) - for edge threshold change the value of edge_thresh in this line

---

## Acknowledgements

This work was partially supported by a Swarnajayanti Fellowship from the Department  of  Science  and  Technology,  India  (DST/SJF/ETA-02/2015-16),  a  Mindtree  Chair  research grant, the Robert Bosch Centre for Cyber Physical Systems, IISc, an IoE research grant from Indian Institute of Science, the Divecha Centre for Climate Change, IISc and the Swedish Research Council (Vetenskapsradet) grant E0531901. 

All datasets used for testing and validation are available at the ECMWF website: https://apps.532ecmwf.int/datasets/data/interim-full-daily/levtype=pl/ 

---

## Copyright

Copyright (c) 2020 Visualization & Graphics Lab (VGL), Indian Institute of Science. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Authors : Karran Pandey

Contact : karran13 [AT] gmail.com
