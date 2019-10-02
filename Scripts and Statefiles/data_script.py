# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
import os
from sys import argv

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

data_dir=argv[1]
output_dir=argv[2]
script_dir=argv[3]

# load state
LoadState(script_dir+'nc_to_vtk.pvsm', LoadStateDataFileOptions='Choose File Names',
    DataDirectory=script_dir,
    a_300_areancFileName=[data_dir])

# find source
transform1 = FindSource('Transform1')

# set active source
SetActiveSource(transform1)

# save data
print(output_dir+data_dir.rsplit('/',1)[-1])
os.mkdir(output_dir+data_dir.rsplit('/',1)[-1])
output_dir=output_dir+data_dir.rsplit('/',1)[-1]
os.mkdir(output_dir+'/scalar field')
SaveData(output_dir+'/scalar field/scalar_field.vtk', proxy=transform1, Writetimestepsasfileseries=1)

