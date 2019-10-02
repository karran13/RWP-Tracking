# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
import os
from sys import argv

data_path=[]

data_dir=argv[1]
script_dir=argv[2]

for root, dirs, files in os.walk(os.path.abspath(data_dir+"scalar field")):
    for file in files:
        data_path.append(os.path.join(root, file))

data_path.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

# load state
LoadState(script_dir+'min_resub.pvsm', LoadStateDataFileOptions='Choose File Names',
    DataDirectory=script_dir,
    scalar_field_vtkFileNames=data_path)



geodesicClustering = FindSource('Geodesic Clustering')

os.mkdir(data_dir+"clustered min")

SaveData(data_dir+'clustered min/clustered_min.vtu', proxy=geodesicClustering, Writetimestepsasfileseries=1)



#### uncomment the following to render all views
#RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
