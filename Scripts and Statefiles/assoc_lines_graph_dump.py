from paraview.simple import *
import os
from sys import argv

clustered_max_path=[]
clustered_min_path=[]
scalar_field_path=[]

data_dir=argv[1]
script_dir=argv[2]


for root, dirs, files in os.walk(os.path.abspath(data_dir+"clustered max")):
    for file in files:
        clustered_max_path.append(os.path.join(root, file))

clustered_max_path.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

for root, dirs, files in os.walk(os.path.abspath(data_dir+"clustered min")):
    for file in files:
        clustered_min_path.append(os.path.join(root, file))

clustered_min_path.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

for root, dirs, files in os.walk(os.path.abspath(data_dir+"scalar field")):
    for file in files:
        scalar_field_path.append(os.path.join(root, file))

scalar_field_path.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))




LoadState(script_dir+'resub_assoc_path_ranking_rep_selection.pvsm', LoadStateDataFileOptions='Choose File Names',
    DataDirectory=script_dir,
    clustered_max_vtuFileName=clustered_max_path,
    clustered_min_vtuFileName=clustered_min_path,
    scalar_field_vtkFileNames=scalar_field_path)


renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [991, 492]

# set active view
SetActiveView(renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [5.584622422587262, -8.05339805714692, 1061.03884213629]
renderView1.CameraFocalPoint = [5.584622422587262, -8.05339805714692, 150.0]
renderView1.CameraParallelScale = 159.57758956061497


#os.mkdir(data_dir+"RWP snapshots")

# save animation
#SaveAnimation(data_dir+'RWP snapshots/rwp.png', renderView1, ImageResolution=[991, 492],
#    FrameWindow=[0, 123])

pathopt = FindSource('Threshold4')

# set active source
SetActiveSource(pathopt)

os.mkdir(data_dir+"graph_data")

# save data
SaveData(data_dir+'graph_data/assoc_lines.vtu', proxy=pathopt, Writetimestepsasfileseries=1)
