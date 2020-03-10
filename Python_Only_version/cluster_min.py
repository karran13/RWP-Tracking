import vtk
import numpy as np
from sklearn import cluster
from scipy import spatial
import math

def addMaxData(scalar_field):
	scalar_values=scalar_field.GetPointData().GetArray("v")
	#print(extract_edges.GetOutput())
	count=0
	check_flag_arr=np.full(scalar_field.GetNumberOfPoints(),0.0)
	max_flag_arr=np.full(scalar_field.GetNumberOfPoints(),0.0)
	isMax=vtk.vtkIntArray()
	isMax.SetNumberOfComponents(0)
	isMax.SetName("is min")
	VertexIdentifiers=vtk.vtkIntArray()
	VertexIdentifiers.SetNumberOfComponents(0)
	VertexIdentifiers.SetName("VertexIdentifiers_")
	pointLocator=vtk.vtkPointLocator()
	pointLocator.SetDataSet(scalar_field)
	pointLocator.BuildLocator()

	for i in range(scalar_field.GetNumberOfPoints()):
	    VertexIdentifiers.InsertNextValue(i)
	    if(check_flag_arr[i]==1):
	        isMax.InsertNextValue(0)
	        continue
	    cellIdList=vtk.vtkIdList()
	    connectedVertices=vtk.vtkIdList()
	    test_ids=vtk.vtkIdList()
	    pointLocator.FindPointsWithinRadius(1.5,scalar_field.GetPoint(i),test_ids)
	    #print("Point Locator:")
	    #print(test_ids)
	    connectedVertices=test_ids
	    max_flag=1
	    for j in range(connectedVertices.GetNumberOfIds()):
	        if(scalar_values.GetTuple1(connectedVertices.GetId(j))<scalar_values.GetTuple1(i)):
	            max_flag=0
	        else: 
	            check_flag_arr[connectedVertices.GetId(j)]=1
	    if(max_flag==1):
	        max_flag_arr[i]=1
	        isMax.InsertNextValue(1)
	        for j in range(connectedVertices.GetNumberOfIds()):
	            check_flag_arr[connectedVertices.GetId(j)]=1

	        if(scalar_values.GetTuple1(i)<=-5):
	            count+=1
	    else:
	        isMax.InsertNextValue(0)
	print(count)
	scalar_field.GetPointData().AddArray(isMax)
	scalar_field.GetPointData().AddArray(VertexIdentifiers)
	return (scalar_field)	

def interpolateCellVals(inputs):
	numCells=inputs.GetNumberOfCells()
	scalar_v=inputs.GetPointData().GetArray("v")
	cell_scalars=vtk.vtkFloatArray()
	cell_scalars.SetNumberOfComponents(1)
	cell_scalars.SetNumberOfTuples(numCells)
	cell_scalars.SetName("Cell V")
	for i in range(numCells):
	    cell=inputs.GetCell(i)
	    num_points=cell.GetNumberOfPoints()
	    func_value=0
	    for j in range(num_points):
	        pid=cell.GetPointId(j)
	        func_value+=(scalar_v.GetTuple1(pid))
	    func_value/=num_points
	    cell_scalars.SetTuple1(i,func_value)    
        inputs.GetCellData().AddArray(cell_scalars)
	return (inputs)

def clipDataset(dataset,scalar_name,scalar_val):
	clip_dataset=vtk.vtkClipDataSet()
	dataset.GetPointData().SetScalars(dataset.GetPointData().GetArray("v"))
	clip_dataset.SetValue(scalar_val)
	clip_dataset.SetInputData(dataset)
	clip_dataset.InsideOutOn()
	clip_dataset.Update()
	return (clip_dataset.GetOutput())	

def addConnectivityData(dataset):
	connectivity_filter=vtk.vtkConnectivityFilter()
	connectivity_filter.SetInputData(dataset)
	connectivity_filter.SetExtractionModeToAllRegions()
	connectivity_filter.ColorRegionsOn()
	connectivity_filter.Update()
	return (connectivity_filter.GetOutput())

def extractPosMaxIds(scalar_field):
	pos_max_ids=vtk.vtkIdTypeArray()
	num_pts=scalar_field.GetNumberOfPoints()
	is_max_arr=scalar_field.GetPointData().GetArray("is min")
	scalar_arr=scalar_field.GetPointData().GetArray("v")
	for i in range(num_pts):
		if(is_max_arr.GetTuple1(i)==1 and scalar_arr.GetTuple1(i)<=-5):
			pos_max_ids.InsertNextValue(i)
	return pos_max_ids

def extractSelectionIds(scalar_field,id_list):
	selectionNode=vtk.vtkSelectionNode()
	selectionNode.SetFieldType(1)
	selectionNode.SetContentType(4)
	selectionNode.SetSelectionList(id_list)
	selection=vtk.vtkSelection()
	selection.AddNode(selectionNode)
	extractSelection=vtk.vtkExtractSelection()
	extractSelection.SetInputData(0,scalar_field)
	extractSelection.SetInputData(1,selection)
	extractSelection.Update()
	return extractSelection.GetOutput()	


def clusterMax(scalar_field,connectivity_clipped_scalar_field,max_points):
	#import scalar field and critical point data objects
	scalar_field=connectivity_clipped_scalar_field
	maxima_points=max_points
	base_field=scalar_field

	geometryFilter=vtk.vtkGeometryFilter()
	geometryFilter.SetInputData(scalar_field)
	geometryFilter.Update()
	scalar_field=geometryFilter.GetOutput()

	print(scalar_field)

	triangleFilter=vtk.vtkTriangleFilter()
	triangleFilter.SetInputData(scalar_field)
	triangleFilter.Update()
	scalar_field=triangleFilter.GetOutput()

	print(scalar_field)

	maxima_point_id=maxima_points.GetPointData().GetArray("vtkOriginalPointIds")
	num_points=maxima_points.GetNumberOfPoints()
	#maxima_point_id=maxima_points.GetPointData().GetArray("VertexIdentifiers")
	#print(maxima_point_id)

	maxima_regions=maxima_points.GetPointData().GetArray("RegionId")

	point_region_id=scalar_field.GetPointData().GetArray("RegionId")
	num_regions = int(np.max(point_region_id)+1)

	dist_matrix=np.full((num_points,num_points),400)

	dijkstra=vtk.vtkDijkstraGraphGeodesicPath()
	dijkstra.SetInputData(scalar_field)

	#region_distance_array=[[[0 for col in range(0)]for row in range(0)]for clusters in range(num_regions)]

	locator=vtk.vtkCellLocator()
	locator.SetDataSet(base_field)
	locator.BuildLocator()
	cellIds=vtk.vtkIdList()

	cell_v=base_field.GetCellData().GetArray("Cell V")

	co_ords=np.empty((0,3))
	for i in range(num_points):
		co_ords=np.append(co_ords,[maxima_points.GetPoint(i)],axis=0)

	for i in range(num_points):
	    for j in range(i+1,num_points):
	        min_v=1000
	        max_v=-1000
	        av_v=0
	        p0 = [0,0,0]
	        p1 = [0,0,0]
	        dist = 0.0
	        region_1=maxima_regions.GetTuple1(i)
	        region_2=maxima_regions.GetTuple1(j)
	        if(region_1!=region_2):
	            continue
	        dijkstra.SetStartVertex(int(maxima_point_id.GetTuple1(i)))
	        dijkstra.SetEndVertex(int(maxima_point_id.GetTuple1(j)))
	        dijkstra.Update()
	        pts = dijkstra.GetOutput().GetPoints()
	        for ptId in range(pts.GetNumberOfPoints()-1):
	            pts.GetPoint(ptId, p0)  
	            pts.GetPoint(ptId+1, p1)
	            dist += math.sqrt(vtk.vtkMath.Distance2BetweenPoints(p0, p1))
	        dist_matrix[i][j]=dist
	        dist_matrix[j][i]=dist
	        locator.FindCellsAlongLine(co_ords[i],co_ords[j],0.001,cellIds)
	        for k in range(cellIds.GetNumberOfIds()):
	            if(cell_v.GetTuple1(cellIds.GetId(k))>max_v):
	                max_v=cell_v.GetTuple1(cellIds.GetId(k))
		dist_matrix[i][j]=dist_matrix[i][j]+max_v
	        dist_matrix[j][i]=dist_matrix[i][j]

	
	region_array=[[0 for col in range(0)]for row in range(num_regions)]
	cluster_assign=np.full(num_points,0)
	
	median_dist=-np.median(dist_matrix)

	print('dist matrix computed')
	
	for i in range(num_points):
	    #print(point_region_id.GetTuple1(int(maxima_point_id.GetTuple1(i))))
	    region_array[int(point_region_id.GetTuple1(int(maxima_point_id.GetTuple1(i))))].append(i)
	    
	
	prev_max=0
	
	for k in range(num_regions):
	    if(len(region_array[k])==1):
	        cluster_assign[region_array[k][0]]=prev_max
	        prev_max+=1 
	        continue
	    if(len(region_array[k])==2):
	        cluster_assign[region_array[k][0]]=prev_max
	        cluster_assign[region_array[k][1]]=prev_max
	        prev_max+=1 
	        continue
	        
	#    print(len(region_array[k]))
	    num_cluster=int(len(region_array[k]))
	    new_dist=np.full((num_cluster,num_cluster),0)
	    #print(new_dist)
	
	    for i in range(num_cluster):
	        for j in range(i+1,num_cluster):
	            new_dist[i][j]=dist_matrix[region_array[k][i]][region_array[k][j]]
	            new_dist[j][i]=new_dist[i][j]
	
	#print(new_dist)
	    
	        
	    if(num_cluster==0):
	        continue
	
	    sim_matrix=np.negative(new_dist)
	
	    #print(sim_matrix)
	
	    af_clustering = cluster.AffinityPropagation(preference=np.full(num_cluster,median_dist/5.0),affinity='precomputed')
	    af_clustering.fit(sim_matrix)
	    clusters=af_clustering.labels_ + prev_max
	    prev_max=np.max(clusters)+1
	    #print(clusters)
	
	    for i in range(num_cluster):
	        cluster_assign[region_array[k][i]]=clusters[i]
	    #print(cluster_assign)
	
	
	
	
	cluster_id=vtk.vtkIntArray()
	cluster_id.SetNumberOfComponents(1)
	cluster_id.SetNumberOfTuples(num_points)
	cluster_id.SetName("Cluster ID")
	
	
	#print(cluster_assign)
	
	
	for i in range(num_points):
	    cluster_id.SetTuple1(i,cluster_assign[i])
	
	
	#clustered_output=self.GetOutput()
	maxima_points.GetPointData().AddArray(cluster_id)
	#clustered_output.ShallowCopy(maxima_points)
	return maxima_points
	#print(dijkstra.GetOutput())
	
	#print(pts)
	
	
	
	
file_reader=vtk.vtkRectilinearGridReader()

file_reader.SetFileName('forecast_bust_0.vtk')

file_reader.Update()

scalar_field = file_reader.GetOutput()

print(scalar_field)

scalar_field = addMaxData(scalar_field)

print(scalar_field)

scalar_field = interpolateCellVals(scalar_field)

print(scalar_field)

clipped_scalar_field = clipDataset(scalar_field,"v",0)

print(clipped_scalar_field)

connectivity_clipped_scalar_field = addConnectivityData(clipped_scalar_field)

print(connectivity_clipped_scalar_field)

max_points = extractSelectionIds(connectivity_clipped_scalar_field, extractPosMaxIds(connectivity_clipped_scalar_field))

print(max_points)

max_points = clusterMax(scalar_field,connectivity_clipped_scalar_field,max_points)

print(max_points)

vtuFileWriter=vtk.vtkXMLUnstructuredGridWriter()
vtuFileWriter.SetInputDataObject(max_points)
vtuFileWriter.SetFileName('clustered_min.vtu')
vtuFileWriter.Update()

