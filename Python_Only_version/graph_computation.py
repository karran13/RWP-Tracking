import vtk
import numpy as np
import math
from collections import defaultdict
from scipy.spatial import distance
import networkx as nx
import itertools

def hav_distance(lat1,lon1,lat2,lon2):
    
    r_earth = 6371.0
    
    circum = 2*np.pi*r_earth*np.cos(np.radians(30))
        
    dlat = math.radians(lat1 - lat2)
    
    dlon = math.radians(lon1 - lon2)
    
    a = (math.sin(dlat/2))**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * (math.sin(dlon/2))**2
    c = 2 * np.arctan2(math.sqrt(a), math.sqrt(1-a))
    distance = r_earth * c           

    #print(wavenumber)
    return distance

def getIsoContour(scalar_field,val):
	contourFilter=vtk.vtkContourFilter()
	scalar_field.GetPointData().SetScalars(scalar_field.GetPointData().GetArray("v"))
	contourFilter.SetValue(0,val)
	contourFilter.SetInputData(scalar_field)
	contourFilter.Update()
	return contourFilter.GetOutput()

def computeGradients(scalar_field):
	gradientFilter=vtk.vtkGradientFilter()
	scalar_field.GetPointData().SetScalars(scalar_field.GetPointData().GetArray("v"))
	gradientFilter.SetInputData(scalar_field)
	gradientFilter.Update()
	return gradientFilter.GetOutput()	

def computeAssocGraph(max_points,min_points,iso_contour):
	num_contour_pts=iso_contour.GetNumberOfPoints()
	point_grad=iso_contour.GetPointData().GetArray("Gradients")
	max_cluster_ids=max_points.GetPointData().GetArray("Cluster ID")
	min_cluster_ids=min_points.GetPointData().GetArray("Cluster ID")
	num_max_pts=max_points.GetNumberOfPoints()
	num_min_pts=min_points.GetNumberOfPoints()
	
	num_max_clusters=int(np.max(max_cluster_ids)+1)
	
	num_min_clusters=int(np.max(min_cluster_ids)+1)
	
	
	cluster_max_arr=np.full(num_max_clusters,0.0)
	
	cluster_min_arr=np.full(num_min_clusters,100.0)
	
	cluster_max_point=np.full((num_max_clusters,2),0.0)
	
	cluster_min_point=np.full((num_min_clusters,2),0.0)

	assoc_index_array=np.full((num_max_clusters,num_min_clusters),0.0)

	line_dir_array=np.full((num_max_clusters,num_min_clusters),0.0)

	assoc_set=set()

	max_scalars=max_points.GetPointData().GetArray("v")
	min_scalars=min_points.GetPointData().GetArray("v")

	for i in range(num_max_pts):
	    if(cluster_max_arr[int(max_cluster_ids.GetTuple1(i))]<max_scalars.GetTuple1(i)):
	        cluster_max_arr[int(max_cluster_ids.GetTuple1(i))]=max_scalars.GetTuple1(i)
	        cluster_max_point[int(max_cluster_ids.GetTuple1(i))][0]=max_points.GetPoint(i)[0]
	        cluster_max_point[int(max_cluster_ids.GetTuple1(i))][1]=max_points.GetPoint(i)[1]
	#print(cluster_point_array)

	for i in range(num_min_pts):
	    if(cluster_min_arr[int(min_cluster_ids.GetTuple1(i))]>min_scalars.GetTuple1(i)):
	        cluster_min_arr[int(min_cluster_ids.GetTuple1(i))]=min_scalars.GetTuple1(i)
	        cluster_min_point[int(min_cluster_ids.GetTuple1(i))][0]=min_points.GetPoint(i)[0]
	        cluster_min_point[int(min_cluster_ids.GetTuple1(i))][1]=min_points.GetPoint(i)[1]
	

	bound_arr=np.full(num_contour_pts,0)

	bound_pts=vtk.vtkIntArray()
	bound_pts.SetNumberOfComponents(0)
	bound_pts.SetName("boundary flag")

	bound_max=vtk.vtkIntArray()
	bound_max.SetNumberOfComponents(0)
	bound_max.SetName("closest max cluster")

	bound_min=vtk.vtkIntArray()
	bound_min.SetNumberOfComponents(0)
	bound_min.SetName("closest min cluster")


	assoc_dict={(-1,-1):0}
	max_num_bd_pts_dict={-1:0}
	min_num_bd_pts_dict={-1:0}
	


	for i in range(num_contour_pts):
	    contour_point = iso_contour.GetPoint(i)
	    max_dist=1000
	    min_dist=1000
	    max_id=-1
	    min_id=-1
	    curr_max_dir_deriv=0
	    curr_min_dir_deriv=0
	    grad_vector=[point_grad.GetTuple3(i)[0],point_grad.GetTuple3(i)[1]]
	    curr_max_scalar=0
	    curr_min_scalar=0

	    for j in range(num_max_pts):
	        max_point=max_points.GetPoint(j)
	        curr_max_id=max_cluster_ids.GetTuple1(j)
	        max_dir_vector=[max_point[0]-contour_point[0],max_point[1]-contour_point[1]]
	        max_dir_deriv=max_dir_vector[0]*grad_vector[0]+max_dir_vector[1]*grad_vector[1]
	        curr_max_dist=(max_dir_vector[0]**2+max_dir_vector[1]**2)**0.5
	        #if(max_dir_deriv>0):
	        if(curr_max_dist<max_dist):
	            max_dist=curr_max_dist
	            max_id=curr_max_id
	            curr_max_dir_deriv=max_dir_deriv
	            curr_max_scalar=max_scalars.GetTuple1(j)
	            curr_max_x=max_point[0]

    
	    #if(curr_max_dir_deriv<0):
	    #    max_id=-1

	    #if(curr_max_scalar<25):
	    #    max_id=-1

	    for j in range(num_min_pts):
	        min_point=min_points.GetPoint(j)
	        curr_min_id=min_cluster_ids.GetTuple1(j)
	        min_dir_vector=[min_point[0]-contour_point[0],min_point[1]-contour_point[1]]
	        min_dir_deriv=min_dir_vector[0]*grad_vector[0]+min_dir_vector[1]*grad_vector[1]
	        curr_min_dist=(min_dir_vector[0]**2+min_dir_vector[1]**2)**0.5
	        #if(min_dir_deriv<0):
	        if(curr_min_dist<min_dist):
	            min_dist=curr_min_dist
	            min_id=curr_min_id
	            curr_min_dir_deriv=min_dir_deriv
	            curr_min_scalar=min_scalars.GetTuple1(j)
	            curr_min_x=min_point[0]


	    #if(curr_min_dir_deriv>0):
	    #    min_id=-1

	    max_id=int(max_id)
	    min_id=int(min_id)

	    if((max_id,min_id) in assoc_dict):
	        assoc_dict[(max_id,min_id)]+=1
	    else:
	        assoc_dict[(max_id,min_id)]=1

	    if(max_id in max_num_bd_pts_dict):
	        max_num_bd_pts_dict[max_id]+=1
	    else:
	        max_num_bd_pts_dict[max_id]=1        

	    if(min_id in min_num_bd_pts_dict):
	        min_num_bd_pts_dict[min_id]+=1
	    else:
	        min_num_bd_pts_dict[min_id]=1        
   

	    if((int(max_id)!=-1) & (int(min_id)!=-1)):
	        assoc_index_array[int(max_id)][int(min_id)]+=1
	        if(curr_max_x<curr_min_x):
	            line_dir_array[int(max_id)][int(min_id)]+=1
	        else:
	            line_dir_array[int(max_id)][int(min_id)]-=1
    
	    max_id=int(max_id)
	    min_id=int(min_id)
	    assoc_set.add((int(max_id),int(min_id)))

	max_count={-1:0}
	min_count={-1:0}


	print(max_count)
	print(min_count)


	#print(assoc_dict)
	#print(max_num_bd_pts_dict)
	#print(min_num_bd_pts_dict)

	#print(assoc_index_array)
	print(assoc_set)

	track_lines=vtk.vtkPolyData()

	iso_contour.GetPointData().AddArray(bound_pts)
	iso_contour.GetPointData().AddArray(bound_max)
	iso_contour.GetPointData().AddArray(bound_min)


	max_index = vtk.vtkIntArray()
	max_index.SetNumberOfComponents(0)
	max_index.SetName("Max Cell")

	min_index = vtk.vtkIntArray()
	min_index.SetNumberOfComponents(0)
	min_index.SetName("Min Cell")

	min_scalar_intensity_index = vtk.vtkFloatArray()
	min_scalar_intensity_index.SetNumberOfComponents(0)
	min_scalar_intensity_index.SetName("Min Scalar Intensity")

	max_scalar_intensity_index = vtk.vtkFloatArray()
	max_scalar_intensity_index.SetNumberOfComponents(0)
	max_scalar_intensity_index.SetName("Max Scalar Intensity")


	association_index=vtk.vtkFloatArray()
	association_index.SetNumberOfComponents(0)
	association_index.SetName("Assoc Index")

	cluster_dist=vtk.vtkFloatArray()
	cluster_dist.SetNumberOfComponents(0)
	cluster_dist.SetName("Distance")

	line_dir=vtk.vtkFloatArray()
	line_dir.SetNumberOfComponents(0)
	line_dir.SetName("Line dir")


	appendFilter=vtk.vtkAppendPolyData()
	
                

	for elem in assoc_set:

	    if(elem[0]==-1):
	        continue
	    if(elem[1]==-1):
	        continue

	    max_index.InsertNextValue(int(elem[0]))
	    min_index.InsertNextValue(int(elem[1]))
        
	    max_centre=cluster_max_point[elem[0]]
	    min_centre=cluster_min_point[elem[1]]

	    #if(max_centre[0]<min_centre[0]):
	    #    line_dir.InsertNextValue(1)
	    #else:
	    #    line_dir.InsertNextValue(0)

	    if(line_dir_array[elem[0]][elem[1]]>=0):
	        line_dir.InsertNextValue(1)
	    else:    
	        line_dir.InsertNextValue(0)
    
	    min_scalar=0
	    max_scalar=0 

	    if(cluster_max_arr[int(elem[0])]<-cluster_min_arr[int(elem[1])]):
	        min_scalar=cluster_max_arr[int(elem[0])]
	        max_scalar=-cluster_min_arr[int(elem[1])]
	    else:
	        min_scalar=-cluster_min_arr[int(elem[1])]
	        max_scalar=cluster_max_arr[int(elem[0])]
	    assoc_weight=assoc_dict[(elem[0],elem[1])]*(1/max_num_bd_pts_dict[elem[0]]+1/min_num_bd_pts_dict[elem[1]])*0.5
	    distance= ((max_centre[1]-min_centre[1])**2 + (max_centre[0]-min_centre[0])**2)**0.5
	    distance=hav_distance(max_centre[0],max_centre[1],min_centre[0],min_centre[1])
	    #wave_num=wavenumber(max_centre[0],max_centre[1],min_centre[0],min_centre[1])
	    #if(math.isnan(wave_num)):
	    #    distance=(max_scalar+min_scalar)/(distance)
	    #else:
	    #    distance=(max_scalar+min_scalar)/(distance)
	    distance=(max_scalar+min_scalar)/(distance)
	    

	    max_scalar_intensity_index.InsertNextValue(max_scalar)
	    min_scalar_intensity_index.InsertNextValue(min_scalar)
	    association_index.InsertNextValue(assoc_weight)
	    cluster_dist.InsertNextValue(distance)
	    
	    track_points=vtk.vtkPoints()
	    line=vtk.vtkLine()
	    lines=vtk.vtkCellArray()
	    track_points.InsertNextPoint(max_centre[0],max_centre[1],0)
	    track_points.InsertNextPoint(min_centre[0],min_centre[1],0)
	    line.GetPointIds().SetId(0,0)
	    line.GetPointIds().SetId(0,1)
	    lines.InsertNextCell(line)
	    linesPolyData=vtk.vtkPolyData()
	    linesPolyData.SetPoints(track_points)
	    linesPolyData.SetLines(lines)
	    appendFilter.AddInputData(track_lines)
	    appendFilter.AddInputData(linesPolyData)
	    appendFilter.Update()
	    track_lines=appendFilter.GetOutput()
	
	track_lines.GetCellData().AddArray(max_index)
	track_lines.GetCellData().AddArray(min_index)          
	track_lines.GetCellData().AddArray(max_scalar_intensity_index)
	track_lines.GetCellData().AddArray(min_scalar_intensity_index)
	track_lines.GetCellData().AddArray(association_index)
	track_lines.GetCellData().AddArray(cluster_dist)
	track_lines.GetCellData().AddArray(line_dir)

	return track_lines

def scalarThresh(assoc_graph,var_name,scalar_min,scalar_max):
	select_ids=vtk.vtkIdTypeArray()
	num_lines=assoc_graph.GetNumberOfCells()
	scalar_arr=assoc_graph.GetCellData().GetArray(var_name)
	for i in range(num_lines):
		if(scalar_arr.GetTuple1(i)>=scalar_min and scalar_arr.GetTuple1(i)<=scalar_max):
			select_ids.InsertNextValue(i)
	selectionNode=vtk.vtkSelectionNode()
	selectionNode.SetFieldType(0)
	selectionNode.SetContentType(4)
	selectionNode.SetSelectionList(select_ids)
	selection=vtk.vtkSelection()
	selection.AddNode(selectionNode)
	extractSelection=vtk.vtkExtractSelection()
	extractSelection.SetInputData(0,assoc_graph)
	extractSelection.SetInputData(1,selection)
	extractSelection.Update()
	return extractSelection.GetOutput()	

def addEdgeWeights(max_points,min_points,graph_lines,scalar_thresh):
	
	scalar_tol=30

	#print(graph_lines)

	max_scalars=max_points.GetPointData().GetArray('v')
	min_scalars=min_points.GetPointData().GetArray('v')

	max_cluster_ids=max_points.GetPointData().GetArray('Cluster ID')
	min_cluster_ids=min_points.GetPointData().GetArray('Cluster ID')

	num_max_pts=max_points.GetNumberOfPoints()
	num_min_pts=min_points.GetNumberOfPoints()

	max_pt_dict=defaultdict(list)
	min_pt_dict=defaultdict(list)

	for i in range(num_max_pts):
	    scalar=max_scalars.GetTuple1(i)
	    if(scalar>15):
	        cluster_id=max_cluster_ids.GetTuple1(i)
	        point_cords=max_points.GetPoint(i)
	        point_tuple=(point_cords,cluster_id,scalar)
	        max_pt_dict[cluster_id].append(point_tuple)
	
	for i in range(num_min_pts):
	    scalar=min_scalars.GetTuple1(i)
	    if(scalar<-15):
	        cluster_id=min_cluster_ids.GetTuple1(i)
	        point_cords=min_points.GetPoint(i)
	        point_tuple=(point_cords,cluster_id,-scalar)
	        min_pt_dict[cluster_id].append(point_tuple)
	
	line_max_ids=graph_lines.GetCellData().GetArray('Max Cell')
	line_min_ids=graph_lines.GetCellData().GetArray('Min Cell')
	line_max_scalars=graph_lines.GetCellData().GetArray('Max Scalar Intensity')
	line_min_scalars=graph_lines.GetCellData().GetArray('Min Scalar Intensity')
	line_dists=graph_lines.GetCellData().GetArray('Distance')
	line_dirs=graph_lines.GetCellData().GetArray('Line dir')
	
	new_line_dists=vtk.vtkFloatArray()
	new_line_dists.SetName('Updated Dist')
	
	new_line_dirs=vtk.vtkFloatArray()
	new_line_dirs.SetName('Updated Line Dir')
	
	#print(min_pt_dict)
	
	num_lines=graph_lines.GetNumberOfCells()
	
	for i in range(num_lines):
	    max_id=line_max_ids.GetTuple1(i)
	    min_id=line_min_ids.GetTuple1(i)
	    line_max_scalar=line_max_scalars.GetTuple1(i)
	    line_min_scalar=line_min_scalars.GetTuple1(i)
	    line_dist=line_dists.GetTuple1(i)
	    line_dir=line_dirs.GetTuple1(i)
	
	    cluster_max_pts=max_pt_dict[max_id]
	    cluster_min_pts=min_pt_dict[min_id]
	
	    print(line_dist)
	    curr_dist=0.0
	    print(line_dir)
	
	    #new_line_dir=line_dir
	    high_val_flag=0
	
	    if(line_max_scalar>50 and line_min_scalar>50):
	        high_val_flag=1
	
	    for max_pt in cluster_max_pts:
	        if(max_pt[2]<scalar_thresh):
	            continue
	        if(max_pt[2]<line_max_scalar-scalar_tol and high_val_flag==0):
	            continue
	
	        for min_pt in cluster_min_pts:
	            if(min_pt[2]<scalar_thresh):
	                continue
	            if(min_pt[2]<line_min_scalar-scalar_tol and high_val_flag==0):
	                continue
	            curr_dist=((max_pt[0][0]-min_pt[0][0])**2 + (max_pt[0][1]-min_pt[0][1])**2)**0.5
	            curr_dist=hav_distance(max_pt[0][0],max_pt[0][1],min_pt[0][0],min_pt[0][1])
	            curr_dist=(max_pt[2]+min_pt[2])/curr_dist
	
	            if(curr_dist>line_dist):
	                print('gain')
	                line_dist=curr_dist
	
	    print(line_dist)
	    new_line_dists.InsertNextValue(line_dist)
	
	    new_line_dir=0
	
	    for max_pt in cluster_max_pts:
	        if(max_pt[2]<30):
	            continue
	        curr_line_dir=0
	        for min_pt in cluster_min_pts:
	            if(min_pt[2]<30):
	                continue
	            if(min_pt[0][0]<max_pt[0][0]):
	                curr_line_dir-=1
	            else:
	                curr_line_dir+=1
	        if(curr_line_dir>=0):
	            new_line_dir+=1
	        else:
	            new_line_dir-=1
	
	    if(new_line_dir>=0):
	        new_line_dirs.InsertNextValue(1)
	    else:
	        new_line_dirs.InsertNextValue(0)

	    #line_dists.SetTuple1(i,line_dist)
	graph_lines.GetCellData().AddArray(new_line_dists)
	graph_lines.GetCellData().AddArray(new_line_dirs)
	return graph_lines
    
	

def get_path_wt(path,order,max_data_dict,min_data_dict):
    path_wt=0
    for i in range(len(path)-1):
        if(order==1):
            if(i%2==0):
                curr_edge_wt=get_edge_wt(path[i],path[i+1],max_data_dict,min_data_dict)
                if(max_data_dict[path[i]][0][0]>min_data_dict[path[i+1]][0][0]):
                    print('issue')
                    #curr_edge_wt=-curr_edge_wt
                    return -1
            else:
                curr_edge_wt=get_edge_wt(path[i+1],path[i],max_data_dict,min_data_dict)
                if(min_data_dict[path[i]][0][0]>max_data_dict[path[i+1]][0][0]):
                    print('issue')
                    #curr_edge_wt=-curr_edge_wt
                    return -1
        else:
            if(i%2==0):
                curr_edge_wt=get_edge_wt(path[i+1],path[i],max_data_dict,min_data_dict)
                if(min_data_dict[path[i]][0][0]>max_data_dict[path[i+1]][0][0]):
                    print('issue')
                    #curr_edge_wt=-curr_edge_wt
                    return -1
            else:
                curr_edge_wt=get_edge_wt(path[i],path[i+1],max_data_dict,min_data_dict)
                if(max_data_dict[path[i]][0][0]>min_data_dict[path[i+1]][0][0]):
                    print('issue')
                    #curr_edge_wt=-curr_edge_wt
                    return -1
        #curr_edge_wt=G.edges()[path[i],path[i+1]]['weight']
        path_wt+=curr_edge_wt
    return path_wt

def get_edge_wt(max_id,min_id,max_data_dict,min_data_dict):
#    min_id=-min_id
#    if(min_id==-100):
#        min_id=0
    max_scalar=max_data_dict[max_id][2]
    max_pt=max_data_dict[max_id][0]
    min_scalar=min_data_dict[min_id][2]
    min_pt=min_data_dict[min_id][0]
    dist=((max_pt[0]-min_pt[0])**2 + 4*(max_pt[1]-min_pt[1])**2)**0.5
    edge_wt=(max_scalar+min_scalar)/dist
    return edge_wt

def getRankedPaths(max_points,min_points,graph_lines):        

	max_scalars=max_points.GetPointData().GetArray('v')
	min_scalars=min_points.GetPointData().GetArray('v')

	max_cluster_ids=max_points.GetPointData().GetArray('Cluster ID')
	min_cluster_ids=min_points.GetPointData().GetArray('Cluster ID')

	num_max_pts=max_points.GetNumberOfPoints()
	num_min_pts=min_points.GetNumberOfPoints()

	max_pt_dict=defaultdict(list)
	min_pt_dict=defaultdict(list)

	max_data_dict={}
	min_data_dict={}

	scalar_thresh=30

	for i in range(num_max_pts):
	    scalar=max_scalars.GetTuple1(i)
	    if(scalar>scalar_thresh):
	        cluster_id=max_cluster_ids.GetTuple1(i)
	        point_cords=max_points.GetPoint(i)
	        point_tuple=(point_cords,cluster_id,scalar)
	        max_pt_dict[cluster_id].append(i)
	        max_data_dict[i]=point_tuple
	
	for i in range(num_min_pts):
	    scalar=min_scalars.GetTuple1(i)
	    if(scalar<-scalar_thresh):
	        cluster_id=min_cluster_ids.GetTuple1(i)
	        point_cords=min_points.GetPoint(i)
	        point_tuple=(point_cords,cluster_id,-scalar)
	        min_pt_dict[cluster_id].append(i)
	        min_data_dict[i]=point_tuple
	
	#print(max_pt_dict)
	#print(min_pt_dict)
	#print(get_edge_wt(11,-66))
	
	max_ids=graph_lines.GetCellData().GetArray("Max Cell")
	min_ids=graph_lines.GetCellData().GetArray("Min Cell")
	dist_wts=graph_lines.GetCellData().GetArray("Updated Dist")
	line_dirs=graph_lines.GetCellData().GetArray("Updated Line Dir")
	
	num_lines=graph_lines.GetNumberOfCells()
	
	edge_list=[]
	
	for i in range(num_lines):
	    max_id=max_ids.GetTuple1(i)
	    min_id=min_ids.GetTuple1(i)
	    if(min_id==0):
	        min_id=100
	    dist=dist_wts.GetTuple1(i)
	    line_dir=line_dirs.GetTuple1(i)
	    if(int(line_dir)==1):
	        edge_list.append((max_id,-min_id,dist))
	    else:
	        edge_list.append((max_id,-min_id,dist)) 
	    
	   
	G = nx.Graph()
	G.add_weighted_edges_from(edge_list)
	path_list=[]
	
	#print(G.edges())
	
	
	start_leaves=[x for x in G.nodes()]
	end_leaves=[x for x in G.nodes()]
	
	#print(start_leaves)
	#print(end_leaves)
	
	for source in start_leaves:
	    for sink in end_leaves:
	        if(nx.has_path(G,source=source,target=sink)):
	            for path in nx.all_simple_paths(G,source=source,target=sink):
	                #print(path)
	                #print(set(path))
	                path_list.append(path)
	
	#print(path_list)
	
	path_pt_dict={}
	path_wt_dict={}
	path_order_dict={}
	
	for path in path_list:
	    cluster_lists=[]
	    order_flag=0
	    if(path[0]>=0):
	        order_flag=1
	    for node in path:
	        if(node>=0):
	            cluster_lists.append(max_pt_dict[node])
	        else:
	            if(node==-100):
	                cluster_lists.append(min_pt_dict[0])
	            else:
	                cluster_lists.append(min_pt_dict[-node])
	    #print(cluster_lists)
	    path_combinations=itertools.product(*cluster_lists)
	    #print(list(path_combinations))
	    max_wt=0
	    for comb in list(path_combinations):
	        curr_wt=(get_path_wt(comb,order_flag,max_data_dict,min_data_dict))
	        if(curr_wt>max_wt):
	            max_wt=curr_wt
	            path_pt_dict[tuple(path)]=comb
	    path_wt_dict[tuple(path)]=max_wt
	    path_order_dict[tuple(path)]=order_flag
	
	
	print(path_pt_dict)
	print(path_wt_dict)
	#    print(get_path_wt(path))    




	top_paths=list(filter(lambda f: not any([(path_wt_dict[tuple(f)]<path_wt_dict[tuple(g)] and len(set(f)&set(g))!=0) for g in 		path_list]),path_list)) 

	#print(top_paths)

	#print(len(path_list))
	#print(len(top_paths))

	relevant=vtk.vtkIntArray()
	relevant.SetNumberOfComponents(0)
	relevant.SetName("Top Path")


	path_points=vtk.vtkPoints()
	path_polydata=vtk.vtkPolyData()
	path_lines=vtk.vtkCellArray()

	edge_wts=vtk.vtkFloatArray()
	edge_wts.SetName('Edge Weights')

	curr_pt_id=0


	for path in top_paths:
	    order_flag=path_order_dict[tuple(path)]            
	    path_point_ids=path_pt_dict[tuple(path)]
	    for i in range(len(path_point_ids)):
	        point_id=path_point_ids[i]
	        if(order_flag==1):
	            if(i%2==0):
	                point=max_data_dict[point_id][0]
	                if(i<len(path_point_ids)-1):
	                    line_wt=get_edge_wt(point_id,path_point_ids[i+1],max_data_dict,min_data_dict)
	            else:
	                point=min_data_dict[point_id][0]
	                if(i<len(path_point_ids)-1):
	                    line_wt=get_edge_wt(path_point_ids[i+1],point_id,max_data_dict,min_data_dict)
	        else:
	            if(i%2==0):
	                point=min_data_dict[point_id][0]
	                if(i<len(path_point_ids)-1):
	                    line_wt=get_edge_wt(path_point_ids[i+1],point_id,max_data_dict,min_data_dict)
	            else:
	                point=max_data_dict[point_id][0]
	                if(i<len(path_point_ids)-1):
	                    line_wt=get_edge_wt(point_id,path_point_ids[i+1],max_data_dict,min_data_dict)
	        path_points.InsertNextPoint([point[0],point[1],0])
	        if(i<len(path_point_ids)-1):
	            line=vtk.vtkLine()
	            line.GetPointIds().SetId(0,curr_pt_id)
	            line.GetPointIds().SetId(1,curr_pt_id+1)
	            path_lines.InsertNextCell(line)
	            edge_wts.InsertNextValue(line_wt)
	            print(curr_pt_id)
	            print(curr_pt_id+1)                      
	        curr_pt_id+=1
	            
	
	#print(path_points)
	path_polydata.SetPoints(path_points)
	path_polydata.SetLines(path_lines)
	path_polydata.GetCellData().AddArray(edge_wts) 
	#print(path_polydata)                       
	
	return path_polydata

def csvOutput(graph_lines):

	csv_output=vtk.vtkPolyData()

	csv_points=vtk.vtkPoints()

	csv_lines=vtk.vtkCellArray()

	num_cells=graph_lines.GetNumberOfCells()

	curr_pt_id=0


	for i in range(num_cells):
	    line=vtk.vtkLine()
	    curr_line=graph_lines.GetCell(i)
	    curr_pts=curr_line.GetPoints()
	    csv_points.InsertNextPoint(curr_pts.GetPoint(0))
	    csv_points.InsertNextPoint(curr_pts.GetPoint(1))
	    line.GetPointIds().SetId(0,curr_pt_id)
	    line.GetPointIds().SetId(1,curr_pt_id+1)
	    csv_lines.InsertNextCell(line)
	    curr_pt_id+=2

	csv_output.SetPoints(csv_points)
	csv_output.SetLines(csv_lines)
	#print(csv_output)
	return (csv_output)
    
#insert points, update curr_pt_id insert lines

file_reader=vtk.vtkXMLUnstructuredGridReader()

file_reader.SetFileName('clustered_max.vtu')

file_reader.Update()

maxima_points = file_reader.GetOutput()

file_reader=vtk.vtkXMLUnstructuredGridReader()

file_reader.SetFileName('clustered_min.vtu')

file_reader.Update()

minima_points = file_reader.GetOutput()

file_reader=vtk.vtkRectilinearGridReader()

file_reader.SetFileName('forecast_bust_0.vtk')

file_reader.Update()

scalar_field = file_reader.GetOutput()

scalar_field = computeGradients(scalar_field)

print(scalar_field)

zero_contour = getIsoContour(scalar_field,0.0)

print (zero_contour)

assoc_graph = computeAssocGraph(maxima_points,minima_points,zero_contour)

print(assoc_graph)

scalar_thresh=30

scalar_pruned_assoc_graph = scalarThresh(assoc_graph,"Min Scalar Intensity",scalar_thresh,100)

print(scalar_pruned_assoc_graph)

edge_wted_graph = addEdgeWeights(maxima_points,minima_points,scalar_pruned_assoc_graph,scalar_thresh)

print(edge_wted_graph)

edge_thresh=0.02

edge_wt_pruned_graph = scalarThresh(edge_wted_graph,"Updated Dist",edge_thresh,1000)

print(edge_wt_pruned_graph)

ranked_paths_graph = getRankedPaths(maxima_points,minima_points,edge_wt_pruned_graph)

print(ranked_paths_graph)

ordered_final_RWP_graph = csvOutput(ranked_paths_graph)

print(ordered_final_RWP_graph)

vtuFileWriter=vtk.vtkXMLPolyDataWriter()
vtuFileWriter.SetInputDataObject(ordered_final_RWP_graph)
vtuFileWriter.SetFileName('RWP_Graph.vtp')
vtuFileWriter.Update()


