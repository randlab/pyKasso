'''Python tools to interact with EPA SWMM.
Chloé Fandel, June 2022.
cfandel@carleton.edu

Add-on module in development for pyKasso.
Creates SWMM input files using .txt or .csv files of nodes and links (designed to work with pyKasso)
Uses a template .inp file to enable editing data in different sections (junctions, conduits, etc.) of the input file. 
Data can be added from .csv files, or can be created as pandas dataframes.
Can also run SWMM directly from Python.

See SWMM documentation for descriptions of parameters.

Notes: 
- Not all sections and not all parameters are enabled for editing yet, and many defaults are not yet modifiable, 
but the basic structure of the code should make it easy to add functions in the same pattern.
- template.inp files must not have any commas in them (otherwise the importer thinks they are column breaks)
- placeholders in template.inp files must not have any other characters (even spaces) in the same line
- order of placeholder list and data filename list must be the same
- common error: the template filename isn't correct, so the insert_data() function can't find the section where it's supposed to insert the data
- MANY SWMM parameters are currently just being set as default, but this can be modified as needed
- To print timeseries to the .rpt text file, run the SWMM command as >> swmm5 projecct.inp project.rpt (leaving off the .out file from the command)
  This will be slower, but more human-friendly. 
'''


#IMPORTS#
########################################################################################
import sys
import pandas as pd
import numpy as np
import subprocess as sp
#import rasterio
#import gdal
#import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import shapely
#from shapely import geometry as shg
#from shapely import ops as sho
#import shapely.speedups
#import pyproj
#from pysheds.grid import Grid
import networkx as nx
import karstnet as kn
#from swmmtoolbox import swmmtoolbox as st
#import mapping



#MODEL SETUP#
#########################################################################################

def import_template(template_filename='template.inp'):
    
    '''Imports the template.inp file into a pandas dataframe.
    Inputs:
    template_filename: filename/path for the template.inp file
                       this must be a text file correctly formatted for SWMM, with placeholder strings in the sections to be edited
                       to create a sample template file, use the tutorial1.inp file from the SWMM documentation, and replace the data
                       in the [JUNCTIONS] section with the word 'junctions'
    Outputs:
    template:          pandas dataframe of the template file, with one column, and a line for each line of text.'''
    
    template = pd.read_csv(template_filename, header=None, skip_blank_lines=False)  #import template .inp file to dataframe
    
    return template


#########################################################################################

def import_data(data_filename):
    
    '''Imports data from a csv file into a pandas dataframe with columns.
    Inputs:
    data_filename:     filename/path for the data.csv file for a section
                       this must be a csv file of data formatted with the first row as a comment row, 
                       and the second row with the correct SWMM column headers for that section 
                       (see SWMM documentation for description of data in each section)
    Outputs:
    template:          pandas dataframe of the data file, with columns, and a line for each item'''
    
    data = pd.read_csv(data_filename, header=1, na_filter=False)    #import data to dataframe that has separate columns
    
    return data


#########################################################################################
def importDEM(filename, show=False):
    '''Import DEM from a tif file using gdal package.
    Return a dem object, and xyz extent and resolution.
    (this can be used to set the model extent)
    NOTE: vertical (z) resolution can't be extracted from the raster!
    
    filename: string indicating the filename (must be a rectangular tif)
    show:     option to show a plot of the DEM or not.
    
    xmin:     minimum x value (same for ymin, zmin)
    xmax:     maximum x value (same for ymax, zmax)
    xres:     x resolution, aka number of columns, aka number of cells along x axis (NOT pixel width)
    etc.
    '''
    
    dem = gdal.Open(filename)    #DEM must be rectangular tif 
    dema = dem.ReadAsArray()     #copy of DEM as a numpy array (defaults to integers)
    dema = dema.astype(float)    #convert integer array to float array
    dema[dema==0] = np.nan       #replace zeros with NaNs (have to convert array to float first)

    ulx, pixelwidthx, xskew, uly, yskew, pixelheighty = dem.GetGeoTransform() #get resolution and coordinate info (for some reason the order of skew and pixel size is flipped for y axis?!)
    ncol = dem.RasterXSize            #number of columns (aka number of cells along x axis)
    nrow = dem.RasterYSize            #number of rows (aka number of cells along y axis)
    lrx = ulx + (ncol * pixelwidthx)  #lower right x coord = upper left x coord + (width of raster cells in x direction * number of raster cells in x direction)
    lry = uly + (nrow * pixelheighty)

    #Get min and max elevations (z):
    #note: gdal's built-in GetRasterBand and GetStatistics return an incorrect zmin (WHY?!)
    zmin = np.nanmin(dema)
    zmax = np.nanmax(dema)
    
    #Assign useful names:
    xmin = ulx
    xmax = lrx
    xres = ncol
    dx =   abs(pixelwidthx)
    ymin = lry
    ymax = uly
    dy =   abs(pixelheighty)
    yres = nrow
    zres = 'na'     #can't be extracted from raster

    #Print results & display raster:
    if show==True:
        print('Raster dimensions: \nxmin: {:<12} xmax: {:<12} xres: {} \nymin: {:<12} ymax: {:<12} yres: {} \nzmin: {:<12} zmax: {:<12} zres: {}'.format(
            xmin,xmax,xres,ymin,ymax,yres,zmin,zmax,zres))
        plt.imshow(dema, extent=(xmin,xmax,ymin,ymax), vmin=zmin, vmax=zmax) #plot raster as image
        #print(gdal.Info(dem))  #for more detailed file info, uncomment this line
        
    return dem,dema,xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres


#########################################################################################

def sks2swmm(nodes_file, links_file, spring_file=None, dim=3, simplify=True, springID=4, elev_file=None):
    
    '''Imports network info (nodes and links) from SKS outputs to pandas dataframes for nodes and links. 
    
    Inputs:
    nodes_file:  text file with two (or three) columns, for X,Y (and Z) coordinates of simplified node locations
    links_file:  text file with two columns, for upstream node and downstream node of each link
    spring_file: text file with three columns, for X,Y,Z coordinates of springs (aka outlets). If None, outlets will not be identified and SWMM cannot be run.
    elev_file:   text file with one columns, for elevation of each node (only if using a 2D network with elevations added after)
    dim:         number of dimensions (2 or 3)
    simplify:    True/False, whether to simplify the SKS network or not (defaults to True)
    springID:    integer indicating spring/outlet identification method 
                 0: outfalls=nodes within a certain tolerance of springs, 1: outfalls=nodes with min distance to springs, 2: outfall=lowest-elev node, 
                 3: outfall=last node, 4: outfalls= new nodes at spring XYZ, with new links to closest junction node, 
                 5: outfalls = shift nodes with min distance to springs to spring coord (current best method)
    
    Outputs:
    nodes:       pandas df with all nodes, and columns: X, Y, Z, type, Name. Type is either 'junction', or 'outfall'
    links:       pandas df with all links, and columns: InNode, OutNode, Name
    '''
    
    #Import & simplify complex network with karstnet:
    if simplify==True:
        print('importing & simplifying network')
        #Import:
        name = nodes_file[:-10]                                 #get base name string for karstnet format (nodes & links must have same base name, with extension '_nodes.dat'and '_links.dat')
        network = kn.from_nodlink_dat(name)                     #import network using karstnet
        nodesSimple = network.graph_simpl.nodes                 #get node indices for simple network (0-based indexing)
        links = network.graph_simpl.edges                       #get edge tuples (aka links or pipes) for simple network

        #Convert data to format used by swmmpy (pandas dfs):
        nodes = pd.DataFrame.from_dict(network.pos3d, orient='index', columns=['X','Y','Z']) #create a df from node coord dic for entire unsimplified network
        nodes = nodes[nodes.index.isin(nodesSimple)]                    #select only the nodes in the simplified network
        nodes['Name'] = np.arange(0,len(nodes))                         #add name column using 0-based indices
        nodes['type'] = 'junction'                                      #add column classifying all nodes as junctions for now
        links = pd.DataFrame(links, columns=['In','Out'])               #convert links to df

        #Replace InNode and OutNode names with re-indexed node names from simplified network:
        ins, outs = [],[]                                                #create empty lists to store info
        for i in range(len(links)):                                      #loop over each link index
            link = links.iloc[i]                                         #get current link
            ins.append(nodes.Name[nodes.index==link.In].values[0])       #get & store new simplified InNode index
            outs.append(nodes.Name[nodes.index==link.Out].values[0])     #same for OutNode
        links['InNode'] = ins                                            #add column with simplified InNode indices
        links['OutNode'] = outs                                          #add column with simplified OutNode indices
        links['Name'] = links.index.values                #add name column with 0-based indices
        links.drop(columns=['In','Out'], inplace=True)    #drop unneeded columns with original complex node indices
        nodes.reset_index(drop=True, inplace=True)        #resart indexing from 0

    #Import an already-simple network:
    else:
        print('importing unsimplified network')
        #Import nodes:
        if dim == 2:
            nodes         = pd.read_csv(nodes_file, header=None, names=['X','Y'], delim_whitespace=True) #import node data from SKS txt output file
            nodes['Z']    = pd.read_csv(elev_file, header=None, names=['Z'])    #import node elev data from SKS and add it as a new column
        if dim == 3:
            nodes         = pd.read_csv(nodes_file, header=None, names=['X','Y','Z'], delim_whitespace=True) #import node data from SKS txt output file        
        nodes['type'] = 'junction'                                              #add a column for node type & fill with junction as default
        nodes['Name'] = nodes.index                                             #add a column for node name (based on index)

        #Import links:
        links         = pd.read_csv(links_file, header=None, names=['InNode','OutNode'], delim_whitespace=True) #import link data from SKS txt output file 
        links         = links - 1                                               #convert node indices to Python 0-based  indices
        links['InNode']  = [int(node) for node in links.InNode]                 #convert to integer
        links['OutNode'] = [int(node) for node in links.OutNode]                #convert to integer
        links['Name'] = links.index                                             #add a column for link name (use row index)

    #Identify springs/outlets:
    #PROBLEM: Currently method 5 is best, but network simplification deletes some springs, so it is still inaccurate. 
    #         Method 4 forces all three springs to be present even if connections are unrealistic
    #         method options: 0: within a certain tolerance of springs, 1: min distance to springs, 2: lowest-elev node, 3: last node, 4: add nodes & links, 5: shift closest nodes to true spring coordinates
    try:
        springs       = pd.read_csv(spring_file, names=['X','Y','Z','system'], header=None, delim_whitespace=True)    #import spring location data
        #print('springs:\n', springs)
        if springID==0:   #select nodes within a certain tolerance
            print('assigning springs by tolerance method')
            for i in range(len(springs)):
                nodes.type[np.isclose(nodes.X, springs.X.iloc[i], atol=.000001, rtol=.001) 
                         & np.isclose(nodes.Y, springs.Y.iloc[i], atol=.000001, rtol=.0001)
                         & np.isclose(nodes.Z, springs.Z.iloc[i], atol=.000001, rtol=.0001)] = 'outfall'     #change node type to outfall based on approx XYZ coord of springs
        if springID==1:   #select nearest nodes 
            print('assigning springs by nearest-node method')
            nodes['dQE'] = ((nodes.X - springs.iloc[0,0])**2 + (nodes.Y - springs.iloc[0,1])**2 + (nodes.Z - springs.iloc[0,2])**2)**0.5
            nodes.type[nodes.dQE==np.min(nodes.dQE)] = 'outfall'
            nodes['dQA'] = ((nodes.X - springs.iloc[1,0])**2 + (nodes.Y - springs.iloc[1,1])**2 + (nodes.Z - springs.iloc[1,2])**2)**0.5
            nodes.type[nodes.dQA==np.min(nodes.dQA)] = 'outfall'
            nodes['dQS'] = ((nodes.X - springs.iloc[2,0])**2 + (nodes.Y - springs.iloc[2,1])**2 + (nodes.Z - springs.iloc[2,2])**2)**0.5
            nodes.type[nodes.dQS==np.min(nodes.dQS)] = 'outfall'
            nodes.drop(labels=['dQE','dQA','dQS'], axis='columns', inplace=True)    #drop  unneeded columns    
        if springID==2:  #select lowest-elev nodes
            print('assigning springs by lowest-elevation method')
            nodes.type[nodes.Z == np.min(nodes.Z)] = 'outfall'                  #change node type to outfall for lowest elev. point
        if springID==3:  #select highest-index nodes
            print('assigning springs by highest-index method')
            nodes.type[nodes.Name==np.max(nodes.index.values)] = 'outfall'      #change node type to outfall for last node
        if springID==4:      #add new nodes at spring XYZ, then link to closest existing node
            print('assigning springs by adding spring nodes and linking to closest junction node')
            springs['type'] = 'outfall'         #add column to assign outfall type to spring nodes
            springs['Name'] = np.max(nodes.index.values) + springs.index + 1    #add column to give springs names (continuing node indexing)
            nodes = nodes.append(springs, ignore_index=True, sort=True)         #append springs to node df
            for i in range(len(springs)):                                       #loop over springs
                spring = nodes[nodes.type=='outfall'].iloc[i]                   #get spring info
                #upnodes = nodes[nodes.Z >= spring.Z]                           #select only upstream nodes (higher elev than spring
                nodes['dist2spring'] = ((nodes.X - spring.X)**2 + (nodes.Y - spring.Y)**2 + (nodes.Z - spring.Z)**2)**0.5  #calculate dist from each point to spring
                junctions = nodes[nodes.type=='junction']                                                                  #get junction nodes only
                nearest = junctions[junctions.dist2spring==np.min(junctions.dist2spring)].index.values[0]                  #find index of closest junction node
                link = pd.DataFrame({'InNode':nearest, 'OutNode':spring.Name, 'Name':np.max(links.index.values)+1}, index=[0])    #identify new link InNode, OutNode, Name
                #print('new link\n', link)
                links = links.append(link, ignore_index=True, sort=True)
        if springID==5:   #select nearest nodes & shift them to true spring coordinates
            print('assigning springs by nearest-node method & shifting them to true coordinates')
            nodes['dQE'] = ((nodes.X - springs.iloc[0,0])**2 + (nodes.Y - springs.iloc[0,1])**2 + (nodes.Z - springs.iloc[0,2])**2)**0.5
            nodes.type[nodes.dQE==np.min(nodes.dQE)] = 'outfall'         #convert to outfall
            nodes.X[nodes.dQE==np.min(nodes.dQE)] = springs.iloc[0,0]    #shift X 
            nodes.Y[nodes.dQE==np.min(nodes.dQE)] = springs.iloc[0,1]    #shift Y
            nodes.Z[nodes.dQE==np.min(nodes.dQE)] = springs.iloc[0,2]    #shift Z
            nodes['dQA'] = ((nodes.X - springs.iloc[1,0])**2 + (nodes.Y - springs.iloc[1,1])**2 + (nodes.Z - springs.iloc[1,2])**2)**0.5
            nodes.type[nodes.dQA==np.min(nodes.dQA)] = 'outfall'
            nodes.X[nodes.dQA==np.min(nodes.dQA)] = springs.iloc[1,0]    #shift X 
            nodes.Y[nodes.dQA==np.min(nodes.dQA)] = springs.iloc[1,1]    #shift Y
            nodes.Z[nodes.dQA==np.min(nodes.dQA)] = springs.iloc[1,2]    #shift Z
            nodes['dQS'] = ((nodes.X - springs.iloc[2,0])**2 + (nodes.Y - springs.iloc[2,1])**2 + (nodes.Z - springs.iloc[2,2])**2)**0.5
            nodes.type[nodes.dQS==np.min(nodes.dQS)] = 'outfall'
            nodes.X[nodes.dQS==np.min(nodes.dQS)] = springs.iloc[2,0]    #shift X 
            nodes.Y[nodes.dQS==np.min(nodes.dQS)] = springs.iloc[2,1]    #shift Y
            nodes.Z[nodes.dQS==np.min(nodes.dQS)] = springs.iloc[2,2]    #shift Z
            nodes.drop(labels=['dQE','dQA','dQS'], axis='columns', inplace=True)    #drop  unneeded columns    
    except:
        print('no outfalls assigned')
        
    #Remove duplicates:
    #NOTE: Currently disabled - not needed unless network has multiple links between the same two nodes
    run=False
    if run==True:
        links['duplicate'] = False
        links['pairs'] = [[links.loc[i].InNode,links.loc[i].OutNode] for i in links.index] #add a column with both in and out nodes for each link as a list
        for i in links.index:
            pairFlip = [links.loc[i].OutNode,links.loc[i].InNode]          #reverse the in and out nodes
            if pairFlip in list(links.pairs.values):                       #check if there is a duplicate conduit between these two nodes
                #print('duplicate flipped link identified', links.loc[i].pairs)
                links.at[i,'duplicate'] = True
                if links.loc[i].InNode > links.loc[i].OutNode:             #if duplicate, remove the link that starts with a higher node number
                    #print('removing flipped link ', links.loc[i].pairs)
                    links.drop(index=i, inplace=True)

        links.drop(labels=['duplicate','pairs'], axis='columns', inplace=True)              #drop  unneeded columns    
        links.reset_index(drop=True,inplace=True)                                           #reset indexing to start from zero
        links['Name'] = links.index                                                         #update link name (use row index)
    
    return nodes, links


##################################################################################

def import_subcatchments(data, nodes=None, find_widths=False, mask_array=None):
    
    '''Import subcatchment information and return inputs for set_subcatchments. 
    Either takes preset polygons from a shapefile and converts to a geopandas polygon object, 
    or calculates subcatchments automatically for each node from a DEM raster.
    
    Inputs:
    data:        shapefile of polygons defining each subwatershed OR the path to a DEM raster of model extent.
    nodes:       if using DEM, either a pandas dataframe of nodes or the path to a nodes txt file, with columns: X,Y,Z
    find_widths: if using DEM, whether or not to calculate the subcatchment widths. If False, a simple width calculation can be done in set_subcatchments()
    mask_array:  if using DEM, provide an optional mask array with nan values in the cells to be masked, and the final subcatchment array will have those cells masked
    
    Outputs:
    either
    polygons:    geopandas polygon object outlining the subcatchments
    
    or
    catchs:      array where the value in each cell indicates the id number of the subcatchment that cell belongs to (ie which node it drains to)
    areas:       list of surface areas of each subcatchment
    widths:      list of widths of each catchment (see SWMM doc for what width parameter means)
    '''
    
    #For pre-defined polygons:
    try:
        polygons      = gpd.read_file(data)                        #import subcatchment polygon shapefile
        print('loading polygon shapefile')
        polygons.sort_values(by='typ',inplace=True)                     #sort by subcatchment number
        return polygons                                                 #return geopandas polygon objects

    #For auto-calculating subwatersheds:
    except:
        print('assigning subcatchments to nodes automatically using pysheds')
        
        #Get & format data:
        try:     #if importing a txt file
            nodes = pd.read_csv(nodes, sep='\t', names=['X','Y','Z'])      #load node point coordinates from txt into a pandas df
        except:  #if already a dataframe
            nodes = nodes
        grid = Grid.from_raster(data, data_name='dem')       #load DEM raster to grid (if data_name='dem', returns an object called grid.dem)
        dem,dema,xmin,xmax,xres,dx,ymin,ymax,yres,dy,zmin,zmax,zres = mapping.importDEM(data,show=False) #get dem spatial info
        xy = np.column_stack((nodes.X.values,nodes.Y.values))       #put X and Y values into a numpy array of dim (n,2)
        dirmap = (64, 128, 1, 2, 4, 8, 16, 32) #specify what value to assign to each compass direction (N,NE,E,SE,S,SW,W,NW). Try (64, 128, 1, 2, 4, 8, 16, 32).
        
        #Clean data:
        grid.fill_depressions('dem',        out_name='demDepFill')  #fill depressions in DEM (smooths out local depressions)
        grid.fill_pits(       'demDepFill', out_name='demPitFill')  #fill pits in DEM (smooths point depressions)
        grid.resolve_flats(   'demPitFill', out_name='demClean')    #give a slight slope to flat regions

        #Find flow directions, accumulation, and pour points:
        grid.flowdir(     data='demClean', out_name='dir', dirmap=dirmap)  #calculate flow direction raster (what direction will water flow at each cell?)     
        grid.accumulation(data='dir',      out_name='acc', dirmap=dirmap)  #calculate flow accumulation (which cells will water flow towards?)     
        
        threshold = 0.005 * np.nanmax(grid.acc)     #set snapping threshold - only cells with acc values higher than this will be considered (20-50 is good)
        #print('threshold:',threshold)
        accThresh = grid.acc.copy()                 #make a copy to preserve original
        accThresh[accThresh<=threshold] = 0         #make all cells with value <= threshold ineligible for snapping
        xyPour = grid.snap_to_mask(accThresh, xy, return_dist=False)     #snap pour points to the high-accumulation cells closest to the desired nodes
        if find_widths==True:
            cell_dist = grid.cell_distances(grid.dir, inplace=False, as_crs=pyproj.Proj('+proj=utm +zone=32')) #get array of distances from each cell to downstream neighbor using local projection (defaults to north zone)

        #Debugging:
        #f2 = plt.figure(figsize=(20,20)) 
        #plt.imshow(grid.acc,extent=[xmin,xmax,ymin,ymax])
        #plt.scatter(nodes.X,nodes.Y,c='k',s=10)
        
        prevent_duplicates=True
        if prevent_duplicates == True:
            print('resolving duplicates')
            #If any nodes get snapped to the same pour point, replace one of the duplicates with the original point location:
            count = 0 #iniialize counter to track how many while loops are happening
            while len(set(map(tuple, xyPour))) < len(xyPour):  #check for duplicates (if length of set is smaller, duplicates are present)
                count = count + 1  #add one to counter for every run-through
                #print('duplicate pour points present, loop', count)
                points = [[0, 0]]               #create empty list to store points checked
                threshold = 0.5*threshold       #lower the threshold by half
                for i in range(len(xyPour)):    #loop over snapped points
                    pt = xyPour[i]              #get current snapped point
                    #plt.scatter(pt[0],pt[1], s=30, c='none', edgecolors='b')
                    if (pt==points).all(1).any():  #if current pt already in the list (meaning another pt has already been snapped to same coordinates)
                        #print('pour point', pt, 'already snapped to!')
                        if threshold > 1:           #if threshold value is still large enough for a decrease to matter
                            #print('decreasing threshold to', threshold)
                            accThreshLow = grid.acc.copy()                  #make a copy to preserve original
                            accThreshLow[accThreshLow<=threshold] = 0       #set all cells with value <= threshold equal to zero 
                            node = np.asarray([nodes.X.iloc[i], nodes.Y.iloc[i]])   #get XY coordinates of original unsnapped node
                            newpt = grid.snap_to_mask(accThreshLow, node, return_dist=False)    #snap pour pt to high-acc. cell closest to original node
                            points.append(newpt)       #add new point to list of snapped points
                            #plt.scatter(newpt[0],newpt[1], c='none', edgecolors='c', s=50)  
                        if threshold <=1:                   #if threshold can't meaningfully be decreased more
                            #print('threshold cannot be decreased more')
                            node = np.asarray([nodes.X.iloc[i], nodes.Y.iloc[i]])  #get XY coordinates of original unsnapped node
                            pointslist = [list(point) for point in points]         #convert list of arrays to list of lists to enable truth testing
                            if (node==pointslist).all(1).any():                    #if original node is already in the list
                                if count<20:   #if haven't looped too many times
                                    #print('original point', node, 'already in list! shifting duplicate southeast by 1') #print notification message
                                    points.append([nodes.X.iloc[i] + 1,nodes.Y.iloc[i] - 1])     #append shifted point 
                                    #plt.scatter(nodes.X, nodes.Y-0.5, c='none', edgecolors='r', s=70)
                                else:          #if looping too many times, shift a different amount
                                    #print('original point', node, 'already in list! shifting duplicate southeast by 1.5') #print notification message
                                    points.append([nodes.X.iloc[i] + 1.5,nodes.Y.iloc[i] - 1.5])     #append shifted point 
                                    #plt.scatter(nodes.X, nodes.Y-0.5, c='none', edgecolors='r', s=70)
                            else: 
                                #print('replacing duplicate with original node:', [nodes.X.iloc[i], nodes.Y.iloc[i]] ) #print notification message
                                points.append([nodes.X.iloc[i], nodes.Y.iloc[i]])  #append original point value to the list
                                #plt.scatter(nodes.X, nodes.Y, c='none', edgecolors='g', s=70)
                    else:
                        #print('point not already in list - appending')
                        points.append(pt)       #if point isn't a duplicate, append to list of snapped points
                xyPour = np.asarray(points[1:]) #convert list to np array (removing initial 0,0 point)                
                #print('set size',len(set(map(tuple, xyPour))), 'list size', len(xyPour)) #debugging
        
        #Debugging:
        #print('xyPour', len(xyPour), xyPour)
        #plt.scatter(xyPour[:,0],xyPour[:,1], c='none', edgecolors='b', s=30) #debugging
        #for i, pt in enumerate(xyPour):
            #plt.annotate(str(i), pt, color='white')
        
        # Delineate subcatchments
        #f3 = plt.figure(figsize=(20,50))            #debugging
        labels = np.zeros(grid.shape, dtype=int)    #create empty array of same shape as grid, to store subcatchment label numbers
        catchments = []                             #create empty catchment list to store subcatchment arrays
        dists = []                                  #create empty distances list to store flow distance arrays
        for i in range(len(xyPour)):                #loop over pour points
            x, y = xyPour[i,:]                      #get x and y coord
            try:
                catchment = grid.catchment(x=x, y=y, data='dir', xytype='label', inplace=False) #calculate catchment for current point
            except:
                catchment = np.zeros(grid.shape)    #if the calculation fails (e.g. for a skipped point with xy=(nan.nan)), the catchment array is all zeros 
            catchments.append(catchment)            #append to list of catchment arrays
            mask = (catchment != 0)                 #create a True/False array where True cells are in the current subcatchment
            labels = labels + mask                  #add 1 to each cell in labels if the value of that cell in the mask = True
            if find_widths==True:
                dist = grid.flow_distance(x, y, data='dir', weights=cell_dist, dirmap=dirmap, out_name='dist', xytype='label', inplace=False) #get flow dist from each cell in the subcatchment to the outlet of that subcatchment (in cells)
                dists.append(dist)                  #append distance array to list of distance arrays
            #Debugging
            #ax = plt.subplot((len(xyPour)//3) +1,3,i+1)
            #ax.imshow(labels, cmap='viridis') 
        #print('labels',np.unique(labels[~np.isnan(labels)])) #debugging
        
        # Find non-overlapping catchments
        for index, catchment in enumerate(catchments):      #loop over catchments (and get the index by using enumerate)
            mask = (catchment > 0)                          #create a True/False array where True cells are in the current catchment            
            try:
                label = labels[mask].min()                      #get the minimum label value for the cells in the current catchment
            except: 
                label = 0                                       #if no cells are in the catchment, the label value = 0
            catchments[index] = np.where(mask & (labels==label), catchment, 0) #assign value of catchment cells to cells with min label value, otherwise assign zero (this removes upstream catchments?)
            if find_widths==True:
                dists[index] = np.where(mask & (labels==label), dists[index], np.nan)  #do the same process for the flow distance arrays
        catchs = sum([np.where(catchment, index + 1, 0)   #assign index label to each subcatchment & zero to cells outside, then sum each cell
                                  for index, catchment in enumerate(catchments)]).astype(float)
        catchs[catchs == 0] = np.nan          #replace zeros with nans
        catchs = catchs - 1                   #convert catchment labels to 0-based indices
        #print('catchs',len(np.unique(catchs[~np.isnan(catchs)])), np.unique(catchs[~np.isnan(catchs)]))  #debugging
        
        remove_mismatches=False
        if remove_mismatches==True:
            #Remove catchments that don't correspond to a node (duplicates)
            #This should be obsolete now that duplicate pour points are resolved
            catchs_notin = np.isin(catchs, nodes.Name, invert=True) #bool array of where catch ID doesn't match a node ID
            catchs[catchs_notin] = np.nan  #convert cells where catchment ID doesn't match a node ID to nan
        
        #Mask cells outside desired area:
        try:
            catchs[np.isnan(mask_array)] = np.nan     #assign nan to cells that have nan value in the mask array
        except:
            pass
        
        #plt.imshow(catchs, extent = [xmin, xmax, ymin, ymax], cmap='viridis') #debugging
        
        #Calculate areas and widths:
        catchIDs, counts = np.unique(catchs[~np.isnan(catchs)], return_counts=True) #get list of unique non-NaN values in array, & count times each occurs
        catchIDs = np.array(catchIDs, dtype=int)        #convert to int
        areas = np.full(len(nodes), np.nan)             #create nan array of correct length to preserve indexing if nodes with no catchment present 
        counts = dict(zip(catchIDs, counts))            #convert to dictionary to preserve indexing
        for node in nodes.index.values:                 #loop over nodes
            try:
                count = counts[node]                    #get count for that node's catchment 
            except:
                count = 0                               #if node has no catchment, assign zero
            areas[node] = count*dx*dy                   #calculate subcatchment areas: ncells * cell width * cell height
        if find_widths==True:
            longestPath = [np.nanmax(dist) for dist in dists] #get list of longest overland flow path for each subcatchment (calculated in previous section)
            avgPath = [np.nanmean(dist) for dist in dists]    #get list of mean flow path length for each subcatch
            widths = [np.round(areas[i] / avgPath[i], 2) for i in range(len(areas))] #method 2: SWMM user manual suggests area/avg flow path length
        areas = areas / 10000                                 #convert from sq.m to hectares for SWMM input file (WHYYY??? this is a ridiculous unit!)

        #Debugging
        #if len(np.unique(catchs[~np.isnan(catchs)])) != len(nodes):    #check that each node has a catchment
               #print('warning:', len(nodes), 'nodes, but', len(np.unique(catchs[~np.isnan(catchs)])), 'subcatchments')

        #Return outputs:
        if find_widths==True:
            return catchs, areas, widths
        else:
            return catchs, areas
        

#########################################################################################

def get_order(nodes, links, areas, alpha=0.1, beta=0.001):
    '''Calculate conduit order and diameters by converting network to a directed graph.
    Use equations provided in Borghi et al. 2016 Section 4.1:
    conduit order = area drained by conduit / total catchment area
    conduit diameter = 2 * alpha * e^(order * beta), where alpha and beta are parameters to be optimized. 
    
    Inputs:
    nodes: pandas dataframe of nodes with columns: X, Y, Z, type, Name   (returned by swmmpy.sks2swmm())
    links: pandas dataframe of links with columns: InNode, OutNode, Name (returned by swmmpy.sks2swmm())
    areas: list of subcatchment areas (returned by swwmpy.import_subcatchments())
    alpha: Adjustment parameter. In Borghi et al. 2016, alpha = 0.2. For Gottesacker, alpha = 0.1
    beta:  Adjustment parameter. In Borghi et al. 2016, beta = 1.5.  For Gottesacker, beta = 0.001
    
    Outputs:
    order:     list of conduit orders (indexed by name). Higher order=drains larger area=lower elevation=closer to springs
    diameters: list of conduit diameters (indexed by name)
    '''
    
    #Create directed graph from links df:
    G = nx.from_pandas_edgelist(links, source='InNode', target='OutNode', edge_attr='Name', create_using=nx.DiGraph())

    #Calculate cumulative area drained by each node (see Borghi et al. 2012, Figure 9):
    for node in G.nodes:
        ancestors = list(nx.algorithms.dag.ancestors(G, node)) #get list of all upstream nodes
        ancestors = np.asarray(ancestors, dtype=int)           #convert to array for use as indexer
        area = sum(areas[ancestors])                           #sum area values for all upstream nodes
        G.nodes[node]['area'] = area + areas[node]             #add subcatch area of current node
        #print(node, ancestors, G.nodes[node]['area'])

    #Convert to conduit diameter (see Borghi et al. 2016 section 4.1)
    #conduit order = area drained by conduit / total catchment area
    #conduit diameter = 2* alpha * e^(order * beta)
    diameters = pd.DataFrame(index=links.Name, columns=['diameter'] )  #set up empty list to store diameters
    order = pd.DataFrame(index=links.Name, columns=['order'] )  #set up empty list to store diameters
    #print('total area', sum(areas))
    
    for link in G.edges:                  #loop over links
        #print(G.edges[link]['Name'])
        fromNode = link[0]                #get from node for each link
        G.edges[link]['area'] = G.nodes[fromNode]['area']  #assign that link the area of the from node
        G.edges[link]['order'] = G.nodes[fromNode]['area'] / sum(areas)  #order is ratio of link area to total catchment area
        order.loc[G.edges[link]['Name'], 'order'] = G.edges[link]['order']       #store in list for use later    
        G.edges[link]['diameter'] = 2*alpha*np.exp(beta*G.nodes[fromNode]['area']) #diameter is based on exponential relationship to order
        diameters.loc[G.edges[link]['Name'], 'diameter'] = G.edges[link]['diameter']       #store in list for use later
        #print('area', G.edges[link]['area'], 'order', G.edges[link]['order'], 'diameter',G.edges[link]['diameter'])
    diameters = diameters.diameter.values  #store as list
    order = np.asarray(order.order.values, dtype='float')   #store as array of floats (defaults to weird 'o dtype')
    
    return order, diameters


#########################################################################################
#OLD VERSION
def old_get_timeseries(rainfallOld, baseflowOld, scaNew, scaOld, timestrings):
    '''Recalculates Rainfall (units: mm) and Baseflow (L/s?) timeseries for automatically-assigned subcatchments, based on original subcatchments assigned by Zhao.
       This will need to be re-written if the original rainfall data is distributed differently.
       
       rainfallOld: df of original rainfall timeseries, with strings for date and time, and original subcatchment indices (units: mm)
                    with columns: Name,      Date,       Time,     Value, datetime,            type,     catch
                                  Rainfall1, 11/15/2013, 01:00:00, 0.0,   11/15/2013 01:00:00, Rainfall, 1.0 
       baseflowOld: df of original baseflow timeseries, with same columns as rainfall (units: L/s?)
       scaNew:      array of new subcatchment indices, of dim (nrow,ncol), with 0-based indexing
       scaOld:      array of old subcatchment indices, of dim (nrow,ncol), with 1-based indexing
       timestrings: list of strings indicating the timesteps, will be same for all timeseries
       
       Returns:
       Pnewcatchs: df of recalculated Rainfall tseries for new subcatchments
       Bnewcatchs: df of recalculated Baseflow tseries for new subcatchments'''

    t = pd.to_datetime(timestrings)                                 #convert timestrings to datetime indices
    catchIDs, counts  = np.unique(scaNew[~np.isnan(scaNew)], return_counts=True) #get list of unique non-NaN values in array, and number of times each value occurs
    newcatchlist = np.array(catchIDs+1, dtype=int)                               #convert to int and to 1-based indexing 
    Pnewcatchs = pd.DataFrame(np.zeros((len(t),len(newcatchlist))), index=t, columns=newcatchlist)   #make an empty df for new rainfall tseries
    Bnewcatchs = pd.DataFrame(np.zeros((len(t),len(newcatchlist))), index=t, columns=newcatchlist)   #make an empty df for new baseflow tseries

    for nsc in newcatchlist:      #loop over new catch indices
        #print(nsc)
        ncells = len(scaNew[scaNew==nsc])                   #count number of cells in current (new) catch
        oldcatchlist = np.unique(scaOld[scaNew==nsc])       #get list of old catch indices in current new catch
        oldcatchlist = oldcatchlist[oldcatchlist!=0]        #remove zero values (zeros are null)
        Pnew = pd.DataFrame(np.zeros(len(t)), index=np.arange(len(t)), columns=['Psum'])    #make a new dataframe of zeros with length (ntimesteps)
        Bnew = pd.DataFrame(np.zeros(len(t)), index=np.arange(len(t)), columns=['Bsum'])    #make a new dataframe of zeros with length (ntimesteps)
        for oldsc in oldcatchlist:                        #loop over old catch indices
            noldcells = np.sum(scaOld[scaNew==nsc]==oldsc)   #count number of cells in new catchment from each old catchment
            #print(oldsc,noldcells)
            Pold = rainfallOld[rainfallOld.catch==oldsc]  #get rainfall timeseries for old catch
            Bold = baseflowOld[baseflowOld.catch==oldsc]  #get baseflow timeseries for old catch
            Pold['flowxncells'] = Pold.Value * noldcells  #multiply flow value by number of cells from old catchment (to weight)
            Bold['flowxncells'] = Bold.Value * noldcells
            #P['flowWeighted'][P.catch==oldsc] = P.flow[P.catch==oldsc] * noldcells
            #Pnew['Psum'] = P.flowWeighted[P.catch==oldsc].values + Pnew.Psum.values
            Pnew['Psum'] = Pold.flowxncells.values + Pnew.Psum.values   #add to new df as sum of previous weighted flow & current weighted flow
            Bnew['Bsum'] = Bold.flowxncells.values + Bnew.Bsum.values   #add to new df as sum of previous weighted flow & current weighted flow
        Pnew['Pweighted'] = Pnew.Psum / ncells                          #divide flow values by ncells in new subcatchment to unweight
        Bnew['Bweighted'] = Bnew.Bsum / ncells
        Pnewcatchs[nsc] = Pnew.Pweighted.values
        Bnewcatchs[nsc] = Bnew.Bweighted.values
        
    return Pnewcatchs, Bnewcatchs

####################################################################################

def get_timeseries(rainfallgrid, baseflowgrid, grid, nodes, catchs, datetimes, dim):
    '''Recalculates Rainfall (units: mm/m2) and Baseflow (L/s?) timeseries for automatically-assigned subcatchments.
       
       rainfallgrid:  dataframe of rainfall (fast-flow) timeseries, one per grid cell (100mx100m cells)  
                      index =  date/time and columns = cell number (units: mm/m2)
       baseflowgrid:  dataframe of baseflow (slow-flow) timeseries, one per grid cell (100mx100m cells)
                      index =  date/time and columns = cell number (units: mm/m2)(units: L/s?)
       grid:          dataframe of cell IDs and XY coordinates
       catchs:        array of 50x50m cells with the value of each cell indicating its subcatchment number
                      returned by pysheds
       datetimes:     list of datetime indices, will be same for all timeseries
       dim:           list of model dimensions: [xres,yres,zres,xmin,ymin,zmin,xmax,ymax,zmax]
       
       Returns:
       rainfall: df of recalculated Rainfall tseries for new subcatchments
       baseflow: df of recalculated Baseflow tseries for new subcatchments'''
        
    #Convert raster array of subcatchments to list of polygons:
    grid[grid==-9999] = np.nan              #replace values outside of catchment with NaN
    catchs = catchs.astype('float32')       #convert catchs to float32 so that rasterio can use it      
    polygons = np.full(len(nodes), np.nan)  #set up array of nans of correct length (in case some nodes don't have catchments) to store polygons
    polygons = list(polygons)               #convert to list to be able to store objects in it
    #print('polygons', len(polygons), polygons)
    for node in nodes.index.values:         #loop over node IDs
        catch = int(node)                   #look at catchment for each node
        #print(catch)
        shapes = rasterio.features.shapes(catchs)  #convert subcatchment array to shapes
        polygon = [shg.Polygon(shape[0]["coordinates"][0]) for shape in shapes if shape[1] == catch] #convert outline of each subcatchment to polygon
        polygons[catch] = polygon           #insert current subcatchment polygon at correct location in array of all subcatchment polygons
        
        #Debugging:
        #print(polygon)
        #for poly in polygon:               #loop over polygons in list 
            #x,y = poly.exterior.xy         #get coordinates of outline
            #plt.plot(x,y)                  #plot outline
            #polygons.append(poly)          #add current subcatchment polygon to list of all subcatchment polygons
    #print('all polygons', len(polygons), polygons)
        
    #Identify which points in the inflow grid are in each subcatchment (SLOW):
    #assuming grid XY coordinates represent center of cell and converting to cell indices
    grid.minorCatch = np.nan                        #reset the minor catchments to all be NaNs
    rows,cols,lays = mapping.xyz2rowcollay(grid.x, grid.y, np.zeros(len(grid)), dim, flip=False) #convert grid points to row, lay, col
    for node in nodes.index.values:                 #loop over nodes
        catch = int(node)
        #print(catch, polygons[catch])
        for i in grid.index:                        #loop over points in grid
            #print(i)
            point = shg.Point(cols[i], rows[i])     #create shapely point object
            for poly in polygons[int(catch)]:       #loop over polygons in list of polygons describing current catchment
                #print('poly', poly)
                if np.isnan(grid.minorCatch[i]):    #check that point hasn't already been assigned a catchment
                    if point.within(poly):          #check if the point is inside the current catchment polygon
                        grid.minorCatch[i] = catch  #set the catchment for current point to be the current catchment
                    #print(point)
                    #plt.scatter(cols[i],rows[i])
    #plt.gca().invert_yaxis()
    #plt.gca().set_aspect('equal')
    
    #Compile inflow data for the new subcatchments (at 100mx100m resolution):
    rainfall = pd.DataFrame(index=datetimes,columns=[])    #create empty df to store data
    #for catch in np.unique(catchs[~np.isnan(catchs)]):     #loop over catchments
    for node in nodes.index.values:                       #loop over nodes (indicating catchments) to preserve indexing
        catch = int(node)
        cells = grid.pointID[grid.minorCatch==catch]       #select all grid points in minor subcatchment 
        catchRainfall = rainfallgrid[cells]                #get all timeseries in subcatchment  [mm/m2]
        catchRainfall['Avg'] = catchRainfall.sum(axis=1) / len(cells)  #avg all cells in subcatchment to get [mm/catchment]
        rainfall[catch] = catchRainfall.Avg                #store catchment rainfall in new column in new df

    baseflow = pd.DataFrame(index=datetimes,columns=[])
    #for catch in np.unique(catchs[~np.isnan(catchs)]):
    for catch in nodes.index.values:                       #loop over nodes (indicating catchments) to preserve indexing
        cells = grid.pointID[grid.minorCatch==catch]      #select all grid points in minor subcatchment #1
        catchBaseflow = baseflowgrid[cells]               #get all timeseries in subcatchment #1 [mm/m2]
        catchBaseflow['Avg'] = catchBaseflow.sum(axis=1) / len(cells)  #avg all cells in subcatchment #1 to get [mm/catchment]
        baseflow[catch] = catchBaseflow.Avg          #store catchment rainfall in new column in new df

    return rainfall, baseflow, polygons


#########################################################################################

def set_raingages(nodes, raintype='VOLUME', tintvl='1:00', snowctch=1.0, datatype='TIMESERIES'):

    '''Create raingages dataframe from raster subcatchment data.
    Inputs:
    catch:    array or raster with an integer in each cell indicating the ID of the subcatchment that cell belongs to 
    (subcatchment ID number is the same as the ID number of the node that subcatchment drains to) (returnable by import_subcatchments())
    (the raingage number must be the same as the number of the subcatchment it represents)
    *SWMMargs:  can be either a single value or a list of same length as number of junctions. See SWMM doc for details.'''
    
    #catchIDs  = np.unique(catch[~np.isnan(catch)])    #get list of unique non-NaN values in array
    #catchIDs = np.array(catchIDs, dtype=int)          #convert to int   
    
    raingages = pd.DataFrame()                        #create empty dataframe to store info
    #raingages['Name'] = [str(s+1) for s in catchIDs]  #set raingage name from catch name (1-based indexing)
    raingages['Name'] = [str(int(s+1)) for s in nodes.index.values]  #set raingage name from node (1-based indexing)
    raingages['RainType']   = raintype
    raingages['tIntvl']     = tintvl
    raingages['SnowCtch']   = snowctch
    #raingages['DataSource'] = [datatype+' Rainfall'+str(s) for s in catchIDs] #assign a rainfall timeseries to each raingage of the same number (0-based indexing)
    raingages['DataSource'] = [datatype+' Rainfall'+str(int(s)) for s in nodes.index.values] #assign a rainfall timeseries to each raingage of the same number (0-based indexing)
    
    return raingages


#########################################################################################

def set_junctions(nodes, maxdepth=0, initdepth=0, surdepth=200, aponded=0):
    
    '''Create junctions dataframe from node data.
    Inputs:
    nodes:      pandas df with all nodes, and columns: X, Y, Z, type, Name. Type is either 'junction', or 'outfall'
    *SWMMargs:  can be either a single value or a list of same length as number of junctions. See SWMM doc for details.'''
    
    junctions  = nodes[nodes.type!='outfall'].copy()                              #split nodes into only junction-type nodes(not outfalls)
    junctions.drop(labels=['X','Y','type'], axis='columns', inplace=True)   #drop unneeded columns
    junctions.rename({'Z':'InvertElev'},    axis='columns', inplace=True)   #rename to SWMM's column names
    junctions['MaxDepth']  = maxdepth                                              #add required columns
    junctions['InitDepth'] = initdepth
    junctions['SurDepth']  = surdepth
    junctions['Aponded']   = aponded
    colnames = ['Name','InvertElev','MaxDepth','InitDepth','SurDepth','Aponded']    #list of column names in the correct order
    junctions = junctions.reindex(columns=colnames)                                 #reorder the column names
    
    return junctions


#########################################################################################

def set_outfalls(nodes, outtype='FREE', stage='', tidegate='NO'):
    
    '''Create outfalls dataframe from node data.
    Inputs:
    nodes:      pandas df with all nodes, and columns: X, Y, Z, type, Name. Type is either 'junction', or 'outfall'
    *SWMMargs:  can be either a single value or a list of same length as number of junctions. See SWMM doc for details.'''
    
    outfalls = nodes[nodes.type=='outfall'].copy()                          #split nodes into only outfall-type nodes(not junctions)
    outfalls.drop(labels=['X','Y','type'], axis='columns', inplace=True)   #drop unneeded columns
    outfalls.rename({'Z':'InvertElev'},    axis='columns', inplace=True)   #rename to SWMM's column names
    outfalls['OutType']  = 'FREE'                                          #add required columns
    outfalls['Stage']    = ''
    outfalls['TideGate'] = 'NO'
    colnames = ['Name','InvertElev','OutType','Stage','TideGate']           #list of column names in the correct order
    outfalls = outfalls.reindex(columns=colnames)                           #reorder the column names
    
    return outfalls


#########################################################################################

def set_conduits(links, nodes, dim=3, manningN=0.01, inoffset='*', outoffset='*', initflow=0, maxflow=0):

    '''Create conduits dataframe from link and node data.
    Inputs:
    links:      pandas df with all links, and columns: InNode, OutNode, Name
    nodes:      pandas df with all nodes, and columns: X, Y, Z, type, Name. Type is either 'junction', or 'outfall'
    dim:        number of dimensions (2 or 3)
    *SWMMargs:  can be either a single value or a list of same length as number of junctions. See SWMM doc for details.'''
    
    conduits = links.copy()                                      #create a dataframe by copying links
    fromX = nodes.X.loc[conduits.InNode]                            #calculate coordinates for link start and end points
    fromY = nodes.Y.loc[conduits.InNode]
    toX   = nodes.X.loc[conduits.OutNode]
    toY   = nodes.Y.loc[conduits.OutNode]
    if dim==2:
        conduits['Length']    = [((toX.iloc[ind] - fromX.iloc[ind])**2 + (toY.iloc[ind] - fromY.iloc[ind])**2)**0.5 for ind in links.index]  #calculate length using distance formula and add it as a new column
    if dim==3:
        fromZ = nodes.Z.loc[conduits.InNode]
        toZ   = nodes.Z.loc[conduits.OutNode]
        conduits['Length']    = [((toX.iloc[ind] - fromX.iloc[ind])**2 + 
                                  (toY.iloc[ind] - fromY.iloc[ind])**2 + 
                                  (toZ.iloc[ind] - fromZ.iloc[ind])**2)**0.5 for ind in links.index] #calculate distance between points and add as new column
    conduits['ManningN']  = manningN                                 #add other columns
    conduits['InOffset']  = inoffset
    conduits['OutOffset'] = outoffset
    conduits['InitFlow']  = initflow
    conduits['MaxFlow']   = maxflow
    colnames = ['Name','InNode','OutNode','Length','ManningN','InOffset','OutOffset','InitFlow','MaxFlow']  #list of column names in the correct order
    conduits = conduits.reindex(columns=colnames)                                                           #reorder the column names
    
    return conduits


#########################################################################################

def set_coordinates(nodes):

    '''Create coordinates dataframe from node data.
    Inputs:
    nodes:      pandas df with all nodes, and columns: X, Y, Z, type, Name. Type is either 'junction', or 'outfall'
    '''
    
    coordinates = nodes.copy()                                                          #create a dataframe by copying nodes
    coordinates.drop(labels=['Z','type'],              axis='columns', inplace=True)    #drop unneeded columns
    coordinates.rename({'X':'X-Coord', 'Y':'Y-Coord'}, axis='columns', inplace=True)    #rename to SWMM's column names
    colnames = ['Name','X-Coord','Y-Coord']                                             #list of column names in the correct order
    coordinates = coordinates.reindex(columns=colnames)                                 #reorder the column names
    
    return coordinates


#########################################################################################

def set_xsections(links, shape='FORCE_MAIN', g1=3.0, g2=1000.0, g3=0, g4=0, barrels=1, culvert=''):

    '''Create xsections dataframe from link data.
    Inputs:
    links:      pandas df with all links, and columns: InNode, OutNode, Name
    *SWMMargs:  can be either a single value or a list of same length as number of junctions. See SWMM doc for details.'''
    
    xsections = links.copy()
    xsections.drop(columns=['InNode','OutNode'], inplace=True)    #drop unneeded columns
    xsections['Shape']   = 'FORCE_MAIN'    #set shape of conduit (see SWMM doc table D-1) - FORCE_MAIN is a circle
    xsections['Geom1']   = g1             #set diameter in m
    xsections['Geom2']   = g2          #set roughness in mm
    xsections['Geom3']   = g3               #not changing remaining parameters
    xsections['Geom4']   = g4
    xsections['Barrels'] = barrels
    xsections['Culvert'] = culvert
    
    return xsections


#########################################################################################

def set_inflows(nodes, tseriesname='Baseflow', catch=None, par='FLOW', partype='FLOW', unitfactor=1.0, scalefactor=1.0, baseval='', basepattern=''):
    
    '''Create inflows dataframe from node data.
    Inputs:
    nodes:  pandas df with all nodes, and columns: X, Y, Z, type, Name. Type is either 'junction', or 'outfall'
    tseriesname: string identifying the name prefix of the timeseries to use. A node number will be appended.
    catch:  optional array or raster with an integer in each cell indicating the ID of the subcatchment that cell belongs to 
    (subcatchment ID number must be same as node number that subcatch drains to, and as inflow timeseries number) 
    (catch is returnable by import_subcatchments())
    *SWMMargs:  can be either a single value or a list of same length as number of junctions. See SWMM doc for details.
    
    '''
    #catchIDs, counts  = np.unique(catch[~np.isnan(catch)], return_counts=True)  #get list of unique non-NaN values in array, & # of times each value occurs
    #catchIDs = np.array(catchIDs, dtype=int)                                    #convert to int
    inflows = pd.DataFrame()                                                     #create empty dataframe to store info
    #inflows['Node'] = [str(s) for s in catchIDs]                                #assign node name to be same as subcatch name
    inflows['Node'] = [str(s) for s in nodes.index.values]                       #assign node name to be same as subcatch name
    inflows['Parameter']   = par                                                 #add needed columns
    if catch:
        inflows['TimeSeries']  = [tseriesname+str(s) for s in catchIDs]          #assign timeseries name from subcatch (0-based indexing)
    inflows['TimeSeries']  = [tseriesname+str(s) for s in nodes.index.values]    #assign timeseries name from node index
    inflows['ParType']     = partype              
    inflows['UnitFactor']  = unitfactor              
    inflows['ScaleFactor'] = scalefactor            
    inflows['BaseVal']     = baseval              
    inflows['BasePattern'] = basepattern             
    
    return inflows


#########################################################################################

def set_subcatchments(nodes, catch, areas, width='auto', pctimpv=0.25, slope=0.25, curblen=0, snowpack='', plot=False):

    '''Create subcatchments dataframe from node and catchment data.
    Inputs:
    nodes:    pandas df with all nodes, and columns: X, Y, Z, type, Name. Type is either 'junction', or 'outfall'
    catch:    either a geopandas polygons object, or an array or raster with a 0-based integer index in each cell indicating the subcatchment that cell belongs to (either returned by import_subcatchments())
    areas:    a list of the area (in hectares) of each subcatchment (can be returned by import_subcatchments()). If None, defaults to 2*sqrt(area) as per Larry Bodnaruk's OpenSWMM recommendation.
    width:    a list of the representative width of each subcatchment (see SWMM doc for what this parameter means). (can be returned by import_subcatchments())
    *SWMMargs:  can be either a single value or a list of same length as number of junctions. See SWMM doc for details.'''
    
    #Calculate widths if not provided:
    if width=='auto':
        print('calculating width')
        width   = np.round(2*(areas*10000)**.5, 2)                    #note: this is VERY ROUGH. This parameter is squishy. See SWMM doc. 10000 is to convert from ha to m
    
    #For array or raster type subcatchment data:
    try:
        #Each subcatchment name, index, and rain gage are the same
        #catchIDs, counts  = np.unique(catch[~np.isnan(catch)], return_counts=True) #get list of unique non-NaN values, & # of times each value occurs
        #catchIDs = np.array(catchIDs, dtype=int)                         #convert to int (still 0-based indexing)      
        outlets = pd.DataFrame()                                         #create empty dataframe to store outlet info
        outlets['subcatchment'] = nodes.index.values+1                   #convert to 1-based & store (doesn't matter for rasters, but compatible with polygons)
        outlets['Name']         = nodes.index.values                     #store 0-based node name that each subcatch drains to 
        #outlets['subcatchment'] = catchIDs+1                             #convert to 1-based & store (doesn't matter for rasters, but compatible with polygons)
        #outlets['Name']         = catchIDs                               #store 0-based node name that each subcatch drains to (see previous line)
        
        
    #For polygon-type subcatchment data:
    except:
        
        #Create junction points using shapely, to use for spatial analysis:
        #jnodes = nodes[nodes.type=='junction']                              #split nodes into only junction-type nodes(not outfalls)
        jnodes=nodes                                                        #if outfalls can receive inflow, then all nodes are available
        points = [shg.Point(jnodes.loc[i].X, jnodes.loc[i].Y, jnodes.loc[i].Z) for i in jnodes.index]    #create list of point objects by looping over points in df
        jpts = jnodes                                                                     #create df for junction points starting from junction data
        jpts['geometry'] = points                                                            #add a column with the points
        
        #Loop over subcatchments to identify draining node:
        outlets = pd.DataFrame()                                    #create an empty dataframe to store outlet info
        nsc  = len(catch.typ)                                    #get number of subcatchments
        if plot==True:
            f,ax = plt.subplots(nsc,1, figsize=(10,nsc*5))          #create figure and axis objects

        for i in range(nsc):                                        #loop over subcatchments
            sc = catch.typ.iloc[i]                               #get subcatchment index (not the same as i)
            #print('calculating for subcatchment ', sc)
            #Get current polygon
            loc = catch['typ']==sc                               #get location index of current subcatchment
            S = catch.loc[loc]                                   #get current subcatchment
            poly = S.iloc[0].geometry                               #get polygon info for desired polygon

            #Check which points are in selected polygon:
            pip = []                                                #create empty list to store points in polygon
            for j in jpts.index:                                  #loop over point indices
                if jpts.loc[j].geometry.within(poly):            #if the point is in the polygon
                    pip.append(jpts.loc[j])                      #add it to the list
            pip = pd.DataFrame(pip)                                 #convert to dataframe for easier manipulation

            #Plot:
            if plot==True:
                catch.plot(ax=ax[i], facecolor='gray')           #plot the polygons
                S.plot(ax=ax[i])                                    #plot selected sub-catchment
                pts = [ax[i].scatter((p.x),(p.y),color='k', s=8) for p in jpts.geometry]          #plot all points
                if len(pip)>0:                                                                      #if there are points in the polygon
                    pts = [ax[i].scatter((p.x),(p.y),color='r', s=8) for p in pip.geometry]         #plot only points in polygon
                ax[i].set_title(str(sc))

            #Adjust if several points in one polygon, or if none:
            if len(pip)==1:                    #if polygon holds exactly one point
                outlet = pip                   #that point is the outlet
            if len(pip)>1:                     #if polygon holds more than one point
                lowp = pip[pip.Z==min(pip.Z)]  #select lowest-elevation point from points in polygon
                outlet = lowp                  #the low point is the outlet
                if plot==True:                
                    ax[i].scatter(lowp.X,lowp.Y, c='r', s=30, marker='X') 
            if len(pip)==0:                    #if polygon holds no points
                c = poly.centroid              #get centroid of polygon
                mpo = shg.MultiPoint(jpts.geometry)  #create a multi-point object for all points under consideration
                nearest_geometry = sho.nearest_points(c,mpo)[1]  #calculate the nearest point (note that nearest_points returns two objects, the centroid and the actual point)
                nearest = jpts[jpts.X==nearest_geometry.x]       #get that point from the detailed points dataframe (by checking the x and y coords)
                nearest = nearest[nearest.Y==nearest_geometry.y]
                outlet = nearest                                 #the nearest point is the outlet
                if plot==True:
                    ax[i].scatter(nearest.X,nearest.Y, c='r', s=30, marker='X')
            #print('number of points in subcatchment: ', len(pip))
            outlet['subcatchment'] = sc                         #add column with subcatchment index
            outlets = outlets.append(outlet)                    #store outlet info in a dataframe
            
    outlets.drop_duplicates(subset=['subcatchment'], keep='first',inplace=True)  #remove duplicate subcatchmets

    subcatchments = pd.DataFrame()                          #create new empty dataframe
    subcatchments['Name']     = ['S'+str(s) for s in outlets.subcatchment] #assign subcatchments (add an S to name to distinguish from nodes)
    subcatchments['Raingage'] = outlets.subcatchment.values #assign raingages (rain gage number = catchment number)
    subcatchments['Outlet']   = outlets.Name.values         #assign newly-calculated outlets
    subcatchments['TotArea']  = areas                       #can this be auto-calculated from polygons?
    subcatchments['PctImpv']  = pctimpv
    subcatchments['Width']    = width                       
    subcatchments['PctSlope'] = slope
    subcatchments['CurbLen']  = curblen
    subcatchments['SnowPack'] = snowpack
    
    return subcatchments


#########################################################################################

def set_subareas(nodes, nimperv=0.2, nperv=0.1, simperv=0.15, sperv=0.05, pctzero=100, routeto='IMPERVIOUS', pctrouted=100):

    '''Create subareas dataframe from subcatchment data.
    Note: polygon catchment data not currently supported - to-do
    Inputs:
    catch:    array or raster with an integer in each cell indicating the ID of the subcatchment that cell belongs to (subcatchment ID number is the same as the ID number of the node that subcatchment drains to) (returnable by import_subcatchments())
    *SWMMargs:  can be either a single value or a list of same length as number of junctions. See SWMM doc for details.'''
    
    #catchIDs, counts  = np.unique(catch[~np.isnan(catch)], return_counts=True) #get list of unique non-NaN values in array, and number of times each value occurs
    #catchIDs = np.array(catchIDs, dtype=int)                                   #convert to int      
    
    subareas = pd.DataFrame()                                                  #create empty dataframe to store info
    #subareas['Subcatchment'] = ['S'+str(s+1) for s in catchIDs] #assign subcatchs (add S to name to distinguish from nodes)
    subareas['Subcatchment'] = ['S'+str(int(s+1)) for s in nodes.index.values] #assign subcatchs (add S to name to distinguish from nodes)
    subareas['nImperv']      = nimperv
    subareas['nPerv']        = nperv
    subareas['sImperv']      = simperv
    subareas['sPerv']        = sperv
    subareas['PctZero']      = pctzero
    subareas['RouteTo']      = routeto
    subareas['PctRouted']    = pctrouted
    
    return subareas

    
#########################################################################################
def set_infiltration(nodes, maxrate=3, minrate=0.5, decay=4, drytime=7, maxinfil=0):

    '''Create infiltration dataframe from subcatchment data.
    Note: polygon catchment data not currently supported - to-do
    Inputs:
    catch:    array or raster with an integer in each cell indicating the ID of the subcatchment that cell belongs to (subcatchment ID number is the same as the ID number of the node that subcatchment drains to) (returnable by import_subcatchments())
    *SWMMargs:  can be either a single value or a list of same length as number of junctions. See SWMM doc for details.'''
    
    #catchIDs, counts  = np.unique(catch[~np.isnan(catch)], return_counts=True) #get list of unique non-NaN values in array, and number of times each value occurs
    #catchIDs = np.array(catchIDs, dtype=int)                                   #convert to int      
    infiltration = pd.DataFrame()                                                  #create empty dataframe to store info
    #infiltration['Subcatchment'] = ['S'+str(s+1) for s in catchIDs]                #assign subcatchments (add an S to name to distinguish from nodes)
    infiltration['Subcatchment'] = ['S'+str(int(s+1)) for s in nodes.index.values] #assign subcatchments (add S to name to distinguish from nodes)
    infiltration['MaxRate']      = maxrate
    infiltration['MinRate']      = minrate
    infiltration['Decay']        = decay
    infiltration['DryTime']      = drytime
    infiltration['MaxInfil']     = maxinfil
    
    return infiltration


#########################################################################################

def set_timeseries(nodes, filenames, seriesnames):
    '''
    Import and format timeseries data so that it can be used in the SWMM input file.
    Concatenates several timeseries into one large dataframe.
    Each timeseries must correspond to the inflows to a specific inlet node.
    nodes:      nodes dataframe
    filenames:  list of strings indicating the paths to the data files with inflow data. 
                data files must have three columns with the date, time, and inflow value in volume/time
                filenames must all have the format "path"+"node number"+"file extension"
                ex: inputs/inflow5.txt
    seriesnames: list of strings indicating just the name of each timeseries to match those used in inflows section
                 must be in same order as filenames
      '''
    
    timeseries = pd.DataFrame(columns=['Name','Date','Time','Value'])
    for i,file in enumerate(filenames):
        data = pd.read_csv(file, delim_whitespace=True, header=None, names=['Date','Time','Value'])
        data.insert(0, 'Name', seriesnames[i])
        timeseries = pd.concat([timeseries,data], ignore_index=True)
        
    return timeseries


#########################################################################################

def set_report(nodes2report, links2report):

    '''Create report dataframe.
    Inputs:
    nodes2report: list of node name strings to save data for
    links2report: list of link name strings to save data for
    
    Note: currently the options for which item types to report can't be changed.
    Add this in the future.
    '''
    
    report = pd.DataFrame()                                                            #create empty dataframe
    report['Type']   = ['INPUT','CONTROLS','SUBCATCHMENTS','NODES',      'LINKS']      #add column for item types to report (see SWMM doc for options)
    report['Option'] = ['NO',    'NO',     'NONE',         nodes2report, links2report] #set options for each item type being reported
    
    return report


#########################################################################################

def insert_data(template, placeholder, data, show=False):
    
    '''Inserts data from a pandas dataframe into a template dataframe by replacing the section placeholder in the template.
    The template file and data file must be read in first. Best to make a copy of the template.
    
    Inputs:
    template:     pandas dataframe with one column, and one line per line of text in the template input file.
                  NOTE: each section that will be edited MUST have a section placeholder string at the desired insert location, 
                  and a section header string with dashes for the column spacing, like so: ;;---- ------ ------
    placeholder:  string corresponding to the placeholder string to be replaced
    data:         pandas dataframe with rows of data split into columns, to be inserted at the placeholder location
                  NOTE: the first column will get converted to an integer
    show:         True/False. If True, prints the updated section
    
    Outputs:
    new:          pandas dataframe with one column, and one line per line of text in the desired output file, 
                  with the new data inserted where the placeholder was in the original file
    '''
    
    #Get row index and location of placeholder:
    ind = template.index[template[0]==placeholder]              #df.index gets the index, df[0] looks in column 0 (returns an index object not just the name)
    print(ind, placeholder)
    loc = template.index.get_loc(ind[0])                        #get the integer position of the index specified (need to select first item (index name) in index object)

    #Special fast method for timeseries section:
    if placeholder=='timeseries':
        data = data.round(6).astype(str) + '\t'                 #round floats, convert all values to strings, and add a tab separator
        data_strings = pd.DataFrame(data.values.sum(axis=1))    #concatenate all columns into one and save as new df with one string per line
    
    #Slower pretty formatting for other sections:
    else:    
        #Convert data to formatted line strings:
        #Create format string based on template file:
        fmt_line = template[0].loc[loc-1]                           #get string to model formatting on (the row of dashes in the template file)
        cols = fmt_line.split()                                     #split string into chunks separated by whitespace
        l = [len(item) for item in cols]                            #get length of each column
        form = ''                                                   #create empty string to fill
        form = [form + '{d['+ str(i) + ']:<' + str(l[i]+1) + '}' for i in range(len(l))]  #concatenate formatting info and return a list of format strings (one per column)
        form = ''.join(form)                                        #join format strings into one for entire row

        #Insert values into format string
        #print('\n', 'data','\n', data)                     #optional print to check if correct
        data = data.round(5)                                #round data to 6 decimal places (to correct for artifacts caused by binary/float conversions)
        data_strings = pd.DataFrame(columns=[0])            #create empty dataframe to fill
        for ind in data.index:                              #loop over lines of data
            dl = data.loc[ind].tolist()                     #get line of data to be formatted
            if type(dl[0])!=type('str'):                    #if name isn't a string
                dl[0] = int(dl[0])                          #make sure name column is an integer not a float 
            line = form.format(d=dl)                        #insert each item in list into format string
            data_strings.loc[ind] = line                    #insert line string into new dataframe

    #Replace placeholder with new data strings:
    #Split original df into two, one for everything before the placeholder, one for everything after:
    dfA = template[template.index < loc]                            #create df for everything above placeholder
    dfB = template[template.index > loc]                            #create df for everything below placeholder

    #Append the three dfs to each other (part above, part to insert, part below):
    new = pd.concat([dfA, data_strings], ignore_index=True)               #append additional part to top part
    new = pd.concat([new, dfB],          ignore_index=True)               #append bottom part to new df
    
    #show=True
    if show==True:
        print(new.iloc[loc-2:loc+len(data_strings)])               #print updated section if desired

    return new


#########################################################################################

def write_input_from_csv(inputfile, placeholders, data_filenames, template_filename='template.inp'):
    
    '''Write a SWMM input file using a template.inp file with placeholders in each sections,
    and data.csv files for each section of data. The placeholders will be replaced with data.
    
    Inputs:
    inputfile:           string for name of input file to write (must end in .inp). Example: project.inp
    placeholders:        list of placeholder strings to be replaced with data. Ex: ['junctions', 'conduits']
    data_filenames:      list of file name strings for data.csv files to insert into the input file. Ex: ['junctions.csv', 'conduits.csv']
    template_filename:   string for name of template file (defaults to 'template.inp')
    
    Outputs:
    project.inp:         SWMM5 input text file
    '''
    
    template = import_template(template_filename)                                 #import template file from txt
    for i in range(len(placeholders)):                                            #loop over list of placeholders
        data = import_data(data_filenames[i])                                     #import data to insert from csv
        template = insert_data(template,  placeholders[i], data)                  #replace placeholder string with data
    template.to_csv(inputfile, header=False, index=False, quoting=3, na_rep=" ")  #write dataframe to .inp text file with specified name
        

#########################################################################################

def write_input(inputfile, placeholders, data, template_filename='template.inp'):
    
    '''Write a SWMM input file using a template.inp file with placeholders in each sections,
    and pandas dataframes holding the data to be inserted for each section. 
    The placeholders will be replaced with data.
    Section order must be the same in placeholders and data.
    
    Inputs:
    inputfile:           string for name of input file to write (must end in .inp). Example: project.inp
    placeholders:        list of placeholder strings to be replaced with data. Ex: ['junctions', 'conduits']
    data:                list of dataframes to be inserted (obtained using import_data())
    template_filename:   string for name of template file (defaults to 'template.inp')
    
    Outputs:
    project.inp:         SWMM5 input text file
    '''
    
    template = import_template(template_filename)                                  #import template file from txt
    for i in range(len(placeholders)):                                             #loop over list of placeholders
        template = insert_data(template,  placeholders[i], data[i])                #replace placeholder string with data
    template.to_csv(inputfile, header=False, index=False, quoting=3, na_rep=" ")   #write dataframe to .inp text file with specified name
        

#########################################################################################

def run(inputfile, reportfile, outputfile, exe='swmm5.exe'):
    
    '''Run SWMM using the specified input file, and create the specified output files.
    Inputs:
    inputfile:  name of .inp file to use with SWMM (must be formatted according to SWMM documentation)
    reportfile: name of .rpt file that SWMM will write to
    outputfile: name of .out file that SWMM will write binary outputs to
    exe:
    
    Outputs:
    project.rpt: SWMM report file
    project.out: SWMM output file
    '''
    
    p = sp.Popen([exe, inputfile, reportfile, outputfile], stdout=sp.PIPE, universal_newlines=True)          #run SWMM (and report process output)
    for line in p.stdout:          #for each line of process output
        if 'hour' not in line:     #if the line doesn't include the string 'hour' (to avoid a huge mass of text for each timestep)
            print(line)  

            
#POST-PROCESSING#            
#########################################################################################

def plot_map(nodes, links, dim=3, grid=None, poly=None, fig=None, ax=None, labels=False, c='k', a=1, lw=1):
    
    '''Plot system map with subcatchments, nodes, and conduits.
    Inputs:
    nodes:    pandas dataframe with columns:  X, Y, Z, Name, type
    links:    pandas dataframe with columns:  InNode, OutNode, Name
    dim:      number of dimensions (2 or 3) (note: choosing 2 for a 3d system will plot a flattened map view)
    grid:     csv file or dataframe containing grid information (xmin, xmax, xres, dx, etc.)
    poly:     optional, geopandas polygons object for subcatchments, returned by import_subcatchments(shapefile)
    fig:      figure to plot on (if not defined, new figure is created)
    ax:       axes to plot on (if not defined, new axes are created)
    labels:   True/False, whether to display node labels
    c:        color (default is black)
    a:        transparency (default is opaque)
    lw:       line width (default is 1)

    '''
    
    fromX = nodes.X.loc[links.InNode]                  #calculate coordinates for link start and end points
    fromY = nodes.Y.loc[links.InNode]
    toX   = nodes.X.loc[links.OutNode]
    toY   = nodes.Y.loc[links.OutNode]
    if dim==3:
        fromZ = nodes.Z.loc[links.InNode]
        toZ   = nodes.Z.loc[links.OutNode]
    
    if ax==None:
        if dim==2:
            fig,ax = plt.subplots(figsize=(10,10))      #create 2d figure & axis objects
        if dim==3:
            fig = plt.figure(figsize=(10,10))           #create empty figure
            ax = fig.add_subplot(111, projection='3d')  #add 3D subplot axes (requires Axes3D module from mpl_toolkits.mplot3d)
            
    if np.size(lw)==1:                                  #if a single value is given for the line width
        lw = np.full(len(links),lw)                     #convert to an array with that value for each link
    
    if np.any(poly):
        poly.plot(ax=ax, cmap='jet')                        #plot subcatchments
    
    if dim==2:
        ax.scatter(nodes.X, nodes.Y, c='k', s=8, alpha=a)   #plot nodes
        ax.scatter(nodes[nodes.type=='outfall'].X, nodes[nodes.type=='outfall'].Y, c='deepskyblue', s=30, alpha=a)  #plot spring nodes in a dif color
        for ind in links.index:                             #loop over link indices
            ax.plot((fromX.iloc[ind],toX.iloc[ind]),(fromY.iloc[ind],toY.iloc[ind]), c=c, alpha=a, linewidth=float(lw[ind])*2)                                 #plot links
        if labels==True:
            for ind in nodes.index:                         #loop over node indices
                ax.annotate(ind,xy=(nodes.X[ind],nodes.Y[ind]+100))     #label nodes with index
        ax.set_aspect('equal')
        
    if dim==3:
        ax.scatter(nodes.X, nodes.Y, nodes.Z, color='k', s=2, alpha=a)  #plot nodes
        ax.scatter(nodes[nodes.type=='outfall'].X, nodes[nodes.type=='outfall'].Y, nodes[nodes.type=='outfall'].Z,
                   c='deepskyblue', s=20, alpha=a)                                #plot spring nodes in a dif color
        for ind in links.index:                                         #loop over link indices
            ax.plot((fromX.iloc[ind],toX.iloc[ind]),(fromY.iloc[ind],toY.iloc[ind]),(fromZ.iloc[ind],toZ.iloc[ind]), c=c, alpha=a, linewidth=float(lw[ind]))    #plot links
        if labels==True:
            for ind in nodes.index:                                                 #loop over node indices
                ax.text(nodes.X[ind],nodes.Y[ind],nodes.Z[ind],  '%s' % (str(ind))) #label nodes with index
        try:
            grid = pd.read_csv(grid)
        except:
            try:
                ax.auto_scale_xyz([grid.xmin, grid.xmax], [grid.ymin, grid.ymax], [grid.zmin, grid.zmax])
            except:
                pass
    
    return fig,ax
  
    
#########################################################################################

def getQ(swmm_output_file, Q=[], columns=[]):           
    '''Get spring discharge timeseries from SWMM binary output file using swmmtoolbox package.
       
       Inputs:
       swmm_output_file: binary SWMM output file path (ex: 'modelname.out')
       Q:                optional: empty dataframe with index set to be the same datetimes as in Qobs timeseries
                         if not provided, a new dataframe will be created
       columns:          optional: list of column name strings to use (one per node with data being retrieved)
                         if not provided, original node name will be used
       
       Outputs:
       Q: dataframe with datetime index and one column per outlet node, filled with discharge values'''
    
    node_info = st.catalog(swmm_output_file, itemtype='node')                    #get node info (i.e. what data is available)
    node_info = pd.DataFrame(node_info, columns=['itemtype','name','datatype'])  #convert to dataframe
    nodenames = node_info.name.unique()                                          #get list of unique node names
    datatypes = node_info.datatype.unique()      
    for i, node in enumerate(nodenames):                                         #loop over available spring nodes
        label = '{},{},{}'.format(node_info.itemtype[0], node, datatypes[4])     #create label to identify desired data series (this is the format required by SWMMtoolbox)
        nodeInflow = st.extract(swmm_output_file, label)                         #get node timeseries data
        try:
            Q = Q.join(nodeInflow)                                               #append timeseries to dataframe
        except:
            Q = nodeInflow                                                       #if df does not already exist, create it
        try:
            Q = Q.rename(columns={nodeInflow.columns[0]: columns[i]})            #if names provided, rename columns
        except:
            pass
    return Q



#########################################################################################

def RMSE(obs, sim):
    '''Calculate RMSE (root mean square error) between simulated and observed timeseries 
    (will auto-trim to min overlapping times and remove nans).
    obs:   observed timeseries (pandas series object)
    sim:   model-predicted timeseries (pandas series object, e.g. one column returned by getQ())'''
    obs = obs[~np.isnan(obs)]       #trim to length of actual data (excluding nan values)
    sim = sim[~np.isnan(sim)]
    obs = obs[:len(sim)].values     #trim to length of sim
    sim = sim[:len(obs)].values     #trim to length of obs
    n = len(sim)
    try:
        rmse = ( sum((obs - sim)**2) / n )**0.5     #calculate RMSE
        rmse = round(rmse, 2)                       #round to two decimal places
        return rmse
    except:
        print('cannot calculate RMSE')
        return np.nan


##########################################################################################

def NSC(obs, sim):
    '''Calculate NSC (Nash-Sutcliffe coefficient) between simulated and observed timeseries
    (will auto-trim to min overlapping times and remove nans).
    obs:   observed timeseries (pandas series object)
    sim:   model-predicted timeseries (pandas series object)'''
    obs = obs[~np.isnan(obs)]       #trim to length of actual data (excluding nan values)
    sim = sim[~np.isnan(sim)]
    obs = obs[:len(sim)].values     #trim to length of sim
    sim = sim[:len(obs)].values     #trim to length of obs
    try:
        nsc = 1 - ( sum((obs - sim)**2) / sum((obs - np.mean(obs))**2) )
        nsc = round(nsc, 2)
        return nsc
    except:
        print('cannot calculate NSC')
        return np.nan
    

###########################################################################################

def VE(obs, sim):
    '''Calculate VE (Volume Error) between simulated and observed timeseries
    (will auto-trim to min overlapping times and remove nans).
    obs:   observed timeseries (pandas series object)
    sim:   model-predicted timeseries (pandas series object)'''
    obs = obs[~np.isnan(obs)]       #trim to length of actual data (excluding nan values)
    sim = sim[~np.isnan(sim)]
    obs = obs[:len(sim)].values     #trim to length of sim
    sim = sim[:len(obs)].values     #trim to length of obs
    try: 
        ve = 100 * ( (sum(sim) - sum(obs)) / (sum(obs)) )
        ve = round(ve, 2)
        return ve
    except:
        print('cannot calculate ve')
        return np.nan

#########################################################################################

def extract_parameters(swmm_inputfile, start, end, col):
    '''Extract the values being used for a desired parameter in a desired section of a SWMM input file.
    swmm_inputfile: path to the input file (should end in .inp)
    start:          string flagging the beginning of the section of interest
    end:            string flagging the start of the following section (used to define the endpoint of the section of interest)
    iPar:           column index of the parameter of interest'''
    
    swmm_inputdf = import_template(swmm_inputfile)
    iStart   = swmm_inputdf.index[swmm_inputdf[0]==start]   #df.index gets the index, df[0] looks in col 0 (returns an index object)
    iEnd     = swmm_inputdf.index[swmm_inputdf[0]==end]     #get index for end of section
    locStart = swmm_inputdf.index.get_loc(iStart[0])        #get the integer position of the index specified (need to select first item (index name) in index object)
    locEnd   = swmm_inputdf.index.get_loc(iEnd[0])          #get position of end
    data  = swmm_inputdf.loc[locStart+3:locEnd-2,0]         #get just the desired section
    lines = [data.iloc[i].split() for i in range(len(data))] #split into an array with each line and column separated out
    parameters = [float(line[col]) for line in lines]       #extract parameter values from line array

    return parameters


#########################################################################################

def report(reportfile, flags=['WARNING','ERROR','Continuity Error']):
    
    '''Print lines with the specified flags from the report file.
    Inputs:
    reportfile:  path to the SWMM project.rpt file
    flags:       list of strings to look for - lines with those strings will be printed'''
    
    report = pd.read_table(reportfile, header=None)         #read report file into dataframe
    for flag in flags:                                      #loop over flags
        ind = report.index[report[0].str.contains(flag)]    #df.index gets the index, df[0] looks in column 0 (returns an index object not just the name)
        if ind.size > 0:                                    #if the flag is found (i.e. the array is not empty)
            loc = report.index.get_loc(ind[0])              #get the location index (i.e. the row number)
            print(report.iloc[loc].values[0])               #print the line containing that index