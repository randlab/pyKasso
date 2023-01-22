#######
# DOC #
#######

PARAMETERS_DOC = """
'data_has_mask' : bool
    Defines if a mask is used or not.
    If true, a mask must be defined with 'mask_data' parameter.
'mask_data' : str || array
    Defines the mask vertices.
    Mask datafile path or list of vertices coordinates.
    Useful only when 'data_has_mask' is true.
'outlets_mode' : str
    Defines the outlets mode.
    'random'    - Full random points
    'import'    - Import points
    'composite' - Add n random points to imported points
'outlets_data' : str || list
    Defines the outlets.
    Outlets datafile path or list of outlets coordinates.
    Useful only when 'outlets_mode' parameter is on 'import' or 'composite'.
'outlets_number' : str
    Defines the number of outlets to generate.
    Useful only when 'outlets_mode' parameter is on 'random' or 'composite'.
'outlets_shuffle' : bool
    Defines whether to shuffle the order of the outlets randomly.
    False - don't shuffle, True - shuffle randomly
    Useful only when iterating over outlets.
'outlets_importance' : list
    Defines the proportion of outlets to be distributed across each iteration.
    Length of array indicates number of outlet iterations,
    each integer indicates number of outlets to run in that iteration,
    sum of integers = total number of outlets
    [1] - a single iteration with all outlets,
    [1,1,1] - three iterations with one outlet in each,
    [1,2,3] - three iterations with one outlet in the first, 2 outlets in the second, and 3 outlets in the third.
    Useful only when iterating over outlets.
'inlets_mode' : str
    Defines the inlets mode.
    'random'    - Full random points
    'import'    - Import points
    'composite' - Add n random points to imported points
'inlets_data' : str || list
    Defines the inlets.
    Inlets datafile path or list of inlets coordinates.
    Useful only when 'inlets_mode' parameter is on 'import' or 'composite'.
'inlets_number' : str
    Defines the number of inlets to generate.
    Useful only when 'inlets_mode' parameter is on 'random' or 'composite'.
'inlets_shuffle' : bool
    Defines whether to shuffle the order of the inlets randomly.
    False - don't shuffle, True - shuffle randomly
    Useful only when iterating over inlets.
'inlets_per_outlet' : list
    Defines the proportion of inlets to be distributed across each outlet.
    Length of array indicates number of outlets,
    each integer indicates number of inlets to assign to that outlet,
    sum of integers = total number of inlets
    [1] - a single iteration with all inlets to one outlet,
    [1,1,1] - three outlets with one inlet in each,
    [1,2,3] - three outlets with one inlet in the first, 2 inlets in the second, and 3 inlets in the third.
    Useful only when iterating over inlets and outlets.
'inlets_importance' : list
    Defines the proportion of inlets to be distributed across each iteration.
    Length of array indicates number of inlet iterations,
    each integer indicates number of inlets to run in that iteration,
    sum of integers = total number of inlets
    [1] - a single iteration with all inlets,
    [1,1,1] - three iterations with one inlet in each,
    [1,2,3] - three iterations with one inlet in the first, 2 inlets in the second, and 3 inlets in the third.
    Useful only when iterating over inlets.
'geology_mode' : str
    Defines the geological mode.
    'null'  - No geology
    'gslib' - Import geology via gslib
    'csv'   - Import geology via csv
    'image' - Import geology via image
'geology_datafile' : str
    Defines the geological datafile path.
    Useful only when 'geology_mode' parameter is not 'null'.
'topography_mode' : str
    Defines the topography mode.
    'null'  - No topography
    'gslib' - Import topography via gslib
    'csv'   - Import topography via csv
'topography_datafile' : str
    Defines the topography datafile path.
    Useful only when 'topography_mode' parameter is not 'null'.
'orientation_mode' : str
    Defines the orientation mode.
    'null'    - No orientation
    'topo'    - Calculate from topography
    'surface' - Calculate from csv file of a surface (useful if using lower surface of karst unit)
'orientation_datafile' : str
    Defines the orientation datafile path.
    Useful only when 'orientation_mode' parameter is not 'null'.
'faults_mode' : str
    Defines the mode for the faults.
    'null'  - No faults
    'gslib' - Import faults via gslib
    'csv'   - Import faults via csv
    'image' - Import faults via image
'faults_datafile' : str
    Defines the faults datafile path.
    Useful only when the 'faults_mode' parameter is on 'gslib', 'csv' or 'image'.
'fractures_mode' : str
    Defines the mode for the fractures.
    'null'   - No fractures
    'gslib'  - Import fractures via gslib
    'csv'    - Import fractures via csv
    'image'  - Import fractures via image
    'random' - Generate fractures randomly
'fracture_datafile' : str
    Defines the fractures datafile path.
    Useful only when the 'fractures_mode' parameter is on 'gslib', 'csv' or 'image'.
'fractures_densities' : list
    Defines the fractures densitiy for each fractures family.
    Useful only when the 'fractures_mode' parameter is on 'random'.
'fractures_min_orientation' : list
    Defines the minimum orientation of the fracturation for each fractures family.
    Useful only when the 'fractures_mode' parameter is on 'random'.
'fractures_max_orientation' : list
    Defines the maximum orientation of the fracturation for each fractures family.
    Useful only when the 'fractures_mode' parameter is on 'random'.
'fractures_min_dip' : list
    Defines the minimum dip of the fracturation for each fractures family.
    Useful only when the 'fractures_mode' parameter is on 'random'.
'fractures_max_dip' : list
    Defines the maximum dip of the fracturation for each fractures family.
    Useful only when the 'fractures_mode' parameter is on 'random'.
'fractures_alpha' : list
    Defines alpha, a parameter in the fracturation law.
    Useful only when the 'fractures_mode' parameter is on 'random'.
'fractures_min_length' : list
    Defines the minimum lenght for all the fractures.
    Useful only when the 'fractures_mode' parameter is on 'random'.
'fractures_max_length' : list
    Defines the maximum lenght for all the fractures.
    Useful only when the 'fractures_mode' parameter is on 'random'.
'algorithm' : str
    Defines the algorithm to use when calculating travel time to spring.
    'Isotropic2' - isotropic 2D
    'Isotropic3' - isotropic 3D
    'Riemann2'   - anisotropic 2D
    'Riemann3'   - anisotropic 3D
    See AGD-HFM documentation for full list of options.
'cost_out' : float, (default: 0.999)
    Defines the fast-marching value for the outside of the study area.
    The value must be between 0 and 1 and should be high to avoid unrealistic conduits.
'cost_aquifer' : float, (default: 0.3)
    Defines the fast-marching value for the aquifer cells.
    Should be between 0 and 1 and lower than aquiclude but higher than conduits.
'cost_aquiclude' : float, (default: 0.8)
    Defines the fast-marching value for the aquiclude cells.
    Should be between 0 and 1 and higher than aquiclude but lower than cost_out
'cost_faults' : float, (default: 0.2)
    Defines the fast-marching value for the faults cells.
    Should be between 0 and 1 and between conduits and cost_out. Higher = conduit will avoid faults, lower = conduits will follow faults
'cost_fractures' : float, (default: 0.2)
    Defines the fast-marching value for the fractures cells.
    Should be between 0 and 1 and between conduits and cost_out. Higher = conduit will avoid fractures, lower = conduits will follow fractures
'cost_conduits' : float, (default: 0.01)
    Defines the fast-marching value for the conduits cells.
    Should be between 0 and 1 but lower than aquifer (for conduits to preferentially follow each other)
'cost_ratio' : float, (default: 0.25)
    Defines the fast-marching ratio of travel cost parallel to gradient / travel cost prependicular to gradient.
    Should be between 0 and 0.5.
'geology_id' : list
    Defines the geology id (from geology datafile) to consider in the simulation.
    Useful only when the 'geology_mode' parameter is on 'gslib' or 'csv'.
'rand_seed' : int
    Defines the random seed.
    May help for reproduicity.
'verbosity' : int
    Define the verbosity (how much output to print during runs).
    0 - print minimal output, 1 - print medium output, 2 - print max output
"""
