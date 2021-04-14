#import math
#import skfmm
#import karstnet as kn
#import matplotlib.pyplot as plt

#import pyvista   as pv
#import pyvistaqt as pvqt0

    ##############################
    # karst simulation functions #
    ##############################

    def compute_karst_network(self, verbose = False):
        """
        Compute the karst network according to the parameters.

        Save the results in the `karst_simulations` list attribute.
        """
        # 1 - Initialize the parameters
        self._initialize_karst_network_parameters()

        # 2 - Compute conduits for each generation
        self._compute_iterations_karst_network()

        # 3 - Gather the conduits points to construct a node network
        edges, nodes = self._compute_nodes_network()

        # 4 - Calculate the karst network statistics indicators with karstnet and save karst network
        k = kn.KGraph(edges, nodes)
        stats = k.characterize_graph(verbose)

        maps = self.maps.copy()

        points = {}
        points['inlets']  = self.inlets
        points['outlets'] = self.outlets

        network = {}
        network['edges'] = edges
        network['nodes'] = nodes

        config = self.settings

        self.karst_simulations.append(KarstNetwork(maps, points, network, stats, config))
        return None

    # 1
    def _initialize_karst_network_parameters(self):
        """
        Initialize the karst network parameters.
        """
        self.nbr_iteration = len(self.settings['importance_factor'])
        self.inlets        = self._set_inlets_repartition()

        # Raster maps
        self.maps = {}
        self.maps['phi']      = np.ones((self.grid.ynum, self.grid.xnum))
        self.maps['velocity'] = np.zeros((self.nbr_iteration, self.grid.ynum, self.grid.xnum))
        self.maps['time']     = np.zeros((self.nbr_iteration, self.grid.ynum, self.grid.xnum))
        self.maps['karst']    = np.zeros((self.nbr_iteration, self.grid.ynum, self.grid.xnum))

        # Vector maps
        self.nodeID      = 0
        self.conduits    = []
        self.outletsNode = []
        return None

    # 1
    def _set_inlets_repartition(self):
        """
        Distribute the inlets points between the different simulating generation according to the `importance_factor` setting.
        """
        inlets_repartition = []
        inlets_proportion  = []
        total = float(sum(self.settings['importance_factor']))

        for iteration in range(self.nbr_iteration):
            inlets_proportion.append(float(self.settings['importance_factor'][iteration])/total)
            inlets_repartition.append(round(inlets_proportion[iteration]*len(self.inlets)))
        inlets_repartition[-1] = len(self.inlets)-sum(inlets_repartition[0:-1])

        inlets = []
        i = 0
        for k,repartition_nbr in enumerate(inlets_repartition):
            for num in range(repartition_nbr):
                inlets.append((self.inlets[i][0],self.inlets[i][1],k))
                i += 1
        return inlets

    # 2
    def _compute_iterations_karst_network(self):
        """
        Compute each generation of karst conduits.
        """
        # Define phi map according to outlets emplacements
        #print('-START-')
        self._compute_phi_map()

        # Compute velocity map and travel time for each iteration and draw network
        for iteration in range(self.nbr_iteration):
            self._compute_velocity_map(iteration)
            self._compute_time_map(iteration)
            self._compute_karst_map(iteration)
#            print('iteration:{}/{}'.format(iteration+1,self.nbr_iteration))
#        print('- END -')
        return None

    # 2
    def _compute_phi_map(self):
        """
        Compute the phi map.
        """
        for (x,y) in self.outlets:
            X = int(math.ceil((x - self.grid.x0 - self.grid.dx/2) / self.grid.dx))
            Y = int(math.ceil((y - self.grid.y0 - self.grid.dy/2) / self.grid.dy))
            self.maps['phi'][Y][X] = -1
        return None

    # 2
    def _compute_velocity_map(self, iteration):
        """
        Compute the velocity map.
        """
        # If it's the first iteration, iniatialize the velocity map according to the geological settings.
        if iteration == 0:
            # Geology
            if self.geology.data['geology']['mode'] == 'null':
                self.maps['velocity'][0] = np.full((self.grid.ynum, self.grid.xnum), self.settings['code_aquifere'])
            elif self.geology.data['geology']['mode'] == 'image':
                self.maps['velocity'][0] = np.where(self.geology.data['geology']['data']==1, self.settings['code_aquiclude'], self.settings['code_aquifere'])
            elif self.geology.data['geology']['mode'] == 'import':

                tableFMM = {}
                if len(self.settings['geology_id']) != len(self.settings['geology_velocity']):
                    print("- _compute_velocity_map() - Error : number of lithologies does not match with number of FMM code.")
                    sys.exit()

                for geology, codeFMM in zip(self.settings['geology_id'], self.settings['geology_velocity']):
                    if geology in self.geology.data['geology']['stats']['ID']:
                        tableFMM[geology] = codeFMM
                    else:
                        print("- initialize_velocityMap() - Warning : no geology n {} found.".format(geology))
                        tableFMM[geology] = self.settings['code_out']

                for y in range(self.grid.ynum):
                    for x in range(self.grid.xnum):
                        geology = self.geology.data['geology']['data'][y][x]
                        self.maps['velocity'][0][y][x] = tableFMM[geology]

            # Faults
            self.maps['velocity'][0] = np.where(self.geology.data['faults']['data'] > 0, self.settings['code_faults'], self.maps['velocity'][0])

            # Fractures
            self.maps['velocity'][0] = np.where(self.geology.data['fractures']['data'] > 0, self.settings['code_fractures'], self.maps['velocity'][0])

            # If out of polygon
            if self.mask is not None:
                self.maps['velocity'][0] = np.where(self.mask==1, self.settings['code_out'], self.maps['velocity'][0])

        else:
            self.maps['velocity'][iteration] = self.maps['velocity'][iteration-1]
            self.maps['velocity'][iteration] = np.where(self.maps['karst'][iteration-1] > 0, self.settings['code_conduits'], self.maps['velocity'][iteration])
        return None

    # 2
    def _compute_time_map(self, iteration):
        """
        Compute the time map.
        """
        try:
            self.maps['time'][iteration] = skfmm.travel_time(self.maps['phi'], self.maps['velocity'][iteration], dx=self.grid.dx, order=2)
        except:
            try:
                self.maps['time'][iteration] = skfmm.travel_time(self.maps['phi'], self.maps['velocity'][iteration], dx=self.grid.dx, order=1)
            except:
                raise
        return None

    # 2
    def _compute_karst_map(self, iteration):
        """
        Compute the karst map.
        """
        factor = 2
        self.step = self.grid.dx/factor

        grad_y, grad_x = np.gradient(self.maps['time'][iteration], self.grid.dx, self.grid.dy)

        get_X = lambda x : int(math.ceil((x - self.grid.x0 - self.grid.dx/2) / self.grid.dx))
        get_Y = lambda y : int(math.ceil((y - self.grid.y0 - self.grid.dy/2) / self.grid.dy))

        # get karst map from previous iteration
        if iteration > 0:
            self.maps['karst'][iteration] = self.maps['karst'][iteration-1]

        for (x,y,i) in self.inlets:
            # check if inlet iteration match with actual iteration
            if i == iteration:
                conduit = Conduit(iteration)
                X = get_X(x)
                Y = get_Y(y)
                while self.maps['time'][iteration][Y][X] != np.amin(self.maps['time'][iteration]):

                    # If X,Y is on an old karstic channel, then stop
                    if self.maps['karst'][iteration-1][Y][X] == 1:
                        conduit.add_node(self.nodeID,x,y,1) # node type 1 (conjonction)
                        self.nodeID += 1
                        break
                    else:
                        self.maps['karst'][iteration][Y][X] = 1

                    # If X,Y is on an outlet, then stop
                    #if iteration == 0 and self.phiMap[Y][X] == -1:
                    if self.maps['phi'][Y][X] == -1:
                        self.maps['karst'][iteration][Y][X] = 1
                        conduit.add_node(self.nodeID,x,y,2) # node type 2 (end)
                        self.nodeID += 1
                        break

                    # else conduit is here
                    conduit.add_node(self.nodeID,x,y,0)
                    self.nodeID += 1

                    # new move
                    alpha = math.sqrt((self.step**2)*(1/(grad_x[Y][X]**2+grad_y[Y][X]**2)))
                    dx = grad_x[Y][X] * alpha
                    dy = grad_y[Y][X] * alpha

                    # check if we are going out boundaries
                    X_ = get_X(x - dx)
                    Y_ = get_Y(y - dy)
                    if (X_ < 0) or (Y_ < 0) or (X_ > self.grid.xnum - 1) or (Y_ > self.grid.ynum - 1):
                        dx_,dy_ = self._check_boundary_conditions(iteration,X,Y)
                        x = x + dx_
                        y = y + dy_
                    else: # otherwise acts normal
                        x = x - dx
                        y = y - dy

                    X = get_X(x)
                    Y = get_Y(y)

                self.conduits.append(conduit)
        return None

    # 2
    def _check_boundary_conditions(self,iteration,X,Y):
        borders_conditions = [     #(X,Y)
            [(0,-1),(0,1),(1,0) ], # X = 0
            [(-1,0),(1,0),(0,1) ], # Y = 0
            [(0,-1),(0,1),(-1,0)], # X = xnum
            [(-1,0),(1,0),(0,-1)]] # Y = ynum

        corners_conditions = [ #(X,Y)
            [(1,0) ,(0,1) ],   # X = 0    && Y = 0
            [(1,0) ,(0,-1)],   # X = 0    && Y = ynum
            [(-1,0),(0,1) ],   # X = xnum && Y = ynum
            [(-1,0),(0,-1)]]   # X = xnum && Y = 0

        time_values = []
        rank = 0
        if (X <= 0) and (Y <= 0):
            for row in corners_conditions[0]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = corners_conditions[0][rank]
        elif (X <= 0) and (Y >= self.grid.ynum-1):
            for row in corners_conditions[1]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = corners_conditions[1][rank]
        elif (X >= self.grid.xnum-1) and (Y >= self.grid.ynum-1):
            for row in corners_conditions[2]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = corners_conditions[2][rank]
        elif (X >= self.grid.xnum-1) and (Y <= 0):
            for row in corners_conditions[3]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = corners_conditions[3][rank]
        elif X <= 0:
            for row in borders_conditions[0]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = borders_conditions[0][rank]
        elif Y <= 0:
            for row in borders_conditions[1]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = borders_conditions[1][rank]
        elif X >= self.grid.xnum-1:
            for row in borders_conditions[2]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = borders_conditions[2][rank]
        elif Y >= self.grid.ynum-1:
            for row in borders_conditions[3]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = borders_conditions[3][rank]

        else:
            print(X,Y)
        dx = moove[0] * self.step
        dy = moove[1] * self.step

        return (dx,dy)

    # 3
    def _compute_nodes_network(self):
        """
        This script computes the nodes network with all the previous calculated conduits.
        """
        # Add outlets
        outlets = self.points.points['outlets']
        for outlet in outlets:
            self.outletsNode.append((self.nodeID,outlet[0],outlet[1]))
            self.nodeID += 1

        # Connect regular nodes
        [conduit.compute_edges() for conduit in self.conduits]

        for conduit in self.conduits:
            dist = []
            last_node = conduit.nodes[-1]
            x = last_node[1]
            y = last_node[2]
            # connect last node on old karstic conduit
            if last_node[-1] == 1:
                older_conduits = [conduit_ for conduit_ in self.conduits if conduit_.iteration < conduit.iteration]
                all_nodes = []
                for older_conduit in older_conduits:
                    for node in older_conduit.nodes:
                        all_nodes.append(node)
                for node in all_nodes:
                    dx   = abs(x - node[1])
                    dy   = abs(y - node[2])
                    dist.append((math.sqrt(dx**2 + dy**2),node[0]))
                min_node = min(dist)
                conduit.edges.append((last_node[0],min_node[1]))

            # connect last node on outlet coordinates
            if last_node[-1] == 2:
                for outlet in self.outletsNode:
                    dx = abs(x - outlet[1])
                    dy = abs(y - outlet[2])
                    dist.append((math.sqrt(dx**2 + dy**2),outlet[0]))
                min_outlet = min(dist)
                conduit.edges.append((last_node[0],min_outlet[1]))

        # Get data for karstnet
        EDGES = []
        NODES = {}
        for conduit in self.conduits:
            for edge in conduit.edges:
                EDGES.append(edge)
            for node in conduit.nodes:
                NODES[node[0]] = (node[1],node[2])
        for outlet in self.outletsNode:
            NODES[outlet[0]] = (outlet[1],outlet[2])
        return (EDGES,NODES)


    ###########################
    # visualization functions #
    ###########################

    def show_catchment(self, data='geology', title=None, mask=False, cmap='binary'):
        """
        Show the entire study domain.
        """
        fig, ax1 = plt.subplots()
        if title is None:
            title = data
        fig.suptitle(title, fontsize=16)

        if   data == 'geology':
            d = self.geology.data['geology']['data']
        elif data == 'faults':
            d = self.geology.data['faults']['data']
        elif data == 'fractures':
            d = self.geology.data['fractures']['data']

        if mask==True:
            if   data == 'geology':
                d = self.geology_masked['geology']
            elif data == 'faults':
                d = self.geology_masked['faults']
            elif data == 'fractures':
                d = self.geology_masked['fractures']

        im1 = ax1.imshow(d, extent=_extents(self.grid.x) + _extents(self.grid.y), origin='lower', cmap=cmap)
        fig.colorbar(im1, ax=ax1)
        if self.settings['data_has_polygon']:
            closed_polygon = self.polygon.polygon[:]
            closed_polygon.append(closed_polygon[0])
            x,y = zip(*closed_polygon)
            ax1.plot(x,y, color='red', label='polygon')
        for key in self.points.points:
            x,y = zip(*self.points.points[key])
            ax1.plot(x,y,'o',label=key)
        ax1.set_aspect('equal', 'box')
        plt.legend(loc='upper right')
        plt.show()
        return None

    def _show_karst_network(self, sim=-1, iteration=-1, cmap='binary'):
        """
        Show the simulated karst network.
        """
        karst_network = self.karst_simulations[sim]

        fig, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2, 2, sharex=True, sharey=True)
        fig.suptitle('Karst Network', fontsize=16)

        ax1.imshow(karst_network.maps['phi'], extent=_extents(self.grid.x) + _extents(self.grid.y), origin='lower', cmap=cmap)
        ax1.set_title('Phi')

        ax2.imshow(karst_network.maps['velocity'][iteration], extent=_extents(self.grid.x) + _extents(self.grid.y), origin='lower', cmap=cmap)
        ax2.set_title('Velocity')

        ax3.imshow(karst_network.maps['time'][iteration], extent=_extents(self.grid.x) + _extents(self.grid.y), origin='lower', cmap=cmap)
        ax3.set_title('Time')

        ax4.imshow(karst_network.maps['karst'][iteration], extent=_extents(self.grid.x) + _extents(self.grid.y), origin='lower', cmap=cmap)
        ax4.set_title('Karst')

        fig.subplots_adjust(hspace=0.5)
        plt.show()
        return None

    def show(self, data=None, title=None, cmap='binary', probability=False):
        """
        Show the entire study domain.
        """
        if data is None:
            data = self.karst_simulations[-1].maps['karst'][-1]

        if probability == True:
            data = self._compute_average_paths()

        fig, ax1 = plt.subplots()
        if title is not None:
            fig.suptitle(title, fontsize=16)

        im1 = ax1.imshow(data, extent=_extents(self.grid.x) + _extents(self.grid.y), origin='lower', cmap=cmap)

        fig.colorbar(im1, ax=ax1)

        if self.settings['data_has_polygon']:
            closed_polygon = self.polygon.polygon[:]
            closed_polygon.append(closed_polygon[0])
            x,y = zip(*closed_polygon)
            ax1.plot(x,y, color='red', label='polygon')
        for key in self.points.points:
            x,y = zip(*self.points.points[key])
            ax1.plot(x,y,'o',label=key)
        ax1.set_aspect('equal', 'box')
        plt.legend(loc='upper right')
        plt.show()
        return None

    def _compute_average_paths(self, show=False):
        """
        Compute the mean of all the simulations.
        """
        karst_maps = []
        for karst_simulation in self.karst_simulations:
            karst_maps.append(karst_simulation.maps['karst'][-1])

        karst_prob = sum(karst_maps)/len(karst_maps)

        if self.mask is not None:
            karst_prob = np.ma.MaskedArray(karst_prob, self.mask)

        return karst_prob

    def compare_stats(self, data=None, mean=False):
        """
        Compare statistics between reference indicators and calculated networks.
        """

        cpd  = []
        cv_d = []
        cv_l = []
        o_e  = []
        l_e  = []
        spl  = []
        m_d  = []
        cvd  = []
        vars = [cpd,cv_d,cv_l,o_e,l_e,spl,m_d,cvd]

        if data == None:
            i = -1
        else:
            i = 0

        for karst_network in self.karst_simulations[i:]:
            stats_list = []
            stats_list.append(karst_network.stats['cpd'])
            stats_list.append(karst_network.stats['cv degree'])
            stats_list.append(karst_network.stats['cv length'])
            stats_list.append(karst_network.stats['orientation entropy'])
            stats_list.append(karst_network.stats['length entropy'])
            stats_list.append(karst_network.stats['aspl'])
            stats_list.append(karst_network.stats['mean degree'])
            stats_list.append(karst_network.stats['correlation vertex degree'])

            cpd.append(karst_network.stats['cpd'])
            cv_d.append(karst_network.stats['cv degree'])
            cv_l.append(karst_network.stats['cv length'])
            o_e.append(karst_network.stats['orientation entropy'])
            l_e.append(karst_network.stats['length entropy'])
            spl.append(karst_network.stats['aspl'])
            m_d.append(karst_network.stats['mean degree'])
            cvd.append(karst_network.stats['correlation vertex degree'])

        if mean == False:
            print('\n')
            print('STATS for modelisation')
            print('%-10s%-12s%-12s%-12s%-12s' % ('Var', 'Value', 'Min ref', 'Max ref', 'Result'))
            print(54*'-')
            for stat_calc, key in zip(stats_list, self.reference_statistics):
                if stat_calc > self.reference_statistics[key][1]:
                    result = 'out +'
                elif stat_calc < self.reference_statistics[key][0]:
                    result = 'out -'
                else:
                    result = 'IN'
                print('%-10s%-12s%-12s%-12s%-12s' % (key, round(stat_calc,3), round(self.reference_statistics[key][0],3), round(self.reference_statistics[key][1],3), result))

        else:
            print('\n')
            print('MEAN STATS')
            print('%-10s%-12s%-12s%-12s%-12s' % ('Var', 'Value', 'Min ref', 'Max ref', 'Result'))
            print(54*'-')
            for var, key in zip(vars, self.reference_statistics):
                mean = sum(var)/len(var)
                if mean > self.reference_statistics[key][1]:
                    result = 'out +'
                elif mean < self.reference_statistics[key][0]:
                    result = 'out -'
                else:
                    result = 'IN'
                print('%-10s%-12s%-12s%-12s%-12s' % (key, round(mean,3), round(self.reference_statistics[key][0],3), round(self.reference_statistics[key][1],3), result))
        return None



#####################
class KarstNetwork():
    """
    A class for stroring a calculated karst network.
    """
    def __init__(self, maps, points, network, stats, settings):
        self.maps     = maps
        self.points   = points
        self.network  = network
        self.stats    = stats
        self.settings = settings
