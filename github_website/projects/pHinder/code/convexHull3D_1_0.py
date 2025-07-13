# The difference between convexhull_v2_1 and convexhull_v2_2 is that the emerging hull
# is now calculated in a step-wise manner that ensures general position of the hull vertices.
global zero
zero = 1e-10 # was 1e-3

from compGeometry import Vertex, Edge, Triangle, gp2D, gp3D, planeCoefficients3D

class convexHull3D:

	#def __init__(self,list_of_vertices,minArea=1e-3,sortBy="x",noChangeID=0):
	def __init__(self,list_of_vertices,minArea=zero,sortBy="x",noChangeID=0): # XXXX

		#A list of Cartesian coordinates in dimension 3.
		self.list_of_vertices = list_of_vertices
		self.non_mutated_vertices = list_of_vertices
		self.edges = {}
		self.horizon_edges = {}
		self.horizon_cone = {}
		self.i = 0
		self.loop_count = 0
		self.number_of_cumulative_search_steps = 0
		#self.number_of_cumulative_traverse_steps = 0

		# Minimum facet area for collinearity testing.
		self.minArea = minArea
		#self.collinear = []

		self.encountered = {}
		#self.iterator_count = 0

		"""The class attribute self.start_edge is critical to the
		   efficiency of the calculation of the convex hull.  It
		   ensures that traversals of the edges of the current
		   convex hull (the hull changes throughout the calculation)
		   begin within and remain within the edges that comprise
		   the current convex hull.  Rather than deleting edges that
		   are not part of the convex hull these edges are simply
		   disconnected in the function self.build_cone.  As a result,
		   the class attribute self.edges contains a subset of old
		   edges that are not part of the current convex hull, and a
		   subset of edges, old and new, that are a part of the
		   convex hull.  The attribute self.start_edge is an edge
		   that is always in the subset of new edges that are part
		   of the current convex hull.  For those class methods that
		   loop over the edges of the current convex hull, e.g., the
		   rotation functions, the attribute self.start_edges ensures
		   that the loop is initialized with an edge that is in the
		   current convex hull.
		"""
		self.start_edge = None

		#Call the class methods that calculate the hull.
		self.sort_vertices(sortBy=sortBy,noChangeID=noChangeID)
		self.make_simplex_3d()
		self.calculate_the_convex_hull()
		self.facets = None
		self.hyperPlanes = None

		#PPPP: Find the right place to put this reminder.
		#!!!!Be careful transforming the data.
		#Translating it by 100 gives the same answer.
		#Translating it by 1000 givea a slightly different answer,
		#because the coordinate become large, the determinant becomes
		#larger, which introduces errors into the mathematics.
		#See O'rourke for an explanation. The main point here is that
		#if the coordinates are large (how large?) the code may not
		#behave properly.  Should I take the time to account for this
		#using integer math?
	
	def sort_vertices(self, sortBy="x",noChangeID=0):
		temp_vertices = {}
		i = 1
		for v in self.list_of_vertices:
			#from compGeometry import Vertex
			nV = Vertex((v.x,v.y,v.z))
			# Assign new vertex ID.
			if not noChangeID:
				nV.id = i
			# Retain original vertex ID.
			else:
				nV.id = v.id
			try:
				nV.data = v.data
			except:
				nV.data = None

			"""Map the vertices by x-coordinate.
			   The y- and z-coordinate are included in the key of the
			   mapping to ensure that the key is unique."""
			if sortBy == "x":
				temp_vertices.update({(nV.x,nV.y,nV.z):nV})
			elif sortBy == "y":
				temp_vertices.update({(nV.y,nV.x,nV.z):nV})
			elif sortBy == "z":
				temp_vertices.update({(nV.z,nV.x,nV.y):nV})
			else:
				# Default is sort by x coordinate.
				temp_vertices.update({(nV.x,nV.y,nV.z):nV})
			i += 1

		#Sort the set of vertices by x-coordinate.
		keys = sorted(temp_vertices)

		#Replace the unsorted vertex list with the sorted vertex list.
		self.list_of_vertices = []
		for k in keys:
			self.list_of_vertices.append(temp_vertices[k])

	def make_simplex_3d(self):

		# XXXX This is where I ensure general position in the initial tetrahedron simplex.

		# 3 Vertices.
		#############
		v1, v2 = self.list_of_vertices[0], self.list_of_vertices[1]
		vGp = 0
		while not vGp:
			v3 = self.list_of_vertices[2]
			vGp = gp2D(v1, v2, v3)
			if not vGp:
				v3.jostle()

		# 4 Vertices.
		#############
		vGp = 0
		while not vGp:
			v4 = self.list_of_vertices[3]
			vGp = gp3D(v1,v2,v3,v4)
			if not vGp:
				v4.jostle()

		#vList = self.list_of_vertices
		# #from compGeometry import testGp3D, testGp4D
		# iterate, v1, v2, v3, v4 = 1, None, None, None, None
		# while iterate:

		# 	# Vertices that need to be jostled into general position.
		# 	#########################################################
		# 	vGp = []

		# 	# 3 Vertices.
		# 	#############
		# 	v1, v2 = vList[0], vList[1]
		# 	i, testGP = 2, 0
		# 	while not testGP and i < len(vList):
		# 		v3 = vList[i]
		# 		#testGP = testGp3D(v1,v2,v3)
		# 		testGp = gp2D(v1, v2, v3)
		# 		if not testGP:
		# 			vGp.append((v3,i))
		# 			v3 = None
		# 		i += 1

		# 	#print("From 3: ", vGp, i )
		# 	#print(v1, v2, v3)

		# 	# 4 Vertices.
		# 	#############
		# 	testGP = 0
		# 	while not testGP and i < len(vList):
		# 		v4 = vList[i]
		# 		#testGP = testGp4D(v1,v2,v3,v4)
		# 		testGp = gp3D(v1,v2,v3,v4)
		# 		if not testGP:
		# 			vGp.append((v4,i))
		# 			v4 = None
		# 		i += 1

		# 	#print("From 4: ", vGp, i)
		# 	#print(v1, v2, v3, v4)

		# 	# Jostle any vertices that were not in general position.
		# 	########################################################
		# 	if vGp:
		# 		for vTuple in vGp:
		# 			#print("jostling iteration", iterate)
		# 			vTuple[0].jostle()

		# 	if None in [v1,v2,v3,v4]:
		# 		iterate += 1
		# 	else:
		# 		iterate = 0

		# 	print(v1, v2, v3, v4)

		# Delete used vertices.
		#######################
		for v in [v1,v2,v3,v4]:
			self.list_of_vertices.remove(v)

		# Now that general position has been established,
		# build the initial hull in dimension 4. 
		#################################################

		# General position vertices.
		############################


		"""There are four vertices (0-faces) in the intial simplex:
		   0, 1, 2, 3.
		   There are six edges (1-faces) in the initial simplex:
		   01, 02, 03, 12, 13, 23.
		   There are four facets (d-1-faces) in the initial simplex:
		   012, 013, 023, 123.
		   There is one volume (d-face) in the intial simplex:
		   0123."""
		#list_of_simplex_vertices = self.list_of_vertices[0:4]
		#v0 = list_of_simplex_vertices[0]
		#v1 = list_of_simplex_vertices[1]
		#v2 = list_of_simplex_vertices[2]
		#v3 = list_of_simplex_vertices[3]

		#Delete the subset of vertices used to make the simplex.
		#del self.list_of_vertices[0:4]

		#Make the facets of the simplex.
		t012 = Triangle(v1, v2 , v3, VertexR=v4)
		t012.make_ltriangle(); t012.set_edge_topology()
		self.edges.update(t012.get_edges())

		t023 = Triangle(v1, v3, v4, VertexR=v2)
		t023.make_ltriangle(); t023.set_edge_topology()
		self.edges.update(t023.get_edges())

		t031 = Triangle(v1, v4, v2, VertexR=v3)
		t031.make_ltriangle(); t031.set_edge_topology()
		self.edges.update(t031.get_edges())

		t123 = Triangle(v2, v3, v4, VertexR=v1)
		t123.make_ltriangle(); t123.set_edge_topology()
		self.edges.update(t123.get_edges())

		t130 = Triangle(v2, v4, v1, VertexR=v3)
		t130.make_ltriangle(); t130.set_edge_topology()
		self.edges.update(t130.get_edges())

		t231 = Triangle(v3, v4, v2, VertexR=v1)
		t231.make_ltriangle(); t231.set_edge_topology()
		self.edges.update(t231.get_edges())

		#Randomly choose an edge from the simplex.
		#This edge is the start point for iterating over the edges
		#of the current convex hull.
		self.start_edge = self.edges[list(self.edges)[0]]

	def calculate_the_convex_hull(self):

		# Get dependencies.
		###################
		#from compGeometry import planeCoefficients3D
		#from compGeometry import Vertex
		
		i = 0
		while i < len(self.list_of_vertices):
			
			# Get the next vertex in the set.
            #################################
			vT, ivGp = self.list_of_vertices[i], 0
			
			number_of_search_steps = 0
			#XXXX: Explain what determines a horizon ridge.  Negative
			#volume, the vertex is below the plane of triangle face;
			#therefore, the face is a member of the hull.  Positive volume,
			#the vertex is above the plane of the triangle face; therefore
			#the face is visible and is not a member of the new hull.
			#PPPP: Move through the qes by topology. Start the traversal of
			#of the edges of the hull with a horizon ridge of the previous
			#hull.  Further explanation. Negative: its outside? Vice Versa?
			"""Traverse the set of horizon edges in a counterclockwise
			   direction. Guarantee that the right face of the Edge
			   object is the visible face of the horizon edge. This
			   ensures that traversal of the horizon edges is
			   counterclockwise in its direction, and it ensures that
			   the horizon cone is comprised of the right face of each
			   edge object."""
			if not self.number_of_cumulative_search_steps:
				edges = self.edge_iterator(self.start_edge)
			else:
				"""The set of vertices is sorted by their x-ccordinate.
				   As a result, at least one edge of the horizon cone
				   built during the last iteration of this while loop
				   will be visible.  Because of this property, the
				   visibility of every edge of the convex hull need not
				   be considered, and the calculation of the hull becomes
				   much more efficient.

				   It is important to note that the horizon cone is
				   comprised of both the subset of horizon edges and the
				   subset of the new horizon cone edges from the last
				   iteration of this while loop.

				   It is also important to note that this added
				   efficiency is based on a key assumption:  if one or
				   more facets are visible from the current vertex, the
				   previous horizon cone contains an edge that is
				   incident upon one of these visible faces.  In other
				   words, it is assumed that the previous horizon cone
				   contains a current horizon edge.  This assumption
				   holds only when the set of vertices has been sorted
				   so that the next vertex is outside of the current
				   hull.
				"""
				edges = self.horizon_cone

			# Vertex3D coordinates.
			#######################
			x , y, z = vT.x, vT.y, vT.z

			iResults = {}
			for e in edges.values():

				number_of_search_steps += 1

				#Facet 1 : left facet.
				######################
				lfAboveBelow = e.lf.beneath_beyond(vT)
				# c = planeCoefficients3D(e.lf.v1, e.lf.v2, e.lf.v3)
				# # Have to divide by 6. to match volume given by beneath_beyond()!
				# #################################################################
				# lfAboveBelow = (c[0]*x + c[1]*y + c[2]*z + c[3])/6.
				# aL = e.lf.beneath_beyond(vT, returnVolume=1)
				# print(lfAboveBelow, aL)

				#Facet 1 : right facet.
				#######################
				rfAboveBelow = e.rf.beneath_beyond(vT)
				# c = planeCoefficients3D(e.rf.v1, e.rf.v2, e.rf.v3)
				# # Have to divide by 6. to match volume given by beneath_beyond()!
				# #################################################################
				# rfAboveBelow = (c[0]*x + c[1]*y + c[2]*z + c[3])/6. # 2017.05.14 XXXX was .6 not 6 !
				# aR = e.rf.beneath_beyond(vT, returnVolume=1)
				# print(rfAboveBelow, aR)

				if (not lfAboveBelow) or (not rfAboveBelow):
					ivGp = 1
					break 

				# # Test for general position.  The zero tolerance of the algorithm is critical.
				# # 1e-3 seems to be working well, but I can't be sure it is a general solution.
				# # It appears to be the best emperical solution.
				# ##############################################################################
				# if abs(lfAboveBelow) <= zero or abs(rfAboveBelow) <= zero: 
				# 	ivGp = 1
				# 	break 

				# # Convert volume results after general position is evaluated.
				# #############################################################
				# if lfAboveBelow > 0.:
				# 	lfAboveBelow = -1
				# else:
				# 	lfAboveBelow = 1
				# if rfAboveBelow > 0.:
				# 	rfAboveBelow = -1
				# else:
				# 	rfAboveBelow = 1

				# Save results for classification in case the vertex is in GP.
				##############################################################
				iResults.update({e:((lfAboveBelow,e.lf),(rfAboveBelow,e.rf))})

			# This vertex is not in general position.  Save to fix....
			###################################################################
			if ivGp:
				print("jostling.....")
				vT.jostle()
				continue

			# vT is in GP so classify iResults.
			####################################
			for e in iResults:

				# Unpack the result.
				####################
				f1AboveBelow = iResults[e][0][0]
				f1 = iResults[e][0][1]
				f2AboveBelow = iResults[e][1][0]
				f2 = iResults[e][1][1]

				# Classify result.
				# This loop finds the first horizon edge, then efficiently 
				# identifies the remaining horizon edges using a function 
				# call to traverse_the_horizon_edges(). The function returns
				# all of the new horizon edges and then "breaks" the loop.
				# Thus, the behavior of this loop is to progress until the
				# first horizon edge is encountered, then it is terminated.
				############################################################
				l, r = f1AboveBelow, f2AboveBelow
				if l != r:
					if l == 1:
						self.traverse_the_horizon_edges(vT,self.edges[e.flipid])
					else:
						self.traverse_the_horizon_edges(vT,e)
					self.number_of_cumulative_search_steps += number_of_search_steps
					break

				
			# Iterate to next vertex index.
			###############################
			#print("Vertex added: ", i)
			i += 1
			#print("breaking here....")
			#raise SystemExit 

	def traverse_the_horizon_edges(self,vertex,start_edge):

		self.horizon_edges = [start_edge,]

		"""Interestingly, counterclockwise rotation over the edges
		   incident on the origin vertex of the current horizon edge
		   to find the next horizon edge is more efficient than
		   clockwise rotation.  In other words, when rotating
		   counterclockwise over the edges incident on the origin
		   vertex, then next horizon edge is almost always
		   encountered in fewer iterations than if the search for
		   the next horizon edge proceeds in the clockwise direction.
		   Ties are not uncommon: the next horizon edge is encountered
		   in the same number of iterations whether the search proceeds
		   in either a counterclockwise or clockwise direction.
		"""
		ccw_origin_edges = self.rotate_ccw_about_the_origin_vertex(start_edge)

		i = 1
		e = ccw_origin_edges[i]
		ePrev = e.id
		while e.id != start_edge.flipid:
			#self.number_of_cumulative_traverse_steps += 1
			"""This is the general solution for traversing the horizon
			   edges in a counterclockwise direction.  Given a horizon
			   edge (start_edge), rotate over the set of edges that are
			   incident upon the origin vertex of start_edge.  The
			   orientation of these edges is not predictable, i.e., each
			   each edge is a directed line that points back to,
			   or away from the origin vertex.  In this algorithm, it is
			   necessary to select for the edge that is a directed line
			   that points back to the origin vertex.  The direction of
			   line is subsequently reversed.  This ensures that
			   as the algorithm moves from horizon edge to horizon edge
			   it proceeds in a counterclockwise direction.  This also
			   ensures that the right face of each horizon edge is the
			   visible face, which simplifies the process of building
			   the new horizon cone"""
			# if e.destination_vertex != start_edge.origin_vertex:
			# 	e = self.edges[e.flipid]
			# 	ePrev = e.id

			if e.destination_vertex != self.edges[ePrev].origin_vertex:
				e = self.edges[e.flipid]
			ePrev = e.id

			# # Have to see if vT is in general position here.
			# #Facet 1 : left facet.
			# ######################
			# #from compGeometry import planeCoefficients3D, testGp4D
			# #x, y, z = vertex.x, vertex.y, vertex.z
			# #c = planeCoefficients3D(e.lf.v1, e.lf.v2, e.lf.v3)
			# #lfAboveBelow = c[0]*x + c[1]*y + c[2]*z + c[3]

			# #test1 = testGp4D(e.lf.v1, e.lf.v2, e.lf.v3, vertex)
			# test1 = gp3D(e.lf.v1, e.lf.v2, e.lf.v3, vertex, returnVolume=1)

			# #Facet 1 : left facet.
			# ######################
			# #c = planeCoefficients3D(e.rf.v1, e.rf.v2, e.rf.v3)
			# #rfAboveBelow = c[0]*x + c[1]*y + c[2]*z + c[3]

			# #test2 = testGp4D(e.lf.v1, e.lf.v2, e.lf.v3, vertex)
			# test2 = gp3D(e.rf.v1, e.rf.v2, e.rf.v3, vertex, returnVolume=1) # XXXX was lf 

			# print(start_edge.id, e.id, test1, test2)

			# #print("while...1")
			# while (not test1 or not test2):
			# 	print("vertex is not in general position relative to previous horizon cone")
			# 	vertex.jostle()
			# 	#test1 = testGp4D(e.lf.v1, e.lf.v2, e.lf.v3, vertex)
			# 	test1 = gp3D(e.lf.v1, e.lf.v2, e.lf.v3, vertex, returnVolume=1)
			# 	#test2 = testGp4D(e.lf.v1, e.lf.v2, e.lf.v3, vertex)
			# 	test2 = gp3D(e.rf.v1, e.rf.v2, e.rf.v3, vertex, returnVolume=1) # XXXX was lf 
			# 	#print(test1, test2)
			# #print("while...2")

			# Test for general position.  The zero tolerance of the algorithm is critical.
			# 1e-3 seems to be working well, but I can't be sure it is a general solution.
			# It appears to be the best emperical solution.
			##############################################################################
			#print(lfAboveBelow, rfAboveBelow)
			#if abs(lfAboveBelow) <= zero or abs(rfAboveBelow) <= zero:
			#	print("should be trying to correct.....")
			#	vertex.jostle()
			#	#ivGp = 1
			#	#break 

			l = self.edges[e.id].lf.beneath_beyond(vertex)
			r = self.edges[e.id].rf.beneath_beyond(vertex)

			if l != r:
				self.horizon_edges.append(self.edges[e.id])
				ccw_origin_edges =  self.rotate_ccw_about_the_origin_vertex(self.edges[e.id])
				i = 1
				e = ccw_origin_edges[i]
			else:
				i += 1
				e = ccw_origin_edges[i]
 
			# print(l, r, i)
			# # e is horizon edge.
			# if l != r:
			# 	if l == 1:
			# 		i = 1 
			# 		self.horizon_edges.append(self.edges[e.id])
			# 		print("start rotate 1...")
			# 		ccw_origin_edges =  self.rotate_ccw_about_the_origin_vertex(self.edges[e.id])
			# 		print("end rotate 1...")
			# 		e = ccw_origin_edges[i]
			# 	else:
			# 		self.horizon_edges.append(self.edges[e.flipid])
			# 		print("start rotate 2...")
			# 		ccw_origin_edges =  self.rotate_ccw_about_the_origin_vertex(self.edges[e.flipid])
			# 		e = ccw_origin_edges[i]
			# 		print("end rotate 2...")
			# else:
			# 	i += 1
			# 	e = ccw_origin_edges[i]

			# # e is horizon edge.
			# ####################
			# if l != r:
			# 	print("L is ", l # XXXX I believe I am expected L to alway be one by my convention)
			# 	if l == 1:
			# 		self.horizon_edges.append(self.edges[e.flipid])
			# 		self.traverse_the_horizon_edges(vT,self.edges[e.flipid])
			# 	else:
			# 		self.horizon_edges.append(self.edges[e.id])
			# 		self.traverse_the_horizon_edges(vT,e)

			# 	#self.number_of_cumulative_search_steps += number_of_search_steps
			# 	#self.horizon_edges.append(self.edges[e.id])
			# 	ccw_origin_edges =  self.rotate_ccw_about_the_origin_vertex(self.edges[e.id])
			# 	i = 1
			# 	e = ccw_origin_edges[i]

			# 	# Found the next horizon edge so move on.
			# 	#########################################
			# 	break
			# 	# self.horizon_edges.append(self.edges[e.id])
			# 	# ccw_origin_edges =  self.rotate_ccw_about_the_origin_vertex(self.edges[e.id])
			# 	# i = 1
			# 	# e = ccw_origin_edges[i]
			# # e is not a horizon edge.
			# ##########################
			# else:
			# 	i += 1
			# 	e = ccw_origin_edges[i]

		#print("made it here...3")
		self.build_the_horizon_cone(vertex)
		self.link_the_horizon_cone(vertex)

	def build_the_horizon_cone(self,vertex):
		updated_horizon_edges = [] 
		for e in self.horizon_edges:
			t = Triangle(e.lf.v1, e.lf.v2, e.lf.v3, VertexR=vertex, skip_orientation=1)
			t.make_ltriangle(); t.set_edge_topology()
			updated_edges = t.get_edges()
			updated_edge = updated_edges[e.id]
			updated_sym_edge = updated_edges[e.flipid]
			updated_horizon_edges.append(updated_edge)
			self.edges.update({updated_edge.id:updated_edge})
			self.edges.update({updated_sym_edge.id:updated_sym_edge})
		self.horizon_edges = updated_horizon_edges

	def link_the_horizon_cone(self,vertex):
		"""Copy the first horizon edge in the list self.horizon_edges
		   and append it to the end of the list.  This ensures that
		   the last two faces of the horizon cone are properly linked.
		"""
		self.horizon_edges.append(self.horizon_edges[0])

		i = 0
		#A dictionary of the new cone edges.
		horizon_cone = {}
		while i < len(self.horizon_edges)-1:
			#Update the new horizon cone with the current horizon edges.
			horizon_cone.update({self.horizon_edges[i].id:
                    			 self.horizon_edges[i]})

			#Origin vertex of the current edge.
			o = self.horizon_edges[i].origin_vertex
			#Destination vertex of the current edge.
			d = self.horizon_edges[i].destination_vertex
			#The new vertex of the hull. 
			n = vertex
			#The origin vertex of the next edge in a ccw direction.
			o_ccw = self.horizon_edges[i+1].origin_vertex

			#Orientation is important.
			#From horizon edges identify new directed edge (n,o).
			#Make the right face of this new edge first; therefore,
			#change the direction of the directed edge to (o.n).
			t = Triangle(o, n, o_ccw, VertexR=d, skip_orientation=1)
			t.make_ltriangle(); t.set_edge_topology()
			self.edges.update(t.get_edges())
			#Update the new horizon cone with the newly built edges of
			#the horizon cone. 
			horizon_cone.update(t.get_edges())
			i += 1

		"""Replace the old horizon cone with the newly built horizon
		   cone.  This update is at the heart of the efficiency of
		   the calculation.  
		"""
		self.horizon_cone = horizon_cone
		"""Relatedly, update the class attribute self.start_edge with
		   the first edge of the horizon edge.  This ensures that the
		   search for an initial horizon edge is limited to the
		   edges of the current hull.  Moreover, this ensures that
		   the search for an intial horizon edge is limited to the
		   edges of the horizon cone.
		"""
		self.start_edge = self.horizon_edges[0]

	def test_count_horizon_edges_and_visible_edges(self,vertex):
		horizon_edge, visible_edge = 0, 0
		for q in self.edges.values():
			l = q.lf.beneath_beyond(vertex)
			r = q.rf.beneath_beyond(vertex)
			if l != r: 
				horizon_edge += 1
			if l == 1 and r == 1: 
				visible_edge += 1
		result = ""
		result += "\n%-14.14s% 4i %59.59s\n" % \
        		  ("TEST OUTCOME:",
        		   self.loop_count,
        		   "test_count_horizon_edges_and_visible_edges")
		result += "%-79.79s\n" % (79*"-")
		result += "%-40.40s%-39.0f\n" % ("NUMBER OF HORIZON EDGES:",
                                		 horizon_edge/2.)
		result += "%-40.40s%-39.0f\n" % ("NUMBER OF VISIBLE EDGES:",
                                		 visible_edge/2.)
		result += 2*"\n"
		return result

	def test_show_me_the_horizon_edges(self,vertex):
		result = ""
		result += "\n%-14.14s% 4i %59.59s\n" % \
        		  ("TEST OUTCOME:",
        		   self.loop_count,
        		   "test_show_me_the_horizon_edges")
		result += "%-79.79s\n" % (79*"-")
		result += "%-25s%20s%34s\n" % ("HORIZON EDGE",
                            		   "LEFT FACE VOLUME",
                            		   "RIGHT FACE VOLUME")
		for q in self.edges.values():
			l = q.lf.beneath_beyond(vertex)
			r = q.rf.beneath_beyond(vertex)
			if l != r:
				result += "%-25s%20i%34i\n" % (q,l,r) 
		result += 2*"\n"
		return result

	def test_efficiency_of_updating_the_hull(self):
		i = len(self.edges)
		efficiency = 1. - (float(self.i)/float(i))
		result = ""
		result += "\n%-14.14s% 4i %59.59s\n" % \
        		  ("TEST OUTCOME:",
        		   self.loop_count,
        		   "test_efficiency_of_updating_the_hull")
		result += "%-79.79s\n" % (79*"-")
		result += "%-40s%-39i\n" % ("NUMBER OF EDGES OF THE HULL:", i)
		result += "%-40s%-39i\n" % ("NUMBER OF HYPERPLANE CALCULATIONS:",
                            		   self.i)
		result += "%-40s%-39.2f\n" % ("EFFICIENCY RATIO:",efficiency)
		result += 2*"\n"
		return result

	def test_eulers_rule(self):
		#Euler's rule for polyhedra (i.e., d=3): V - E + F = 2.
		edges = self.edge_iterator(self.start_edge)

		nVertices = {}
		nEdges = {}
		nTriangles = {}
		for e in edges.values():
			nVertices.update({e.origin_vertex.id:e.origin_vertex})
			nVertices.update({e.destination_vertex.id:e.destination_vertex})
			nEdges.update({e.id:e})
			nTriangles.update({e.lf.id:e.lf})
			nTriangles.update({e.rf.id:e.rf})

		v = len(nVertices)
		e = len(nEdges)/2.
		t = len(nTriangles)

		# A successful test evaluates to 2.
		euler = v - e + t

		result = ""
		result += "\n%-14.14s% 4i %59.59s\n" % ("TEST OUTCOME:",
                                        		self.loop_count,
                                    		   "test_eulers_rule")
		result += "%-79.79s\n" % (79*"-")
		result += "%-40.40s%-39.0f\n" % ("NUMBER OF VERTICES OF THE HULL:",v)
		result += "%-40.40s%-39.0f\n" % ("NUMBER OF EDGES OF THE HULL:",e)
		result += "%-40.40s%-39.0f\n" % ("NUMBER OF TRIANGLES OF THE HULL:",t)
		result += "%-40.40s%-39.0f\n" % ("EULER\'S RULE:",v - e + t)
		result += 2*"\n"
		return result

	def rotate_ccw_over_the_right_face(self,start_edge):
		rotate = [start_edge,]
		next_edge = start_edge.Rnext
		n = 0
		while next_edge != start_edge.id:
			if n%2:
				next_edge = self.edges[next_edge].flipid
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Dprev
			else:
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Rnext
			n += 1
		return rotate

	def rotate_ccw_over_the_left_face(self,start_edge):
		rotate = [start_edge,]
		next_edge = start_edge.Lnext
		n = 0
		while next_edge != start_edge.id:
			if n%2:
				next_edge = self.edges[next_edge].flipid
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Oprev
			else:
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Lnext
			n += 1
		return rotate

	def rotate_cw_over_the_right_face(self,start_edge):
		rotate = [start_edge,]
		next_edge = start_edge.Rprev
		n = 0
		while next_edge != start_edge.id:
			if n%2:
				next_edge = self.edges[next_edge].flipid
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Onext
			else:
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Rprev
			n += 1
		return rotate

	def rotate_cw_over_the_left_face(self,start_edge):
		rotate = [start_edge,]
		next_edge = start_edge.Lprev
		n = 0
		while next_edge != start_edge.id:
			if n%2:
				next_edge = self.edges[next_edge].flipid
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Dnext
			else:
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Lprev
			n += 1
		return rotate

	def rotate_ccw_about_the_origin_vertex(self,start_edge):
		#Rotate ccw over the edges that share the same origin vertex.
		rotate = [start_edge,]
		next_edge = start_edge.Onext
		n = 0
		while next_edge != start_edge.id:
			rotate.append(self.edges[next_edge])
			next_edge = self.edges[next_edge].Onext
		rotate.append(start_edge)
		return rotate

	def rotate_ccw_about_the_destination_vertex(self,start_edge):
		#Rotate ccw over the edges that share the same destination vertex.
		rotate = [start_edge,]
		next_edge = start_edge.Dnext
		n = 0
		while next_edge != start_edge.id:
			rotate.append(self.edges[next_edge])
			next_edge = self.edges[next_edge].Dnext
		return rotate

	def rotate_cw_about_the_origin_vertex(self,start_edge):
		#Rotate cw over the edges that share the same origin vertex.
		rotate = [start_edge,]
		next_edge = start_edge.Oprev
		n = 0
		while next_edge != start_edge.id:
			rotate.append(self.edges[next_edge])
			next_edge = self.edges[next_edge].Oprev
		rotate.append(start_edge)
		return rotate

	def rotate_cw_about_the_destination_vertex(self,start_edge):
		#Rotate cw over the edges that share the same destination vertex.
		rotate = [start_edge,]
		next_edge = start_edge.Dprev
		n = 0
		while next_edge != start_edge.id:
			rotate.append(self.edges[next_edge])
			next_edge = self.edges[next_edge].Dprev
		return rotate

	def edge_iterator(self,start_edge): # XXXX resolve with getEdges
		"""If the unique_ids of two vertices match, this code will
		   enter an infinite loop."""
		self.encountered = {}
		self.encountered.update({start_edge.id:self.edges[start_edge.id]})
		self.encountered.update({start_edge.flipid:self.edges[start_edge.flipid]})
		self.get_the_edges_of_the_hull(start_edge)
		edges = self.encountered
		self.encountered = {}
		#print("NUMBER OF ACCESSED EDGES: ", self.iterator_count)
		return edges

	# 2012.10.09: Legacy function for collecting the edges of the hull.
	#             I will keep for now because of dependencies elsewhere. 
	#             I will eventually replace this call with getEdges.
	def get_the_edges_of_the_hull(self,start_edge): # XXXX
		e = start_edge.Rnext
		while e != start_edge.id:
			#self.iterator_count += 1
			e     = self.edges[e].Rnext
			esym  = self.edges[e].flipid
			if e not in self.encountered.keys():
				self.encountered.update({e:self.edges[e]})
				self.encountered.update({esym:self.edges[esym]})               
				self.get_the_edges_of_the_hull(self.edges[e])
				self.get_the_edges_of_the_hull(self.edges[esym])

	def get_the_vertices_of_the_hull(self):

		# Map all of the vertices used to calculate the convex hull.
		all_verts = {}
		for v in self.non_mutated_vertices:
			all_verts.update({v.id:v})

		# Map the vertices of the convex hull.
		temp_verts_of_hull = {}
		edges = self.edge_iterator(self.start_edge)
		for e in edges.values():
			ov = e.origin_vertex
			dv = e.destination_vertex
			temp_verts_of_hull.update({ov.id:ov, dv.id:dv})
			if ov.id in all_verts: del all_verts[ov.id]
			if dv.id in all_verts: del all_verts[dv.id]

		keys = sorted(temp_verts_of_hull)
		verts_of_hull = []
		for key in keys:
			verts_of_hull.append(temp_verts_of_hull[key])
		return verts_of_hull 

	def get_the_vertices_inside_of_the_hull(self):

		# hull_verts = self.get_the_vertices_of_the_hull()
		# inside_verts = {}
		# #for v in self.list_of_vertices:
		# for v in self.non_mutated_vertices:
		# 	if v not in hull_verts:
		# 		inside_verts.update({v.id:v})

		inside_verts = {}
		for v in self.list_of_vertices:
			result = self.test_if_this_vertex_is_in_the_hull(v)
			# print(result, v.data)
			if result:
				inside_verts.update({v.id:v})
				# print(v.data)
		# for v in self.non_mutated_vertices:
		# 	if v not in hull_verts:
		# 		inside_verts.update({v.id:v})

		keys = sorted(inside_verts)
		verts_inside_of_hull = []
		for key in keys:
			verts_inside_of_hull.append(inside_verts[key])
		return verts_inside_of_hull

	def empty_hull(self,vertex_list):
		for vertex in vertex_list:
			try:
				if self.test_if_this_vertex_is_in_the_hull(vertex):
					"""Return false if the hull is empty"""
					return 0
			except:
				print("\n1--The test vertex is not in general position",)
				print("relative to the hull. " )
				print("2--Attempt to jostle the point.")
				vertex.jostle()
				#Try to jostle the point.
				#from random import random
				#vertex.x += (random()/5.)
				#vertex.y += (random()/5.)
				#vertex.z += (random()/5.)
				try:
					print("3--Attempt to jostle was a success.")
					if self.test_if_this_vertex_is_in_the_hull(vertex):
						return 0
				except:
					print("3--The test vertex is STILL not in general position",)
					print("relative to the hull." )
					print("4--Returning false.\n")
					return 0
		"""Return true if the hull is empty"""
		return 1

	# 2012.10.09: Legacy function for collecting the facets of the hull.
	#             I will keep for now because of dependencies elsewhere. 
	#             I will eventually replace this call with getFacets.
	def get_the_facets_of_the_hull(self):
		facets = {}
		edges = self.edge_iterator(self.start_edge)
		for e in edges.values():
			leftFacetId = list(e.lf.id)
			leftFacetId.sort()
			rightFacetId = list(e.rf.id)
			rightFacetId.sort()
			facets.update({tuple(leftFacetId):e.lf})
			facets.update({tuple(rightFacetId):e.rf})
		return facets.values()

	# 2012.10.09: New, more concise class method call. 
	def getEdges(self):
		return self.edge_iterator(self.start_edge).values()

	# 2012.10.09: New, more concise class method call.
	def getFacets(self):
		if self.facets:
			return self.facets
		else:
			facets = {}
			edges = self.edge_iterator(self.start_edge)
			for e in edges.values():
				leftFacetId = list(e.lf.id)
				leftFacetId.sort()
				rightFacetId = list(e.rf.id)
				rightFacetId.sort()
				facets.update({tuple(leftFacetId):e.lf})
				facets.update({tuple(rightFacetId):e.rf})
			self.facets = facets.values()
			return self.facets

	def getHullAsString(self):
		facets, facetString = self.getFacets(), ""
		for facet in facets:
			fString = "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n"
			facetString += fString % (facet.v1.x, facet.v1.y, facet.v1.z,facet.v2.x, facet.v2.y, facet.v2.z, facet.v3.x, facet.v3.y, facet.v3.z) 
		return facetString


	def getHyperPlanes(self):
		if self.hyperPlanes:
			return self.hyperPlanes
		else:
			#from compGeometry import planeCoefficients3D
			hyperPlanes, facets = {}, self.getFacets()
			for facet in facets:
				c = planeCoefficients3D(facet.v1,facet.v2,facet.v3)
				hyperPlanes.update({facet:c})
			self.hyperPlanes = hyperPlanes
			return self.hyperPlanes

	# def test_if_this_vertex_is_in_the_hull(self,v):
	# 	inside, hyperPlanes = 1, self.getHyperPlanes()
	# 	for facet in hyperPlanes:
	# 		hP = hyperPlanes[facet]
	# 		sign = hP[0]*v.x + hP[1]*v.y + hP[2]*v.z + hP[3]
	# 		# The sign convention is the opposite of what it seems 
	# 		# because of the manner in which the topology of the hull
	# 		# facets has been arranged. 
	# 		if sign > 0: # if sign < 0: # original, but likely wrong topology
	# 			inside = 0
	# 			break
	# 	return inside

	def test_if_this_vertex_is_in_the_hull(self, v, returnVolume=0):

		# Only need to test beneath_beyond for the left face because all edges (i.e., the e.flipid edges) are accessed in the loop.
		# In other words, in the quad edge, what is left face for e.id is right face for e.flipid, and vice versa. 
		###########################################################################################################################
		minVolume = 100000.
		edges = self.getEdges()
		for e in edges:
			l = e.lf.beneath_beyond(v)
			if l > 0: 
				return 0
			elif l < 0 and returnVolume:
				volume = e.lf.beneath_beyond(v, returnVolume=returnVolume)
				if volume < minVolume:
					minVolume = volume

		if returnVolume: 
			return minVolume
		else:
			return 1

	def modulateHull(self, multiplier=1, swell=1, shrink=0):

		from pdbFile import PseudoAtom
		from compGeometry import centroid3D, crossProductCoordinates

		edges = self.getEdges()
		swellFacetsString = ""
		vs = {}
		for edge in edges:

			e = edge
			psa = PseudoAtom()
			v1, v2, v3 = e.rf.v1, e.rf.v2, e.rf.v3
			c = centroid3D([v1, v2, v3], psa)

			psa1 = PseudoAtom()
			if swell:
				psaResult = crossProductCoordinates(v2, psa.v, v1, psa1, multiplier=multiplier, norm=1)
				psaResult.atom_serial_number = v2.data.atom_serial
				psaResult.v.id = v2.data.atom_serial
			else:
				psaResult = crossProductCoordinates(v1, psa.v, v2, psa1, multiplier=multiplier, norm=1)
				psaResult.atom_serial_number = v1.data.atom_serial
				psaResult.v.id = v1.data.atom_serial
			vs[psaResult.v.id] = psaResult.v

		# Calculates and returns a new convex hull instance...
		vs = list(vs.values())
		return convexHull3D(vs)

	def get_LO(self):
		edges = self.edge_iterator(self.start_edge)
		result_class = ConvexHull3D_LO(edges)
		result_class.of_hull = self.get_the_vertices_of_the_hull()
		result_class.in_hull = self.get_the_vertices_inside_of_the_hull()
		return result_class

class ConvexHull3D_LO:

    #THE SIMPLEST VERSION OF A CONVEX HULL INSTANCE THAT REDUCES OVERHEAD.
    
    def __init__(self,edges):
        self.edges   = edges
        self.of_hull = None
        self.in_hull = None

    def test_if_this_vertex_is_in_the_hull(self,vertex):
        for e in self.edges.values():
            l = e.lf.beneath_beyond(vertex)
            r = e.rf.beneath_beyond(vertex)
            if l+r != -2 and l+r != 0:
                return 0
        return 1



        
