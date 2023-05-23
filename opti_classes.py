# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 22:23:39 2023

@author: aangre
"""

import numpy as np
import pandas as pd
import itertools

# Vertex class
class Vertex(object):

    def __init__(self, vertex_id, coordinates):
        self.id = vertex_id         
        self.coordinates = coordinates   
        self.adjacent_vertices = []
        self.connected_edges = []
        
    def __repr__(self):
        rep = 'Vertex(id=' + str(self.id) + ',coords=' + str(self.coords) +')'
        return rep
    
# Edge class
class Edge(object):
    
    def __init__(self, edge_id, vertices, edge_width):
        self.id = edge_id
        self.edge_vertices = vertices
        self.parent_loops = set()
        self.target = edge_width
        
        # Set edge vertices in an increasing order of their ids
        if (vertices[0].id > vertices[1].id):
            self.edge_vertices = [vertices[1], vertices[0]]
        
    def __repr__(self):
        rep = 'Edge(id=' + str(self.id) + ',vertices_ids=(' + str(self.edge_vertices[0].id) + ',' \
                + str(self.edge_vertices[1].id) + '))'
        return rep

# Connection class
class Connection(object):
    
    def __init__(self, connection_id, edges):
        self.id = connection_id
        self.connection_edges = edges
        self.connection_edge_ids = [edges[0].id, edges[1].id]
        self.target = max([x.target for x in self.connection_edges])
              
    def __repr__(self):
        rep  = 'Connection(id=' + str(self.id) + ',edges_ids=(' + str(self.connection_edges[0].id) + ',' \
                + str(self.connection_edges[1].id) + '), target='+ str(self.target) +'))'
        return rep
        
# Loop class
class Loop(object):
    
    def __init__(self, loop_id, edges):
        self.id = loop_id
        self.loop_vertex_ids = []
        self.loop_edges = edges
        self.loop_edge_ids = [self.loop_edges[i].id for i in range(len(self.loop_edges))]
        self.target = max([x.target for x in self.loop_edges])
        
        # Set loop vertices ids
        loop_vertex_ids = set()
        for i in range(len(self.loop_edges)):
            loop_vertex_ids.add(self.loop_edges[i].edge_vertices[0].id)
            loop_vertex_ids.add(self.loop_edges[i].edge_vertices[1].id)
        self.loop_vertex_ids = list(loop_vertex_ids)
        
    def __repr__(self):
        rep = 'Loop(id=' + str(self.id) + ',edges_ids=('
        for edge_id in self.loop_edge_ids:
            rep += str(edge_id) + ','
        rep = rep[:-1] + ')'
        return rep    

def createNewSheet(xlFilePath, sheetName):
    g = Graph()
    g.reset()
    g.fill_vertices(xlFilePath)
    g.fill_edges(xlFilePath)
    g.fill_loops(xlFilePath, sheet=sheetName)
    return g

class Graph(object):
    """ Graph implements undirected graph with edge and vertex extra data"""

    def __init__(self):
        self.vertices = {}
        self.edges = {}
        self.loops = {}
        self.connections = {}
        
    def reset(self):
        self.__init__()
                
    def add_vertex_to_graph(self, vertex_id, coords):
        if vertex_id not in self.vertices:
            x_coord = coords[0]
            y_coord = coords[1]
            vertex = Vertex(vertex_id, [x_coord, y_coord])
            self.vertices[vertex_id] = vertex
            return vertex
        else:
            return self.vertices[vertex_id]
        
    def add_edge_to_graph(self, edge_id, vertices_ids, width_initial):
        if edge_id not in self.edges:
            if (vertices_ids[0] in self.vertices.keys()) and (vertices_ids[1] in self.vertices.keys()):
                edge = Edge(edge_id, [self.vertices[vertices_ids[0]], self.vertices[vertices_ids[1]]], width_initial)
                edge.edge_vertices[0].adjacent_vertices.append(edge.edge_vertices[1])
                edge.edge_vertices[0].connected_edges.append(edge)
                edge.edge_vertices[1].adjacent_vertices.append(edge.edge_vertices[0])
                edge.edge_vertices[1].connected_edges.append(edge)
                self.edges[edge_id] = edge
                return edge
            else:
                raise ValueError("Missing vertices! Please add missing vertices to the graph.")
        else:
            return self.edges[edge_id]
                
    def add_loop_to_graph(self, loop_id, edges_ids):
        if loop_id not in list(self.loops.keys()):
            if(set(edges_ids).issubset(set(self.edges.keys()))):
                loop_edges = [self.edges[edge_id] for edge_id in edges_ids]
                loop = Loop(loop_id, loop_edges)
                self.loops[loop_id] = loop
                # Add connections to the graph
                for i in range(len(edges_ids)-1):
                    self.edges[edges_ids[i]].parent_loops.add(loop_id) 
                    self.edges[edges_ids[i+1]].parent_loops.add(loop_id) 
                return loop
            else:
                raise ValueError("Missing edges! Please add missing edges to the graph.")
        else:
            return self.loops[loop_id]   
        
    def define_connections(self):
        # Check if vertices and edges exist in the graph
        if (bool(self.vertices) and bool(self.edges)):
            connection_edge_ids_list = []
            for vertex in list(self.vertices.values()):
                vertex_connected_edges = vertex.connected_edges
                vertex_connected_edge_ids = list(map(lambda x : x.id, list(vertex_connected_edges)))
                connection_combs = list(set(list(itertools.combinations(vertex_connected_edge_ids, 2))))
                connection_edge_ids_list.extend(connection_combs)
            
            # Setting up all connections
            connection_edge_ids_list = sorted(connection_edge_ids_list)
            for i in range(len(connection_edge_ids_list)):
                connection_edges = [self.edges[connection_edge_ids_list[i][0]], self.edges[connection_edge_ids_list[i][1]]]
                connection = Connection(i, connection_edges)
                self.connections[i] = connection
        else:
            raise ValueError("Vertices or edges do not exist in the graph.")
        
    def get_connection_loop_relations(self, Fw):
        # Extract connections and loops list
        loops_list = [self.loops[x] for x in sorted(self.loops.keys())]
        connections_list = [self.connections[x] for x in sorted(self.connections.keys())]

        # Initialize coefficient(2D) and target(1D) array
        A =  np.zeros((len(connections_list), len(loops_list)))
        Y1 =  np.zeros(len(connections_list))
                
        # store the coefficients and targets in the arrays
        for i in range(len(connections_list)):
            connection_edges = connections_list[i].connection_edge_ids
            for j in range(len(loops_list)):
                loop_edge_ids = loops_list[j].loop_edge_ids
                loop_conns = [sorted([loop_edge_ids[i],loop_edge_ids[i+1]]) for i in range(len(loop_edge_ids)-1)]
                if sorted(connection_edges) in loop_conns:
                    A[i,j] = 1
                Y1[i] = connections_list[i].target/Fw
                
        # Return connection-loops relation matrix and connection target vector
        return A, Y1
    
    def get_edge_loop_relations(self, Fw):
        # Extract edges and loops list
        loops_list = [self.loops[x] for x in sorted(self.loops.keys())]
        edges_list = [self.edges[x] for x in sorted(self.edges.keys())]

        # Initialize coefficient(2D) and target(1D) array
        B =  np.zeros((len(edges_list), len(loops_list)))
        Y2 =  np.zeros(len(edges_list))
                
        # store the coefficients and targets in the arrays
        for i in range(len(edges_list)):
            for j in range(len(loops_list)):
                loop_edge_ids = loops_list[j].loop_edge_ids
                if edges_list[i].id in loop_edge_ids:
                    B[i,j] = 1
                Y2[i] = edges_list[i].target/Fw
                
        # Return edge-loops relation matrix and edge target vector
        return B, Y2
        
    def get_vertex_loop_relations(self, Fw):
        # Extract edges and loops list
        loops_list = [self.loops[x] for x in sorted(self.loops.keys())]
        vertices_list = [self.vertices[x] for x in sorted(self.vertices.keys())]
        
        # Initialize coefficient(2D) and target(1D) array
        V =  np.zeros((len(vertices_list), len(loops_list)))
        Y3 =  np.zeros(len(vertices_list))
                
        # store the coefficients and targets in the arrays
        for i in range(len(vertices_list)):
            for j in range(len(loops_list)):
                loop_vertex_ids = loops_list[j].loop_vertex_ids
                if vertices_list[i].id in loop_vertex_ids:
                    V[i,j] = 1
                #Should be some criteria constraining no of loops through vertex
                Y3[i] = 1.5*max([edge.target for edge in vertices_list[i].connected_edges])/Fw
                
        # Return vertex-loops relation matrix and vertex target vector
        return V, Y3
    
    def fill_vertices(self, xlFilePath):
        df = pd.read_excel(xlFilePath, sheet_name="Vertex")
        vertex_ids = np.array([df.loc[:,"id"].at[i] for i in range(df.shape[0])])
        xcoords = np.array([df.loc[:,"xcoord"].at[i] for i in range(df.shape[0])])
        ycoords = np.array([df.loc[:,"ycoord"].at[i] for i in range(df.shape[0])])
        
        # Fill the graph with vertices
        for i in range(len(vertex_ids)):
            coords = [xcoords[i], ycoords[i]]
            self.add_vertex_to_graph(vertex_ids[i], coords)
        
        # return vertex_ids, xcoords, ycoords
        
    def fill_edges(self, xlFilePath):
        df = pd.read_excel (xlFilePath, sheet_name="Edge")
        edge_ids = np.array([df.loc[:, "id"].at[i] for i in range(df.shape[0])])
        vertex1_ids = np.array([df.loc[:, "vertex1"].at[i] for i in range(df.shape[0])])
        vertex2_ids = np.array([df.loc[:, "vertex2"].at[i] for i in range(df.shape[0])])
        edge_targets = np.array([df.loc[:, "target"].at[i] for i in range(df.shape[0])])
        
        # Fill the graph with edges
        for i in range(len(edge_ids)):
            vertices_ids = [vertex1_ids[i], vertex2_ids[i]]
            self.add_edge_to_graph(edge_ids[i], vertices_ids, edge_targets[i])
        
        # return edge_ids, vertex1_ids, vertex2_ids, edge_targets

    def fill_loops(self, xlFilePath, sheet):
        df = pd.read_excel (xlFilePath, sheet_name=sheet)
        loop_ids = np.array([df.loc[:, "id"].at[i] for i in range(df.shape[0])])
        loop_edge_ids = []
        for i in loop_ids:
            string_edge_ids = str(df.loc[:,"edge_ids"].at[i])
            string_edge_ids = string_edge_ids.replace(" ", "").split(",")
            string_edge_ids = [int(x) for x in string_edge_ids]
            loop_edge_ids.append(string_edge_ids)
        maxLen = max(map(len, loop_edge_ids))
        loop_edge_ids = np.array([i+[np.NaN]*(maxLen-len(i)) for i in loop_edge_ids])
        
        # Define loops in the graph
        for i in range(len(loop_ids)):
            edge_ids = list(loop_edge_ids[i,:np.count_nonzero(~np.isnan(loop_edge_ids[i,:]))])
            self.add_loop_to_graph(loop_ids[i], edge_ids)
            
        # Define connections in the graph
        self.define_connections()
        
        # return loop_edge_ids, loop_targets 
    