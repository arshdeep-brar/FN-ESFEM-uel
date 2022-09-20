
"""

@author: arshdeep
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import math

class ESFEM_Mesh:

    def __init__(self, NodesCoord, EdgeNodeConn, EdgeNodes):

        self.NodesCoord = NodesCoord
        self.EdgeNodeConn = EdgeNodeConn
        self.EdgeNodes = EdgeNodes
        self.totNodes = max(list(self.NodesCoord.keys()))
        self.Nelem = max(list(self.EdgeNodeConn.keys()))
        self.totCentroids = None
        self.Centroids = {}
        self.Points = {}
        self.Edges = {}
        self.SmoothedDomain = {}
        self.normals = {}
        self.IntegrationPoint = {}

    def createSmoothedDomain(self, edge):

        elem = self.EdgeNodeConn[edge]

        if len(elem) == 3:
            Centroid = (1/3) * (np.array(self.NodesCoord[elem[0]]) + np.array(self.NodesCoord[elem[1]]) + \
                                np.array(self.NodesCoord[elem[2]]))
            
            smootheddomain = [self.NodesCoord[elem[0]], list(Centroid), self.NodesCoord[elem[1]]] 
        
        if len(elem) == 4:
            Centroid1 = (1/3) * (np.array(self.NodesCoord[elem[0]]) + np.array(self.NodesCoord[elem[1]]) + \
                                np.array(self.NodesCoord[elem[2]]))
            Centroid2 = (1/3) * (np.array(self.NodesCoord[elem[0]]) + np.array(self.NodesCoord[elem[2]]) + \
                                np.array(self.NodesCoord[elem[3]]))

            smootheddomain = [self.NodesCoord[elem[0]], list(Centroid1), self.NodesCoord[elem[2]], list(Centroid2)]

        return np.array(smootheddomain)

    def Normal(self, line, rot):

        if rot == 'clk':
            length = np.sqrt((line[0,0] - line[1,0])**2 + (line[0,1] - line[1,1])**2) 
            nx = -(1/length) * (line[1,1] - line[0,1])
            ny = (1/length) * (line[1,0] - line[0,0])
        if rot == 'aclk':
            length = np.sqrt((line[0,0] - line[1,0])**2 + (line[0,1] - line[1,1])**2) 
            nx = (1/length) * (line[1,1] - line[0,1])
            ny = -(1/length) * (line[1,0] - line[0,0])

        return np.array([nx, ny])

    def MidPt(self, line):

        return 0.5 * (line[0,:] + line[1,:])

        

    def SD_segments(self, smootheddomain):

        if len(smootheddomain) == 3:
            rot = 'clk'
            line1 = smootheddomain[0:2, :]
            line2 = smootheddomain[1:3, :]
            line3 = smootheddomain[[2,0], :]
            linesegments = (rot, line1, line2, line3)
        
        if len(smootheddomain) == 4:
            rot = 'aclk'
            line1 = smootheddomain[0:2, :]
            line2 = smootheddomain[1:3, :]
            line3 = smootheddomain[2:4, :]
            line4 = smootheddomain[[3,0], :]
            linesegments = (rot, line1, line2, line3, line4)
        
        return linesegments

    def Triangles(self, edge):

        if len(self.EdgeNodeConn[edge]) == 3:
            Nodes = self.EdgeNodeConn[edge]
            Triangles = (np.array([self.NodesCoord[Nodes[0]],
                                 self.NodesCoord[Nodes[1]],
                                 self.NodesCoord[Nodes[2]]]))

        if len(self.EdgeNodeConn[edge]) == 4:
            Nodes = self.EdgeNodeConn[edge]
            Triangles = (np.array([self.NodesCoord[Nodes[0]],
                                 self.NodesCoord[Nodes[1]],
                                 self.NodesCoord[Nodes[2]]]),
                         np.array([self.NodesCoord[Nodes[0]],
                                 self.NodesCoord[Nodes[2]],
                                 self.NodesCoord[Nodes[3]]]))
        return Triangles
    
    def ShapeFunctions(self, X, Triangle):

        Node1 = Triangle[0,:]
        Node2 = Triangle[1,:]
        Node3 = Triangle[2,:]

        A = (1/2) * np.abs(np.linalg.det(np.array([[1, Node1[0], Node1[1]],
                                                   [1, Node2[0], Node2[1]],
                                                   [1, Node3[0], Node3[1]]])))
        
        N1 = 1/(2*A) * ((Node2[0]*Node3[1] - Node3[0]*Node2[1]) + (Node2[1] - Node3[1]) * X[0] + \
                                    (Node3[0] - Node2[0]) * X[1]) 
    
        N2 = 1/(2*A) * ((Node3[0]*Node1[1] - Node1[0]*Node3[1]) + (Node3[1] - Node1[1]) * X[0] + \
                                    (Node1[0] - Node3[0]) * X[1])

        N3 = 1/(2*A) * ((Node1[0]*Node2[1] - Node2[0]*Node1[1]) + (Node1[1] - Node2[1]) * X[0] + \
                                    (Node2[0] - Node1[0]) * X[1])

        return np.array([N1, N2, N3])
    
    def createmesh(self):

        nodeID = max(list(self.NodesCoord.keys())) + 1

        for key in self.EdgeNodes.keys():
            
            self.Edges[key] = np.array([self.NodesCoord[self.EdgeNodes[key][0]], 
                                        self.NodesCoord[self.EdgeNodes[key][1]]])

            thisdomain = self.createSmoothedDomain(edge=key)
            
            self.SmoothedDomain[key] = thisdomain

            segments = self.SD_segments(thisdomain)

            if len(segments) == 4:
                rot, line1, line2, line3 = segments 
                self.normals[key] = (self.Normal(line1, rot), self.Normal(line2, rot),
                                        self.Normal(line3,rot))
                self.IntegrationPoint[key] = (self.MidPt(line1), self.MidPt(line2), 
                                            self.MidPt(line3)) 
                
                if isValue(self.Centroids, thisdomain[1,:]):
                    pass
                else:
                    self.Centroids[nodeID] = list(thisdomain[1,:])
                    nodeID = nodeID + 1  
                
            
            if len(segments) == 5:
                rot, line1, line2, line3, line4 = segments 
                self.normals[key] = (self.Normal(line1, rot), self.Normal(line2, rot),
                                      self.Normal(line3,rot), self.Normal(line4, rot))

                self.IntegrationPoint[key] = (self.MidPt(line1), self.MidPt(line2), 
                                            self.MidPt(line3), self.MidPt(line4))
                
                if isValue(self.Centroids, thisdomain[1,:]):
                    pass
                else: 
                    self.Centroids[nodeID] = list(thisdomain[1,:])
                    nodeID = nodeID + 1

                if isValue(self.Centroids, thisdomain[3,:]):
                    pass
                else: 
                    self.Centroids[nodeID] = list(thisdomain[3,:])
                    nodeID = nodeID + 1

            self.Points = {**self.NodesCoord, **self.Centroids}

            self.totCentroids = len(list(self.Centroids.keys()))
        return

    def ElementTopology(self):

        Elem_topo = {}

        for key in self.SmoothedDomain.keys():

            thesepoints = []
            thisdomain = self.SmoothedDomain[key]

            for i in range(len(thisdomain)):
                point = findkey(self.Points, list(thisdomain[i,:]))
                thesepoints.append(point)
            
            Elem_topo[key] = thesepoints
        
        return Elem_topo
    

    def PlotMesh(self, SD=[], label_nodes=False, label_edge=False, label_centroid=False, save=False, zoom=np.array([])):

        plt.axes()    

        for key in self.Edges.keys():
            Edge = plt.Line2D(self.Edges[key][:,0], self.Edges[key][:,1], lw=2)
            if label_edge:
                textpos = 0.5 * np.sum(self.Edges[key], axis=0)
                plt.text(textpos[0], textpos[1], '{:d}'.format(key), fontsize='small', ha='center', va='center', 
                            color='red', bbox=dict(facecolor='white', edgecolor='red', pad=0.5))
            plt.gca().add_line(Edge)

            if key in SD:
                
                Domain = plt.Polygon(self.SmoothedDomain[key], fc='springgreen', ec='black', ls='--')
                plt.gca().add_patch(Domain)
                
                textpos = (1/len(self.SmoothedDomain[key])) * np.sum(self.SmoothedDomain[key], axis=0) 
                
                plt.text(textpos[0], textpos[1], '{:d}'.format(key), fontsize='small', ha='center', va='center')

                for id in range(len(self.SmoothedDomain[key])):

                    x = self.IntegrationPoint[key][id][0]
                    y = self.IntegrationPoint[key][id][1]

                    dx = self.normals[key][id][0]
                    dy = self.normals[key][id][1]

                    plt.arrow(x, y, dx, dy, color='red', head_width=0.4) 

                Nodes_associated = np.zeros((len(self.EdgeNodeConn[key]), 2))
                
                i = 0
                for node in self.EdgeNodeConn[key]:
                    Nodes_associated[i, :] = self.NodesCoord[node]
                    i = i + 1
                
                plt.scatter(Nodes_associated[:,0], Nodes_associated[:,1], marker='o', color='red')

                
            else:
                Domain = plt.Polygon(self.SmoothedDomain[key], fill=None, ls='--')
                plt.gca().add_patch(Domain)

        if label_nodes:
            for node in self.NodesCoord.keys():
                plt.text(self.NodesCoord[node][0], self.NodesCoord[node][1], '{:d}'.format(node), 
                            ha='center', va='center', fontsize='small', 
                            bbox=dict(facecolor='white', edgecolor='black', pad=0.5))

        if label_centroid:
            for point in self.Centroids.keys():
                plt.text(self.Centroids[point][0], self.Centroids[point][1], '{:d}'.format(point), 
                            ha='center', va='center', fontsize='small', 
                            bbox=dict(facecolor='white', edgecolor='red', pad=0.5))    
        
        xmin, xmax, ymin, ymax = plt.axis('scaled')
        
        if len(zoom) == 0:
            plt.xlim(xmin-5, xmax+5)
            plt.ylim(ymin-5, ymax+5)
        else:
            plt.xlim(zoom[0,0], zoom[0,1])
            plt.ylim(zoom[1,0], zoom[1,1])

        if save:
            if len(SD) == 0:
                filename = 'ESFEM_mesh.png'
            else:
                filename = 'ESFEM_SM_'
                for elem in SD:
                    filename = filename + '{:d}'.format(elem)

            plt.savefig(filename, dpi=300)
        
        plt.show()

        return

##################### Required Functions #########################

def Edgedef(ElemNodes):
    '''
    The function generates the Edges for a given element

    Parameters
    ----------
    ElemNodes : ARRAY
        Connectivity array for the an element.

    Returns
    -------
    Edges : ARRAY OF TUPLES
        An Array containing all the edges with its nodes in Tuple

    '''
    Edges = [] #Initiating an empty list for edges
    
    for i in range(len(ElemNodes)):
        Node1 = ElemNodes[i]
        
        if (i+1) == len(ElemNodes):
            Node2 = ElemNodes[0]    
        else:
            Node2 = ElemNodes[i+1]
        
        Edges.append((Node1, Node2))
    
    return Edges

def EdgeConnectivity(ElemNodeConn):
    '''
    Generates a dictionay containing edge number as keys and Element number
    connected to that edge as a list 

    Parameters
    ----------
    ElemNodeConn : DICT
        Element connectivity

    Returns
    -------
    EdgeConn : DICT
        Edge number: Element number.
    EdgeNumber : DICT
        Edge number: Node Definition.

    '''

    ElemEdges = {}

    for i in ElemNodeConn.keys():
        ElemNodes = ElemNodeConn[i]
        ElemEdges[i] = Edgedef(ElemNodes)

    EdgeConn = {}  #Edge connectivity Dictionary
    
    EdgeNumber = {} #Edge difinition with nodes
    
    AlternateEdge = {} #Edge definition in opposite node direction
    
    ElemEdgeConn = {}

    i = 1
    
    for key in ElemEdges.keys():
        Edges = ElemEdges[key]
        Edgelst = []
        for Edge in Edges:
            
            print(Edge)
            #Check if the edge already exist in Edge Connectivity 
            if Edge in AlternateEdge.values():
                Edgeindex = findkey(AlternateEdge, Edge)
                EdgeConn[Edgeindex].append(key)
                Edgelst.append(-Edgeindex)
            
            else:
                EdgeConn[i] = []
                EdgeConn[i].append(key)
                AlternateEdge[i] = (Edge[1], Edge[0])
                EdgeNumber[i] = Edge
                Edgelst.append(i)
                i += 1
        ElemEdgeConn[key] = Edgelst

    return EdgeConn, ElemEdgeConn, EdgeNumber

def ESFEM_Nodes(ElemNodeConn, EdgeConn, EdgeNumber):
    '''
    Calculates the Edge Node Connectivity for the smoothed domain 
    associated to the edge.

    Parameters
    ----------
    ElemNodeConn: DICT
                Element Node Connectivity of a triangular mesh
    EdgeConn: DICT
            Edge connectivity - all the elements associated to the
            edge
    
    Returns
    -------
    EdgeNodeConn: DICT
                Edge Node Connectivity of the smoothed domain
    '''

    EdgeNodeConn = {}

    for edge in EdgeConn.keys():
        elems = EdgeConn[edge]
        rawNodelist = []
        for elem in elems:
        
            nodes = ElemNodeConn[elem]
            
            for node in nodes:
                
                if node in rawNodelist:
                    pass

                else:
                    rawNodelist.append(node)

        Node1, Node2 = EdgeNumber[edge]
        Nodelist = []
        Nodelist.append(Node1)
        Nodelist.append(Node2)

        if len(rawNodelist) == 3:
            for node in rawNodelist:
                if node in Nodelist:
                    pass
                else:
                    Nodelist.append(node)

        else:
            Nodelist.insert(1, rawNodelist[-1])
            for node in rawNodelist:
                if node in Nodelist:
                    pass
                else:
                    Nodelist.append(node)

        EdgeNodeConn[edge] = Nodelist

    return EdgeNodeConn  
        
def findkey(dict, val):
    for key, value in dict.items():
        if np.isclose(val, value).all() :
            return key

def isValue(dict, val):
    for value in dict.values():
        if np.isclose(val, value).all():
            return True
    
    return False


def findArea(NodesCoord, ElemNodeConn):
    '''
    Calculates the area of the triangle from the coordinates of the
    nodes

    Parameters
    ----------
    NodesCoord : DICT
            Stores the coordinates of the node with node number as the
            key

    ElemNodeConn : DICT
            Connectivity matrix for an element    

    Returns
    -------
    Area : DICT
        Area of the Element
    '''
    
    Area = {}

    for key in ElemNodeConn.keys():
        Nodes = ElemNodeConn[key] 
        temp = np.array([[NodesCoord[Nodes[0]][0], NodesCoord[Nodes[0]][1], 1],
                         [NodesCoord[Nodes[1]][0], NodesCoord[Nodes[1]][1], 1],
                         [NodesCoord[Nodes[2]][0], NodesCoord[Nodes[2]][1], 1]])
        Area[key] = (1/2) * abs(np.linalg.det(temp))
    
    return Area
    
def ReadInpfile(filename):
    '''
    Function to read the ABAQUS input file

    Parametres:
    ----------
    filename: STRING 
        Name of the ABAQUS input file in '*.inp' 
    
    '''

    #Opening the input file 
    with open(filename, 'r') as Inpfile:
        #Raw data collected as string
        ElemNodeConn_raw = []

        #Flag to be raised when the in data for Element
        # node connectivity starts
        ElemFlag = False
        
        for line in Inpfile:        
        
            if '*' in line and ElemFlag:
                break
            
            if ElemFlag:
                ElemNodeConn_raw.append(line[:-1])
            
            if '*Element' in line:
                ElemFlag = True
        
    ElemNodeConn = {}

    #Converting the raw data to a appropriate data type and storing it
    # in a dictionary
    for i in range(len(ElemNodeConn_raw)):
        temp = [int(s) for s in ElemNodeConn_raw[i].split(',')]
        ElemNodeConn[temp[0]] = temp[1:]

    # Opening the input file for Node Coordinates
    with open(filename, 'r') as Inpfile:
        #Raw data collected as string
        NodeCoord_raw = []
        
        #Flag to be raised when the in data for Node 
        # Coordinates starts
        NodeFlag = False
        
        for line in Inpfile:        
            if '*' in line and NodeFlag:
                break
            
            if NodeFlag:
                NodeCoord_raw.append(line[:-1])
            
            if '*Node' in line:
                NodeFlag = True

    NodeCoord = {}

    #Converting the raw data to a appropriate data type and storing it
    # in a dictionary
    for i in range(len(NodeCoord_raw)):
        temp = [float(s) for s in NodeCoord_raw[i].split(',')]
        NodeCoord[int(temp[0])] = temp[1:]
    
    return ElemNodeConn, NodeCoord

def writeInputfile(EdgeNodeConn, Inputfile, filename, UELprop, stressType):
    '''
    Writes the input file for ESFEM element

    Parameters
    ----------
    EdgeNodeConn: DICT
                Dictionary containing the information for Nodes associated to 
                smoothed domain of an Edge
    Inputfile: STR
              Name of the input file
    filename: STR
             Name of the newcreated file
    UELprop: LIST
            LIST of all the properties of the UEL.

    '''

    fp = int(math.log(2*(list(EdgeNodeConn.keys())[-1]))) + 1

    with open(Inputfile) as Inpfile:
        InpfileContent = Inpfile.readlines()
    
    IsNodeline = False
    for index, line in enumerate(InpfileContent):

        if IsNodeline and '*' in line:

            if stressType == '2D':
                UEL4Nodesdef = '*USER ELEMENT, TYPE=U2004, NODES=4, COORDINATES=2, PROPERTIES=2, VARIABLES=10, UNSYMM\n'
                UEL3Nodesdef = '*USER ELEMENT, TYPE=U2003, NODES=3, COORDINATES=2, PROPERTIES=2, VARIABLES=10, UNSYMM\n'
            elif stressType == '3D':
                UEL4Nodesdef = '*USER ELEMENT, TYPE=U2004, NODES=4, COORDINATES=2, PROPERTIES=2, VARIABLES=19, UNSYMM\n'
                UEL3Nodesdef = '*USER ELEMENT, TYPE=U2003, NODES=3, COORDINATES=2, PROPERTIES=2, VARIABLES=19, UNSYMM\n'
            else: 
                print('Incorrect stress type : ', stressType)
                return
                 
            UELdof = '{:d},{:d}\n'.format(1,2)
            
            InpfileContent.insert(index, UEL3Nodesdef)
            InpfileContent.insert(index+1, UELdof)
            InpfileContent.insert(index+2, UEL4Nodesdef)
            InpfileContent.insert(index+3, UELdof)

            Elemdeclare = '*Element, type=U2003, elset=esfemfreeedge\n'
            InpfileContent.insert(index+4, Elemdeclare)

            idx = index + 5
            for key in EdgeNodeConn.keys():
                if len(EdgeNodeConn[key]) == 3:
                    dataline = '{:{width}d}, {:{width}d}, {:{width}d}, {:{width}d}\n'.format(
                        key, EdgeNodeConn[key][0], EdgeNodeConn[key][1], EdgeNodeConn[key][2],
                        width=fp) 
                    InpfileContent.insert(idx, dataline)
                    idx = idx + 1

            Elemdeclare = '*Element, type=U2004, elset=esfemedge\n'
            InpfileContent.insert(idx, Elemdeclare)
            
            idx = idx + 1

            for key in EdgeNodeConn.keys():
                if len(EdgeNodeConn[key]) == 4:
                    dataline = '{:{width}d}, {:{width}d}, {:{width}d}, {:{width}d}, {:{width}}\n'.format(
                        key, EdgeNodeConn[key][0], EdgeNodeConn[key][1], EdgeNodeConn[key][2], EdgeNodeConn[key][3],
                        width=fp) 
                    InpfileContent.insert(idx, dataline)
                    idx = idx + 1

            break
        
        if '*Node' in line:
            IsNodeline = True

    lastelem = list(EdgeNodeConn.keys())[-1]
    isghostelem = False
    # stop = False
    elemno = 1
    for index, line in enumerate(InpfileContent):

        if isghostelem and '*' in line:
            break
        
        if isghostelem:
            replacethis = '{:d}'.format(elemno)
            withthis = '{:{width}d}'.format(lastelem+elemno, width=fp)
            InpfileContent[index] = InpfileContent[index].replace(replacethis, withthis, 1)
            elemno = elemno + 1

        if 'CPS3' in line:
            isghostelem = True
    
    for index, line in enumerate(InpfileContent):
        if '*End Part' in line:
            UEL3propline = '*UEL PROPERTY, ELSET=esfemfreeedge\n'
            UEL4propline = '*UEL PROPERTY, ELSET=esfemedge\n'
            UELpropdataline = ''
            for i in range(len(UELprop)):
                UELpropdataline = UELpropdataline + '{:.{prec}f}, '.format(UELprop[i], prec=4)
            UELpropdataline = UELpropdataline + '\n'
            
            InpfileContent.insert(index, UEL3propline)
            InpfileContent.insert(index+1, UELpropdataline)
            InpfileContent.insert(index+2, UEL4propline)
            InpfileContent.insert(index+3, UELpropdataline)
            break

    with open(filename, 'w') as newInpfile:
        newInpfile.seek(0)
        newInpfile.writelines(InpfileContent)

    return

def write_MeshInfo(ESFEM_Mesh, filename):

    ESFEM_Mesh.createmesh()

    ElemTopo = ESFEM_Mesh.ElementTopology()

    fileContent = []

    fileContent.append('{:13} {:5d}\n'.format('Nelem', ESFEM_Mesh.Nelem))
    fileContent.append('{:13} {:5d}\n'.format('totNodes', ESFEM_Mesh.totNodes))
    fileContent.append('{:13} {:5d}\n'.format('totCentroids', ESFEM_Mesh.totCentroids))

    fileContent.append('Points data starts here : \n')

    for key in ESFEM_Mesh.Points.keys():
        line = '{:5d}, {:{width}.{prec}f}, {:{width}.{prec}f}\n'.format(key, ESFEM_Mesh.Points[key][0], 
                                                    ESFEM_Mesh.Points[key][1], prec=16, width=24)
        
        fileContent.append(line)

    fileContent.append('Element topology data starts here : \n')

    for key in ElemTopo.keys():
        thiselem = ElemTopo[key]
        line = '{:5d}'.format(key)

        for i in range(len(thiselem)):
            line  =  line + ', {:5d}'.format(thiselem[i])

        if len(thiselem) == 3:
            line = line + ', {:5d}'.format(0)
        line = line + '\n'

        fileContent.append(line)

    with open(filename, 'w') as meshInfo:
        meshInfo.writelines(fileContent)

    return 



######################### End of functions #######################


######################### Code Runs here #########################

###  Reading innput files and stroing it in dictionaries  ####

import argparse

def main(Inputfile, meshInfofile, stressType):
    
    if '.inp' not in Inputfile:
        print('Incorrect file name')

    if os.path.isfile(Inputfile):
        print('Extracting data from input file ...')
        ElemNodeConn, NodesCoord = ReadInpfile(Inputfile)
        print('Data Extraction complete ...')

        print('Overlaying ESFEM mesh ... (It can take few minutes)')
        EdgeConn, ElemEdgeConn, EdgeNumber = EdgeConnectivity(ElemNodeConn)
        EdgeNodeConn = ESFEM_Nodes(ElemNodeConn, EdgeConn, EdgeNumber)
        print('Overlaying Completed ...')

        name = Inputfile.split('.')

        newInpfile = name[0] + '_ESFEM.inp'

        print('Writing ESFEM input file ...')
        writeInputfile(EdgeNodeConn, Inputfile, newInpfile, UELprop=[70000., 0.3], stressType=stressType)
        print('Completed writing input file ...')

        print('Generating Mesh Connectivity data ...')
        Mesh = ESFEM_Mesh(NodesCoord, EdgeNodeConn, EdgeNodeConn)
        print('Completed ...')

        print('Writing Mesh Info file ... (It can take few minutes)')
        write_MeshInfo(Mesh, meshInfofile)
        print('Completed ...')

        print('Done !!!')


    else:
        print('Error file does not exist ')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Arguments for the preprocessing file of an ESFEM element')
    parser.add_argument('-i', '--inpfile', type=str, metavar='', help='The name of the input file with .inp extension')
    parser.add_argument('-s', '--stress', type=str, metavar='', help='Type of stress used for analysis(2D or 3D)')
    args = parser.parse_args()

    main(args.inpfile, 'MeshInfo.txt', args.stress)

