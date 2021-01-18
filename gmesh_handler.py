import sys
import os
import numpy as np
import time
import pdb

class read:

    def __init__(self, filepath):
        self.path = filepath
        self.type = 1
        self.triangule = 2
        self.element_start = 0
        self.physical_names = 0

        if not os.path.isfile('./mesh/'+self.path):
            print('The file {} does not exist!'.format(self.path))
            sys.exit()
        else:
            print('File found!')
            print('Setting path for BST')
            #with open("path", 'w') as path_file:
            #    path_file.write(self.path)
            #    path_file.write(str(9))
              #  path_file.write(self.nodes)
              #  path_file.write(self.elements)
            print('Done.')
            with open('./mesh/'+self.path, 'r') as file:
                while True:
                    line = file.readline().rstrip()
                    if line == '$MeshFormat':
                        form = file.readline().rstrip()
                        version = form[0]
                        if int(version) == 2:
                            print('Mesh format: {}'.format(form))
                        else:
                            print('Version not supported!')
                            #terminar programa
                            sys.exit()
                    if line == "$EndMeshFormat":
                        file.close()
                        break

    def read_nodes(self):

        with open('./mesh/'+self.path, 'r') as file:
            while True:
                line = file.readline().rstrip()
                if line == '$PhysicalNames':
                    physical = file.readline().rstrip()
                    self.physical_names = int(physical)
                if line == '$Nodes':
                    nodes_read = file.readline().rstrip()
                    self.nodes = int(nodes_read)
                    print('Nodes: {}'.format(self.nodes))
                    self.nodes_mat = np.zeros([self.nodes, 3])
                    for i in range(self.nodes):
                        node = np.fromstring(file.readline().rstrip(), sep=" ")
                        self.nodes_mat[i] = node[1 : ]
                elif line == '$EndNodes':
                    file.close()
                    break
        if self.nodes < 1:
            print('No nodes were found. Exiting...')
            #no encontro nodos. salir del programa
        else:
            print('Ok!')
        return self.nodes_mat
    
    def find_elements(self, ElementType):
        self.toSearch = ElementType
        self.element_start = 0
        print('Looking for type: {} elements'.format(ElementType))
        with open('./mesh/'+self.path, 'r') as file:
            while True:
                line = file.readline().rstrip()
                if line == '$Elements':
                    elements = file.readline().rstrip()
                    self.elements = int(elements)
                    if self.toSearch == 4:
                        self.element_mat = np.empty([0,4], dtype='int')
                    else:
                        self.element_mat = np.empty([0,3], dtype='int')
                    for i in range(self.elements):
                        element = np.fromstring(file.readline().rstrip(), sep=" ", dtype='int')
                        if element[self.type] == self.toSearch:
                            if self.toSearch == 4:
                                self.element_mat = np.append(self.element_mat, [element[-4 :]], axis=0)
                            else:
                                self.element_mat = np.append(self.element_mat, [element[-3 :]], axis=0)
                        else:
                            self.element_start = self.element_start + 1
                    print('Elements: {}'.format(np.size(self.element_mat,0)))
                elif line == '$EndElements':
                    break
        if np.size(self.element_mat, 0) < 1:
            print('No elements of the type were found. Retry')
        else:
            print('Done.')
        return self.element_mat
    
    def write(self, D, D_name, F, F_name, sigma, sigma_name, dof_per_node):
        nodes = int(len(D)/dof_per_node)
        elements = int(len(sigma))
        flag = 0
        with open('./mesh/'+self.path, 'r') as file:
            while True:
                line = file.readline().rstrip()
                if line == "$EndElements":

                    for i in range(10):
                        line = file.readline().rstrip()
                        if line == '$NodeData':
                            flag = 1
                            print('The file already has nodes and elements data!')
                    break
        if flag == 0:
            file.close()
            with open('./mesh/'+self.path, 'a') as file:
            # Escribe los desplazamientos de cada nodo
                if dof_per_node == 2:
                    file.write('\n')
                    file.write('$NodeData\n')
                    file.write('1\n')
                    file.write('"'+D_name+'"')
                    file.write('\n')
                    file.write('1\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('{}\n'.format(nodes))
                    index = 0
                    for i in range(len(D)):
                        if i % 2 ==0:
                            index+=1
                            file.write(str(index))
                            file.write(' ')
                            file.write(str(D[i].item()))
                            file.write(' ')
                            file.write(str(D[i+1].item()))
                            file.write(' ' + str(0))
                            file.write('\n')
                    file.write('$EndNodeData')
                #Escribe las fuerzas en cada nodo
                    file.write('\n')
                    file.write('$NodeData\n')
                    file.write('1\n')
                    file.write('"'+F_name+'"')
                    file.write('\n')
                    file.write('1\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('{}\n'.format(nodes))
                    index = 0
                    for i in range(len(F)):
                        if i % 2 ==0:
                            index+=1
                            file.write(str(index))
                            file.write(' ')
                            file.write(str(F[i].item()))
                            file.write(' ')
                            file.write(str(F[i+1].item()))
                            file.write(' ' + str(0))
                            file.write('\n')
                    file.write('$EndNodeData')
                # Escribe la propiedad de cada elemento
                    file.write('\n')
                    file.write('$ElementData\n')
                    file.write('1\n')
                    file.write('"'+sigma_name+'"')
                    file.write('\n')
                    file.write('1\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('0\n')
                    file.write('1\n')
                    file.write('{}\n'.format(elements))
                    index = 0
                    for i in range(len(sigma)):
                        file.write(str(i+1+self.element_start))
                        file.write(' ')
                        file.write(str(sigma[i].item()))
                        file.write('\n')
                    file.write('$EndElementData')
                
                elif dof_per_node == 1:
                    file.write('\n')
                    file.write('$NodeData\n')
                    file.write('1\n')
                    file.write('"'+D_name+'"')
                    file.write('\n')
                    file.write('1\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('{}\n'.format(nodes))
                    index = 0
                    for i in range(len(D)):
                        #if i % 2 ==0:
                        index+=1
                        file.write(str(index))
                        file.write(' ')
                        file.write(str(D[i].item()))
                        file.write(' ' + str(0))
                        #file.write(str(D[i+1].item()))
                        file.write(' ' + str(0))
                        file.write('\n')
                    file.write('$EndNodeData')
                #Escribe las fuerzas en cada nodo
                    file.write('\n')
                    file.write('$NodeData\n')
                    file.write('1\n')
                    file.write('"'+F_name+'"')
                    file.write('\n')
                    file.write('1\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('{}\n'.format(nodes))
                    index = 0
                    for i in range(len(F)):
                        #if i % 2 ==0:
                        index+=1
                        file.write(str(index))
                        file.write(' ')
                        file.write(str(F[i].item()))
                        file.write(' ' + str(0))
                        #file.write(str(F[i+1].item()))
                        file.write(' ' + str(0))
                        file.write('\n')
                    file.write('$EndNodeData')
                # Escribe la propiedad de cada elemento
                if dof_per_node == 3:
                    file.write('\n')
                    file.write('$NodeData\n')
                    file.write('1\n')
                    file.write('"'+D_name+'"')
                    file.write('\n')
                    file.write('1\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('{}\n'.format(nodes))
                    index = 0
                    for i in range(len(D)):
                        if i % 3 ==0:
                            index+=1
                            file.write(str(index))
                            file.write(' ')
                            file.write(str(D[i].item()))
                            file.write(' ')
                            file.write(str(D[i+1].item()))
                            file.write(' ')
                            file.write(str(D[i+2].item())+'\n')
                    file.write('$EndNodeData')
                #Escribe las fuerzas en cada nodo
                    file.write('\n')
                    file.write('$NodeData\n')
                    file.write('1\n')
                    file.write('"'+F_name+'"')
                    file.write('\n')
                    file.write('1\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('{}\n'.format(nodes))
                    index = 0
                    for i in range(len(F)):
                        if i % 3 ==0:
                            index+=1
                            file.write(str(index))
                            file.write(' ')
                            file.write(str(F[i].item()))
                            file.write(' ')
                            file.write(str(F[i+1].item()))
                            file.write(' ')
                            file.write(str(F[i+2].item())+'\n')
                    file.write('$EndNodeData')
                # Escribe la propiedad de cada elemento
                    file.write('\n')
                    file.write('$ElementData\n')
                    file.write('1\n')
                    file.write('"'+sigma_name+'"')
                    file.write('\n')
                    file.write('1\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('0\n')
                    file.write('1\n')
                    file.write('{}\n'.format(elements))
                    index = 0
                    for i in range(len(sigma)):
                        file.write(str(i+1+self.element_start))
                        file.write(' ')
                        file.write(str(sigma[i].item()))
                        file.write('\n')
                    file.write('$EndElementData')
                '''
                    file.write('\n')
                    file.write('$ElementData\n')
                    file.write('1\n')
                    file.write('"'+sigma_name+'"')
                    file.write('\n')
                    file.write('1\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('0\n')
                    file.write('1\n')
                    file.write('{}\n'.format(elements))
                    index = 0
                    for i in range(len(sigma)):
                        file.write(str(i+1+self.element_start))
                        file.write(' ')
                        file.write(str(sigma[i].item()))
                        file.write('\n')
                    file.write('$EndElementData')
                '''
                print('File succesfully generated!')
        else:
            print('No changes were made!')
    
    def clean(self):
        buffer = []
        with open('./mesh/'+self.path, 'r') as file:
            while True:
                line = file.readline().rstrip()
                buffer.append(line)
                if line == '$EndElements':
                    break
        file.close()
        with open('./mesh/'+self.path, 'w') as file:
            for lines in buffer:
                if lines == '$EndElements':
                    file.write(lines)
                else:
                    file.write(lines+'\n')

        print('File cleaned!')
    
    def BST_config(self):
            print('Editing BST configuration file')
            if self.physical_names == 0:
                print('No physical names were found!')
                print('Exiting...')
                exit()
            else:
                with open("config", 'w') as path_file:
                    path_file.write(self.path)
                    path_file.write('\n')
                    path_file.write(str(self.physical_names))
                    path_file.write('\n')
                    path_file.write(str(self.nodes))
                    path_file.write('\n')
                    path_file.write(str(self.elements))
                    path_file.write('\n')
                print('Done.')

    def find_physical_elements(self, ElementType):
        self.toSearch = ElementType
        self.element_start = 0
        print('Looking for type: {} elements'.format(ElementType))
        with open('./mesh/'+self.path, 'r') as file:
            while True:
                line = file.readline().rstrip()
                if line == '$Elements':
                    elements = file.readline().rstrip()
                    self.elements = int(elements)
                    if self.toSearch == 1:
                        self.element_mat = np.empty([0,4], dtype='int')
                    elif self.toSearch == 2:
                        self.element_mat = np.empty([0,5], dtype='int')

                    for i in range(self.elements):
                        element = np.fromstring(file.readline().rstrip(), sep=" ", dtype='int')
                        if element[self.type] == self.toSearch:
                            if self.toSearch == 1:
                                self.element_mat = np.append(self.element_mat, [element[-4 :]], axis=0)
                            elif self.toSearch == 2:
                                self.element_mat = np.append(self.element_mat, [element[-5 :]], axis=0)
                                # Obtener el physical group de los elementos tipo 2
                                #print(element)
                                #self.element_mat = np.append(self.element_mat, [element[0]], axis=0)
                                #self.element_mat = np.append(self.element_mat, [element[-3:]], axis=0)
                        else:
                            self.element_start = self.element_start + 1
                    print('Elements: {}'.format(np.size(self.element_mat,0)))
                elif line == '$EndElements':
                    break
        if np.size(self.element_mat, 0) < 1:
            print('No elements of the type were found. Retry')
        else:
            print('Done.')
        return self.element_mat

    def simple_write(self, data, data_name, dof_per_node):

        nodes = int(len(data)/dof_per_node)
        with open('./mesh/'+self.path, 'a') as file:
        # Escribe los desplazamientos de cada nodo
            if dof_per_node == 2:
                file.write('\n')
                file.write('$NodeData\n')
                file.write('1\n')
                file.write('"'+data_name+'"')
                file.write('\n')
                file.write('1\n')
                file.write('0\n')
                file.write('3\n')
                file.write('0\n')
                file.write('3\n')
                file.write('{}\n'.format(nodes))
                index = 0
                for i in range(len(data)):
                    if i % 2 ==0:
                        index+=1
                        file.write(str(index))
                        file.write(' ')
                        file.write(str(data[i].item()))
                        file.write(' ')
                        file.write(str(data[i+1].item()))
                        file.write(' ' + str(0))
                        file.write('\n')
                file.write('$EndNodeData')

            elif dof_per_node == 3:
                file.write('\n')
                file.write('$NodeData\n')
                file.write('1\n')
                file.write('"'+data_name+'"')
                file.write('\n')
                file.write('1\n')
                file.write('0\n')
                file.write('3\n')
                file.write('0\n')
                file.write('3\n')
                file.write('{}\n'.format(nodes))
                index = 0
                for i in range(len(data)):
                    if i % dof_per_node ==0:
                        index+=1
                        file.write(str(index))
                        file.write(' ')
                        file.write(str(data[i].item()))
                        file.write(' ')
                        file.write(str(data[i+1].item()))
                        file.write(' ')
                        file.write(str(data[i+2].item()))
                        file.write('\n')
                file.write('$EndNodeData')
    
    def clip_write(self, data, data_name, dof_per_node):

        nodes = int(len(data)/dof_per_node)
        tags = 62
        for timestamp in range(tags):
            with open('./mesh/'+self.path, 'a') as file:
            # Escribe los desplazamientos de cada nodo
                if dof_per_node == 2:
                    file.write('\n')
                    file.write('$NodeData\n')
                    file.write('1\n')
                    file.write('"'+data_name+'"')
                    file.write('\n')
                    file.write('1\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('0\n')
                    file.write('3\n')
                    file.write('{}\n'.format(nodes))
                    index = 0
                    for i in range(len(data)):
                        if i % 2 ==0:
                            index+=1
                            file.write(str(index))
                            file.write(' ')
                            file.write(str(data[i].item()))
                            file.write(' ')
                            file.write(str(data[i+1].item()))
                            file.write(' ' + str(0))
                            file.write('\n')
                    file.write('$EndNodeData')

                if dof_per_node == 3:
                    file.write('\n')
                    file.write('$NodeData\n')
                    file.write('1\n')
                    file.write('"'+data_name+'"')
                    file.write('\n')
                    file.write('1\n')
                    file.write(str(timestamp)+'\n')
                    file.write('3\n')
                    file.write(str(timestamp)+'\n')
                    file.write('3\n')
                    file.write('{}\n'.format(nodes))
                    index = 0
                    for i in range(len(data)):
                        if i % dof_per_node ==0:
                            index+=1
                            file.write(str(index))
                            file.write(' ')
                            file.write(str(np.cos(timestamp/10)*data[i].item()))
                            file.write(' ')
                            file.write(str(np.cos(timestamp/10)*data[i+1].item()))
                            file.write(' ')
                            file.write(str(np.cos(timestamp/10)*data[i+2].item()))
                            file.write('\n')
                    file.write('$EndNodeData')

    def dynamic_clip_write(self, data, data_name, dof_per_node, tags):

        nodes = int(data.shape[1]/dof_per_node)
        for timestamp in range(tags):
            with open('./mesh/'+self.path, 'a') as file:
            # Escribe los desplazamientos de cada nodo
                if dof_per_node == 2:
                    file.write('\n')
                    file.write('$NodeData\n')
                    file.write('1\n')
                    file.write('"'+data_name+'"')
                    file.write('\n')
                    file.write('1\n')
                    file.write(str(timestamp)+'\n')
                    file.write('3\n')
                    file.write(str(timestamp)+'\n')
                    file.write('3\n')
                    file.write('{}\n'.format(nodes))
                    index = 0
                    for i in range(data.shape[1]):
                        if i % 2 ==0:
                            index+=1
                            file.write(str(index))
                            file.write(' ')
                            file.write(str(data[timestamp, i].item()))
                            file.write(' ')
                            file.write(str(data[timestamp, i+1].item()))
                            file.write(' ' + str(0))
                            file.write('\n')
                    file.write('$EndNodeData')

                if dof_per_node == 3:
                    file.write('\n')
                    file.write('$NodeData\n')
                    file.write('1\n')
                    file.write('"'+data_name+'"')
                    file.write('\n')
                    file.write('1\n')
                    file.write(str(timestamp)+'\n')
                    file.write('3\n')
                    file.write(str(timestamp)+'\n')
                    file.write('3\n')
                    file.write('{}\n'.format(nodes))
                    index = 0
                    for i in range(data.shape[1]):
                        if i % dof_per_node == 0:
                            index+=1
                            file.write(str(index))
                            file.write(' ')
                            file.write(str(data[timestamp, i].item()))
                            file.write(' ')
                            file.write(str(data[timestamp, i+1].item()))
                            file.write(' ')
                            file.write(str(data[timestamp, i+2].item()))
                            file.write('\n')
                    file.write('$EndNodeData')

    
