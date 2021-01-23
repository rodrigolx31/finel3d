import numpy as np 
import time
from alive_progress import alive_bar
import pdb

class finite_elements:
    
    def __init__(self, nodes, elements, dof_per_node):
        self.nodes = nodes
        self.elements = elements
        self.dof_per_node = dof_per_node
        self.areas = np.empty([np.size(elements,0),0])
        self.long = np.empty([np.size(elements, 0), 0])
    
    def get_areas(self, ElementType, areas):
        self.areas = np.empty([np.size(self.elements,0),0])
        if ElementType == 2:
            for i in range(np.size(self.elements,0)):
                i, j, m = self.elements[i, :] - 1 #le resto 1 a cada elemento para que coincida con los indices 
                x, y = self.nodes[:,0], self.nodes[:,1]
                A = np.array([[1, x[i], y[i]],
                              [1, x[j], y[j]],
                              [1, x[m], y[m]]])
                area = np.linalg.det(A) / 2
                self.areas = np.append(self.areas, area)

        if ElementType == 1:
            self.areas = areas * np.ones([np.size(self.elements, 0)])
        return self.areas
    def get_stiff_mat(self, ElementType):

        if ElementType == 2:
            D = np.array([[1,0.3,0], [0.3, 1, 0],[0,0,0.35]])
            self.D = (1/0.91)*D
            #self.D = (30e6/0.91)*D

            stiff_matrix = []
            self.B = []
            x, y = self.nodes[:,0], self.nodes[:,1]
            for index in range(np.size(self.elements,0)):
                i, j, m = self.elements[index, :] - 1 #le resto 1 a cada elemento para que coincida con los indices 
                #x, y = self.nodes[:,0], self.nodes[:,1] movido a un lugar donde se ejecute una sola vez

                beta_i = y[j] - y[m]
                beta_j = y[m] - y[i]
                beta_m = y[i] - y[j]

                gamma_i = x[m] - x[j]
                gamma_j = x[i] - x[m]
                gamma_m = x[j] - x[i]
            
                B0 = np.zeros([3,6])
                B0[0,[0,2,4]] = beta_i, beta_j, beta_m
                B0[1,[1,3,5]] = gamma_i, gamma_j, gamma_m
                B0[2] = gamma_i, beta_i, gamma_j, beta_j, gamma_m, beta_m
                B0 = 1/(2*self.areas[index]) * B0
                self.B.append(B0)
                #calcula la matriz de rigidez
                Bt = np.transpose(B0)
                k = np.matmul(Bt, self.D)
                k2 = np.matmul(k,B0)
                k2 = abs(self.areas[index])*k2

                stiff_matrix.append(k2)
            
            self.stiff_mat = np.array(stiff_matrix)
        
        elif ElementType == 1:

            stiff_matrix = []
            #Este metodo obtiene la matriz de rigidez
            self.nodes_per_element = 2
            k_size = self.dof_per_node*self.nodes_per_element
            k_matrix = np.zeros([k_size, k_size])
            for index, element in enumerate(self.elements, start=0):
                #Obtiene el angulo de cada elemento
                #print(element[1:3])
                x = self.nodes[element[1:3]-1]
                #print(x)
                X = x[1,0] - x[0, 0]
                Y = x[1,1] - x[0, 1] 
                theta = np.arctan2(Y,X) #Modificado!!! ojo
                
                #Obtiene la longitud de cada elemento
                a = np.array([x[0,0], x[0,1]]) 
                b = np.array([x[1,0], x[1,1]])
                L = np.linalg.norm(a-b)
                '''
                print('Elemento: {}, Angulo: {}°, Long: {}, Area: {}'.format(index+1,
                                                                np.around(theta*180/np.pi, 2),
                                                                np.around(L,2),
                                                                self.areas[index])) # Muestra angulos de c/elemento
                '''
                Kel = 30e6*self.areas[index]/L
                #Calcula la matriz de rigidez de cada elemento para 2 grados de libertad por nodo
                comp_k_matrix = np.array([[
                np.power(np.cos(theta), 2),
                np.cos(theta)*np.sin(theta),
                -np.power(np.cos(theta), 2),
                -np.cos(theta)*np.sin(theta)],
                [
                np.cos(theta)*np.sin(theta),
                np.power(np.sin(theta), 2),
                -np.cos(theta)*np.sin(theta),
                -np.power(np.sin(theta), 2)],
                [
                -np.power(np.cos(theta), 2),
                -np.cos(theta)*np.sin(theta),
                np.power(np.cos(theta), 2),
                np.cos(theta)*np.sin(theta)],
                [
                -np.cos(theta)*np.sin(theta),
                -np.power(np.sin(theta), 2),
                np.cos(theta)*np.sin(theta),
                np.power(np.sin(theta), 2)]],
                dtype=np.single)
                for i in range(k_size):
                    for j in range(k_size):
                        if self.dof_per_node is 1:
                            k_matrix[i,j] = comp_k_matrix[i*self.nodes_per_element, j*self.nodes_per_element]
                        else: 
                            k_matrix[i,j] = comp_k_matrix[i, j]
                        
                stiff_matrix.append(np.around(Kel*k_matrix,9))
            self.stiff_mat = np.array(stiff_matrix)
        
        elif ElementType == 5:
            stiff_matrix = []
            #Este metodo obtiene la matriz de rigidez
            self.nodes_per_element = 2
            k_size = self.dof_per_node*self.nodes_per_element
            k_matrix = np.zeros([k_size, k_size])
            longitude = [] 
            for index, element in enumerate(self.elements, start=0):
                #Obtiene el angulo de cada elemento
                #print(element[1:3])
                x = self.nodes[element[1:3]-1]
                #print(x)
                X = x[1,0] - x[0, 0]
                Y = x[1,1] - x[0, 1] 
                theta = np.arctan2(Y,X) #Modificado!!! ojo
                
                #Obtiene la longitud de cada elemento
                a = np.array([x[0,0], x[0,1]]) 
                b = np.array([x[1,0], x[1,1]])
                L = np.linalg.norm(a-b)
                #Matriz de rigidez para una viga
                k = np.array([ [12, 6*L, -12, 6*L],
                [6*L, 4*np.power(L,2), -6*L, 2*np.power(L,2)],
                [-12, -6*L, 12, -6*L],
                [6*L, 2*np.power(L,2), -6*L, 4*np.power(L,2)]])
 
                #      E  *  I    / L**3
                kel = 210e9*10e-8/np.power(L, 3)
                stiff_matrix.append(kel * k)
                longitude.append(L)
            self.stiff_mat = np.array(stiff_matrix)
            self.long = np.array(longitude)
        elif ElementType == 4:
            E = 210e3
            nu = 0.3
            D = np.zeros([6, 6])
            D[0,0], D[1,1], D[2,2] = (1 - nu) * np.array([1, 1, 1])
            D[3,3], D[4,4], D[5,5] = ((1 - 2*nu)/2) * np.array([1, 1, 1])
            D[0, 1], D[0,2], D[1, 2] = nu, nu, nu
            D[1, 0], D[2,0], D[2, 1] = nu, nu, nu
            self.D = D * ( E/((1+nu)*(1-(2*nu))) )
            self.B = []
            stiff_matrix = []
            volume = []
            x, y, z = self.nodes[:,0], self.nodes[:,1], self.nodes[:,2]
            for index in range(np.size(self.elements,0)):
                i1, i2, i3, i4 = self.elements[index, :] - 1 #le resto 1 a cada elemento para que coincida con los indices
                # 1
                A = np.array([  [x[i2], y[i2], z[i2]],
                                [x[i3], y[i3], z[i3]],
                                [x[i4], y[i4], z[i4]]])
                alpha_1 = np.linalg.det(A)

                A = np.array([  [1, y[i2], z[i2]],
                                [1, y[i3], z[i3]],
                                [1, y[i4], z[i4]]])
                beta_1 = -1 * np.linalg.det(A)

                A = np.array([  [1, x[i2], z[i2]],
                                [1, x[i3], z[i3]],
                                [1, x[i4], z[i4]]])
                gamma_1 = np.linalg.det(A)

                A = np.array([  [1, x[i2], y[i2]],
                                [1, x[i3], y[i3]],
                                [1, x[i4], y[i4]]])
                delta_1 = -1 * np.linalg.det(A)
                
                # 2
                A = np.array([  [x[i1], y[i1], z[i1]],
                                [x[i3], y[i3], z[i3]],
                                [x[i4], y[i4], z[i4]]])
                alpha_2 = -1 * np.linalg.det(A)

                A = np.array([  [1, y[i1], z[i1]],
                                [1, y[i3], z[i3]],
                                [1, y[i4], z[i4]]])
                beta_2 = np.linalg.det(A)

                A = np.array([  [1, x[i1], z[i1]],
                                [1, x[i3], z[i3]],
                                [1, x[i4], z[i4]]])
                gamma_2 = -1 * np.linalg.det(A)

                A = np.array([  [1, x[i1], y[i1]],
                                [1, x[i3], y[i3]],
                                [1, x[i4], y[i4]]])
                delta_2 = np.linalg.det(A)
                # 3
                A = np.array([  [x[i1], y[i1], z[i1]],
                                [x[i2], y[i2], z[i2]],
                                [x[i4], y[i4], z[i4]]])
                alpha_3 = np.linalg.det(A)

                A = np.array([  [1, y[i1], z[i1]],
                                [1, y[i2], z[i2]],
                                [1, y[i4], z[i4]]])
                beta_3 = -1 * np.linalg.det(A)

                A = np.array([  [1, x[i1], z[i1]],
                                [1, x[i2], z[i2]],
                                [1, x[i4], z[i4]]])
                gamma_3 = np.linalg.det(A)

                A = np.array([  [1, x[i1], y[i1]],
                                [1, x[i2], y[i2]],
                                [1, x[i4], y[i4]]])
                delta_3 = -1 * np.linalg.det(A)
                # 4
                A = np.array([  [x[i1], y[i1], z[i1]],
                                [x[i2], y[i2], z[i2]],
                                [x[i3], y[i3], z[i3]]])
                alpha_4 = -1 * np.linalg.det(A)

                A = np.array([  [1, y[i1], z[i1]],
                                [1, y[i2], z[i2]],
                                [1, y[i3], z[i3]]])
                beta_4 = np.linalg.det(A)

                A = np.array([  [1, x[i1], z[i1]],
                                [1, x[i2], z[i2]],
                                [1, x[i3], z[i3]]])
                gamma_4 = -1 * np.linalg.det(A)

                A = np.array([  [1, x[i1], y[i1]],
                                [1, x[i2], y[i2]],
                                [1, x[i3], y[i3]]])
                delta_4 = np.linalg.det(A)
                
                volume_6 = np.linalg.det([[1, x[i1], y[i1], z[i1]],
                                          [1, x[i2], y[i2], z[i2]],
                                          [1, x[i3], y[i3], z[i3]],
                                          [1, x[i4], y[i4], z[i4]]])

                B1= np.array([[beta_1, 0, 0],
                             [0, gamma_1, 0],
                             [0, 0, delta_1],
                             [gamma_1, beta_1, 0],
                             [0, delta_1, gamma_1],
                             [delta_1, 0, beta_1]])

                B2= np.array([[beta_2, 0, 0],
                              [0, gamma_2, 0],
                              [0, 0, delta_2],
                              [gamma_2, beta_2, 0],
                              [0, delta_2, gamma_2],
                              [delta_2, 0, beta_2]])

                B3= np.array([[beta_3, 0, 0],
                              [0, gamma_3, 0],
                              [0, 0, delta_3],
                              [gamma_3, beta_3, 0],
                              [0, delta_3, gamma_3],
                              [delta_3, 0, beta_3]])
                
                B4= np.array([[beta_4, 0, 0],
                              [0, gamma_4, 0],
                              [0, 0, delta_4],
                              [gamma_4, beta_4, 0],
                              [0, delta_4, gamma_4],
                              [delta_4, 0, beta_4]])

                B = np.concatenate((B1, B2, B3, B4), axis=1)
                B = B / volume_6
                self.B.append(B)
                # Obtengo la matriz de rigidez
                Bt = np.transpose(B)
                k = np.matmul(Bt, self.D)
                k2 = np.matmul(k,B)
                k2 = k2 * (volume_6/6)
                stiff_matrix.append(k2)
                volume.append(volume_6/6)
            self.stiff_mat = np.array(stiff_matrix)
            self.volume_array = np.array(volume)
            #pdb.set_trace()

        return self.stiff_mat
    
    def get_global_mat(self, ElementType, matrix_to_assemble):
        
        if ElementType == 1:
            nodes_per_element = 2
            index_correction = 1
        
        elif ElementType == 2:
            nodes_per_element = 3
            index_correction = 0
        
        elif ElementType ==5:
            nodes_per_element = 2
            index_correction = 1

        elif ElementType ==4:
            nodes_per_element = 4
            index_correction = 0

        #numero de nodos*grados de libertad por nodo
        size = np.size(self.nodes, 0) * self.dof_per_node
        k_global = np.zeros([size,size])
        k = matrix_to_assemble
        assembled = np.zeros([size,size])
        print()
        print('Getting global matrix...')
        with alive_bar(len(self.elements), bar = 'filling', spinner = 'dots_reverse') as bar:
            for element in range(len(self.elements)):
                for i in range(nodes_per_element):
                    for j in range(nodes_per_element):
                        ind_u = self.dof_per_node* (self.elements[element,i + index_correction] - 1)
                        ind_v = self.dof_per_node* (self.elements[element,j + index_correction] - 1)
                        assembled[ind_u:ind_u+self.dof_per_node, ind_v:ind_v+self.dof_per_node] += k[element, self.dof_per_node*i:self.dof_per_node*i+self.dof_per_node, self.dof_per_node*j:self.dof_per_node*j+self.dof_per_node]
                bar()
        return assembled

    def solve_for(self, K, F_vect, D_vect, r, s):
        print('Calculating solution...')
        with alive_bar(6, bar = 'filling', spinner = 'dots_reverse') as bar:
            start_time = time.time()
            end_time = 0
            self.Forces = F_vect
            self.Displacements = D_vect
            self.r = r
            self.s = s
            #K = self.k_global
            bar()

            Kred = K[np.ix_(self.r, self.r)]
            self.Kred = Kred
            bar()

            Kinv = np.linalg.inv(Kred)
            bar()

            Kvin = K[np.ix_(self.r,self.s)]
            bar()

           # Posiciones finales
            self.Displacements[self.r] = np.matmul(Kinv,(self.Forces[self.r] - np.matmul(Kvin,self.Displacements[self.s])))
            bar()

            # Fuerzas
            self.Forces[self.s] = np.matmul(K[np.ix_(self.s,)],self.Displacements)
            end_time = time.time()
            bar()

        print()
        print("Solved in {} seconds!".format((end_time - start_time)))
        
        return self.Forces, self.Displacements

    def get_stress(self, elem_t):
        if elem_t == 2:
            self.stress = np.zeros([len(self.elements),3])

            for i in range(len(self.elements)):
                disp_elem = np.zeros([6,1])
                disp_elem[0::2] = self.Displacements[2*(self.elements[i]-1)]
                disp_elem[1::2] = self.Displacements[2*(self.elements[i]-1)+1]

                b_dot_disp = np.matmul(self.B[i], disp_elem)
                self.stress[i,:] = np.matmul(self.D, b_dot_disp ).T

        elif elem_t == 4:
            self.stress = np.zeros([len(self.elements),6])
            for i in range(len(self.elements)):
                
                # Numero de nodos por elemento por grados de libertad
                disp_elem = np.zeros([12,1])
                disp_elem[0::3] = self.Displacements[3*(self.elements[i]-1)]
                disp_elem[1::3] = self.Displacements[3*(self.elements[i]-1)+1]
                disp_elem[2::3] = self.Displacements[3*(self.elements[i]-1)+2]

                b_dot_disp = np.matmul(self.B[i], disp_elem)
                self.stress[i,:] = np.matmul(self.D, b_dot_disp ).T

        return self.stress
    
    def get_principal_stress(self, type, ElementType):
        if ElementType == 2:
            sigma = self.get_stress(2)
            value = []
            if type=="max":
                for i in range(len(sigma)):
                    value.append(((sigma[i][0]+sigma[i][1])/2) + np.sqrt((((sigma[i][0]-sigma[i][1])/2)**2)+sigma[i][2]**2))
            elif type=="min":
                for i in range(len(sigma)):
                    value.append(((sigma[i][0]+sigma[i][1])/2) - np.sqrt((((sigma[i][0]-sigma[i][1])/2)**2)+sigma[i][2]**2))
        elif ElementType == 4:
            sigma = self.get_stress(4)
            value = []
            for i in range(len(sigma)):
                stress_tens = np.zeros([3, 3])
                stress_tens[0, 0] = sigma[i][0]
                stress_tens[1, 1] = sigma[i][1]
                stress_tens[2, 2] = sigma[i][2]
                stress_tens[1, 2] = sigma[i][3]
                stress_tens[0, 2] = sigma[i][4]
                stress_tens[0, 1] = sigma[i][5]
                stress_tens[2, 1] = sigma[i][3]
                stress_tens[2, 0] = sigma[i][4]
                stress_tens[1, 0] = sigma[i][5]
                stress, autov = np.linalg.eig(stress_tens)
                value.append(np.amax(stress))

        return value
    
    def mass_matrix(self, ElementType, matrix_type, density):
        if ElementType == 5:
            matrix = []
            if matrix_type == 'reduced':
                for elem in range(len(self.elements)):
                    mass_mat = np.zeros([4,4])
                    L2 = np.power(self.long[elem], 2)
                    mass_mat[0,0], mass_mat[1,1], mass_mat[2,2], mass_mat[3,3] = 12, L2, 12, L2
                    matrix.append((density*self.areas[elem]/24)*self.long[elem]*mass_mat)
            elif matrix_type == 'consistent':
                for elem in range(len(self.elements)):
                    A = self.areas[elem]
                    L = self.long[elem]
                    L2 = np.power(self.long[elem], 2)
                    mass_mat = [[156, 22*L, 54, -13*L], 
                                [22*L, 4*L2, 13*L, -3*L2],
                                [54, 13*L, 156, -22*L],
                                [-13*L, -3*L2, -22*L, 4*L2]]
                    mass_mat *= np.asarray(A*L*density/420)
                    matrix.append(mass_mat)
        elif ElementType == 4:
            matrix = []
            if matrix_type == 'consistent':
                mass_mat = np.array([[2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0],
                                     [0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0],
                                     [0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1],
                                     [1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0],
                                     [0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0],
                                     [0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1],
                                     [1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0],
                                     [0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0],
                                     [0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1],
                                     [1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0],
                                     [0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0],
                                     [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2]])

                for elem in range(len(self.elements)):
                   
                    matrix.append(mass_mat * self.volume_array[elem] * density / 20)
            
            else:
                print("Can't obtain mass matrix. Change to 'consistent'")
                exit()
        elif ElementType == 2:
            matrix = []
            if matrix_type == 'consistent':
                mass_mat = np.array([[2, 0, 1, 0, 1, 0],
                                     [0, 2, 0, 1, 0, 1],
                                     [1, 0, 2, 0, 1, 0],
                                     [0, 1, 0, 2, 0, 1],
                                     [1, 0, 1, 0, 2, 0],
                                     [0, 1, 0, 1, 0, 2]])

                for elem in range(len(self.elements)):
                   
                    matrix.append(mass_mat * self.areas[elem] * density / 12)

        else:
            print("Can't obtain mass matrix.Element not configured.")
            exit()

        np_matrix = np.array(matrix)
        return np_matrix
 
