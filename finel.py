import gmesh_handler as gmesh
import finite_elements as fe
import numpy as np
from scipy.linalg import eigh
import boundaries_handler as boundaries
from pyBST import BST
import os
import subprocess
import sys
import platform
from consolemenu import *
from consolemenu.items import *
import pdb
import matplotlib.pyplot as plt
import time
from alive_progress import alive_bar


with open('gmsh_path', 'r') as file:
    gmsh_path = file.readline().rstrip()

if 'Linux' in platform.system():
    os.system('clear')

elif 'Darwin' in platform.system():
    os.system('clear')

else:
    os.system('cls')

if sys.argv[1] == "-f":
    print("Opening {}".format(sys.argv[2]))
    file_name = sys.argv[2]

elif sys.argv[1] == "-v":
    #file = gmesh.read(sys.argv[2])
    #file.read_nodes()
    #elem = file.find_elements(ElementType=2)
    #file.BST_config()
    #subprocess.Popen('./bst')
    print('Feature not available')
    exit()

elif sys.argv[1] == "-g":
    if 'Linux' in platform.system():
        os.system('/usr/local/bin/gmsh -open ./mesh/'+sys.argv[2])
    elif 'Darwin' in platform.system():
        os.system('/Applications/Gmsh.app/Contents/MacOS/gmsh -open ./mesh/'+sys.argv[2])
    else:
        subprocess.Popen(gmsh_path + "/gmsh.exe -open ./mesh/"+sys.argv[2])
    exit()

elif sys.argv[1] == "-c":
    file_name = sys.argv[2]
    while True:
        answer = input("Clean file? (y/n): ")
        if answer == 'y':
            file = gmesh.read(sys.argv[2])
            file.clean()
            print("Done.")
            exit()
        if answer =='n':
            exit()

else:
    file_name = "mesh"

menu_list=["1D Element (Line)", "2D Element (Triangle)", "3D Element (Tetrahedron)", "Beams"]
menu = SelectionMenu(menu_list,"Finite Elements 1.1 - Select element type")

menu.show()
menu.join()

selection = menu.selected_option

if selection == 0:
    #DOF_nodes = 1
    ELEM_type = 1
    menu_list=["1D Solution", "2D Solution", "3D Solution", "--------", "Dynamics"]

elif selection == 1:
    #DOF_nodes = 2
    ELEM_type = 2
    menu_list=["--------", "2D Solution", "3D Solution", "--------", "Dynamics"]

elif selection == 2:
    ELEM_type = 4
    menu_list=["--------", "--------", "3D Solution", "Modal Analysis", "Dynamics"]

elif selection == 3:
    #DOF_nodes = 2
    ELEM_type = 4
    menu_list=["--------", "2D Solution", "3D Solution", "Modal Analysis", "Dynamics"]

elif selection == 4:
    exit()


#menu_list=["1D Solution (1 DOF)", "2D Solution (2 DOF)", "3D Solution (3 DOF)"]
menu = SelectionMenu(menu_list,"Finite Elements 1.1 - Select the degree of freedom per node")

menu.show()
menu.join()

selection2 = menu.selected_option


if selection2 == 0:
    DOF_nodes = 1
    #ELEM_type = 1
    solve = 'Stress-Strain'

elif selection2 == 1:
    DOF_nodes = 2
    #ELEM_type = 2
    solve = 'Stress-Strain'

elif selection2 == 2:
    DOF_nodes = 3
    #ELEM_type = 4
    solve = 'Stress-Strain'

elif selection2 == 3:

    if selection == 2:
        DOF_nodes = 3
    if selection == 3:
        DOF_nodes = 2
    solve = 'Normal-Modes'

elif selection2 == 4:
    
    if selection == 1:
        DOF_nodes = 2
    if selection == 2:
        DOF_nodes = 3
    solve = 'dynamics'

elif selection2 == 5:
    exit()

if solve == 'Stress-Strain':
##########################################################################################################
#                               Solucion tension-deformacion
##########################################################################################################
    file = gmesh.read(file_name)
    nodes_mat = file.read_nodes()
    if ELEM_type == 5:
        # Las vigas se buscan como elementos tipo 1 en gmsh
        elem = file.find_elements(ElementType=1)
    else:
        # El resto de los elementos se buscan de manera normal    
        elem = file.find_elements(ElementType=ELEM_type)

    #file.BST_config()
    #print('Opening Boundaries Selection Tool...')
    #print('BST disabled!')
    #subprocess.call('./bst')
    print('Opening pyBST...')
    print("Please generate boundary condition file!")
    print()
    print('physical type, to set boundaries conditions')
    print('1: line      2: surf     4: volume')
    PHY_type = input('> ')
    BST(file_name=file_name, num_of_nodes = len(nodes_mat), elem_type=int(PHY_type))
    while True:
        answer = input("Continue? (q to quit)")
        if answer == 'y':
            break
        if answer =='q':
            exit()
        else:
            break
    print('Reading boundaries file')
    boundaries = boundaries.read('nodes_selection_file.dat')
    
    f_elem = fe.finite_elements(nodes=nodes_mat, elements=elem, dof_per_node=DOF_nodes)
    if ELEM_type == 5:
        print('Getting beams areas')
        a = f_elem.get_areas(ElementType=1, areas=1e-3)
    else:
        pass

    #boundaries = boundaries.read('nodes_selection_file.dat')
    #lines = file.find_elements(ElementType=1)

    #Lineas con el physical group
    physical_lines = file.find_physical_elements(ElementType=int(PHY_type))
    nodes = boundaries.get_boundaries()
    r, s = boundaries.get_r_s(nodes, dof_per_node=DOF_nodes)

    print()
    forces = {}
    PS = physical_lines[0][0]
    new_force = input('Forces in Physical Group {}: '.format(physical_lines[0][0]))
    forces[physical_lines[0][0]] = float(new_force)

    for i in range(1,len(physical_lines)):
        if physical_lines[i][0] in forces:
            pass
        else:
            new_force = input('Forces in Physical Group {}: '.format(physical_lines[i][0]))
            forces[physical_lines[i][0]] = float(new_force)
            PS = physical_lines[i][0]
    print()
    print('Physical Group : |Force|')
    print(forces)
    print()

    print('Solving...')
    F = boundaries.get_stress(nodes_mat, physical_lines, forces, nodes, ElementType=ELEM_type, dof_per_node=DOF_nodes)
    f_elem = fe.finite_elements(nodes=nodes_mat, elements=elem, dof_per_node=DOF_nodes)
    a = f_elem.get_areas(ElementType=ELEM_type, areas=1)
    mat_rig = f_elem.get_stiff_mat(ElementType=ELEM_type)
    K = f_elem.get_global_mat(ElementType=ELEM_type, matrix_to_assemble=mat_rig)
    D = np.zeros(f_elem.dof_per_node*len(nodes)).reshape(-1,1)
    F, D = f_elem.solve_for(K, F, D, r, s)

    file.clean()
    if ELEM_type == 2:
        #f_elem.get_stress(ELEM_type)
        tensiones = f_elem.get_principal_stress('max', ElementType=2)
        #Cuenta los elementos nuevamente para que las tensiones se escriban correctamente
        countElements = file.find_elements(ElementType=2)
        file.write(D,'Desplazamientos(m)',F,'Fuerzas(N)', tensiones, 'Tensiones(psi)', dof_per_node=2)

    elif ELEM_type == 4:
        #f_elem.get_stress()
        tensiones = f_elem.get_principal_stress('max', ElementType=4)
        #Cuenta los elementos nuevamente para que las tensiones se escriban correctamente
        countElements = file.find_elements(ElementType=4)
        file.write(D,'Desplazamientos(m)',F,'Fuerzas(N)', tensiones, 'Tensiones(psi)', dof_per_node=3)

    elif ELEM_type == 1:

        file.write(D,'Desplazamientos(m)',F,'Fuerzas(N)', F, 'Fuerzas(N)', DOF_nodes)
    
    elif ELEM_type == 5:

        file.write(D,'Desplazamientos(m)',F,'Fuerzas(N)', F, 'Fuerzas(N)', DOF_nodes)
    while True:
        answer = input("Open gmsh? (y/n): ")
        if answer == 'y':
            if 'Linux' in platform.system():
                os.system('/usr/local/bin/gmsh -open ./mesh/'+ sys.argv[2])
            elif 'Darwin' in platform.system():
                os.system('/Applications/Gmsh.app/Contents/MacOS/gmsh -open ./mesh/'+sys.argv[2])
            else:
                subprocess.Popen(gmsh_path + "/gmsh.exe -open ./mesh/"+sys.argv[2])
            break
        if answer =='n':
            exit()
##########################################################################################################
#                               Solucion Modos Normales de vibracion
##########################################################################################################
elif solve == 'Normal-Modes':
    # ElementType == 5 ! ojo que para gmsh son elementos tipo 1
    file = gmesh.read(file_name)
    nodes_mat = file.read_nodes()

    if ELEM_type == 5:
        print('Using element 1 as Beam')
        elem = file.find_elements(ElementType=1)
    else:
        elem = file.find_elements(ElementType=ELEM_type)

    # BST! para los gl empotrados
    #file.BST_config()
    #print('Opening Boundaries Selection Tool...')
    #subprocess.call('./bst')
    
    print('Opening pyBST...')
    print("Please generate boundary condition file!")
    print()
    print('physical type, to set boundaries conditions')
    print('1: line      2: surf     4: volume')
    PHY_type = input('> ')
    BST(file_name=file_name, num_of_nodes = len(nodes_mat), elem_type=int(PHY_type))
    while True:
        answer = input("Continue? (q to quit)")
        if answer == 'y':
            break
        if answer =='q':
            exit()
        else:
            break
    print('Reading boundaries file')
    boundaries = boundaries.read('nodes_selection_file.dat')

    f_elem = fe.finite_elements(nodes=nodes_mat, elements=elem, dof_per_node=DOF_nodes)
    
    if ELEM_type == 5:
        print('Getting beams areas')
        a = f_elem.get_areas(ElementType=1, areas=1e-3)
    else:
        pass
    
    # matriz de rigidez
    mat_rig = f_elem.get_stiff_mat(ElementType=ELEM_type)
    K_glob = f_elem.get_global_mat(ElementType=ELEM_type, matrix_to_assemble=mat_rig)
    #matriz de masa reducida
    #np.savetxt('rigidez_global.txt',K_glob,fmt='%.2f')
    if ELEM_type == 5:
        red_mass = f_elem.mass_matrix(ElementType=ELEM_type, matrix_type='consistent', density=7850)
        m_glob = f_elem.get_global_mat(ElementType=ELEM_type, matrix_to_assemble=red_mass)

    if ELEM_type == 4:
        mass_mat = f_elem.mass_matrix(ElementType=ELEM_type, matrix_type='consistent', density=7850/1e9)
        m_glob = f_elem.get_global_mat(ElementType=ELEM_type, matrix_to_assemble=mass_mat)
    #condiciones de contorno del problema
    nodes = boundaries.get_boundaries()
    r, s = boundaries.get_r_s(nodes, dof_per_node=DOF_nodes)
    
    #calculo de autovalores
    print('solving...')
    start_time = time.time()
    autoval, autov = eigh(K_glob[np.ix_(r,r)], m_glob[np.ix_(r,r)])
    print("Solved in %s seconds" % (time.time() - start_time))
    print('Done')
    #ans= input('Modes to graph: ')
    freq = np.sqrt(autoval) / (2*np.pi)
    disp = np.zeros(DOF_nodes*len(nodes))
    # Imprime las frecuencias de los modos normales
    #MDF-COMMENT    supongo que aca elegis el modo que escribis
    #MDF-COMMENT    num_freq = 2
    max_frec = min(8, autov.shape[1])
    print("Saving {:d} first modes...".format(max_frec))
    #MDF-COMMENT    print('Frequency: {}'.format(freq[0:6]))
    print('Frequency: {}'.format(freq[0:max_frec]))
    #print(str(np.round(freq[num_freq], 1)) + ' Hz')
    #print()
    #pdb.set_trace()
    file.clean()
    for num_freq in range(max_frec):
        print('saving mode {:d}'.format(num_freq))
        for i in range(len(r)):
            disp[r[i]] = autov[i, num_freq]

        if ELEM_type == 4:
            file.clip_write(
                    disp,
                    'Desplazam. MAX {:d}'.format(num_freq),
                    dof_per_node=DOF_nodes
                    )
    
    #Guardar los valores en un archivo para leerlo en dinamica
    disp_n = np.zeros(DOF_nodes*len(nodes))
    for i in range(len(r)):
        disp_n[r[i]] = autov[i, 0] #guarda el modo normal 0 !! cambiar aca      
    with open('modal_analysis.f3d', 'wb') as stored_data:
        np.save(stored_data, disp_n)
        print()
        print('Numpy array saved! can be used in dynamic analysis.')
        print()
####################
    while True:
        answer = input("Open gmsh? (y/n): ")
        if answer == 'y':
            if 'Linux' in platform.system():
                os.system('/usr/local/bin/gmsh -open ./mesh/'+ sys.argv[2])
            elif 'Darwin' in platform.system():
                os.system('/Applications/Gmsh.app/Contents/MacOS/gmsh -open ./mesh/'+sys.argv[2])
            else:
                subprocess.Popen(gmsh_path + "/gmsh.exe -open ./mesh/"+sys.argv[2])
            break
        if answer =='n':
            exit()
##########################################################################################################
#                               Dinamica
##########################################################################################################
elif solve == 'dynamics':
    # ElementType == 5 ! ojo que para gmsh son elementos tipo 1
    file = gmesh.read(file_name)
    nodes_mat = file.read_nodes()

    if ELEM_type == 5:
        print('Using element 1 as Beam')
        elem = file.find_elements(ElementType=1)
    else:
        elem = file.find_elements(ElementType=ELEM_type)

    # BST! para los gl empotrados
    #file.BST_config()
    #print('Opening Boundaries Selection Tool...')
    #subprocess.call('./bst')
    
    print('Opening pyBST...')
    print("Please generate boundary condition file!")
    print()
    print('physical type, to set boundaries conditions')
    print('1: line      2: surf     4: volume')
    PHY_type = input('> ')
    BST(file_name=file_name, num_of_nodes = len(nodes_mat), elem_type=int(PHY_type))
    while True:
        answer = input("Continue? (q to quit)")
        if answer == 'y':
            break
        if answer =='q':
            exit()
        else:
            break
    print('Reading boundaries file')
    boundaries = boundaries.read('nodes_selection_file.dat')

    f_elem = fe.finite_elements(nodes=nodes_mat, elements=elem, dof_per_node=DOF_nodes)
    
    if ELEM_type == 5:
        print('Getting beams areas')
        a = f_elem.get_areas(ElementType=1, areas=1e-3)
    else:
        pass
    
    physical_lines = file.find_physical_elements(ElementType=int(PHY_type))
    nodes = boundaries.get_boundaries()
    r, s = boundaries.get_r_s(nodes, dof_per_node=DOF_nodes)

    print()
    forces = {}
    PS = physical_lines[0][0]
    new_force = input('Forces in Physical Group {}: '.format(physical_lines[0][0]))
    forces[physical_lines[0][0]] = float(new_force)

    for i in range(1,len(physical_lines)):
        if physical_lines[i][0] in forces:
            pass
        else:
            new_force = input('Forces in Physical Group {}: '.format(physical_lines[i][0]))
            forces[physical_lines[i][0]] = float(new_force)
            PS = physical_lines[i][0]
    print()
    print('Physical Group : |Force|')
    print(forces)
    print()
    # matriz de rigidez
    f_elem = fe.finite_elements(nodes=nodes_mat, elements=elem, dof_per_node=DOF_nodes)
    a = f_elem.get_areas(ElementType=ELEM_type, areas=1)
    mat_rig = f_elem.get_stiff_mat(ElementType=ELEM_type)
    K_glob = f_elem.get_global_mat(ElementType=ELEM_type, matrix_to_assemble=mat_rig)
    #matriz de masa reducida
    #np.savetxt('rigidez_global.txt',K_glob,fmt='%.2f')
    if ELEM_type == 5:
        red_mass = f_elem.mass_matrix(ElementType=ELEM_type, matrix_type='consistent', density=7850)
        m_glob = f_elem.get_global_mat(ElementType=ELEM_type, matrix_to_assemble=red_mass)

    if ELEM_type == 4:
        mass_mat = f_elem.mass_matrix(ElementType=ELEM_type, matrix_type='consistent', density=7850/1e9)
        m_glob = f_elem.get_global_mat(ElementType=ELEM_type, matrix_to_assemble=mass_mat)
    
    if ELEM_type == 2:
        mass_mat = f_elem.mass_matrix(ElementType=ELEM_type, matrix_type='consistent', density=7850)
        m_glob = f_elem.get_global_mat(ElementType=ELEM_type, matrix_to_assemble=mass_mat)
    #condiciones de contorno del problema
    nodes = boundaries.get_boundaries()
    r, s = boundaries.get_r_s(nodes, dof_per_node=DOF_nodes)
   
   ################################### 
   # DEFINIR ALGUNOS HIPERPARAMETROS #
   ################################### 0.0025
    steps = 500
    delta_t = 1e-5
    delta_t2 = delta_t**2

   # Calculos
    print()
    print('Solving...')
    start_time = time.time()
    F = boundaries.get_stress(nodes_mat, physical_lines, forces, nodes, ElementType=ELEM_type, dof_per_node=DOF_nodes)
    
    f = np.zeros([steps, len(F)])
    a = np.zeros([steps, len(F)])
    v = np.zeros([steps, len(F)])
    d = np.zeros([steps, len(F)])

    ####################### Fuerza constante
    #for i in range(steps):
    #    f[i, ] = (1)*F.T
    ####################### Desplazamientos de modos normales de vibracion
    
    with open('modal_analysis.f3d', 'rb') as stored_data:
        disp = np.load(stored_data)
        for i in range(len(F)):
            d[0,i] = disp[i]
    

    #calculos
    Minv = np.linalg.inv(m_glob[np.ix_(r, r)])
    # Calculo de la aceleracion inicial
    #a[0, r] = np.matmul(Minv, (f[0, r]-np.matmul(K_glob[np.ix_(r,)], d[0])))
    a[0, r] = np.matmul(Minv, (f[0, r]-np.matmul(K_glob[np.ix_(r,r)], d[0,r])))
    #print('a0')
    #print(a)
    #Este seria el d-1
    d_1 = ( d[0] - delta_t*v[0] + 0.5*(delta_t2)*a[0] )
    #pdb.set_trace()
    #d1
    #d[1, r] = np.matmul(Minv, 2*np.matmul(m_glob[np.ix_(r, r)], d[0, r])-np.matmul(m_glob[np.ix_(r, r)], d_1[r])+
    #   delta_t2*(f[0, r]-np.matmul(K_glob[np.ix_(r,)], d[0,])-np.matmul(m_glob[np.ix_(r, s)], a[0,s])))
    
    #term = delta_t2*f[0, r] + np.matmul(2*m_glob[np.ix_(r, r)] - delta_t2*K_glob[np.ix_(r,r)], d[0,r]) - np.matmul(m_glob[np.ix_(r, r)], d_1[r])
    #d[1, r] = np.matmul(Minv,term)
    print('Calculating dynamic position')
    with alive_bar(steps, bar = 'filling', spinner = 'dots_reverse') as bar:
        '''
        #d2 en adelante
        for i in range(2, steps):
            d[i, r] = np.matmul(Minv, 2*np.matmul(m_glob[np.ix_(r, r)], d[i-1, r])-np.matmul(m_glob[np.ix_(r, r)], d[i-2, r])+
                      delta_t2*(f[i-1, r]-np.matmul(K_glob[np.ix_(r,)], d[i-1,])-np.matmul(m_glob[np.ix_(r, s)], a[i-1, s])))
            bar()
        '''

        
        for i in range(2, steps):
            #term = delta_t2*f[i-1, r] + np.matmul(2*m_glob[np.ix_(r, r)] - delta_t2*K_glob[np.ix_(r,r)], d[i-1,r]) - np.matmul(m_glob[np.ix_(r, r)], d[i-2, r])
            d[i, r] = np.matmul(Minv, delta_t2*f[i-1, r] + np.matmul(2*m_glob[np.ix_(r, r)] - delta_t2*K_glob[np.ix_(r,r)] , d[i-1,r]) - np.matmul(m_glob[np.ix_(r, r)], d[i-2, r]))
            bar()
    print("Solved in %s seconds" % (time.time() - start_time))
    print('Done')
    print()
    file.clean()
    print()
    print('Saving clip...')

    if ELEM_type == 4:
        file.dynamic_clip_write(
                d,
                'Displacement',
                dof_per_node=DOF_nodes,
                tags = steps
                )
    if ELEM_type == 2:
        file.dynamic_clip_write(
                d,
                'Displacement',
                dof_per_node=DOF_nodes,
                tags=steps
                )
    while True:
        answer = input("Open gmsh? (y/n): ")
        if answer == 'y':
            if 'Linux' in platform.system():
                os.system('/usr/local/bin/gmsh -open ./mesh/'+ sys.argv[2])
            elif 'Darwin' in platform.system():
                os.system('/Applications/Gmsh.app/Contents/MacOS/gmsh -open ./mesh/'+sys.argv[2])
            else:
                subprocess.Popen(gmsh_path + "/gmsh.exe -open ./mesh/"+sys.argv[2])
            break
        if answer =='n':
            exit()
