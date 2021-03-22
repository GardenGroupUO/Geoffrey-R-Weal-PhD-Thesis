print 'Loading ase Atom and Atoms'
from ase import Atom, Atoms
from ase.io import read as ASE_read
from ase.io import write as ASE_write
from ase.visualize import view
print 'Loading ase.cluster Packages'
from ase.cluster.icosahedron import Icosahedron
from ase.cluster.decahedron import Decahedron
from ase.cluster.octahedron import Octahedron
print 'Loading numpy'
import numpy as np
print 'Loading view'
from ase.visualize import view
print 'Loading copy, sys, os, shutil'
import copy, sys, os, shutil
print 'Loading matplotlib'
import matplotlib
matplotlib.use('Agg')
print 'Loading sin and cos'
from math import cos, sin, radians
print 'Loading Structural_Recognition_Program from: '
'''
add_to_import_list = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir, os.path.pardir, os.path.pardir, 'Main_Program')
print add_to_import_list
sys.path.append(add_to_import_list)
'''
from Structural_Recognition_Program import Structural_Recognition_Program
print 'Loading RunMinimisation'
from RunMinimisation import RunMinimisation
print 'Running program'

def delete():
    for item in os.listdir('.'):
        if not item == os.path.basename(__file__) and not os.path.isdir(item) and not item.startswith('submit.sl'):# and not item.endswith('.traj'):
            os.remove(item)
            print 'Deleting ' + str(item)
        elif item == 'Diagrams_Signature_Analysis':
            shutil.rmtree(item)
            print 'Deleting ' + str(item)
        elif item == 'Diagrams':
            shutil.rmtree(item)
            print 'Deleting ' + str(item)
        elif item == 'Text_Information':
            shutil.rmtree(item)
            print 'Deleting ' + str(item)
        elif item == 'Text_Information_Signature_Analysis':
            shutil.rmtree(item)
            print 'Deleting ' + str(item)
        else:
            print 'Did not delete ' + str(item)

print '---------------------------'
delete()
print '---------------------------'

def optimise(cluster_name):
    cluster = ASE_read(cluster_name+'.traj')
    cluster_optimised = RunMinimisation(cluster).get_cluster()
    cluster_name_opt = cluster_name + '_Opt'
    ASE_write(cluster_name_opt+'.traj',cluster_optimised)
    return cluster_name_opt, cluster_optimised

def shift_x(normalised_length,length_1,length_2,angle,layer='1st_layer',Use_optimised_structure=False):
    name = 'Au38_Octa_'+str(length_1)+'_'+str(length_2)+'_'+str(angle)
    angle = float(radians(float(angle)))
    if Use_optimised_structure:
        cluster = copy.deepcopy(cluster_main_Opt)
        cluster.set_calculator(None)
        name += '_Opt'
    else:
        cluster = copy.deepcopy(cluster_main)
    x_neighbour_length = cluster[16].x - cluster[32].x
    z_neighbour_length = cluster[37].z - cluster[23].z
    #neighbour_length = (x_neighbour_length**2.0 + z_neighbour_length**2.0)**0.5
    neighbour_length = normalised_length
    layer_1 = [24,19,21,35]
    layer_2 = [23,22,3,18,17,34,36,37,20]
    for index in layer_1:
    	cluster[index].x += length_1*cos(angle)*neighbour_length
    	cluster[index].z += length_1*sin(angle)*neighbour_length
    #if layer == '2nd_layer':
    for index in layer_2:
        cluster[index].x += length_2*cos(angle)*neighbour_length
        cluster[index].z += length_2*sin(angle)*neighbour_length
    ASE_write(name+'.traj',cluster,'traj')
    return name

rCut_low = 2.8840
rCut_high = round(rCut_low * (2.0)**0.5,4)
rCut_resolution = 0.0005
mode = 'total'

def run_CNA_program(Cluster_1_name,Cluster_2_name):
    Cluster_1 = ASE_read(Cluster_1_name+'.traj')
    print 'Made ' + str(Cluster_1_name)
    Cluster_2 = ASE_read(Cluster_2_name+'.traj')
    print 'Made ' + str(Cluster_2_name)
    Structural_Recognition_Program(Cluster_1,Cluster_2,rCut_low,rCut_high,rCut_resolution,mode,name_1=Cluster_1_name,name_2=Cluster_2_name,recognise_multimetallic=False,print_plots=True)

symbol = 'Au'
cluster_main = Octahedron(symbol,4,1)
length_1 = 0; length_2 = 0; angle = 0;
name_original = shift_x(0,0,0,0)
name_original_Opt, cluster_main_Opt = optimise(name_original)

############################################################################################################
############################################################################################################
############################################################################################################
layer = '1st_layer'
############################################################################################################
############################################################################################################
############################################################################################################
normal_length = 2.040
length_1s = np.arange(0.1,1.01,0.1)
length_2 = 0; angle = 0;
for length_1 in length_1s:
    name_comparing = shift_x(normal_length,length_1,length_2,angle,layer=layer,Use_optimised_structure=False)
    run_CNA_program(name_original,name_comparing)
    name_comparing_Opt = shift_x(normal_length,length_1,length_2,angle,layer=layer,Use_optimised_structure=True)
    run_CNA_program(name_original_Opt,name_comparing_Opt)

######################################################
######################################################
######################################################
######################################################
normal_length = 2.885
length_1s = np.arange(0.1,1.01,0.1)
length_2 = 0; angle = 45;
for length_1 in length_1s:
    name_comparing = shift_x(normal_length,length_1,length_2,angle,layer=layer,Use_optimised_structure=False)
    run_CNA_program(name_original,name_comparing)
    name_comparing_Opt = shift_x(normal_length,length_1,length_2,angle,layer=layer,Use_optimised_structure=True)
    run_CNA_program(name_original_Opt,name_comparing_Opt)

############################################################################################################
############################################################################################################
############################################################################################################
layer = '2nd_layer'
############################################################################################################
############################################################################################################
############################################################################################################
normal_length = 2.040
length_1s = np.arange(0.1,1.01,0.1)
length_2 = 0; angle = 0;
for length_1 in length_1s:
    length_2 = length_1
    name_comparing = shift_x(normal_length,length_1,length_2,angle,layer=layer,Use_optimised_structure=False)
    run_CNA_program(name_original,name_comparing)
    name_comparing_Opt = shift_x(normal_length,length_1,length_2,angle,layer=layer,Use_optimised_structure=True)
    run_CNA_program(name_original_Opt,name_comparing_Opt)

######################################################
######################################################
######################################################
######################################################
normal_length = 2.885
length_1s = np.arange(0.1,1.01,0.1)
length_2 = 0; angle = 45;
for length_1 in length_1s:
    length_2 = length_1
    name_comparing = shift_x(normal_length,length_1,length_2,angle,layer=layer,Use_optimised_structure=False)
    run_CNA_program(name_original,name_comparing)
    name_comparing_Opt = shift_x(normal_length,length_1,length_2,angle,layer=layer,Use_optimised_structure=True)
    run_CNA_program(name_original_Opt,name_comparing_Opt)

############################################################################################################
############################################################################################################
############################################################################################################