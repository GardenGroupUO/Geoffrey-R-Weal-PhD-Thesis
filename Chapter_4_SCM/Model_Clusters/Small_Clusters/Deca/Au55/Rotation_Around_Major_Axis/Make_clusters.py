print 'Loading ase Atom and Atoms'
from ase import Atom, Atoms
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
print 'Loading copy'
import copy
print 'Loading sin and cos'
from math import cos, sin, radians
print "Loading Other Things"
import sys, os, shutil
print 'Loading matplotlib'
import matplotlib
matplotlib.use('Agg')

'''
add_to_import_list = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir, os.path.pardir, os.path.pardir, 'Main_Program')
print add_to_import_list
sys.path.append(add_to_import_list)
'''
from Structural_Recognition_Program import Structural_Recognition_Program
from ase.cluster.decahedron import Decahedron
from RunMinimisation import RunMinimisation
#from asap3.Internal.BuiltinPotentials import Gupta
#from ase.optimize import FIRE
from ase.io import read as ASE_read
from ase.visualize import view
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
	return cluster_name_opt

symbol = 'Au'
r0 = 2.8840
latticeconstant = r0 * (2.0)**0.5
cluster_main = Decahedron(symbol,3,3,0,latticeconstant=latticeconstant)
cluster_main_optimised = RunMinimisation(cluster_main).get_cluster()
cluster_main_optimised.set_calculator(None)

'''
def rotate_top(angle,atom_to_remove,name):
    angle = radians(float(angle))
    angle = float(angle)
    cluster = copy.deepcopy(cluster_main)
    rotate_around_x = cluster[3].x; rotate_around_y = cluster[3].y; 
    first_layer_to_rotate = [19,23,7,11,15]
    second_layer_to_rotate = [27,8,54,24,51,48,20,45,42,16,39,36,12,33,30]
    for index in (first_layer_to_rotate+second_layer_to_rotate):
	   old_x = cluster[index].x - rotate_around_x; old_y = cluster[index].y - rotate_around_y
	   new_x = old_x*cos(angle) - old_y*sin(angle) + (rotate_around_x)
	   new_y = old_x*sin(angle) + old_y*cos(angle) + (rotate_around_y)
	   cluster[index].x = new_x; cluster[index].y = new_y
    for index in atom_to_remove:
        del cluster[index]
    #view(cluster)
    #import pdb; pdb.set_trace()
    ASE_write(name+'.traj',cluster,'traj')
'''

def rotate_top(angle,atom_to_remove,name,only_rotate_second_later=False,Optimise=False):
    angle = radians(float(angle))
    angle = float(angle)
    if not Optimise:
        cluster = copy.deepcopy(cluster_main)
    else:
        cluster = copy.deepcopy(cluster_main_optimised)
    rotate_around_x = cluster[3].x; rotate_around_y = cluster[3].y; 
    first_layer_to_rotate = [19,23,7,11,15]
    second_layer_to_rotate = [27,8,54,24,51,48,20,45,42,16,39,36,12,33,30]
    rotate_atoms = second_layer_to_rotate
    if not only_rotate_second_later:
        rotate_atoms += first_layer_to_rotate
    for atom in atom_to_remove:
        del cluster[atom]
        for index2 in range(len(rotate_atoms)):
            if atom < rotate_atoms[index2]:
                rotate_atoms[index2] -= 1
    if Optimise:
        cluster = RunMinimisation(cluster).get_cluster()
        cluster.set_calculator(None)
    for index in rotate_atoms:
       old_x = cluster[index].x - rotate_around_x; old_y = cluster[index].y - rotate_around_y
       new_x = old_x*cos(angle) - old_y*sin(angle) + (rotate_around_x)
       new_y = old_x*sin(angle) + old_y*cos(angle) + (rotate_around_y)
       cluster[index].x = new_x; cluster[index].y = new_y
    #view(cluster)
    #import pdb; pdb.set_trace()
    if Optimise:
        name += '_Opt'
    ASE_write(name+'.traj',cluster,'traj')

rCut_low = 2.8840
rCut_high = round(rCut_low * (2.0)**0.5,4)
rCut_resolution = 0.0005
mode = 'total'


def run_CNA_program(Cluster_1_name,Cluster_2_name):
    Cluster_1 = ASE_read(Cluster_1_name+'.traj')
    #view(Cluster_1)
    print 'Made ' + str(Cluster_1_name)
    Cluster_2 = ASE_read(Cluster_2_name+'.traj')
    #view(Cluster_2)
    print 'Made ' + str(Cluster_2_name)
    Structural_Recognition_Program(Cluster_1,Cluster_2,rCut_low,rCut_high,rCut_resolution,mode,name_1=Cluster_1_name,name_2=Cluster_2_name,recognise_multimetallic=False,print_plots=True)

deca_types = []
#########################################################################
atom_to_remove = []
angle = (0.0/6.0)*36
start_name = 'Au55_rotate_'
name = start_name+str(angle)
#deca_types.append(name)
rotate_top(angle,atom_to_remove,name)
name_opt = optimise(name)

#########################################################################
start_name = 'Au54_rotate_'
deca_ar_type1 = '_100-111_corner'
atom_to_remove = [39]
angle = (0.0/6.0)*36
name = start_name+str(angle)+deca_ar_type1
deca_types.append(name)
rotate_top(angle,atom_to_remove,name)
name_opt = optimise(name)

#########################################################################
start_name = 'Au54_rotate_'
deca_ar_type2 = '_100-111_edge'
atom_to_remove = [42]
angle = (0.0/6.0)*36
name = start_name+str(angle)+deca_ar_type2
deca_types.append(name)
rotate_top(angle,atom_to_remove,name)
name_opt = optimise(name)

#########################################################################
start_name = 'Au54_rotate_'
deca_ar_type3 = '_111_corner'
atom_to_remove = [4]
angle = (0.0/6.0)*36
name = start_name+str(angle)+deca_ar_type3
deca_types.append(name)
rotate_top(angle,atom_to_remove,name)
name_opt = optimise(name)

#########################################################################
start_name = 'Au54_rotate_'
deca_ar_type4 = '_111_edge'
atom_to_remove = [16]
angle = (0.0/6.0)*36
name = start_name+str(angle)+deca_ar_type4
deca_types.append(name)
rotate_top(angle,atom_to_remove,name)
name_opt = optimise(name)

#########################################################################
start_name = 'Au54_rotate_'
deca_ar_type5 = '_100_edge'
atom_to_remove = [38]
angle = (0.0/6.0)*36
name = start_name+str(angle)+deca_ar_type5
deca_types.append(name)
rotate_top(angle,atom_to_remove,name)
name_opt = optimise(name)

#########################################################################
start_name = 'Au54_rotate_'
deca_ar_type6 = '_100_face'
atom_to_remove = [41]
angle = (0.0/6.0)*36
name = start_name+str(angle)+deca_ar_type6
deca_types.append(name)
rotate_top(angle,atom_to_remove,name)
name_opt = optimise(name)

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
# Compare structural features of clusters
start_name = 'Au54_rotate_'
angle = 0.0

for index1 in range(len(deca_types)):
    for index2 in range(index1+1,len(deca_types)):
        print deca_types[index1]
        print deca_types[index2]
        run_CNA_program(deca_types[index1],deca_types[index2])
        run_CNA_program(deca_types[index1]+'_Opt',deca_types[index2]+'_Opt')



