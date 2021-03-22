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
	cluster_optimised.set_calculator(None)
	cluster_name_opt = cluster_name + '_Opt'
	ASE_write(cluster_name_opt+'.traj',cluster_optimised)
	return cluster_name_opt, cluster_optimised
	
symbol = 'Au'
r0 = 2.8840
latticeconstant = r0 * (2.0)**0.5
cluster_main = Decahedron(symbol,3,2,0,latticeconstant=latticeconstant)
cluster_main_optimised = RunMinimisation(cluster_main).get_cluster()
cluster_main_optimised.set_calculator(None)

def rotate_top(angle,atom_to_remove,name,only_rotate_second_later=False,Optimise=False):
    angle = radians(float(angle))
    angle = float(angle)
    if not Optimise:
    	cluster = copy.deepcopy(cluster_main)
    else:
    	cluster = copy.deepcopy(cluster_main_optimised)
    rotate_around_x = cluster[3].x; rotate_around_y = cluster[3].y; 
    first_layer_to_rotate = [9,6,18,15,12]
    second_layer_to_rotate = [24,22,20,38,36,34,32,30,28,26]
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
	
max_rotation = 18

#########################################################################
atom_to_remove = []
angle = (0.0/6.0)*max_rotation
start_name = 'Au39_rotate_'
name_0 = start_name+str(angle)
rotate_top(angle,atom_to_remove,name_0)
name_0_Opt, cluster_main_Opt = optimise(name_0)
angles = np.arange(1.0, 6.1, 1.0)
angles *= max_rotation/6.0
for angle in angles:
	name = start_name+str(angle)
	##########################################
	rotate_top(angle,atom_to_remove,name)
	run_CNA_program(start_name+'0.0',name)
	#name_opt = optimise(name)
	#run_CNA_program(name_0_opt,name_opt)
	#### Optimised Structure ####
	rotate_top(angle,atom_to_remove,name,Optimise=True)
	run_CNA_program(name_0_Opt,name+'_Opt')

#########################################################################
atom_to_remove = []
angle = (0.0/6.0)*max_rotation
start_name = 'Au39_rotate_'
suffex = '_only_surface_atoms'
name_0 = start_name+str(angle)+suffex
rotate_top(angle,atom_to_remove,name_0,only_rotate_second_later=False)
name_0_Opt, cluster_main_Opt = optimise(name_0)
angles = np.arange(1.0, 6.1, 1.0)
angles *= max_rotation/6.0
for angle in angles:
	name = start_name+str(angle)+suffex
	##########################################
	rotate_top(angle,atom_to_remove,name,only_rotate_second_later=False)
	run_CNA_program(name_0,name)
	#name_opt = optimise(name)
	#run_CNA_program(name_0_opt,name_opt)
	#### Optimised Structure ####
	rotate_top(angle,atom_to_remove,name,only_rotate_second_later=False,Optimise=True)
	run_CNA_program(name_0_Opt,name+'_Opt')
	
#########################################################################
atom_to_remove = [32]
angle = (0.0/6.0)*max_rotation
start_name = 'Au38_rotate_'
deca_ar_type = '_100_corner'
name_0 = start_name+str(angle)+deca_ar_type
rotate_top(angle,atom_to_remove,name_0)
name_0_Opt, cluster_main_Opt = optimise(name_0)
angles = np.arange(1.0, 6.1, 1.0)
angles *= max_rotation/6.0
for angle in angles:
	name = start_name+str(angle)+deca_ar_type
	##########################################
	rotate_top(angle,atom_to_remove,name)
	run_CNA_program(name_0,name)
	#name_opt = optimise(name)
	#run_CNA_program(name_0_opt,name_opt)
	#### Optimised Structure ####
	rotate_top(angle,atom_to_remove,name,Optimise=True)
	run_CNA_program(name_0_Opt,name+'_Opt')
	
#########################################################################
atom_to_remove = [34]
angle = (0.0/6.0)*max_rotation
start_name = 'Au38_rotate_'
deca_ar_type = '_100_edge'
name_0 = start_name+str(angle)+deca_ar_type
rotate_top(angle,atom_to_remove,name_0)
name_0_Opt, cluster_main_Opt = optimise(name_0)
angles = np.arange(1.0, 6.1, 1.0)
angles *= max_rotation/6.0
for angle in angles:
	name = start_name+str(angle)+deca_ar_type
	##########################################
	rotate_top(angle,atom_to_remove,name)
	run_CNA_program(name_0,name)
	#name_opt = optimise(name)
	#run_CNA_program(name_0_opt,name_opt)
	#### Optimised Structure ####
	rotate_top(angle,atom_to_remove,name,Optimise=True)
	run_CNA_program(name_0_Opt,name+'_Opt')


#########################################################################
atom_to_remove = [15]
angle = (0.0/6.0)*max_rotation
start_name = 'Au38_rotate_'
deca_ar_type = '_111_edge'
name_0 = start_name+str(angle)+deca_ar_type
rotate_top(angle,atom_to_remove,name_0)
name_0_Opt, cluster_main_Opt = optimise(name_0)
angles = np.arange(1.0, 6.1, 1.0)
angles *= max_rotation/6.0
for angle in angles:
	name = start_name+str(angle)+deca_ar_type
	##########################################
	rotate_top(angle,atom_to_remove,name)
	run_CNA_program(name_0,name)
	#name_opt = optimise(name)
	#run_CNA_program(name_0_opt,name_opt)
	#### Optimised Structure ####
	rotate_top(angle,atom_to_remove,name,Optimise=True)
	run_CNA_program(name_0_Opt,name+'_Opt')

#########################################################################
atom_to_remove = [3]
angle = (0.0/6.0)*max_rotation
start_name = 'Au38_rotate_'
deca_ar_type = '_111_corner'
name_0 = start_name+str(angle)+deca_ar_type
rotate_top(angle,atom_to_remove,name_0)
name_0_Opt, cluster_main_Opt = optimise(name_0)
angles = np.arange(1.0, 6.1, 1.0)
angles *= max_rotation/6.0
for angle in angles:
	name = start_name+str(angle)+deca_ar_type
	##########################################
	rotate_top(angle,atom_to_remove,name)
	run_CNA_program(name_0,name)
	#name_opt = optimise(name)
	#run_CNA_program(name_0_opt,name_opt)
	#### Optimised Structure ####
	rotate_top(angle,atom_to_remove,name,Optimise=True)
	run_CNA_program(name_0_Opt,name+'_Opt')

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
# Compare structural features of clusters

start_name = 'Au38_rotate_'
angle = 0.0

deca_ar_type1 = '_100_corner'
deca_ar_type2 = '_100_edge'
deca_ar_type3 = '_111_edge'
deca_ar_type4 = '_111_corner'

deca_types = [deca_ar_type1,deca_ar_type2,deca_ar_type3,deca_ar_type4]
for index1 in range(len(deca_types)):
    for index2 in range(index1+1,len(deca_types)):
        name1 = start_name+str(angle)+deca_types[index1]
        name2 = start_name+str(angle)+deca_types[index2]
        print name1
        print name2
        run_CNA_program(name1,name2)
        run_CNA_program(name1+'_Opt',name2+'_Opt')

#import pdb; pdb.set_trace()



