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
from numpy import array, add, subtract, multiply, divide, matmul, dot, cross
print 'Loading view'
from ase.visualize import view
print 'Loading copy, sys, os'
import copy, sys, os
print 'Loading matplotlib'
import matplotlib
matplotlib.use('Agg')
print 'Loading sin and cos'
from math import cos, sin, radians, pi
print 'Loading Structural_Recognition_Program.'
from Structural_Recognition_Program import Structural_Recognition_Program
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

def move_cluster_to_COM(cluster):
    COM = cluster.get_center_of_mass()
    for atom in cluster:
        atom.position -= COM
    return cluster

symbol = 'Au'
cluster_main = Octahedron(symbol,4,1)
cluster_main = move_cluster_to_COM(cluster_main)
original_cluster_name = 'Au38_original'
ASE_write(original_cluster_name+'.traj',cluster_main)

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

######################################################
######################################################
cluster_name = 'Au38_Face_Sheered'
run_CNA_program(original_cluster_name,cluster_name)
######################################################