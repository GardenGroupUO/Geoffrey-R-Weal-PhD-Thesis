print 'Loading ase Atom and Atoms'
from ase import Atom, Atoms
from ase.io import read as ASE_read
from ase.io import write as ASE_write
print 'Loading ase view'
from ase.visualize import view
print 'Loading ase.cluster Packages'
from ase.cluster.icosahedron import Icosahedron
from ase.cluster.decahedron import Decahedron
from ase.cluster.octahedron import Octahedron
print 'Loading numpy'
import numpy as np
from numpy import array, pi, subtract, multiply, matmul, divide, dot, cross
print 'Loading view'
from ase.visualize import view
print 'Loading copy, sys, os, shutil'
import copy, sys, os, shutil
print 'Loading matplotlib'
import matplotlib
matplotlib.use('Agg')
print 'Loading sin and cos'
from math import cos, sin, radians
print 'Loading Structural_Recognition_Program'
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


def move_cluster_to_COM(cluster):
    COM = cluster.get_center_of_mass()
    for atom in cluster:
        atom.position -= COM
    return cluster

def move_cluster(cluster,centre_point):
    cluster = copy.deepcopy(cluster)
    for atom in cluster:
        atom.position -= centre_point
    return cluster

def Twisting(input_cluster,Twisting_Axis,maximum_twisting_angle,atoms_to_twist=[]):
	output_cluster = copy.deepcopy(input_cluster)
	if atoms_to_twist == []:
		atoms_to_twist = range(len(output_cluster))
	def get_unit_vector(vector):
		norm = (vector.transpose().dot(vector))[0][0]**0.5
		unit_vector = divide(vector,norm)
		return unit_vector
	Twisting_Axis = array(Twisting_Axis)
	Twisting_Axis.shape = (3,1)
	Twisting_Axis_unit_vector = get_unit_vector(Twisting_Axis)
	COM = array(output_cluster.get_center_of_mass())
	COM.shape = (3,1)
	
	Identity = array([[1,0,0],[0,1,0],[0,0,1]])
	Tx = Twisting_Axis_unit_vector[0][0]; Ty = Twisting_Axis_unit_vector[1][0]; Tz = Twisting_Axis_unit_vector[2][0]
	ucross = array([[0,-Tz,Ty],[Tz,0,-Tx],[-Ty,Tx,0]])
	ucrossu = matmul(Twisting_Axis_unit_vector,Twisting_Axis_unit_vector.transpose())
	
	t_values = []
	for index in range(len(output_cluster)):
		atom_position = output_cluster[index].position
		atom_position.shape = (3,1)
		XminusC = subtract(atom_position,COM)
		t_value = float(dot(Twisting_Axis_unit_vector.transpose(),(XminusC))[0][0])
		t_values.append(t_value)
	
	t_value_max = max(t_values)
	t_value_min = min(t_values)
	
	for index in atoms_to_twist:
		atom_position = array(output_cluster[index].position)
		atom_position.shape = (3,1)
		t_value = float(t_values[index])
		angle = (float(maximum_twisting_angle)/180.0)*pi
		if t_value > 0:
			angle *= t_value/t_value_max
		else:
			angle *= -abs(t_value/t_value_min)
		Rotation_Matrix = multiply(cos(angle),Identity) + multiply(sin(angle),ucross) + multiply((1.0 - cos(angle)),ucrossu)
		new_atom_position = matmul(Rotation_Matrix,atom_position)
		new_atom_position.shape = (1,3)
		new_atom_position = np.ndarray.tolist(new_atom_position)[0]
		output_cluster[index].position = new_atom_position
	
	return output_cluster
	
rCut_low = 2.8840
rCut_high = round(rCut_low * (2.0)**0.5,4)
rCut_resolution = 0.0005
mode = 'total'

symbol = 'Au'
latticeconstant = rCut_high
cluster_main = Octahedron(symbol,4,1,latticeconstant=latticeconstant)
cluster_main = move_cluster_to_COM(cluster_main)
original_cluster_name = 'Au38_original'
ASE_write(original_cluster_name+'.traj',cluster_main)
original_cluster_name_Opt, cluster_main_Opt = optimise(original_cluster_name)

def run_CNA_program(Cluster_1_name,Cluster_2_name):
    Cluster_1 = ASE_read(Cluster_1_name+'.traj')
    print 'Made ' + str(Cluster_1_name)
    Cluster_2 = ASE_read(Cluster_2_name+'.traj')
    print 'Made ' + str(Cluster_2_name)
    Structural_Recognition_Program(Cluster_1,Cluster_2,rCut_low,rCut_high,rCut_resolution,mode,name_1=Cluster_1_name,name_2=Cluster_2_name,recognise_multimetallic=False,print_plots=True)

######################################################
#import pdb; pdb.set_trace()
######################################################
maximum_twisting_angles = range(5,31,5)
Twisting_Axis = [0,0,1]
for maximum_twisting_angle in maximum_twisting_angles:
	cluster_name = 'Au38_twist_z_axis_'+str(maximum_twisting_angle)+'_degrees'
	Twisted_cluster = Twisting(cluster_main,Twisting_Axis,maximum_twisting_angle)
	ASE_write(cluster_name+'.traj',Twisted_cluster)
	run_CNA_program(original_cluster_name,cluster_name)
	#cluster_name_Opt = optimise(cluster_name)
	#run_CNA_program(original_cluster_name_Opt,cluster_name_Opt)
	##### Optimised Structures #####
	cluster_name_Opt = 'Au38_twist_z_axis_'+str(maximum_twisting_angle)+'_degrees_Opt'
	Twisted_cluster_Opt = Twisting(cluster_main_Opt,Twisting_Axis,maximum_twisting_angle)
	Twisted_cluster_Opt.set_calculator(None)
	ASE_write(cluster_name_Opt+'.traj',Twisted_cluster_Opt)
	run_CNA_program(original_cluster_name_Opt,cluster_name_Opt)
############################################################################################################
############################################################################################################
############################################################################################################
maximum_twisting_angles = range(5,31,5)
Twisting_Axis = [1,1,1]
for maximum_twisting_angle in maximum_twisting_angles:
	cluster_name = 'Au38_twist_[1,1,1]_axis_'+str(maximum_twisting_angle)+'_degrees'
	Twisted_cluster = Twisting(cluster_main,Twisting_Axis,maximum_twisting_angle)
	ASE_write(cluster_name+'.traj',Twisted_cluster)
	run_CNA_program(original_cluster_name,cluster_name)
	#cluster_name_Opt = optimise(cluster_name)
	#run_CNA_program(original_cluster_name_Opt,cluster_name_Opt)
	##### Optimised Structures #####
	cluster_name_Opt = 'Au38_twist_[1,1,1]_axis_'+str(maximum_twisting_angle)+'_degrees_Opt'
	Twisted_cluster_Opt = Twisting(cluster_main_Opt,Twisting_Axis,maximum_twisting_angle)
	Twisted_cluster_Opt.set_calculator(None)
	ASE_write(cluster_name_Opt+'.traj',Twisted_cluster_Opt)
	run_CNA_program(original_cluster_name_Opt,cluster_name_Opt)
############################################################################################################
############################################################################################################
############################################################################################################
maximum_twisting_angles = range(5,31,5)
Twisting_Axis = [1,1,1]
atoms_to_twist = [24,35,23,37,36,32,33]
for maximum_twisting_angle in maximum_twisting_angles:
	cluster_name = 'Au38_twist_[1,1,1]_axis_'+str(maximum_twisting_angle)+'_degrees_oneFace'
	Twisted_cluster = Twisting(cluster_main,Twisting_Axis,maximum_twisting_angle,atoms_to_twist=atoms_to_twist)
	ASE_write(cluster_name+'.traj',Twisted_cluster)
	run_CNA_program(original_cluster_name,cluster_name)
	#cluster_name_Opt = optimise(cluster_name)
	#run_CNA_program(original_cluster_name_Opt,cluster_name_Opt)
	##### Optimised Structures #####
	cluster_name_Opt = 'Au38_twist_[1,1,1]_axis_'+str(maximum_twisting_angle)+'_degrees_oneFace_Opt'
	Twisted_cluster_Opt = Twisting(cluster_main_Opt,Twisting_Axis,maximum_twisting_angle,atoms_to_twist=atoms_to_twist)
	Twisted_cluster_Opt.set_calculator(None)
	ASE_write(cluster_name_Opt+'.traj',Twisted_cluster_Opt)
	run_CNA_program(original_cluster_name_Opt,cluster_name_Opt)
############################################################################################################
############################################################################################################
############################################################################################################