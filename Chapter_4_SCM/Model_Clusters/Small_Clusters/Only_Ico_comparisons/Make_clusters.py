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
from numpy import array, pi, add, subtract, multiply, matmul, divide, dot, cross
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
from Structural_Recognition_Program import Structural_Recognition_Program, Cluster_CNA_Profile
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

symbol = 'Au'

rCut_low = 2.8840
rCut_high = round(rCut_low * (2.0)**0.5,4)
rCut_resolution = 0.0005
mode = 'total'

latticeconstant = rCut_high

def makeMotif(motif_inputs,latticeconstant,remove_atoms=[]):
	if isinstance(motif_inputs,int) or len(motif_inputs) == 1:
		if isinstance(motif_inputs,int):
			motif_inputs = [motif_inputs]
		cluster = Icosahedron(symbol,motif_inputs[0],latticeconstant=latticeconstant)
		name = 'Ico'
	elif len(motif_inputs) == 2:
		cluster = Octahedron(symbol,motif_inputs[0],motif_inputs[1],latticeconstant=latticeconstant)
		name = 'Octa'
	elif len(motif_inputs) == 3:
		cluster = Decahedron(symbol,motif_inputs[0],motif_inputs[1],motif_inputs[2],latticeconstant=latticeconstant)
		name = 'Deca'
	else:
		print "Error"
		import pdb; pdb.set_trace()
		exit()
	remove_atoms.sort(reverse=True)
	if any(remove_atoms.count(x) > 1 for x in remove_atoms):
		print 'You have got two of the same atom entered in to this list, check this.'
		print 'remove_atoms = ' + str(remove_atoms)
		import pdb; pdb.set_trace()
		exit()
	for remove_atom_index in remove_atoms:
		del cluster[remove_atom_index]
	name += '_' + symbol + str(len(cluster))
	motif_details = ''
	for motif_input in motif_inputs:
		motif_details += str(motif_input) + '_'
	motif_details = motif_details[:-1]
	name += '_' + motif_details
	if not remove_atoms == []:
		name += '_atoms_removed_' + str(len(remove_atoms)) + '_cluster_'
		counter = 1
		while True:
			if not name + str(counter) + '.traj' in os.listdir('.'):
				name += str(counter)
				break
			counter += 1
	return cluster, name

def makeRossetteReconstruction(cluster,vertex_atom,atoms_around_vertex):
	cluster = copy.deepcopy(cluster)
	sum_x = 0; sum_y = 0; sum_z = 0;
	for atom_around_vertex in atoms_around_vertex:
		sum_x += cluster[atom_around_vertex].x
		sum_y += cluster[atom_around_vertex].y
		sum_z += cluster[atom_around_vertex].z
	def getAverage(sum_position):
		return sum_position/float(len(atoms_around_vertex))
	average_x = getAverage(sum_x)
	average_y = getAverage(sum_y)
	average_z = getAverage(sum_z)
	vector_centre = np.array([average_x, average_y, average_z])
	vector_diff = np.array(cluster[atoms_around_vertex[0]].position) - vector_centre
	vector_2 = np.array(cluster[atoms_around_vertex[1]].position) - vector_centre
	axis_vector = np.cross(np.ndarray.tolist(vector_diff),np.ndarray.tolist(vector_2))
	vector_centre.shape = (3,1)
	vector_diff.shape = (3,1)
	vector_2.shape = (3,1)
	axis_vector.shape = (3,1)
	def get_unit_vector(vector):
		norm = (vector.transpose().dot(vector))[0][0]**0.5
		unit_vector = divide(vector,norm)
		return unit_vector
	axis_unit_vector = get_unit_vector(axis_vector)
	Identity = array([[1,0,0],[0,1,0],[0,0,1]])
	Tx = axis_unit_vector[0][0]; Ty = axis_unit_vector[1][0]; Tz = axis_unit_vector[2][0]
	ucross = array([[0,-Tz,Ty],[Tz,0,-Tx],[-Ty,Tx,0]])
	ucrossu = matmul(axis_unit_vector,axis_unit_vector.transpose())
	def rotate(angle):
		angle = radians(angle)
		Rotation_Matrix = multiply(cos(angle),Identity) + multiply(sin(angle),ucross) + multiply((1.0 - cos(angle)),ucrossu)
		return add(vector_centre,matmul(Rotation_Matrix,vector_diff))
	angle = 0.0
	for atom_index in (atoms_around_vertex + [vertex_atom]):
		#import pdb; pdb.set_trace()
		atom_position = rotate(angle)
		atom_position.shape = (1,3)
		atom_position = np.ndarray.tolist(atom_position)[0]
		cluster[atom_index].position = atom_position
		angle += 60.0
	return cluster

def run_CNA_program(Cluster_1_name,Cluster_2_name):
	Cluster_1 = ASE_read(Cluster_1_name+'.traj')
	print 'Made ' + str(Cluster_1_name)
	Cluster_2 = ASE_read(Cluster_2_name+'.traj')
	print 'Made ' + str(Cluster_2_name)
	Structural_Recognition_Program(Cluster_1,Cluster_2,rCut_low,rCut_high,rCut_resolution,mode,name_1=Cluster_1_name,name_2=Cluster_2_name,recognise_multimetallic=False,print_plots=True)

def optimise(cluster, cluster_name):
	cluster_optimised = RunMinimisation(cluster).get_cluster()
	cluster_name_opt = cluster_name + '_Opt'
	ASE_write(cluster_name_opt+'.traj',cluster_optimised)
	return cluster_name_opt

def Cluster_CNA_profile(Cluster_name,CNA_signature_to_monitor):
	Cluster = ASE_read(Cluster_name+'.traj')
	print 'Made ' + str(Cluster_name)
	return Cluster_CNA_Profile(Cluster,rCut_low,rCut_high,rCut_resolution,CNA_signature_to_monitor,name=Cluster_name,recognise_multimetallic=False,print_plots=True)

#######################################################################
#######################################################################
#######################################################################
# Ico Structures Vs lots of other Ico structures
#######################################################################
#######################################################################
#######################################################################
Main_CNA_signatures_to_monitor = [[5,5,5],[4,2,1],[4,2,2]]
CNA_signatures_for_defect_icos = [[5,4,4],[4,3,3]]
CNA_signatures_for_111_surfaces = [[3,2,2],[3,1,1],[3,0,0]]
CNA_signatures_for_100_surfaces = [[2,1,1]]
CNA_signatures_for_110_surfaces = [[2,0,0]]
CNA_signatures_to_monitor = Main_CNA_signatures_to_monitor + CNA_signatures_for_defect_icos + CNA_signatures_for_111_surfaces + CNA_signatures_for_100_surfaces + CNA_signatures_for_110_surfaces
# 55 ico vs ico cluster - Rosette structure
motif_inputs = 3
cluster_main_ico, original_cluster_name_ico = makeMotif(motif_inputs,latticeconstant)
ASE_write(original_cluster_name_ico+'.traj',cluster_main_ico)
Cluster_CNA_profile(original_cluster_name_ico,CNA_signatures_to_monitor)
vertex_atom = 18
atoms_around_vertex = [42,36,35,46,17]
rosette_cluster = makeRossetteReconstruction(cluster_main_ico,vertex_atom,atoms_around_vertex)
rosette_cluster_name_ico = original_cluster_name_ico + '_rosette'
ASE_write(rosette_cluster_name_ico+'.traj',rosette_cluster)
run_CNA_program(original_cluster_name_ico,rosette_cluster_name_ico)
Cluster_CNA_profile(rosette_cluster_name_ico,CNA_signatures_to_monitor)
### Optimized
original_cluster_name_ico_opt = optimise(cluster_main_ico, original_cluster_name_ico)
Cluster_CNA_profile(original_cluster_name_ico_opt,CNA_signatures_to_monitor)
rosette_cluster_name_ico_opt = optimise(rosette_cluster, rosette_cluster_name_ico)
Cluster_CNA_profile(rosette_cluster_name_ico_opt,CNA_signatures_to_monitor)
run_CNA_program(original_cluster_name_ico_opt,rosette_cluster_name_ico_opt)
#######################################################################
# 54 ico vs ico cluster
motif_inputs = 3
remove_atoms = [27]
cluster_main_ico_1, cluster_name_ico_1 = makeMotif(motif_inputs,latticeconstant,remove_atoms=remove_atoms)
ASE_write(cluster_name_ico_1+'.traj',cluster_main_ico_1)
Cluster_CNA_profile(cluster_name_ico_1,CNA_signatures_to_monitor)
remove_atoms = [26]
cluster_main_ico_2, cluster_name_ico_2 = makeMotif(motif_inputs,latticeconstant,remove_atoms=remove_atoms)
ASE_write(cluster_name_ico_2+'.traj',cluster_main_ico_2)
Cluster_CNA_profile(cluster_name_ico_2,CNA_signatures_to_monitor)
remove_atoms = [0]
cluster_main_ico_3, cluster_name_ico_3 = makeMotif(motif_inputs,latticeconstant,remove_atoms=remove_atoms)
ASE_write(cluster_name_ico_3+'.traj',cluster_main_ico_3)
Cluster_CNA_profile(cluster_name_ico_3,CNA_signatures_to_monitor)
## run CNA programs
run_CNA_program(cluster_name_ico_1,cluster_name_ico_2)
run_CNA_program(cluster_name_ico_1,cluster_name_ico_3)
run_CNA_program(cluster_name_ico_2,cluster_name_ico_3)
### Optimized
cluster_name_ico_1_opt = optimise(cluster_main_ico_1, cluster_name_ico_1)
Cluster_CNA_profile(cluster_name_ico_1_opt,CNA_signatures_to_monitor)
cluster_name_ico_2_opt = optimise(cluster_main_ico_2, cluster_name_ico_2)
Cluster_CNA_profile(cluster_name_ico_2_opt,CNA_signatures_to_monitor)
cluster_name_ico_3_opt = optimise(cluster_main_ico_3, cluster_name_ico_3)
Cluster_CNA_profile(cluster_name_ico_3_opt,CNA_signatures_to_monitor)
run_CNA_program(cluster_name_ico_1_opt,cluster_name_ico_2_opt)
run_CNA_program(cluster_name_ico_1_opt,cluster_name_ico_3_opt)
run_CNA_program(cluster_name_ico_2_opt,cluster_name_ico_3_opt)
#######################################################################
# 147 ico cluster - Rosette structure
motif_inputs = 4
cluster_main, original_cluster_name_ico = makeMotif(motif_inputs,latticeconstant)
ASE_write(original_cluster_name_ico+'.traj',cluster_main)
Cluster_CNA_profile(original_cluster_name_ico,CNA_signatures_to_monitor)
vertex_atom = 74
atoms_around_vertex = [119,121,83,73,93]
rosette_cluster = makeRossetteReconstruction(cluster_main,vertex_atom,atoms_around_vertex)
rosette_cluster_name_ico = original_cluster_name_ico + '_rosette'
ASE_write(rosette_cluster_name_ico+'.traj',rosette_cluster)
Cluster_CNA_profile(rosette_cluster_name_ico,CNA_signatures_to_monitor)
run_CNA_program(original_cluster_name_ico,rosette_cluster_name_ico)
### Optimized
original_cluster_name_ico_opt = optimise(cluster_main, original_cluster_name_ico)
Cluster_CNA_profile(original_cluster_name_ico_opt,CNA_signatures_to_monitor)
rosette_cluster_name_ico_opt = optimise(rosette_cluster, rosette_cluster_name_ico)
Cluster_CNA_profile(rosette_cluster_name_ico_opt,CNA_signatures_to_monitor)
run_CNA_program(original_cluster_name_ico_opt,rosette_cluster_name_ico_opt)
#######################################################################
# 146 ico vs ico cluster
motif_inputs = 4
remove_atoms = [74]
cluster_main_ico_1, cluster_name_ico_1 = makeMotif(motif_inputs,latticeconstant,remove_atoms=remove_atoms)
ASE_write(cluster_name_ico_1+'.traj',cluster_main_ico_1)
Cluster_CNA_profile(cluster_name_ico_1,CNA_signatures_to_monitor)
remove_atoms = [121]
cluster_main_ico_2, cluster_name_ico_2 = makeMotif(motif_inputs,latticeconstant,remove_atoms=remove_atoms)
ASE_write(cluster_name_ico_2+'.traj',cluster_main_ico_2)
Cluster_CNA_profile(cluster_name_ico_2,CNA_signatures_to_monitor)
remove_atoms = [122]
cluster_main_ico_3, cluster_name_ico_3 = makeMotif(motif_inputs,latticeconstant,remove_atoms=remove_atoms)
ASE_write(cluster_name_ico_3+'.traj',cluster_main_ico_3)
Cluster_CNA_profile(cluster_name_ico_3,CNA_signatures_to_monitor)
remove_atoms = [0]
cluster_main_ico_4, cluster_name_ico_4 = makeMotif(motif_inputs,latticeconstant,remove_atoms=remove_atoms)
ASE_write(cluster_name_ico_4+'.traj',cluster_main_ico_4)
Cluster_CNA_profile(cluster_name_ico_4,CNA_signatures_to_monitor)
## run CNA programs
run_CNA_program(cluster_name_ico_1,cluster_name_ico_2)
run_CNA_program(cluster_name_ico_1,cluster_name_ico_3)
run_CNA_program(cluster_name_ico_1,cluster_name_ico_4)
run_CNA_program(cluster_name_ico_2,cluster_name_ico_3)
run_CNA_program(cluster_name_ico_2,cluster_name_ico_4)
run_CNA_program(cluster_name_ico_3,cluster_name_ico_4)
### Optimized
cluster_name_ico_1_opt = optimise(cluster_main_ico_1, cluster_name_ico_1)
Cluster_CNA_profile(cluster_name_ico_1_opt,CNA_signatures_to_monitor)
cluster_name_ico_2_opt = optimise(cluster_main_ico_2, cluster_name_ico_2)
Cluster_CNA_profile(cluster_name_ico_2_opt,CNA_signatures_to_monitor)
cluster_name_ico_3_opt = optimise(cluster_main_ico_3, cluster_name_ico_3)
Cluster_CNA_profile(cluster_name_ico_3_opt,CNA_signatures_to_monitor)
cluster_name_ico_4_opt = optimise(cluster_main_ico_4, cluster_name_ico_4)
Cluster_CNA_profile(cluster_name_ico_4_opt,CNA_signatures_to_monitor)
run_CNA_program(cluster_name_ico_1_opt,cluster_name_ico_2_opt)
run_CNA_program(cluster_name_ico_1_opt,cluster_name_ico_3_opt)
run_CNA_program(cluster_name_ico_1_opt,cluster_name_ico_4_opt)
run_CNA_program(cluster_name_ico_2_opt,cluster_name_ico_3_opt)
run_CNA_program(cluster_name_ico_2_opt,cluster_name_ico_4_opt)
run_CNA_program(cluster_name_ico_3_opt,cluster_name_ico_4_opt)
#######################################################################





