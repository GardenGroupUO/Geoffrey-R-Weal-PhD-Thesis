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

#######################################################################
#######################################################################
#######################################################################
# Ico Structures Vs lots of other structures
#######################################################################
#######################################################################
#######################################################################
# 55 ico vs ico cluster - Rosette structure
motif_inputs = 3
cluster_main_ico, original_cluster_name_ico = makeMotif(motif_inputs,latticeconstant)
ASE_write(original_cluster_name_ico+'.traj',cluster_main_ico)
vertex_atom = 18
atoms_around_vertex = [42,36,35,46,17]
rosette_cluster = makeRossetteReconstruction(cluster_main_ico,vertex_atom,atoms_around_vertex)
rosette_cluster_name_ico = original_cluster_name_ico + '_rosette'
ASE_write(rosette_cluster_name_ico+'.traj',rosette_cluster)
run_CNA_program(original_cluster_name_ico,rosette_cluster_name_ico)
### Optimized
original_cluster_name_ico_opt = optimise(cluster_main_ico, original_cluster_name_ico)
rosette_cluster_name_ico_opt = optimise(rosette_cluster, rosette_cluster_name_ico)
run_CNA_program(original_cluster_name_ico_opt,rosette_cluster_name_ico_opt)
#######################################################################
# 55 deca vs ico cluster - Rosette structure
cluster_main_deca, original_cluster_name_deca = makeMotif([3,3,0],latticeconstant)
ASE_write(original_cluster_name_deca+'.traj',cluster_main_deca)
run_CNA_program(original_cluster_name_deca, rosette_cluster_name_ico)
### Optimized
original_cluster_name_deca_opt = optimise(cluster_main_deca, original_cluster_name_deca)
run_CNA_program(original_cluster_name_deca_opt,rosette_cluster_name_ico_opt)
#######################################################################
# 55 octa vs ico cluster - Rosette structure
cluster_main_octa, original_cluster_name_octa = makeMotif([5,2],latticeconstant)
ASE_write(original_cluster_name_octa+'.traj',cluster_main_octa)
run_CNA_program(original_cluster_name_octa, rosette_cluster_name_ico)
### Optimized
original_cluster_name_octa_opt = optimise(cluster_main_octa, original_cluster_name_octa)
run_CNA_program(original_cluster_name_octa_opt,rosette_cluster_name_ico_opt)
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
# 54 cluster - Deca Vs Ico structure minus
cluster_main_deca, original_cluster_name_deca = makeMotif([4,1,0],latticeconstant)
ASE_write(original_cluster_name_deca+'.traj',cluster_main_deca)
original_cluster_name_deca_opt = optimise(cluster_main_deca, original_cluster_name_deca)
motif_inputs = 3
remove_atoms = [18]
cluster_main_ico, original_cluster_name_ico_vertex = makeMotif(motif_inputs,latticeconstant,remove_atoms)
ASE_write(original_cluster_name_ico_vertex+'.traj',cluster_main_ico)
original_cluster_name_ico_vertex_opt = optimise(cluster_main_ico, original_cluster_name_ico_vertex)
vertex_atom = 13
atoms_around_vertex = [39,34,33,43,14]
for remove_atom in remove_atoms:
	for index in range(len(atoms_around_vertex)):
		if remove_atom < atoms_around_vertex[index]:
			atoms_around_vertex[index] -= 1
print atoms_around_vertex
rosette_cluster = makeRossetteReconstruction(cluster_main_ico,vertex_atom,atoms_around_vertex)
rosette_cluster_name_ico = original_cluster_name_ico_vertex + '_rosette'
ASE_write(rosette_cluster_name_ico+'.traj',rosette_cluster)
rosette_cluster_name_ico_opt = optimise(rosette_cluster, rosette_cluster_name_ico)
run_CNA_program(original_cluster_name_deca,original_cluster_name_ico_vertex)
run_CNA_program(original_cluster_name_deca, rosette_cluster_name_ico)
### Optimized
run_CNA_program(original_cluster_name_deca_opt,original_cluster_name_ico_vertex_opt)
run_CNA_program(original_cluster_name_deca_opt,rosette_cluster_name_ico_opt)
###
motif_inputs = 3
remove_atoms = [0]
cluster_main_ico_no_centre, original_cluster_name_ico_no_centre = makeMotif(motif_inputs,latticeconstant,remove_atoms)
ASE_write(original_cluster_name_ico_no_centre+'.traj',cluster_main_ico_no_centre)
original_cluster_name_ico_no_centre_opt = optimise(cluster_main_ico_no_centre, original_cluster_name_ico_no_centre)
run_CNA_program(original_cluster_name_deca,original_cluster_name_ico_no_centre)
run_CNA_program(rosette_cluster_name_ico,original_cluster_name_ico_no_centre)
### Optimized
run_CNA_program(original_cluster_name_deca_opt,original_cluster_name_ico_no_centre_opt)
run_CNA_program(rosette_cluster_name_ico_opt,original_cluster_name_ico_no_centre_opt)
# Other compare ico Structures
motif_inputs = 3
cluster_main_ico_edge, original_cluster_name_ico_edge   = makeMotif(motif_inputs,latticeconstant,remove_atoms=[51])
ASE_write(original_cluster_name_ico_edge+'.traj',cluster_main_ico_edge)
original_cluster_name_ico_edge_opt = optimise(cluster_main_ico_edge, original_cluster_name_ico_edge)
run_CNA_program(rosette_cluster_name_ico           ,original_cluster_name_ico_vertex)
run_CNA_program(original_cluster_name_ico_no_centre,original_cluster_name_ico_vertex)
run_CNA_program(rosette_cluster_name_ico,original_cluster_name_ico_no_centre)
### Optimized
run_CNA_program(rosette_cluster_name_ico_opt           ,original_cluster_name_ico_vertex_opt)
run_CNA_program(original_cluster_name_ico_no_centre_opt,original_cluster_name_ico_vertex_opt)
run_CNA_program(rosette_cluster_name_ico_opt           ,original_cluster_name_ico_no_centre_opt)
#######################################################################
#######################################################################
#######################################################################
# Ico 147
#######################################################################
#######################################################################
#######################################################################
# 147 ico cluster - Rosette structure
motif_inputs = 4
cluster_main, original_cluster_name_ico = makeMotif(motif_inputs,latticeconstant)
ASE_write(original_cluster_name_ico+'.traj',cluster_main)
original_cluster_name_ico_opt = optimise(cluster_main, original_cluster_name_ico)
vertex_atom = 74
atoms_around_vertex = [119,121,83,73,93]
rosette_cluster = makeRossetteReconstruction(cluster_main,vertex_atom,atoms_around_vertex)
rosette_cluster_name_ico = original_cluster_name_ico + '_rosette'
ASE_write(rosette_cluster_name_ico+'.traj',rosette_cluster)
rosette_cluster_name_ico_opt = optimise(rosette_cluster, rosette_cluster_name_ico)
run_CNA_program(original_cluster_name_ico,rosette_cluster_name_ico)
### Optimized
run_CNA_program(original_cluster_name_ico_opt,rosette_cluster_name_ico_opt)
#######################################################################
# 147 deca vs ico cluster - Rosette structure
cluster_main_deca, original_cluster_name_deca = makeMotif([4,4,0],latticeconstant)
ASE_write(original_cluster_name_deca+'.traj',cluster_main_deca)
original_cluster_name_deca_opt = optimise(cluster_main_deca, original_cluster_name_deca)
run_CNA_program(original_cluster_name_deca,rosette_cluster_name_ico)
### Optimized
run_CNA_program(original_cluster_name_deca_opt,rosette_cluster_name_ico_opt)
#######################################################################
# 147 octa vs ico cluster - Rosette structure
cluster_main_octa, original_cluster_name_octa = makeMotif([7,3],latticeconstant)
ASE_write(original_cluster_name_octa+'.traj',cluster_main_octa)
original_cluster_name_octa_opt = optimise(cluster_main_octa, original_cluster_name_octa)
run_CNA_program(original_cluster_name_octa,rosette_cluster_name_ico)
### Optimized
run_CNA_program(original_cluster_name_octa_opt,rosette_cluster_name_ico_opt)
#######################################################################
#######################################################################
#######################################################################
# 146 cluster - Deca Vs Ico structure minus
cluster_main_deca, original_cluster_name_deca = makeMotif([3,2,1],latticeconstant)
ASE_write(original_cluster_name_deca+'.traj',cluster_main_deca)
original_cluster_name_deca_opt = optimise(cluster_main_deca, original_cluster_name_deca)
motif_inputs = 4
remove_atoms = [75]
cluster_main_ico, original_cluster_name_ico = makeMotif(motif_inputs,latticeconstant,remove_atoms)
ASE_write(original_cluster_name_ico+'.traj',cluster_main_ico)
original_cluster_name_ico_opt = optimise(cluster_main_ico, original_cluster_name_ico)
vertex_atom = 74
atoms_around_vertex = [119,121,83,73,93]
for remove_atom in remove_atoms:
	for index in range(len(atoms_around_vertex)):
		if remove_atom < atoms_around_vertex[index]:
			atoms_around_vertex[index] -= 1
rosette_cluster = makeRossetteReconstruction(cluster_main_ico,vertex_atom,atoms_around_vertex)
rosette_cluster_name_ico = original_cluster_name_ico + '_rosette'
ASE_write(rosette_cluster_name_ico+'.traj',rosette_cluster)
rosette_cluster_name_ico_opt = optimise(rosette_cluster, rosette_cluster_name_ico)
run_CNA_program(original_cluster_name_ico , rosette_cluster_name_ico)
run_CNA_program(original_cluster_name_deca,original_cluster_name_ico)
run_CNA_program(original_cluster_name_deca, rosette_cluster_name_ico)
### Optimized
run_CNA_program(original_cluster_name_ico_opt , rosette_cluster_name_ico_opt)
run_CNA_program(original_cluster_name_deca_opt,original_cluster_name_ico_opt)
run_CNA_program(original_cluster_name_deca_opt, rosette_cluster_name_ico_opt)
###
motif_inputs = 4
remove_atoms = [0]
cluster_main_ico_no_centre, original_cluster_name_ico_no_centre = makeMotif(motif_inputs,latticeconstant,remove_atoms)
ASE_write(original_cluster_name_ico_no_centre+'.traj',cluster_main_ico_no_centre)
original_cluster_name_ico_no_centre_opt = optimise(cluster_main_ico_no_centre, original_cluster_name_ico_no_centre)
run_CNA_program(original_cluster_name_deca, original_cluster_name_ico_no_centre)
run_CNA_program(rosette_cluster_name_ico  , original_cluster_name_ico_no_centre)
### Optimized
run_CNA_program(original_cluster_name_deca_opt, original_cluster_name_ico_no_centre_opt)
run_CNA_program(rosette_cluster_name_ico_opt  , original_cluster_name_ico_no_centre_opt)
#######################################################################
# 146 cluster - Octa Vs Ico structure minus
cluster_main_octa, original_cluster_name_octa_vertex = makeMotif([6,0],latticeconstant)
ASE_write(original_cluster_name_octa_vertex+'.traj',cluster_main_octa)
original_cluster_name_octa_vertex_opt = optimise(cluster_main_octa, original_cluster_name_octa_vertex)
run_CNA_program(original_cluster_name_octa_vertex,original_cluster_name_ico)
run_CNA_program(original_cluster_name_octa_vertex, rosette_cluster_name_ico)
run_CNA_program(original_cluster_name_octa_vertex, original_cluster_name_ico_no_centre)
### Optimized
run_CNA_program(original_cluster_name_octa_vertex_opt,original_cluster_name_ico_opt)
run_CNA_program(original_cluster_name_octa_vertex_opt, rosette_cluster_name_ico_opt)
run_CNA_program(original_cluster_name_octa_vertex_opt, original_cluster_name_ico_no_centre_opt)
# Other compare ico Structures
motif_inputs = 4
cluster_main_ico_edge, original_cluster_name_ico_edge   = makeMotif(motif_inputs,latticeconstant,remove_atoms=[119])
ASE_write(original_cluster_name_ico_edge+'.traj',cluster_main_ico_edge)
original_cluster_name_ico_edge_opt = optimise(cluster_main_ico_edge, original_cluster_name_ico_edge)
run_CNA_program(original_cluster_name_octa_vertex  ,original_cluster_name_ico_edge)
run_CNA_program(original_cluster_name_ico_no_centre,original_cluster_name_ico_edge)
### Optimized
run_CNA_program(original_cluster_name_octa_vertex_opt  ,original_cluster_name_ico_edge_opt)
run_CNA_program(original_cluster_name_ico_no_centre_opt,original_cluster_name_ico_edge_opt)
###
cluster_main_ico_face, original_cluster_name_ico_face  = makeMotif(motif_inputs,latticeconstant,remove_atoms=[122])
ASE_write(original_cluster_name_ico_face+'.traj',cluster_main_ico_face)
original_cluster_name_ico_face_opt = optimise(cluster_main_ico_face, original_cluster_name_ico_face)
run_CNA_program(original_cluster_name_octa_vertex  ,original_cluster_name_ico_face)
run_CNA_program(original_cluster_name_ico_no_centre,original_cluster_name_ico_face)
run_CNA_program(original_cluster_name_ico_edge     ,original_cluster_name_ico_face)
### Optimized
run_CNA_program(original_cluster_name_octa_vertex_opt  ,original_cluster_name_ico_face_opt)
run_CNA_program(original_cluster_name_ico_no_centre_opt,original_cluster_name_ico_face_opt)
run_CNA_program(original_cluster_name_ico_edge_opt     ,original_cluster_name_ico_face_opt)
###
#import pdb; pdb.set_trace()
