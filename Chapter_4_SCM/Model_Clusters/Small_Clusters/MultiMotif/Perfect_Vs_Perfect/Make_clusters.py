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
from Structural_Recognition_Program import Structural_Recognition_Program, Cluster_CNA_Profile
from RunMinimisation import RunMinimisation
print 'Running program'

Main_CNA_signatures_to_monitor = [[5,5,5],[4,2,1],[4,2,2]]
CNA_signatures_for_defect_icos = [[5,4,4],[4,3,3]]
CNA_signatures_for_111_surfaces = [[3,2,2],[3,1,1],[3,0,0]]
CNA_signatures_for_100_surfaces = [[2,1,1]]
CNA_signatures_for_110_surfaces = [[2,0,0]]
CNA_signatures_to_monitor = Main_CNA_signatures_to_monitor + CNA_signatures_for_defect_icos + CNA_signatures_for_111_surfaces + CNA_signatures_for_100_surfaces + CNA_signatures_for_110_surfaces
def Cluster_CNA_profile(Cluster_name,CNA_signature_to_monitor):
	Cluster = ASE_read(Cluster_name+'.traj')
	print 'Made ' + str(Cluster_name)
	return Cluster_CNA_Profile(Cluster,rCut_low,rCut_high,rCut_resolution,CNA_signature_to_monitor,name=Cluster_name,recognise_multimetallic=False,print_plots=True)

def delete():
	for item in os.listdir('.'):
		if not item == os.path.basename(__file__) and not os.path.isdir(item) and not item.startswith('submit.sl'):
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
	ASE_write(name+'.traj',cluster,'traj')
	return name
	
def run_CNA_program(Cluster_1_name,Cluster_2_name):
	Cluster_1 = ASE_read(Cluster_1_name+'.traj')
	print 'Made ' + str(Cluster_1_name)
	Cluster_2 = ASE_read(Cluster_2_name+'.traj')
	print 'Made ' + str(Cluster_2_name)
	Structural_Recognition_Program(Cluster_1,Cluster_2,rCut_low,rCut_high,rCut_resolution,mode,name_1=Cluster_1_name,name_2=Cluster_2_name,recognise_multimetallic=False,print_plots=True)

def optimise(cluster_name):
	cluster = ASE_read(cluster_name+'.traj')
	cluster_optimised = RunMinimisation(cluster).get_cluster()
	cluster_name_opt = cluster_name + '_Opt'
	ASE_write(cluster_name_opt+'.traj',cluster_optimised)
	return cluster_name_opt

#######################################################################
#######################################################################
#######################################################################
# Comparing perfect structures
#######################################################################
#######################################################################
#######################################################################
# The 55 range, ico, deca and octa

Cluster_1_name = makeMotif(3,latticeconstant)
Cluster_2_name = makeMotif([5,2],latticeconstant)
run_CNA_program(Cluster_1_name,Cluster_2_name)
Cluster_1_name_opt = optimise(Cluster_1_name)
Cluster_2_name_opt = optimise(Cluster_2_name)
run_CNA_program(Cluster_1_name_opt,Cluster_2_name_opt)
Cluster_CNA_profile(Cluster_1_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_1_name_opt,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name_opt,CNA_signatures_to_monitor)
#
Cluster_1_name = makeMotif(3,latticeconstant)
Cluster_2_name = makeMotif([3,3,0],latticeconstant)
run_CNA_program(Cluster_1_name,Cluster_2_name)
Cluster_1_name_opt = optimise(Cluster_1_name)
Cluster_2_name_opt = optimise(Cluster_2_name)
run_CNA_program(Cluster_1_name_opt,Cluster_2_name_opt)
Cluster_CNA_profile(Cluster_2_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name_opt,CNA_signatures_to_monitor)
#
Cluster_1_name = makeMotif([5,2],latticeconstant)
Cluster_2_name = makeMotif([3,3,0],latticeconstant)
run_CNA_program(Cluster_1_name,Cluster_2_name)
Cluster_1_name_opt = optimise(Cluster_1_name)
Cluster_2_name_opt = optimise(Cluster_2_name)
run_CNA_program(Cluster_1_name_opt,Cluster_2_name_opt)
#######################################################################
# 85, deca and octa
Cluster_1_name = makeMotif([5,0],latticeconstant)
Cluster_2_name = makeMotif([4,2,0],latticeconstant)
run_CNA_program(Cluster_1_name,Cluster_2_name)
Cluster_1_name_opt = optimise(Cluster_1_name)
Cluster_2_name_opt = optimise(Cluster_2_name)
run_CNA_program(Cluster_1_name_opt,Cluster_2_name_opt)
Cluster_CNA_profile(Cluster_1_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_1_name_opt,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name_opt,CNA_signatures_to_monitor)
#######################################################################
# 116, deca and octa
Cluster_1_name = makeMotif([6,2],latticeconstant)
Cluster_2_name = makeMotif([4,3,0],latticeconstant)
run_CNA_program(Cluster_1_name,Cluster_2_name)
Cluster_1_name_opt = optimise(Cluster_1_name)
Cluster_2_name_opt = optimise(Cluster_2_name)
run_CNA_program(Cluster_1_name_opt,Cluster_2_name_opt)
Cluster_CNA_profile(Cluster_1_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_1_name_opt,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name_opt,CNA_signatures_to_monitor)
#######################################################################
# 146, deca and octa
Cluster_1_name = makeMotif([6,0],latticeconstant)
Cluster_2_name = makeMotif([3,2,1],latticeconstant)
run_CNA_program(Cluster_1_name,Cluster_2_name)
Cluster_1_name_opt = optimise(Cluster_1_name)
Cluster_2_name_opt = optimise(Cluster_2_name)
run_CNA_program(Cluster_1_name_opt,Cluster_2_name_opt)
Cluster_CNA_profile(Cluster_1_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_1_name_opt,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name_opt,CNA_signatures_to_monitor)
#######################################################################
# The 147 range, ico, deca and octa
Cluster_1_name = makeMotif(4,latticeconstant)
Cluster_2_name = makeMotif([7,3],latticeconstant)
run_CNA_program(Cluster_1_name,Cluster_2_name)
Cluster_1_name_opt = optimise(Cluster_1_name)
Cluster_2_name_opt = optimise(Cluster_2_name)
run_CNA_program(Cluster_1_name_opt,Cluster_2_name_opt)
Cluster_CNA_profile(Cluster_1_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_1_name_opt,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name_opt,CNA_signatures_to_monitor)
#
Cluster_1_name = makeMotif(4,latticeconstant)
Cluster_2_name = makeMotif([4,4,0],latticeconstant)
run_CNA_program(Cluster_1_name,Cluster_2_name)
Cluster_1_name_opt = optimise(Cluster_1_name)
Cluster_2_name_opt = optimise(Cluster_2_name)
run_CNA_program(Cluster_1_name_opt,Cluster_2_name_opt)
Cluster_CNA_profile(Cluster_2_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name_opt,CNA_signatures_to_monitor)
#
Cluster_1_name = makeMotif([7,3],latticeconstant)
Cluster_2_name = makeMotif([4,4,0],latticeconstant)
run_CNA_program(Cluster_1_name,Cluster_2_name)
Cluster_1_name_opt = optimise(Cluster_1_name)
Cluster_2_name_opt = optimise(Cluster_2_name)
run_CNA_program(Cluster_1_name_opt,Cluster_2_name_opt)
#######################################################################
# 156, two decas

Cluster_1_name = makeMotif([2,1,2],latticeconstant)
Cluster_2_name = makeMotif([5,2,0],latticeconstant)
run_CNA_program(Cluster_1_name,Cluster_2_name)
Cluster_1_name_opt = optimise(Cluster_1_name)
Cluster_2_name_opt = optimise(Cluster_2_name)
run_CNA_program(Cluster_1_name_opt,Cluster_2_name_opt)
Cluster_CNA_profile(Cluster_1_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_1_name_opt,CNA_signatures_to_monitor)
Cluster_CNA_profile(Cluster_2_name_opt,CNA_signatures_to_monitor)
#######################################################################

#######################################################################
#######################################################################
#######################################################################
# Comparing inperfect structures to perfect structures.
#######################################################################
#######################################################################
#######################################################################
