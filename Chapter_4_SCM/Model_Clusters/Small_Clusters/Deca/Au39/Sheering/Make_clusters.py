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
from Structural_Recognition_Program import Structural_Recognition_Program, Print_Multiple_Plots
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
	cluster_optimised.set_calculator(None)
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

def Shearing(input_cluster,Shearing_Matrix,atoms_to_shear=[]):
    if atoms_to_shear == []:
        atoms_to_shear = range(len(input_cluster))
    output_cluster = copy.deepcopy(input_cluster)
    for index in atoms_to_shear:
        output_cluster[index].x = Shearing_Matrix[0][0]*output_cluster[index].x + Shearing_Matrix[0][1]*output_cluster[index].y + Shearing_Matrix[0][2]*output_cluster[index].z
        output_cluster[index].y = Shearing_Matrix[1][0]*output_cluster[index].x + Shearing_Matrix[1][1]*output_cluster[index].y + Shearing_Matrix[1][2]*output_cluster[index].z
        output_cluster[index].z = Shearing_Matrix[2][0]*output_cluster[index].x + Shearing_Matrix[2][1]*output_cluster[index].y + Shearing_Matrix[2][2]*output_cluster[index].z
    return output_cluster

rCut_low = 2.8840
rCut_high = round(rCut_low * (2.0)**0.5,4)
rCut_resolution = 0.0005
mode = 'total'

symbol = 'Au'
latticeconstant = rCut_high
cluster_main = Decahedron(symbol,3,2,0,latticeconstant=latticeconstant)
cluster_main = move_cluster_to_COM(cluster_main)
original_cluster_name = 'Au39_original'
ASE_write(original_cluster_name+'.traj',cluster_main)
original_cluster_name_Opt, cluster_main_Opt = optimise(original_cluster_name)

def run_CNA_program(Cluster_1_name,Cluster_2_name):
    Cluster_1 = ASE_read(Cluster_1_name+'.traj')
    print 'Made ' + str(Cluster_1_name)
    Cluster_2 = ASE_read(Cluster_2_name+'.traj')
    print 'Made ' + str(Cluster_2_name)
    CNA_run = Structural_Recognition_Program(Cluster_1,Cluster_2,rCut_low,rCut_high,rCut_resolution,mode,name_1=Cluster_1_name,name_2=Cluster_2_name,recognise_multimetallic=False,print_plots=True)
    return CNA_run

######################################################
#import pdb; pdb.set_trace()
######################################################
'''
cluster_name = 'Au38_shearing_test'
Shearing_Matrix = [[1,0,0],[0,1,0],[0,0,1]]
Sheared_Cluster = Shearing(cluster_main,Shearing_Matrix)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
run_CNA_program(original_cluster_name,cluster_name)
'''
############################################################################################################
CNA_run_shearing_z_in_x_direction = []
############################################################################################################
############################################################################################################
atoms_to_shear = [32,34,36]
displacements = np.arange(0.2, 1.01, 0.2)
for displacement in displacements:
	Shearing_Matrix = [[1,0,displacement],[0,1,0],[0,0,1]]
	cluster_name = 'Au39_Deca_shearing_100-111_face_by_'+str(displacement)
	Sheared_Cluster = Shearing(cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
	ASE_write(cluster_name+'.traj',Sheared_Cluster)
	CNA_run = run_CNA_program(original_cluster_name,cluster_name)
	rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
	CNA_run_shearing_z_in_x_direction.append({'Name':str(displacement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
	#cluster_name_opt = optimise(cluster_name)
	#run_CNA_program(original_cluster_name_opt,cluster_name_opt)
	##### Optimised Structure #####
	cluster_name_Opt = 'Au39_Deca_shearing_100-111_face_by_'+str(displacement)+'_Opt'
	Sheared_Cluster_Opt = Shearing(cluster_main_Opt,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
	ASE_write(cluster_name_Opt+'.traj',Sheared_Cluster_Opt)
	run_CNA_program(original_cluster_name_Opt,cluster_name_Opt)

############################################################################################################
CNA_run_shearing_z_in_x_direction = []
############################################################################################################
############################################################################################################
atoms_to_shear = [15,18]
displacements = np.arange(0.1, 0.51, 0.1)
for displacement in displacements:
	Shearing_Matrix = [[1,0,displacement],[0,1,0],[0,0,1]]
	cluster_name = 'Au39_Deca_shearing_111_face_by_'+str(displacement)
	Sheared_Cluster = Shearing(cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
	ASE_write(cluster_name+'.traj',Sheared_Cluster)
	CNA_run = run_CNA_program(original_cluster_name,cluster_name)
	rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
	CNA_run_shearing_z_in_x_direction.append({'Name':str(displacement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
	#cluster_name_opt = optimise(cluster_name)
	#run_CNA_program(original_cluster_name_opt,cluster_name_opt)
	##### Optimised Structure #####
	cluster_name_Opt = 'Au39_Deca_shearing_111_face_by_'+str(displacement)+'_Opt'
	Sheared_Cluster_Opt = Shearing(cluster_main_Opt,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
	ASE_write(cluster_name_Opt+'.traj',Sheared_Cluster_Opt)
	run_CNA_program(original_cluster_name_Opt,cluster_name_Opt)
	

#import pdb; pdb.set_trace()
















'''
######################################################
displacement = 0.1
cluster_name = 'Au38_shearing_z_'+str(displacement)+'_in_x'
Shearing_Matrix = [[1,0,displacement],[0,1,0],[0,0,1]]
Sheared_Cluster = Shearing(cluster_main,Shearing_Matrix)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_shearing_z_in_x_direction.append({'Name':str(displacement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
displacement = 0.15
cluster_name = 'Au38_shearing_z_'+str(displacement)+'_in_x'
Shearing_Matrix = [[1,0,displacement],[0,1,0],[0,0,1]]
Sheared_Cluster = Shearing(cluster_main,Shearing_Matrix)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_shearing_z_in_x_direction.append({'Name':str(displacement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
displacement = 0.2
cluster_name = 'Au38_shearing_z_'+str(displacement)+'_in_x'
Shearing_Matrix = [[1,0,displacement],[0,1,0],[0,0,1]]
Sheared_Cluster = Shearing(cluster_main,Shearing_Matrix)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_shearing_z_in_x_direction.append({'Name':str(displacement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
displacement = 0.25
cluster_name = 'Au38_shearing_z_'+str(displacement)+'_in_x'
Shearing_Matrix = [[1,0,displacement],[0,1,0],[0,0,1]]
Sheared_Cluster = Shearing(cluster_main,Shearing_Matrix)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_shearing_z_in_x_direction.append({'Name':str(displacement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
displacement = 0.3
cluster_name = 'Au38_shearing_z_'+str(displacement)+'_in_x'
Shearing_Matrix = [[1,0,displacement],[0,1,0],[0,0,1]]
Sheared_Cluster = Shearing(cluster_main,Shearing_Matrix)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_shearing_z_in_x_direction.append({'Name':str(displacement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
displacement = 0.4
cluster_name = 'Au38_shearing_z_'+str(displacement)+'_in_x'
Shearing_Matrix = [[1,0,displacement],[0,1,0],[0,0,1]]
Sheared_Cluster = Shearing(cluster_main,Shearing_Matrix)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_shearing_z_in_x_direction.append({'Name':str(displacement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
displacement = 0.5
cluster_name = 'Au38_shearing_z_'+str(displacement)+'_in_x'
Shearing_Matrix = [[1,0,displacement],[0,1,0],[0,0,1]]
Sheared_Cluster = Shearing(cluster_main,Shearing_Matrix)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_shearing_z_in_x_direction.append({'Name':str(displacement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
############################################################################################################
Print_Multiple_Plots(CNA_run_shearing_z_in_x_direction,'Au38_shearing_z_disp_in_x',mode)
############################################################################################################
############################################################################################################
############################################################################################################
CNA_run_Face_sheared_in_y_direction = []
############################################################################################################
movement = 0.25
cluster_name = 'Au38_111_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,movement,0],[0,1,0],[0,-movement,1]]
atoms_to_shear = [24,35,23,37,36,32,33]
centre_point = cluster_main[37].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_in_y_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
movement = 0.225
cluster_name = 'Au38_111_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,movement,0],[0,1,0],[0,-movement,1]]
atoms_to_shear = [24,35,23,37,36,32,33]
centre_point = cluster_main[37].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_in_y_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
movement = 0.2
cluster_name = 'Au38_111_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,movement,0],[0,1,0],[0,-movement,1]]
atoms_to_shear = [24,35,23,37,36,32,33]
centre_point = cluster_main[37].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_in_y_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
movement = 0.175
cluster_name = 'Au38_111_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,movement,0],[0,1,0],[0,-movement,1]]
atoms_to_shear = [24,35,23,37,36,32,33]
centre_point = cluster_main[37].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_in_y_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
movement = 0.15
cluster_name = 'Au38_111_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,movement,0],[0,1,0],[0,-movement,1]]
atoms_to_shear = [24,35,23,37,36,32,33]
centre_point = cluster_main[37].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_in_y_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
movement = 0.125
cluster_name = 'Au38_111_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,movement,0],[0,1,0],[0,-movement,1]]
atoms_to_shear = [24,35,23,37,36,32,33]
centre_point = cluster_main[37].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_in_y_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
movement = 0.10
cluster_name = 'Au38_111_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,movement,0],[0,1,0],[0,-movement,1]]
atoms_to_shear = [24,35,23,37,36,32,33]
centre_point = cluster_main[37].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_in_y_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
movement = 0.075
cluster_name = 'Au38_111_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,movement,0],[0,1,0],[0,-movement,1]]
atoms_to_shear = [24,35,23,37,36,32,33]
centre_point = cluster_main[37].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_in_y_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
movement = 0.05
cluster_name = 'Au38_111_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,movement,0],[0,1,0],[0,-movement,1]]
atoms_to_shear = [24,35,23,37,36,32,33]
centre_point = cluster_main[37].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_in_y_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
movement = 0.025
cluster_name = 'Au38_111_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,movement,0],[0,1,0],[0,-movement,1]]
atoms_to_shear = [24,35,23,37,36,32,33]
centre_point = cluster_main[37].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_in_y_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
Print_Multiple_Plots(CNA_run_Face_sheared_in_y_direction,'Au38_Face_sheared_in_y_direction',mode)
############################################################################################################
############################################################################################################
############################################################################################################
CNA_run_Face_sheared_100_in_xy_direction = []
############################################################################################################
movement = 0.1
cluster_name = 'Au38_100_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,0,movement],[0,1,movement],[0,0,1]]
atoms_to_shear = [14,23]
centre_point = cluster_main[16].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_100_in_xy_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
movement = 0.2
cluster_name = 'Au38_100_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,0,movement],[0,1,movement],[0,0,1]]
atoms_to_shear = [14,23]
centre_point = cluster_main[16].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_100_in_xy_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
movement = 0.3
cluster_name = 'Au38_100_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,0,movement],[0,1,movement],[0,0,1]]
atoms_to_shear = [14,23]
centre_point = cluster_main[16].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_100_in_xy_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
movement = 0.4
cluster_name = 'Au38_100_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,0,movement],[0,1,movement],[0,0,1]]
atoms_to_shear = [14,23]
centre_point = cluster_main[16].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
view(Sheared_Cluster)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_100_in_xy_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################
movement = 0.5
cluster_name = 'Au38_100_Face_Sheered_by_'+str(movement)+'y'
Shearing_Matrix = [[1,0,movement],[0,1,movement],[0,0,1]]
atoms_to_shear = [14,23]
centre_point = cluster_main[16].position
print centre_point
shifted_cluster_main = move_cluster(cluster_main,centre_point)
Sheared_Cluster = Shearing(shifted_cluster_main,Shearing_Matrix,atoms_to_shear=atoms_to_shear)
view(Sheared_Cluster)
ASE_write(cluster_name+'.traj',Sheared_Cluster)
CNA_run = run_CNA_program(original_cluster_name,cluster_name)
rCuts, percentage_similarities = CNA_run.get_rCut_and_percentage_similarities()
CNA_run_Face_sheared_100_in_xy_direction.append({'Name':str(movement),'rCuts':rCuts,'percentage_similarities':percentage_similarities})
######################################################

'''

