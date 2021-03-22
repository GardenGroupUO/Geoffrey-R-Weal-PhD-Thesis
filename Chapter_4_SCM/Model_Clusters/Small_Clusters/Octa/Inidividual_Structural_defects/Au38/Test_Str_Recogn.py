import sys, os
'''
add_to_import_list = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir, os.path.pardir, os.path.pardir, os.path.pardir, 'Main_Program')
print add_to_import_list
sys.path.append(add_to_import_list)
'''
from Structural_Recognition_Program import Structural_Recognition_Program
from ase.cluster.decahedron import Decahedron
from asap3.Internal.BuiltinPotentials import Gupta
from ase.optimize import FIRE
from ase.io import read as ASE_read
from ase.visualize import view
print 'Loading matplotlib'
import matplotlib
matplotlib.use('Agg')

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

print 'Running Program'

def optimse_Cu(cluster):
	Cu_parameters = {'Cu': [10.960, 2.2780, 0.0855, 1.2240, 2.556]}
	cluster.set_calculator(Gupta(Cu_parameters, cutoff=1000, debug=True))
	dyn = FIRE(cluster)
	try:
		import time
		startTime = time.time()
		dyn.run(fmax=0.01,steps=5000)
		endTime = time.time()
		if not dyn.converged():
			import os
			name = os.path.basename(os.getcwd())
			errorMessage = 'The optimisation of cluster ' + name + ' did not optimise completely.'
			print >> sys.stderr, errorMessage
			print errorMessage
	except:
		print('Local Optimiser Failed for some reason.')

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

run_CNA_program('Au38_cluster_1','Au38_cluster_2')

#########################################################################
start_name = 'Au38_cluster_'
for first_number in range(1,6):
	for second_number in range(1,6):
		if second_number > first_number:
			print str(first_number) + ' Vs ' + str(second_number)
			run_CNA_program(start_name + str(first_number),start_name + str(second_number))
