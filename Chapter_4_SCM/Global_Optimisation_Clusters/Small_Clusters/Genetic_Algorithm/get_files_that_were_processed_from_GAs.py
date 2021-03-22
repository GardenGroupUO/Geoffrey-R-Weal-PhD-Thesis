import os
from shutil import rmtree, copyfile

def get_clusters_from_text_files(clusters_text_path):
	cluster_names = []
	with open(clusters_text_path,'r') as clusters_TXT:
		for line in clusters_TXT:
			if line.startswith('#') or line.startswith('\n'):
				continue
			line = line.rstrip().split()
			cluster1 = int(line[0])
			cluster2 = int(line[1])
			for cluster in [cluster1, cluster2]:
				if not cluster in cluster_names:
					cluster_names.append(cluster)
	return cluster_names

def grab_clusters_from_folder(cluster_type, cluster_names):
	new_folder = cluster_type
	if os.path.exists(new_folder):
		rmtree(new_folder)
	os.mkdir(new_folder)
	original_folder = str(cluster_type)+'_all_structures'
	for cluster in os.listdir(original_folder):
		for cluster_name_to_accept in cluster_names:
			if cluster.startswith(str(cluster_name_to_accept)):
				print('Moving: '+str((original_folder+'/'+cluster,new_folder+'/'+cluster)))
				copyfile(original_folder+'/'+cluster,new_folder+'/'+cluster)
				break

def get_clusters(cluster_type):
	print('============================================================================')
	print(cluster_type)
	what_I_think_file = 'What_I_Think_Strict/What_I_Think_'+str(cluster_type)+'.txt'
	cluster_names = get_clusters_from_text_files(what_I_think_file)
	print(len(cluster_names))
	grab_clusters_from_folder(cluster_type, cluster_names)
	print('============================================================================')



for cluster_type in ['Au38','Cu37']:
	get_clusters(cluster_type)