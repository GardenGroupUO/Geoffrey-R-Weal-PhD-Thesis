import os

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def split_Cluster(cluster_name):
	for index in range(len(cluster_name)):
		if RepresentsInt(cluster_name[index]):
			break
	cluster_symbol = cluster_name[:index]
	cluster_atom_number = int(cluster_name[index:])
	return cluster_symbol, cluster_atom_number

class Make_Lyx_File:
	def __init__(self,lyx_filename):
		self.lyx_filename = lyx_filename+'.lyx'
		self.Lyx_File = open(self.lyx_filename,'w')
		self.Start_of_file()
		self.Collect_CNA_signature_information()
		self.Lyx_File.close()

	def Start_of_file(self):
		self.Lyx_File.write('%% LyX 2.2.3 created this file.  For more info, see http://www.lyx.org/.\n')
		self.Lyx_File.write('%% Do not edit unless you really know what you are doing.\n')
		self.Lyx_File.write('\documentclass[english]{scrreprt}\n')
		self.Lyx_File.write('\usepackage[T1]{fontenc}\n')
		self.Lyx_File.write('\usepackage[latin9]{inputenc}\n')
		self.Lyx_File.write('\usepackage{geometry}\n')
		self.Lyx_File.write('\geometry{verbose,tmargin=1cm,bmargin=1cm,lmargin=1cm,rmargin=1cm}\n')
		self.Lyx_File.write('\setcounter{secnumdepth}{3}\n')
		self.Lyx_File.write('\setcounter{tocdepth}{3}\n')
		self.Lyx_File.write('\usepackage{color}\n')
		self.Lyx_File.write('\usepackage{babel}\n')
		self.Lyx_File.write('\usepackage{float}\n')
		self.Lyx_File.write('\usepackage{graphicx}\n')
		self.Lyx_File.write('\usepackage{setspace}\n')
		self.Lyx_File.write('\usepackage[unicode=true,pdfusetitle,\n')
		self.Lyx_File.write(' bookmarks=true,bookmarksnumbered=false,bookmarksopen=false,\n')
		self.Lyx_File.write(' breaklinks=true,pdfborder={0 0 0},pdfborderstyle={},backref=false,colorlinks=true]\n')
		self.Lyx_File.write(' {hyperref}\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\makeatletter\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LyX specific LaTeX commands.\n')
		self.Lyx_File.write('%% A simple dot to overcome graphicx limitations\n')
		self.Lyx_File.write('\\newcommand{\lyxdot}{.}\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User specified LaTeX commands.\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\usepackage{fancyhdr}\n')
		self.Lyx_File.write('\usepackage{lastpage}\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\pagestyle{fancy}\n')
		self.Lyx_File.write('\\fancyhf{}\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\\rfoot{Page \\thepage \\hspace{1pt} of \pageref{LastPage}}\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\pagenumbering{arabic}\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\@ifundefined{showcaptionsetup}{}{%\n')
		self.Lyx_File.write(' \PassOptionsToPackage{caption=false}{subfig}}\n')
		self.Lyx_File.write('\usepackage{subfig}\n')
		self.Lyx_File.write('\makeatother\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\\begin{document}\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\\title{CNA Data}\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\\author{Geoffrey Weal}\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\date{\date{\\today}}\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\maketitle\n')
		self.Lyx_File.write('\\tableofcontents{}\n')
		self.Lyx_File.write('\n')
		self.Lyx_File.write('\setcounter{page}{1}\n')
		self.Lyx_File.write('\n')

	def Collect_CNA_signature_information(self):
		files = os.listdir('Diagrams_Signature_Analysis')
		for index in range(len(files)-1,-1,-1):
			if not files[index].endswith('.png'):
				del files[index]
		CNA_signature_profiles = CNA_Signature_Profiles()
		for signature_profile in files:
			print 'Processing ' + str(signature_profile)
			CNA_signature_profiles.add(signature_profile)
		import pdb; pdb.set_trace()


	#def Collect_CNA_signature_information_for_one_cluster(self,signature_profile):

class CNA_Signature_Profiles:
	def __init__(self):
		self.cluster_types = []

	def add(self,signature_profile_filename):
		signature_profile_filename = signature_profile_filename
		signature_profile_split = signature_profile_filename.split('_')
		cluster_name = signature_profile_split[0]
		cluster_symbol, cluster_atom_number = split_Cluster(cluster_name)
		motif_type = signature_profile_split[1]
		temp_string = '_'.join(signature_profile_split[3:])
		temp_string = temp_string.split('(')
		name = temp_string[0][:-1].replace('_',' ')
		signature = '('+temp_string[1].split(')')[0]+')'
		signature = eval(signature)
		#import pdb; pdb.set_trace()

		if self.cluster_types == []:
			self.cluster_types.append(CNA_Signature_Profile(cluster_symbol,cluster_atom_number,motif_type,name))
			self.cluster_types[0].add(signature,signature_profile_filename)
		elif self.cluster_types[-1].same_cluster(cluster_symbol,cluster_atom_number,motif_type,name) == 'Higher':
			self.cluster_types.append(CNA_Signature_Profile(cluster_symbol,cluster_atom_number,motif_type,name))
			self.cluster_types[-1].add(signature,signature_profile_filename)
		else:
			for index in range(len(self.cluster_types)):
				same_cluster = self.cluster_types[index].same_cluster(cluster_symbol,cluster_atom_number,motif_type,name)
				if same_cluster == 'Lower':
					self.cluster_types.insert(index,CNA_Signature_Profile(cluster_symbol,cluster_atom_number,motif_type,name))
					self.cluster_types[index].add(signature,signature_profile_filename)
					break
				elif same_cluster == 'Same':
					self.cluster_types[index].add(signature,signature_profile_filename)
					break

	def __str__(self):
		toString = '------------------------------------------------\n'
		toString += self.cluster_name+str(self.cluster_name)+' '+self.motif_type+': '+self.name+'\n'
		toString += 
		toString += '------------------------------------------------\n'

class CNA_Signature_Profile:
	def __init__(self,cluster_symbol,cluster_atom_number,motif_type,name):
		self.cluster_symbol = cluster_symbol
		self.cluster_atom_number = cluster_atom_number
		self.motif_type = motif_type
		if name.endswith('Opt'):
			name = name.replace('Opt','')
		self.name = name
		self.signature_profiles = []
		self.signature_profiles_Opt = []

	def add(self,signature,signature_profile_filename):
		if 'Opt' in signature_profile_filename:
			signature_profiles_to_add_to = self.signature_profiles_Opt
		else:
			signature_profiles_to_add_to = self.signature_profiles

		if signature_profiles_to_add_to == []:
			signature_profiles_to_add_to.append([signature,signature_profile_filename])
		elif signature < signature_profiles_to_add_to[-1][0]:
			signature_profiles_to_add_to.append([signature,signature_profile_filename])
		else:
			for index in range(len(signature_profiles_to_add_to)):
				if signature > signature_profiles_to_add_to[index][0]:
					signature_profiles_to_add_to.insert(index,[signature,signature_profile_filename])
					break
				elif signature == signature_profiles_to_add_to[index][0] and signature_profile_filename == self.signature_profiles[index][1]:
					print 'ERROR HERE!!!!'
					import pdb; pdb.set_trace()
					exit()

	def same_cluster(self,cluster_symbol,cluster_atom_number,motif_type,name):
		if   cluster_symbol < self.cluster_symbol:
			return 'Lower'
		elif cluster_symbol == self.cluster_symbol:
			if   cluster_atom_number < self.cluster_atom_number:
				return 'Lower'
			elif cluster_atom_number == self.cluster_atom_number:
				if   motif_type < self.motif_type:
					return 'Lower'
				elif motif_type == self.motif_type:
					#import pdb; pdb.set_trace()
					if name.endswith(' Opt'):
						name = name.replace(' Opt','')
					if   name < self.name:
						return 'Lower'
					elif name == self.name:
						return "Same"
					elif name > self.name:
						return 'Higher'
				elif motif_type > self.motif_type:
					return 'Higher'
			elif cluster_atom_number > self.cluster_atom_number:
				return 'Higher'
		elif cluster_symbol > self.cluster_symbol:
			return 'Higher'
lyx_filename = 'text'

Make_Lyx_File(lyx_filename)







