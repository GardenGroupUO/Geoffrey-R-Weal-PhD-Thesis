
structures = []
comparisons = []
no_of_structures = 0
no_of_comparisons = 0

filenames = ['What_I_Think_Au54.txt','What_I_Think_Au55.txt','What_I_Think_Au56.txt']

for filename in filenames:
	print '##################################'
	with open(filename,'r') as file:
		for line in file:
			if line.startswith('#'):
				continue
			line = line.rstrip('\n') 
			if line == '':
				continue
			if '(R' in line:
				continue
			structure1 = line.split()[0]
			if structure1 not in structures:
				structures.append(structure1)
			structure2 = line.split()[1]
			if structure2 not in structures:
				structures.append(structure2)

			if structure1 in structures and structure2 in structures:
				comparisons.append((structure1,structure2))

	print len(structures)
	print len(comparisons)
	no_of_structures += len(structures)
	no_of_comparisons += len(comparisons)

print '##################################'
print '##################################'
print '##################################'
print no_of_structures
print no_of_comparisons