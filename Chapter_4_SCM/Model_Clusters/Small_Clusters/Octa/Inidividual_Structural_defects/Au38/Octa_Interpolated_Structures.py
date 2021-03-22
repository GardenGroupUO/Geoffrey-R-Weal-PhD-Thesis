import sys,os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir, os.path.pardir)))
from Make_Interpolated_Structures import Make_Interpolated_Structures

symbol = 'Au'
designed_atom_number = 38

# Octa
########################################################################
# Corner Adatom

# 100 Corner
starting_details = [4,0]
ending_details = [4,1]
atom_indices_to_remove = [0,20,6,11,29,43]
atom_indices_to_remove.remove(43)
atom_indices_to_remove.append(1)
Make_Interpolated_Structures(symbol,designed_atom_number,starting_details,ending_details,atom_indices_to_remove)

# 100 Edge
starting_details = [4,0]
ending_details = [4,1]
atom_indices_to_remove = [0,20,6,11,29,43]
atom_indices_to_remove.remove(43)
atom_indices_to_remove.append(16)
Make_Interpolated_Structures(symbol,designed_atom_number,starting_details,ending_details,atom_indices_to_remove)

########################################################################
# edge Adatom

# 100 Corner
starting_details = [6,0]
ending_details = [4,1]
atom_indices_to_remove =  [  0,  1,  2,  3,  4]
atom_indices_to_remove += [ 23, 26, 17, 16,  5, 12,  9, 22]
atom_indices_to_remove += [ 6, 7, 10, 11, 21, 25, 29, 28, 27, 20, 18, 8]
atom_indices_to_remove += [30,37,34, 49, 46, 67, 68, 83, 84, 87, 76, 79, 62, 61,42,41]
atom_indices_to_remove += [31,32,35, 36, 47, 48, 66, 70, 82, 86, 90, 89, 88, 81,80,65,63,45,43,33]
atom_indices_to_remove += [91,93,92,101,100,116,117,128,129,131,124,127,114,113,98,97]
atom_indices_to_remove += [95,96,103,104,118,122,130,126,125,115,111,99]
atom_indices_to_remove += [132,134,133,141,142,144,139,138]
atom_indices_to_remove += [137,143,140,136,145]
atom_indices_to_remove += [135,64,85,39,13,52]
atom_indices_to_remove.remove(49)
atom_indices_to_remove.append(59)
Make_Interpolated_Structures(symbol,designed_atom_number,starting_details,ending_details,atom_indices_to_remove)

# 100 Edge
starting_details = [6,0]
ending_details = [4,1]
atom_indices_to_remove =  [  0,  1,  2,  3,  4]
atom_indices_to_remove += [ 23, 26, 17, 16,  5, 12,  9, 22]
atom_indices_to_remove += [ 6, 7, 10, 11, 21, 25, 29, 28, 27, 20, 18, 8]
atom_indices_to_remove += [30,37,34, 49, 46, 67, 68, 83, 84, 87, 76, 79, 62, 61,42,41]
atom_indices_to_remove += [31,32,35, 36, 47, 48, 66, 70, 82, 86, 90, 89, 88, 81,80,65,63,45,43,33]
atom_indices_to_remove += [91,93,92,101,100,116,117,128,129,131,124,127,114,113,98,97]
atom_indices_to_remove += [95,96,103,104,118,122,130,126,125,115,111,99]
atom_indices_to_remove += [132,134,133,141,142,144,139,138]
atom_indices_to_remove += [137,143,140,136,145]
atom_indices_to_remove += [135,64,85,39,13,52]
atom_indices_to_remove.remove(49)
atom_indices_to_remove.append(57)
Make_Interpolated_Structures(symbol,designed_atom_number,starting_details,ending_details,atom_indices_to_remove)

########################################################################
# face Adatom

# 100 Corner
starting_details = [6,0]
ending_details = [4,1]
atom_indices_to_remove =  [  0,  1,  2,  3,  4]
atom_indices_to_remove += [ 23, 26, 17, 16,  5, 12,  9, 22]
atom_indices_to_remove += [ 6, 7, 10, 11, 21, 25, 29, 28, 27, 20, 18, 8]
atom_indices_to_remove += [30,37,34, 49, 46, 67, 68, 83, 84, 87, 76, 79, 62, 61,42,41]
atom_indices_to_remove += [31,32,35, 36, 47, 48, 66, 70, 82, 86, 90, 89, 88, 81,80,65,63,45,43,33]
atom_indices_to_remove += [91,93,92,101,100,116,117,128,129,131,124,127,114,113,98,97]
atom_indices_to_remove += [95,96,103,104,118,122,130,126,125,115,111,99]
atom_indices_to_remove += [132,134,133,141,142,144,139,138]
atom_indices_to_remove += [137,143,140,136,145]
atom_indices_to_remove += [135,64,85,39,13,52]
atom_indices_to_remove.remove(10)
atom_indices_to_remove.append(59)
Make_Interpolated_Structures(symbol,designed_atom_number,starting_details,ending_details,atom_indices_to_remove)

# 100 Edge
starting_details = [6,0]
ending_details = [4,1]
atom_indices_to_remove =  [  0,  1,  2,  3,  4]
atom_indices_to_remove += [ 23, 26, 17, 16,  5, 12,  9, 22]
atom_indices_to_remove += [ 6, 7, 10, 11, 21, 25, 29, 28, 27, 20, 18, 8]
atom_indices_to_remove += [30,37,34, 49, 46, 67, 68, 83, 84, 87, 76, 79, 62, 61,42,41]
atom_indices_to_remove += [31,32,35, 36, 47, 48, 66, 70, 82, 86, 90, 89, 88, 81,80,65,63,45,43,33]
atom_indices_to_remove += [91,93,92,101,100,116,117,128,129,131,124,127,114,113,98,97]
atom_indices_to_remove += [95,96,103,104,118,122,130,126,125,115,111,99]
atom_indices_to_remove += [132,134,133,141,142,144,139,138]
atom_indices_to_remove += [137,143,140,136,145]
atom_indices_to_remove += [135,64,85,39,13,52]
atom_indices_to_remove.remove(10)
atom_indices_to_remove.append(57)
Make_Interpolated_Structures(symbol,designed_atom_number,starting_details,ending_details,atom_indices_to_remove)





