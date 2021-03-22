import sys,os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir, os.path.pardir)))
from Make_Interpolated_Structures import Make_Interpolated_Structures

symbol = 'Au'
designed_atom_number = 38

# Deca
# 111 Corner
starting_details = [3,2,0]
ending_details = []
atom_indices_to_remove = [3]
Make_Interpolated_Structures(symbol,designed_atom_number,starting_details,ending_details,atom_indices_to_remove)

# 111 edge
starting_details = [3,2,0]
ending_details = []
atom_indices_to_remove = [15]
Make_Interpolated_Structures(symbol,designed_atom_number,starting_details,ending_details,atom_indices_to_remove)

# 100 Corner
starting_details = [3,2,0]
ending_details = []
atom_indices_to_remove = [32]
Make_Interpolated_Structures(symbol,designed_atom_number,starting_details,ending_details,atom_indices_to_remove)

# 100 Corner
starting_details = [3,2,0]
ending_details = []
atom_indices_to_remove = [34]
Make_Interpolated_Structures(symbol,designed_atom_number,starting_details,ending_details,atom_indices_to_remove)