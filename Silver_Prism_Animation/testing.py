from ase.io import read
from ase.visualize import view
from silver_nanoprism_growing_model import silver_nanoprism_growing_model

prefix = 'small' # 'mid3'  'small' 'large'
path_to_input = prefix+'_initial_seed.xyz'

change_of_creating_new_100_surface = 0.9
max_no_of_atoms_added_in_simulation = 1000

traj_path = silver_nanoprism_growing_model(path_to_input,change_of_creating_new_100_surface=change_of_creating_new_100_surface,max_no_of_atoms_added_in_simulation=max_no_of_atoms_added_in_simulation)

animation = read(traj_path,index=':')
view(animation)