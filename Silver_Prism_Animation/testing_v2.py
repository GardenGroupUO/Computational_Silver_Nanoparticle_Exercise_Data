from ase.io import read
from ase.visualize import view
from silver_nanoprism_growing_model_with_capping import silver_nanoprism_growing_model

prefix = 'small' # 'mid3'  'small' 'large'
path_to_file = prefix+'_initial_seed.xyz'

chance_of_creating_new_100_surface_111_surface_bromine_capping = [0.8,0.02,0.18]
max_no_of_atoms_added_in_simulation = 1000

traj_path = silver_nanoprism_growing_model(path_to_file,chance_of_creating_new_100_surface_111_surface_bromine_capping=chance_of_creating_new_100_surface_111_surface_bromine_capping,max_no_of_atoms_added_in_simulation=max_no_of_atoms_added_in_simulation)

animation = read(traj_path,index=':')
view(animation)

