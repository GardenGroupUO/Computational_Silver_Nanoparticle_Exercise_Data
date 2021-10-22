from ase.io import read
from ase.visualize import view
from silver_nanoprism_growing_model import silver_nanoprism_growing_model

prefix = 'mid3' # 'mid3'  'small' 'large'
path_to_input = prefix+'_initial_seed.xyz'

change_of_creating_new_100_surface = 0.9

traj_path = silver_nanoprism_growing_model(path_to_input,change_of_creating_new_100_surface=change_of_creating_new_100_surface)

animation = read(traj_path,index=':')
view(animation)