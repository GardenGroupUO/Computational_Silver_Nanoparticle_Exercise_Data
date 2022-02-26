import os

from ase import Atom
from ase.io import read, write, Trajectory
from ase.visualize import view

try:
	from surface_finder import doubledict, neighbourlist
	from surface_finder import get_surface_atoms, get_distance

	from three_and_four_fold_sites import get_three_fold_sites, get_four_fold_sites
	from three_and_four_fold_sites import get_applied_three_fold_sites, get_applied_four_fold_sites
	from three_and_four_fold_sites import get_positions_for_new_atoms, update_positions_for_new_atoms, same_position
except Exception as ee:
	from Computational_Silver_Nanoparticle_Exercise_Data.Silver_Prism_Animation.surface_finder import doubledict, neighbourlist
	from Computational_Silver_Nanoparticle_Exercise_Data.Silver_Prism_Animation.surface_finder import get_surface_atoms, get_distance

	from Computational_Silver_Nanoparticle_Exercise_Data.Silver_Prism_Animation.three_and_four_fold_sites import get_three_fold_sites, get_four_fold_sites
	from Computational_Silver_Nanoparticle_Exercise_Data.Silver_Prism_Animation.three_and_four_fold_sites import get_applied_three_fold_sites, get_applied_four_fold_sites
	from Computational_Silver_Nanoparticle_Exercise_Data.Silver_Prism_Animation.three_and_four_fold_sites import get_positions_for_new_atoms, update_positions_for_new_atoms, same_position

from ase.calculators.emt import EMT
from ase.optimize import FIRE

from random import uniform, randrange

def silver_nanoprism_growing_model(path_to_input,chance_of_creating_new_100_surface_111_surface_bromine_capping,max_no_of_atoms_added_in_simulation=1000):

	if not (sum(chance_of_creating_new_100_surface_111_surface_bromine_capping) == 1.0):
		raise Exception('chance_of_creating_new_100_surface_111_surface_bromine_capping must be between 0.0 and 1.0. You gave chance_of_creating_new_100_surface_111_surface_bromine_capping='+str(chance_of_creating_new_100_surface_111_surface_bromine_capping))

	barrier_111_100 = chance_of_creating_new_100_surface_111_surface_bromine_capping[0]
	barrier_100_cap = chance_of_creating_new_100_surface_111_surface_bromine_capping[0] + chance_of_creating_new_100_surface_111_surface_bromine_capping[1]

	system = read(path_to_input)
	system.set_tags(0)
	system.set_calculator(EMT())

	traj_path = '.'.join(path_to_input.split('.')[:-1:])+'_animation.traj'

	nanoparticle_symbol = system[0].symbol

	cutoff = 3.0

	cluster_positions = system.get_positions()

	print('making initial distance matrix')
	print('making initial full neighbours matrix')
	distances_between_atoms = doubledict()
	full_neighbourlist = neighbourlist()
	for index1 in range(len(cluster_positions)):
		for index2 in range(index1+1,len(cluster_positions)):
			distance = get_distance(cluster_positions[index1],cluster_positions[index2])
			distances_between_atoms[(index1,index2)] = distance
			if distance <= cutoff:
				full_neighbourlist.set(index1,index2)

	print('Getting surface neighbour lists')
	surface_neighbourlist = get_surface_atoms(system,distances_between_atoms,full_neighbourlist,cutoff,last_index=True)
	print('getting triangle surfaces')
	triangles = get_three_fold_sites(surface_neighbourlist)
	print('getting square surfaces')
	squares, nearly_squares = get_four_fold_sites(surface_neighbourlist,system,cutoff)
	print('getting new possible positions')
	tri_pos_new_atoms, tri_pos_new_atoms_indices, nearly_squ_pos_new_atoms, nearly_squ_pos_new_atoms_indices, squ_pos_new_atoms, squ_pos_new_atoms_indices = get_positions_for_new_atoms(system,triangles,squares,nearly_squares)

	tags = system.get_tags() #get_chemical_symbols()
	for index in range(len(tags)):
		tags[index] = 0 # 'Ag'
	for indices in squares:
		for index in indices:
			tags[index] = 1 # 'Fe'
	system.set_tags(tags) #set_chemical_symbols(tags)

	if os.path.exists(traj_path):
		os.remove(traj_path)
	with Trajectory(traj_path,'a') as traj_file:
		traj_file.write(system.copy())

	surface_data = []
	surface_data.append(surface_neighbourlist.copy())
	all_squares = []
	all_squares.append(list(squares))
	all_squ_pos_new_atoms_indices = []
	all_squ_pos_new_atoms_indices.append(list(squ_pos_new_atoms_indices)) 

	counter = 0
	while counter < max_no_of_atoms_added_in_simulation:
		counter += 1
		print('----------------------------------')
		print('Adding atom '+str(counter))
		square_triangle_or_cap = uniform(0, 1)
		if square_triangle_or_cap <= barrier_111_100:
			positions_to_add = squ_pos_new_atoms
			positions_to_add_index = squ_pos_new_atoms_indices
			symbol = nanoparticle_symbol
		elif barrier_111_100 < square_triangle_or_cap <= barrier_100_cap:
			positions_to_add = tri_pos_new_atoms
			positions_to_add_index = tri_pos_new_atoms_indices
			symbol = nanoparticle_symbol
		else:
			positions_to_add = squ_pos_new_atoms + tri_pos_new_atoms
			positions_to_add_index = squ_pos_new_atoms_indices + tri_pos_new_atoms_indices	
			symbol = 'Br'	

		print('squares: '+str(len(squ_pos_new_atoms)))
		print('triangles: '+str(len(tri_pos_new_atoms)))
		if len(squ_pos_new_atoms) == 0:
			break

		random_number = randrange(0, len(positions_to_add), 1)
		random_position = positions_to_add[random_number].copy()
		index_set_to_check = tuple(positions_to_add_index[random_number])

		del positions_to_add[random_number]
		del positions_to_add_index[random_number]
		if any([same_position(random_position,atom.position) for atom in system]):
			for index in range(len(squ_pos_new_atoms)-1,-1,-1):
				check_position = squ_pos_new_atoms[index]
				if same_position(random_position,check_position):
					del squ_pos_new_atoms[index]
					del squ_pos_new_atoms_indices[index]
			for index in range(len(tri_pos_new_atoms)-1,-1,-1):
				check_position = tri_pos_new_atoms[index]
				if same_position(random_position,check_position):
					del tri_pos_new_atoms[index]
					del tri_pos_new_atoms_indices[index]

		atom = Atom(symbol=symbol,position=random_position,tag=counter)

		if len(index_set_to_check) == 3:
			if index_set_to_check in triangles:
				triangles.remove(index_set_to_check)
			if index_set_to_check in nearly_squares:
				nearly_squares.remove(index_set_to_check)
		elif len(index_set_to_check) == 4:
			if index_set_to_check in squares:
				squares.remove(index_set_to_check)
		
		system.append(atom)

		cluster_positions = system.get_positions()

		#print('making initial full neighbours matrix')
		end_of_system = len(cluster_positions)-1
		for index in range(end_of_system):
			distance = get_distance(cluster_positions[index],cluster_positions[end_of_system])
			distances_between_atoms[(index,end_of_system)] = distance
			if distance <= cutoff:
				full_neighbourlist.set(index,end_of_system)

		try:
			len(full_neighbourlist.get(end_of_system))
		except:
			print('check')
			dists = []
			end_of_system = len(cluster_positions)-1
			for index in range(end_of_system):
				distance = get_distance(cluster_positions[index],cluster_positions[end_of_system])
				dists.append(distance)
			print(min(dists))
			import pdb; pdb.set_trace()
			return traj_path
			

		#print('Getting surface neighbour lists')
		surface_atoms_turned_bulk = []
		for index in full_neighbourlist.get(end_of_system):
			if len(full_neighbourlist.get(index)) == 12:
				surface_neighbourlist.remove(index)
				surface_atoms_turned_bulk.append(index)
			else:
				surface_neighbourlist.set(index,end_of_system)
		surface_data.append(surface_neighbourlist.copy())

		indices_to_explore = surface_neighbourlist[end_of_system] + [end_of_system]
		#print('getting triangle surfaces')
		triangles = get_applied_three_fold_sites(surface_neighbourlist,triangles,surface_atoms_turned_bulk,indices_to_explore)
		#print('getting square surfaces')
		squares, nearly_squares = get_applied_four_fold_sites(surface_neighbourlist,system,cutoff,   squares,nearly_squares,   surface_atoms_turned_bulk,indices_to_explore)
		all_squares.append(list(squares))
		#print('getting new possible positions')
		tri_pos_new_atoms, tri_pos_new_atoms_indices, nearly_squ_pos_new_atoms, nearly_squ_pos_new_atoms_indices, squ_pos_new_atoms, squ_pos_new_atoms_indices = update_positions_for_new_atoms(system,triangles,squares,nearly_squares,      tri_pos_new_atoms,tri_pos_new_atoms_indices,nearly_squ_pos_new_atoms,nearly_squ_pos_new_atoms_indices,squ_pos_new_atoms,squ_pos_new_atoms_indices,      surface_atoms_turned_bulk,indices_to_explore)
		all_squ_pos_new_atoms_indices.append(list(squ_pos_new_atoms_indices))

		tags = system.get_tags() #get_chemical_symbols()
		for index in range(len(tags)):
			tags[index] = 0 #'Ag'
		for indices in squares: #squ_pos_new_atoms_indices:
			for index in indices:
				tags[index] = 1 #'Fe'
		system.set_tags(tags) #set_chemical_symbols(tags)

		#print('Adding image to Traj')
		with Trajectory(traj_path,'a') as traj_file:
			traj_file.write(system.copy())

	print('----------------------------------')
	print('The simulation has now finished')
	print('----------------------------------')
	return traj_path