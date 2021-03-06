import os

from ase import Atom
from ase.io import read, write, Trajectory
from ase.visualize import view

from surface_finder import doubledict, neighbourlist
from surface_finder import get_surface_atoms, get_distance

from three_and_four_fold_sites import get_three_fold_sites, get_four_fold_sites
from three_and_four_fold_sites import get_applied_three_fold_sites, get_applied_four_fold_sites
from three_and_four_fold_sites import get_positions_for_new_atoms, update_positions_for_new_atoms, same_position

from ase.calculators.emt import EMT
from ase.optimize import FIRE


from random import randrange

prefix = 'mid3' # 'mid3'  'small' 'large'

system = read(prefix+'_initial_seed.xyz')
system.set_tags(0)
system.set_calculator(EMT())
dyn = FIRE(system)
print('locally optimising')
#dyn.run(fmax=0.001,steps=5000)


traj_name = prefix+'_animation.traj'

symbol = system[0].symbol

cutoff=3.0

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

tags = system.get_chemical_symbols()
for index in range(len(tags)):
	tags[index] = 'Ag'
for indices in squares:
	for index in indices:
		tags[index] = 'Fe'
system.set_chemical_symbols(tags)

if os.path.exists(traj_name):
	os.remove(traj_name)
with Trajectory(traj_name,'a') as traj_file:
	traj_file.write(system.copy())

tags = system.get_chemical_symbols()
for index in range(len(tags)):
	tags[index] = 'Ag'
for indices in squ_pos_new_atoms_indices:
	for index in indices:
		tags[index] = 'Ni'
system.set_chemical_symbols(tags)

if os.path.exists('check_'+traj_name):
	os.remove('check_'+traj_name)
with Trajectory('check_'+traj_name,'a') as check_traj_file:
	check_traj_file.write(system.copy())

surface_data = []
surface_data.append(surface_neighbourlist.copy())
all_squares = []
all_squares.append(list(squares))
all_squ_pos_new_atoms_indices = []
all_squ_pos_new_atoms_indices.append(list(squ_pos_new_atoms_indices))

counter = 0
while True:
	counter += 1
	print(counter)

	square_or_triangle = randrange(0, 100, 1)
	positions_to_add = tri_pos_new_atoms if (square_or_triangle == 0) else squ_pos_new_atoms # nearly_squ_pos_new_atoms
	positions_to_add_index = tri_pos_new_atoms_indices if (square_or_triangle == 0) else squ_pos_new_atoms_indices

	'''
	for random_position in (nearly_squ_pos_new_atoms + squ_pos_new_atoms + tri_pos_new_atoms):
		if any([same_position(random_position,atom.position) for atom in system]):
			print('huh???')
			import pdb; pdb.set_trace()
			exit()
	'''
	print('squares: '+str(len(squ_pos_new_atoms)))
	print('triangles: '+str(len(tri_pos_new_atoms)))
	print('----------------------------------')
	if len(squ_pos_new_atoms) == 0:
		print('finishing')
		break

	counter_2 = 0
	while True:
		random_number = randrange(0, len(positions_to_add), 1)
		random_position = positions_to_add[random_number]
		if any([same_position(random_position,atom.position) for atom in system]):
			print('huh')
			import pdb; pdb.set_trace()
			exit()	
		else:
			break
		#counter_2 += 1
		#if counter_2 > 2000:
		#	exit('got to 2000')

	atom = Atom(symbol=symbol,position=random_position,tag=counter)

	index_set_to_check = positions_to_add_index[random_number]
	#import pdb; pdb.set_trace()
	if len(index_set_to_check) == 3:
		if index_set_to_check in triangles:
			triangles.remove(index_set_to_check)
		if index_set_to_check in nearly_squares:
			nearly_squares.remove(index_set_to_check)
	elif len(index_set_to_check) == 4:
		if index_set_to_check in squares:
			squares.remove(index_set_to_check)
	#import pdb; pdb.set_trace()

	del positions_to_add[random_number]
	del positions_to_add_index[random_number]
	
	system.append(atom)
	#dyn = FIRE(system)#,logfile=None)
	#print('locally optimising')
	#dyn.run(fmax=0.001,steps=5000)

	cluster_positions = system.get_positions()

	#print('making initial distance matrix')
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

	tags = system.get_chemical_symbols()
	for index in range(len(tags)):
		tags[index] = 'Ag'
	for indices in squares: #squ_pos_new_atoms_indices:
		for index in indices:
			tags[index] = 'Fe'
	system.set_chemical_symbols(tags)
	#import pdb; pdb.set_trace()

	with Trajectory(traj_name,'a') as traj_file:
		traj_file.write(system.copy())

	tags = system.get_chemical_symbols()
	for index in range(len(tags)):
		tags[index] = 'Ag'
	for indices in squ_pos_new_atoms_indices:
		for index in indices:
			tags[index] = 'Ni'
	system.set_chemical_symbols(tags)

	with Trajectory('check_'+traj_name,'a') as check_traj_file:
		check_traj_file.write(system.copy())

import pdb; pdb.set_trace()