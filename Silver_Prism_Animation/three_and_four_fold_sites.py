import numpy as np
from ase.visualize import view
from itertools import compress

def get_length(vector):
	return np.linalg.norm(vector)

def get_unit_vector(vector):
	unit_vector = vector / get_length(vector)
	return unit_vector

def get_distance(cluster,index1,index2):
	#xx_dist = atom.x - atom2.x
	#yy_dist = atom1.y - atom2.y
	#zz_dist = atom1.z - atom2.z
	#distance = (xx_dist**2.0 + yy_dist**2.0 + zz_dist**2.0)**0.5
	distance = cluster.get_distance(index1,index2)
	return distance

def get_distance_from_numpy(pos1,pos2):
	xx_dist = pos1[0] - pos2[0]
	yy_dist = pos1[1] - pos2[1]
	zz_dist = pos1[2] - pos2[2]
	distance = (xx_dist**2.0 + yy_dist**2.0 + zz_dist**2.0)**0.5
	return distance

def total_outer_bond_angles(cluster,index1,index2,index3,index4):
	indices = [(index1,index2,index3),(index2,index3,index4),(index3,index4,index1),(index4,index1,index2)]
	try:
		bond_angles = cluster.get_angles(indices)
	except Exception as ee:
		print(indices)
		import pdb; pdb.set_trace()
		exit()
	return all([(80.0 < angle < 100.0) for angle in bond_angles])

def get_three_fold_sites(surface_neighbour_list):
	# find neighbouring sets of triangles
	#import pdb; pdb.set_trace()
	triangles = []
	for index1, indices2 in surface_neighbour_list.items():
		for index2 in indices2:
			neighbours_of_1_as_set = set(indices2)
			neighbours_of_2_as_set = set(surface_neighbour_list[index2])
			intersections = neighbours_of_1_as_set & neighbours_of_2_as_set
			for index3 in intersections:
				triangles.append(tuple(sorted([index1,index2,index3])))
	triangles = sorted(list(set(triangles)),key=lambda x:sorted(x,reverse=True))
	return triangles

def get_applied_three_fold_sites(surface_neighbour_list,triangles,surface_atoms_turned_bulk,indices_to_explore):
	# find neighbouring sets of triangles
	for index in surface_atoms_turned_bulk:
		for ii in range(len(triangles)-1,-1,-1):
			if index in triangles[ii]:
				del triangles[ii]

	for index1 in indices_to_explore:
		indices2 = surface_neighbour_list[index1]
		for index2 in indices2:
			neighbours_of_1_as_set = set(indices2)
			neighbours_of_2_as_set = set(surface_neighbour_list[index2])
			intersections = neighbours_of_1_as_set & neighbours_of_2_as_set
			for index3 in intersections:
				triangles.append(tuple(sorted([index1,index2,index3])))

	triangles = sorted(list(set(triangles)),key=lambda x:sorted(x,reverse=True))
	return triangles

#diff = 1 + (2.0**0.5)
def get_four_fold_sites(surface_neighbour_list,cluster,cutoff):
	extended_distance_across_square = cutoff*(2.0**0.5)
	# find neighbouring sets of squares
	squares = []
	nearly_squares = []
	for index1, indices2 in sorted(surface_neighbour_list.items(),key=lambda x:x[0]):
		for ii in range(len(indices2)):
			first_neighbour_to_index1 = indices2[ii]
			neighbours_of_2_as_set = set(surface_neighbour_list[first_neighbour_to_index1])
			for jj in range(ii+1,len(indices2)):
				second_neighbour_to_index1 = indices2[jj]
				neighbours_of_3_as_set = set(surface_neighbour_list[second_neighbour_to_index1])
				intersections = neighbours_of_2_as_set & neighbours_of_3_as_set
				if len(intersections) > 0:
					for index4 in intersections:
						square = tuple(sorted([index1,first_neighbour_to_index1,index4,second_neighbour_to_index1]))
						if not len(set(square)) == len(square):
							continue
						#print(square)
						#atom1 = cluster[index1]
						#atom4 = cluster[index4]
						#are_1_and_4_neighbours = (index1 in indices_4) and (index4 in indices_1)
						distance14 = get_distance(cluster,index1,index4)
						#atom2 = cluster[first_neighbour_to_index1]
						#atom3 = cluster[second_neighbour_to_index1]
						#are_2_and_3_neighbours = (index2 in indices_3) and (index3 in indices_2)
						distance23 = get_distance(cluster,first_neighbour_to_index1,second_neighbour_to_index1)
						angles_of_a_square = total_outer_bond_angles(cluster,index1,first_neighbour_to_index1,index4,second_neighbour_to_index1)
						#if square == tuple(sorted([1,5,14,15])):
						#	import pdb; pdb.set_trace()
						if (distance14 <= extended_distance_across_square) and (distance23 <= extended_distance_across_square) and angles_of_a_square:	
							#if are_1_and_4_neighbours and are_2_and_3_neighbours:
							squares.append(square)
				else:
					within_distance = get_distance(cluster,first_neighbour_to_index1,second_neighbour_to_index1) < extended_distance_across_square
					indices = [(first_neighbour_to_index1,index1,second_neighbour_to_index1)]
					indices.append((index1,first_neighbour_to_index1,second_neighbour_to_index1))
					indices.append((index1,second_neighbour_to_index1,first_neighbour_to_index1))
					bond_angles = cluster.get_angles(indices)
					ok_angle_1 = (80.0 < bond_angles[0] < 100.0) 
					ok_angle_2 = (80.0 < bond_angles[1] < 100.0) 
					ok_angle_3 = (80.0 < bond_angles[2] < 100.0) 
					if within_distance and ok_angle_1 and ok_angle_2 and ok_angle_3:
						nearly_squares = tuple(index1,first_neighbour_to_index1,second_neighbour_to_index1)
						nearly_squares.append(nearly_square)

	squares = sorted(list(set(squares)),key=lambda x:sorted(x,reverse=True))
	nearly_squares = sorted(list(set(nearly_squares)),key=lambda x:sorted(x,reverse=True))
	return squares, nearly_squares

def get_applied_four_fold_sites(surface_neighbour_list,cluster,cutoff, squares,nearly_squares, surface_atoms_turned_bulk,indices_to_explore):
	extended_distance_across_square = cutoff*(2.0**0.5)

	#import pdb; pdb.set_trace()

	# find neighbouring sets of triangles
	for index in surface_atoms_turned_bulk:
		for ii in range(len(nearly_squares)-1,-1,-1):
			if index in nearly_squares[ii]:
				del nearly_squares[ii]
	for index in surface_atoms_turned_bulk:
		for ii in range(len(squares)-1,-1,-1):
			if index in squares[ii]:
				del squares[ii]
	
	for index1 in indices_to_explore:
		indices2 = surface_neighbour_list[index1]
		for ii in range(len(indices2)):
			first_neighbour_to_index1 = indices2[ii]
			neighbours_of_2_as_set = set(surface_neighbour_list[first_neighbour_to_index1])
			for jj in range(ii+1,len(indices2)):
				second_neighbour_to_index1 = indices2[jj]
				neighbours_of_3_as_set = set(surface_neighbour_list[second_neighbour_to_index1])
				intersections = neighbours_of_2_as_set & neighbours_of_3_as_set
				if len(intersections) > 0:
					for index4 in intersections:
						square = tuple(sorted([index1,first_neighbour_to_index1,index4,second_neighbour_to_index1]))
						if not len(set(square)) == len(square):
							continue
						#print(square)
						#atom1 = cluster[index1]
						#atom4 = cluster[index4]
						#are_1_and_4_neighbours = (index1 in indices_4) and (index4 in indices_1)
						distance14 = get_distance(cluster,index1,index4)
						#atom2 = cluster[first_neighbour_to_index1]
						#atom3 = cluster[second_neighbour_to_index1]
						#are_2_and_3_neighbours = (index2 in indices_3) and (index3 in indices_2)
						distance23 = get_distance(cluster,first_neighbour_to_index1,second_neighbour_to_index1)
						angles_of_a_square = total_outer_bond_angles(cluster,index1,first_neighbour_to_index1,index4,second_neighbour_to_index1)
						#if square == tuple(sorted([1,5,14,15])):
						#	import pdb; pdb.set_trace()
						if (distance14 <= extended_distance_across_square) and (distance23 <= extended_distance_across_square) and angles_of_a_square:	
							#if are_1_and_4_neighbours and are_2_and_3_neighbours:
							squares.append(square)
				else:
					within_distance = get_distance(cluster,first_neighbour_to_index1,second_neighbour_to_index1) < extended_distance_across_square
					indices = [(first_neighbour_to_index1,index1,second_neighbour_to_index1)]
					indices.append((index1,first_neighbour_to_index1,second_neighbour_to_index1))
					indices.append((index1,second_neighbour_to_index1,first_neighbour_to_index1))
					bond_angles = cluster.get_angles(indices)
					ok_angle_1 = (80.0 < bond_angles[0] < 100.0) 
					ok_angle_2 = (80.0 < bond_angles[1] < 100.0) 
					ok_angle_3 = (80.0 < bond_angles[2] < 100.0) 
					if within_distance and ok_angle_1 and ok_angle_2 and ok_angle_3:
						nearly_squares = tuple(index1,first_neighbour_to_index1,second_neighbour_to_index1)
						nearly_squares.append(nearly_square)

	squares = sorted(list(set(squares)),key=lambda x:sorted(x,reverse=True))
	nearly_squares = sorted(list(set(nearly_squares)),key=lambda x:sorted(x,reverse=True))

	return squares, nearly_squares


def get_norm_direction_length(cutoff,distance_1_to_centre):
	return (cutoff**2.0 - distance_1_to_centre**2.0)**0.5

round_to_dp = 2
def same_position(add_position,pos_new_atom):
	#same_x = round(add_position[0],round_to_dp) == round(pos_new_atom[0],round_to_dp)
	#same_y = round(add_position[1],round_to_dp) == round(pos_new_atom[1],round_to_dp)
	#same_z = round(add_position[2],round_to_dp) == round(pos_new_atom[2],round_to_dp)
	#return (same_x and same_y and same_z)
	distance = get_distance_from_numpy(add_position,pos_new_atom)
	return (distance <= 2.0)


def get_position_of_atom_above_square(index1, index2, index3, index4, cluster_positions, centre_of_mass):
	atom1_position = cluster_positions[index1]
	atom2_position = cluster_positions[index2]
	atom3_position = cluster_positions[index3]
	atom4_position = cluster_positions[index4]
	centre_of_square = (atom1_position+atom2_position+atom3_position+atom4_position)/4.0

	adjacent_index = []
	for index, pos in ((index2,atom2_position),(index3,atom3_position),(index4,atom4_position)):
		dist = get_length(pos-atom1_position)
		adjacent_index.append((dist,index,pos))
	adjacent_distance, adjacent_index, adjacent_position = min(adjacent_index)
	
	distance_1_to_centre = get_length(centre_of_square-atom1_position)
	vector_1_to_centre = get_unit_vector(centre_of_square-atom1_position)
	vector_2_to_centre = get_unit_vector(centre_of_square-adjacent_position)
	normal_vector = get_unit_vector(np.cross(vector_1_to_centre,vector_2_to_centre))
	#norm_direction_length = distance_1_to_centre
	norm_direction_length = get_norm_direction_length(adjacent_distance,distance_1_to_centre)
	place_atom_position1 = centre_of_square + norm_direction_length*normal_vector
	place_atom_position2 = centre_of_square - norm_direction_length*normal_vector
	if get_length(place_atom_position1-centre_of_mass) > get_length(place_atom_position2-centre_of_mass):
		add_position = place_atom_position1
	else:
		add_position = place_atom_position2
	#return place_atom_position1, place_atom_position2
	return [add_position]

def get_position_of_atom_above_triangle(index1, index2, index3, cluster_positions, centre_of_mass):
	atom1_position = cluster_positions[index1]
	atom2_position = cluster_positions[index2]
	atom3_position = cluster_positions[index3]
	centre_of_triangle = (atom1_position+atom2_position+atom3_position)/3.0
	adjacent_distance = get_length(atom2_position-atom1_position)
	distance_1_to_centre = get_length(centre_of_triangle-atom1_position)
	vector_1_to_centre = get_unit_vector(centre_of_triangle-atom1_position)
	vector_2_to_centre = get_unit_vector(centre_of_triangle-atom2_position)
	normal_vector = get_unit_vector(np.cross(vector_1_to_centre,vector_2_to_centre))
	#norm_direction_length = distance_1_to_centre
	norm_direction_length = get_norm_direction_length(adjacent_distance,distance_1_to_centre)
	place_atom_position1 = centre_of_triangle + norm_direction_length*normal_vector
	place_atom_position2 = centre_of_triangle - norm_direction_length*normal_vector
	if get_length(place_atom_position1-centre_of_mass) > get_length(place_atom_position2-centre_of_mass):
		add_position = place_atom_position1
	else:
		add_position = place_atom_position2
	#return place_atom_position1, place_atom_position2
	return [add_position]

def get_positions_for_new_atoms(system,triangles,squares,nearly_squares):

	centre_of_mass = system.get_center_of_mass()
	cluster_positions = system.get_positions()

	nearly_squ_pos_new_atoms = []
	nearly_squ_pos_new_atoms_indices = []
	'''
	for index1, index2, index3 in nearly_squares:
		atom1_position = cluster_positions[index1]
		atom2_position = cluster_positions[index2]
		atom3_position = cluster_positions[index3]
		distances = []
		distances.append((index1,index2,index3,get_distance_from_numpy(atom1_position,atom2_position)))
		distances.append((index1,index3,index2,get_distance_from_numpy(atom1_position,atom3_position)))
		distances.append((index2,index3,index1,get_distance_from_numpy(atom2_position,atom3_position)))
		first_index, third_index, second_index, _ = max(distances, key=lambda x:x[3])
		direction1 = cluster_positions[first_index] - cluster_positions[second_index]
		direction2 = cluster_positions[third_index] - cluster_positions[second_index]
		fourth_position = cluster_positions[second_index] + direction1 + direction2
		if any([same_position(add_position,pos_new_atom) for pos_new_atom in nearly_squ_pos_new_atoms]):
			continue
		elif any([same_position(add_position,atom.position) for atom in system]):
			continue
		else:
			nearly_squ_pos_new_atoms.append(add_position)
	'''

	squ_pos_new_atoms = []
	squ_pos_new_atoms_indices = []
	for index1, index2, index3, index4 in squares:
		add_positions = get_position_of_atom_above_square(index1, index2, index3, index4, cluster_positions, centre_of_mass)
		for add_position in add_positions:
			if any([same_position(add_position,pos_new_atom) for pos_new_atom in squ_pos_new_atoms]):
				continue
			elif any([same_position(add_position,pos_new_atom) for pos_new_atom in nearly_squ_pos_new_atoms]):
				continue
			elif any([same_position(add_position,atom_position) for atom_position in system.get_positions()]):
				continue
			else:
				squ_pos_new_atoms.append(add_position)
				squ_pos_new_atoms_indices.append((index1, index2, index3, index4))

	tri_pos_new_atoms = []
	tri_pos_new_atoms_indices = []
	for index1, index2, index3 in triangles:
		add_positions = get_position_of_atom_above_triangle(index1, index2, index3, cluster_positions, centre_of_mass)
		for add_position in add_positions:
			if any([same_position(add_position,pos_new_atom) for pos_new_atom in tri_pos_new_atoms]):
				continue
			elif any([same_position(add_position,pos_new_atom) for pos_new_atom in squ_pos_new_atoms]):
				continue
			elif any([same_position(add_position,pos_new_atom) for pos_new_atom in nearly_squ_pos_new_atoms]):
				continue
			elif any([same_position(add_position,atom_position) for atom_position in system.get_positions()]):
				continue
			else:
				tri_pos_new_atoms.append(add_position)
				tri_pos_new_atoms_indices.append((index1, index2, index3))

	return tri_pos_new_atoms, tri_pos_new_atoms_indices, nearly_squ_pos_new_atoms, nearly_squ_pos_new_atoms_indices, squ_pos_new_atoms, squ_pos_new_atoms_indices

def update_positions_for_new_atoms(system,triangles,squares,nearly_squares,       tri_pos_new_atoms,tri_pos_new_atoms_indices,nearly_squ_pos_new_atoms,nearly_squ_pos_new_atoms_indices,squ_pos_new_atoms,squ_pos_new_atoms_indices,      surface_atoms_turned_bulk,indices_to_explore):
	for index in surface_atoms_turned_bulk:
		for ii in range(len(tri_pos_new_atoms_indices)-1,-1,-1):
			if index in tri_pos_new_atoms_indices[ii]:
				del tri_pos_new_atoms[ii]
				del tri_pos_new_atoms_indices[ii]
		for ii in range(len(nearly_squ_pos_new_atoms_indices)-1,-1,-1):
			if index in nearly_squ_pos_new_atoms_indices[ii]:
				del nearly_squ_pos_new_atoms[ii]
				del nearly_squ_pos_new_atoms_indices[ii]
		for ii in range(len(squ_pos_new_atoms_indices)-1,-1,-1):
			if index in squ_pos_new_atoms_indices[ii]:
				del squ_pos_new_atoms[ii]
				del squ_pos_new_atoms_indices[ii]

	centre_of_mass = system.get_center_of_mass()
	cluster_positions = system.get_positions()

	for index1, index2, index3, index4 in squares[::-1]:
		#if (index1, index2, index3, index4) in squ_pos_new_atoms_indices:
		#	continue
		for index in indices_to_explore:
			if index in (index1, index2, index3, index4):
				break
		else:
			continue
		add_positions = get_position_of_atom_above_square(index1, index2, index3, index4, cluster_positions, centre_of_mass)
		for add_position in add_positions:
			pos_in_triangle = [same_position(add_position,pos_new_atom) for pos_new_atom in tri_pos_new_atoms]
			if any(pos_in_triangle):
				indices_to_remove = list(compress(range(len(pos_in_triangle)), pos_in_triangle))
				#if len(indices_to_remove) > 1:
				#	print('check this out')
				#	import pdb; pdb.set_trace()
				for index_to_remove in indices_to_remove:
					del tri_pos_new_atoms[index_to_remove]
					del tri_pos_new_atoms_indices[index_to_remove]
				squ_pos_new_atoms.append(add_position)
				squ_pos_new_atoms_indices.append((index1, index2, index3, index4))	
				continue
			elif any([same_position(add_position,pos_new_atom) for pos_new_atom in squ_pos_new_atoms]):
				continue
			elif any([same_position(add_position,pos_new_atom) for pos_new_atom in nearly_squ_pos_new_atoms]):
				continue
			elif any([same_position(add_position,atom_position) for atom_position in cluster_positions]):
				continue
			else:
				squ_pos_new_atoms.append(add_position)
				squ_pos_new_atoms_indices.append((index1, index2, index3, index4))

	for index1, index2, index3 in triangles[::-1]:
		#if (index1, index2, index3) in tri_pos_new_atoms_indices:
		#	continue
		for index in indices_to_explore:
			if index in (index1, index2, index3):
				break
		else:
			continue
		add_positions = get_position_of_atom_above_triangle(index1, index2, index3, cluster_positions, centre_of_mass)
		for add_position in add_positions:
			if any([same_position(add_position,pos_new_atom) for pos_new_atom in tri_pos_new_atoms]):
				continue
			elif any([same_position(add_position,pos_new_atom) for pos_new_atom in squ_pos_new_atoms]):
				continue
			elif any([same_position(add_position,pos_new_atom) for pos_new_atom in nearly_squ_pos_new_atoms]):
				continue
			elif any([same_position(add_position,atom_position) for atom_position in cluster_positions]):
				continue
			else:
				tri_pos_new_atoms.append(add_position)
				tri_pos_new_atoms_indices.append((index1, index2, index3))

	return tri_pos_new_atoms, tri_pos_new_atoms_indices, nearly_squ_pos_new_atoms, nearly_squ_pos_new_atoms_indices, squ_pos_new_atoms, squ_pos_new_atoms_indices

def get_distance_numpy(positions_1,positions_2):
	diff = positions_1-positions_2
	#import pdb; pdb.set_trace()
	distance = (np.dot(diff,diff)) ** 0.5
	return distance
