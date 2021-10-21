from ase.cluster import Octahedron
from ase.io import Trajectory, read
from ase.visualize import view

clusters = []
clusters = Trajectory('example.traj',mode='w')

def make_image(length,cutoff):
	atoms = Octahedron('Ag',length,cutoff)
	atoms.set_cell([50,50,50])
	atoms.rotate(45,'x')
	atoms.rotate(30,'y')
	atoms.center()
	#clusters.append(atoms)
	clusters.write(atoms)

for cutoff in range(1,8):
	length = 2*cutoff+3
	print((length, cutoff))
	make_image(length,cutoff)

#clusters = Trajectory('example.traj',mode='r')
clusters = read('example.traj',index=':')
view(clusters,viewer='x3d')