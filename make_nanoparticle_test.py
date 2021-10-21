from size import make_nanoparticle
nanoparticle = make_nanoparticle()

atom1 = nanoparticle[475]
atom2 = nanoparticle[9955]

def get_distance(atom1,atom2):
	x1, y1, z1 = atom1.position
	x2, y2, z2 = atom2.position
	x_dist = x1-x2
	y_dist = y1-y2
	z_dist = z1-z2
	distance = (x_dist**2.0 + y_dist**2.0 + z_dist**2.0) ** 0.5
	print(distance)

get_distance(atom1,atom2)