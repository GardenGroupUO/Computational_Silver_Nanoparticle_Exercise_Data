from ase.cluster import Octahedron
from ase.visualize import view

clusters = []

def make_image(length,cutoff):
    atoms = Octahedron('Ag',length,cutoff)
    atoms.set_cell([50,50,50])
    atoms.rotate(45,'x')
    atoms.rotate(30,'y')
    atoms.center()
    clusters.append(atoms)

for cutoff in range(1,8):
    length = 2*cutoff+3
    print((length, cutoff))
    make_image(length,cutoff)



from x3d_viewer import view_x3d
view_x3d(clusters,'testing')