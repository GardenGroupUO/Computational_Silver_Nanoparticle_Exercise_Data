from ase.build import fcc111
from ase.io import write

symbol = 'Ag'
surf_e_100 = 1.1808
                    
surf_e_111 = 1.0
esurf = [surf_e_100,surf_e_111]
surfaces = [(1,0,0),(1,1,1)]

no = 20
z_no = 8
hexagon = fcc111(symbol=symbol, size=(no,no,z_no), a=None, vacuum=None, orthogonal=False, periodic=False)
hexagon.rotate(-30, 'z', center=hexagon[0].position, rotate_cell=False)
center_position = 50
all_x_positions = tuple((atom.x,atom.index) for atom in hexagon)
max_x, max_x_index = max(all_x_positions)
min_x, min_x_index = min(all_x_positions)
diameter = max_x - min_x
hexagon.set_cell((diameter,diameter,diameter))
hexagon.center()
for index in range(len(hexagon)-1,-1,-1):
	if not (diameter*0.245 <= hexagon[index].x <= diameter*0.745):
		del hexagon[index]

write('hexagon.xyz',hexagon)