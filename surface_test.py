from surface_energy_of_square_and_triangle import get_energy_of_bulk_model, get_square_surface_model
from ase.visualize import view
from ase.io import write

element = 'Ag'
E_BM, lattice_constant = get_energy_of_bulk_model(element)
slab = get_square_surface_model(element,lattice_constant)

view(slab)

write('square_slab.xyz',slab)