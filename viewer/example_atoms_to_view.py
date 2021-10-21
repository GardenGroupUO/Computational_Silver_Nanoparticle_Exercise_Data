from ase.cluster import wulff_construction
from ase.data.colors import jmol_colors
from ase.data import atomic_numbers

symbol = 'Ag'
surfaces = [(1,0,0),(1,1,1)]
esurf = [1.0,0.9]
size = 200

atoms = wulff_construction(symbol,surfaces,esurf,size,'fcc',rounding='above')
atoms.rotate(20,'y')
atoms.set_cell([0,0,0])

colours = {}
transparencies = {}

atoms_100 = [124]
atoms_100b = [index for index in range(len(atoms)) if ((atoms.get_distance(index,atoms_100[0]) < 3.0) and not (index == atoms_100[0]))]
for n in atoms_100:
    colours[n] = jmol_colors[atomic_numbers['Mn']] 
for n in atoms_100b:
    colours[n] = jmol_colors[atomic_numbers['Fe']]
for index in range(len(atoms)):
    if index not in atoms_100+atoms_100b:
        transparencies[index] = 0.9

from ase.build import molecule

molecule_name = 'CH3CH2OCH3' #'H2O','NaCl','CH3CH2OCH3','CH3SiH3'
#atoms = molecule(molecule_name)


from ase.build import fcc100, fcc111
from ase.calculators.emt import EMT
from ase.optimize import FIRE
def get_square_surface_model(element,lattice_constant):
    x_number_of_atoms = 6
    y_number_of_atoms = 6
    slab = fcc100(element,size=(x_number_of_atoms,y_number_of_atoms,4),a=lattice_constant)
    slab.center(vacuum=10.0,axis=2,about=None)
    slab.set_calculator(EMT())
    dyn = FIRE(slab,logfile=None)
    dyn.run(fmax=0.01)
    return slab

def get_triangle_surface_model(element,lattice_constant):
    x_number_of_atoms = 6
    y_number_of_atoms = 6
    slab = fcc111(element,size=(x_number_of_atoms,y_number_of_atoms,4),a=lattice_constant)
    slab.center(vacuum=10.0,axis=2,about=None)
    slab.set_calculator(EMT())
    dyn = FIRE(slab,logfile=None)
    dyn.run(fmax=0.01)
    return slab

atoms = get_triangle_surface_model('Ag',3.7)


from x3d_viewer import view_x3d
view_x3d(atoms,'testing',show_unit_cell=True) #,colours=colours,transparencies=transparencies)

from ase.visualize import view
#view(atoms)
#import pdb; pdb.set_trace()