import numpy as np
from tqdm import tqdm

from ase.build import fcc100, fcc111
from ase.lattice.cubic import FaceCenteredCubic

from ase.calculators.emt import EMT
from ase.optimize import FIRE

def get_energy_of_bulk_model(element):
    energies = []
    lattice_params = np.arange(3.2,4.5,0.05)
    print('Getting Bulk Energy. This may take a minute.')
    print('Scanning across lattice constants')
    pbar = tqdm(lattice_params,unit='LC')
    for lc in pbar:
        pbar.set_description("Processing Lattice Constant (LC) = %.2f A" % lc)
        bulk = FaceCenteredCubic(size=(6,6,6),symbol=element,pbc=(1,1,1),latticeconstant=lc)
        bulk.set_calculator(EMT())
        e_atoms = bulk.get_potential_energy()
        energies.append(e_atoms)
    print()
    lowest_energy, lattice_constant = min(zip(energies,lattice_params))
    E_BM = lowest_energy/len(bulk)
    return E_BM, lattice_constant

def get_square_surface_model(element,lattice_constant):
    x_number_of_atoms = 6
    y_number_of_atoms = 6
    z_number_of_atoms = 4
    slab = fcc100(element,size=(x_number_of_atoms,y_number_of_atoms,z_number_of_atoms),a=lattice_constant)
    slab.center(vacuum=10.0,axis=2,about=None)
    slab.set_calculator(EMT())
    dyn = FIRE(slab,logfile=None)
    dyn.run(fmax=0.01)
    return slab, x_number_of_atoms, y_number_of_atoms, z_number_of_atoms

def get_energy_of_square_surface_model(element,lattice_constant):
    slab, x_number_of_atoms, y_number_of_atoms, z_number_of_atoms = get_square_surface_model(element,lattice_constant)
    E_SM = slab.get_potential_energy()
    return E_SM, x_number_of_atoms*y_number_of_atoms*2, x_number_of_atoms*y_number_of_atoms*z_number_of_atoms, slab

def get_triangle_surface_model(element,lattice_constant):
    x_number_of_atoms = 6
    y_number_of_atoms = 6
    z_number_of_atoms = 4
    slab = fcc111(element,size=(x_number_of_atoms,y_number_of_atoms,4),a=lattice_constant)
    slab.center(vacuum=10.0,axis=2,about=None)
    slab.set_calculator(EMT())
    dyn = FIRE(slab,logfile=None)
    dyn.run(fmax=0.01)
    return slab, x_number_of_atoms, y_number_of_atoms, z_number_of_atoms

def get_energy_of_triangle_surface_model(element,lattice_constant):
    slab, x_number_of_atoms, y_number_of_atoms, z_number_of_atoms = get_triangle_surface_model(element,lattice_constant)
    E_SM = slab.get_potential_energy()
    return E_SM, x_number_of_atoms*y_number_of_atoms*2, x_number_of_atoms*y_number_of_atoms*z_number_of_atoms, slab

# -------------------------------------------------------------------------

def eV_per_atom_into_J_per_m_squared(eV_per_atom,slab,surface_atoms):

    eV_to_J = 1.60218e-19

    #Calculating surface area of non-orthogonal surface unit cell

    cell_vecs=slab.get_cell()
    vec1=cell_vecs[0]
    vec2=cell_vecs[1]

    b=vec1[0]
    h=vec2[1]

    surface_area = (b * h /surface_atoms) * 1e-20 ## in m^2

    J_per_m_squared = eV_per_atom * (eV_to_J / surface_area)  ### J m^-2
    return J_per_m_squared