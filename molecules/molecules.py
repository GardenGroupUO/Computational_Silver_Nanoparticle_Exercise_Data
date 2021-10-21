from ase.build import molecule
from ase.io import read, write

molecules = ['NaCl','CH3CH2OCH3','CH3SiH3']

for a_molecule in molecules:
  atoms = molecule(a_molecule)
  write(a_molecule+'.xyz',atoms)