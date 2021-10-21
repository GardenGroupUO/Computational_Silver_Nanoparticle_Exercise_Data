from ase.cluster import wulff_construction

size = 21127

surfaces = [(1,0,0),(1,1,1)]
surf_e_100 = 0.9
surf_e_111 = 1.0

esurf = [surf_e_100,surf_e_111]
atoms = wulff_construction('Ag',surfaces,esurf,size,'fcc',rounding='above')
atoms.cell=[0,0,0]
print('This nanoparticle contains '+str(len(atoms))+' atoms.')

index1 = 20937
index2 = 432

distance = atoms.get_distance(index1,index2)
print('The diameter of this nanoparticle is '+str(round(distance,2))+' Angstroms')