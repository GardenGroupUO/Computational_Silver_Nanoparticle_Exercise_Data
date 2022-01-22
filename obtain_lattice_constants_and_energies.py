from surface_energy_of_square_and_triangle import get_lattice_contant_and_energy_of_bulk_model_from_lattice_curve, make_LatticeCurve_plot

elements = ['Al', 'Cu', 'Ag', 'Au', 'Ni', 'Pd', 'Pt']
EMT_data = {}
for element in elements:
    E_BM, lattice_constant = get_lattice_contant_and_energy_of_bulk_model_from_lattice_curve(element,lc_low=3.0,lc_high=5.0)
    EMT_data[element] = (E_BM, lattice_constant)
    print(str(element)+': '+str((E_BM, lattice_constant)))
    make_LatticeCurve_plot(element)

print(EMT_data)
with open('EMT_data.txt','w') as EMT_dataTXT:
    EMT_dataTXT.write(str(EMT_data))