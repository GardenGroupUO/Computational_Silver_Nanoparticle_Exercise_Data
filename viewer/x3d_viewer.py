import os

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from IPython.display import HTML

if 'Computational_Silver_Nanoparticle_Exercise_Data' in os.listdir('.'):
    from Computational_Silver_Nanoparticle_Exercise_Data.viewer.x3d_HTML_writer import write_x3d_html
else:
    try:
        from viewer.x3d_HTML_writer import write_x3d_html
    except Exception:
        from x3d_HTML_writer import write_x3d_html

def view_x3d(atoms,notebook_name=None,colours={},transparencies={},show_unit_cell=False):
    """View atoms inline in a jupyter notbook. This command
    should only be used within a jupyter/ipython notebook.
    
    Args:
        atoms - ase.Atoms, atoms to be rendered"""
    
    output = StringIO()
    write_x3d_html(atoms, output, notebook_name=notebook_name, colours=colours, transparencies=transparencies, show_unit_cell=show_unit_cell)
    data = output.getvalue()
    
    # The following will make the physical html page that can be used for debugging the chemistry viewer
    with open(os.path.realpath(__file__).replace('/x3d_viewer.py','')+'/'+'html_file.html','w') as write_to:
        write_to.write(data) 
    
    output.close()
    return HTML(data)