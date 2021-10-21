import os

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from IPython.display import HTML

if 'Computational_Silver_Nanoparticle_Exercise_Data' in os.listdir('.'):
    from Computational_Silver_Nanoparticle_Exercise_Data.movie_viewer.x3d_HTML_movie_writer import write_x3d_movie_html
else:
    try:
        from movie_viewer.x3d_HTML_movie_writer import write_x3d_movie_html
    except Exception:
        from x3d_HTML_movie_writer import write_x3d_movie_html

def view_x3d_movie(images,notebook_name=None,colours=[]):
    """View images inline in a jupyter notbook. This command
    should only be used within a jupyter/ipython notebook.
    """
    
    output = StringIO()
    write_x3d_movie_html(images, output, notebook_name=notebook_name, colours=colours)
    data = output.getvalue()
    
    # The following will make the physical html page that can be used for debugging the chemistry viewer
    with open(os.path.realpath(__file__).replace('/x3d_movie_viewer.py','')+'/'+'html_file.html','w') as write_to:
        write_to.write(data) 
    
    output.close()
    return HTML(data)