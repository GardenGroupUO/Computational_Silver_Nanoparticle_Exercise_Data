from ase.io import read
systems = read('mid3_animation.traj',index=':')

from x3d_movie_viewer import view_x3d_movie
view_x3d_movie(systems,'testing')