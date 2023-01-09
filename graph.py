import os
from pathlib import Path

from sage.misc.persist import save
from sage.plot.plot import plot

from symbolic_matrix import *
from precision import *

root_folder = Path(__file__).parent
image_folder = root_folder / 'images'

os.makedirs(image_folder, exist_ok=True)


def save_and_show_graph(graph, filename):
    save(graph, str(image_folder / filename))
    os.system(f'display "{image_folder / filename}" &')
