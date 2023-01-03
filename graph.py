import os
from pathlib import Path

from sage.misc.persist import save
from sage.plot.plot import plot

from hankel_matrix import *
from precision import *

root_folder = Path(__file__).parent
image_folder = root_folder / 'images'

os.makedirs(image_folder, exist_ok=True)


def save_and_show_graph(graph, filename):
    save(graph, str(image_folder / filename))
    os.system(f'display "{image_folder / filename}" &')


def show_interested_area(depth, weight_hankel, x_range):
    truncked_x_range = [Real20(x) for x in x_range]
    g = plot(
        lambda v: evaluate_det(weight_hankel, v),
        x_range,
        plot_points=200,
        title=f'depth:{depth} {truncked_x_range}',
    )
    save_and_show_graph(
        g, f'depth_{depth}-{truncked_x_range[0]}_{truncked_x_range[1]}.png')
