import os
from pathlib import Path

from sage.misc.persist import save
from sage.misc.persist import load

root_folder = Path(__file__).parent
dump_foler = root_folder / 'dumps'
os.makedirs(dump_foler, exist_ok=True)


def save_obj(obj, filename):
    save(obj, str(dump_foler / filename))


def save_current_ranges(l, depth, current_ranges):
    filename = f'current_ranges-l_{l}-depth_{depth}.sobj'
    save_obj(current_ranges, filename)


def load_obj(filename):
    try:
        return load(str(dump_foler / filename))
    except FileNotFoundError:
        return None


def load_current_ranges(l, depth):
    filename = f'current_ranges-l_{l}-depth_{depth}.sobj'
    return load_obj(filename)
