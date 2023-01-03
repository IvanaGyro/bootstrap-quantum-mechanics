import time

_time_table = {}


def begin(name):
    _time_table[name] = time.time()


def end(name):
    print(f'{name}: {time.time() - _time_table[name]}')
    del _time_table[name]
