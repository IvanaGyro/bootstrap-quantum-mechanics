from functools import cache
from sage.all_cmdline import *
import operator
import time
from pprint import pprint
from pathlib import Path
import os

import sage.libs.ecl

sage.libs.ecl.ecl_eval("(ext:set-limit 'ext:heap-size 0)")

time_table = {}


def begin(name):
    time_table[name] = time.time()


def end(name):
    print(f'{name}: {time.time() - time_table[name]}')
    del time_table[name]


l = 1
energy = var('E', domain='real')
Real300 = RealField(10000)
Real20 = RealField(20)


@cache
def expected_distance(n, l):
    if n <= -2:
        return 0
    if n == -1:
        return -2 * energy
    if n == 0:
        return 1
    return ((
        - 4 * (2 * (n + 1) - 1) * expected_distance(n - 1, l) \
        - n * ((n + 1) * (n - 1) - 4 * l * (l + 1)) * expected_distance(n - 2, l)
        ) / (8 * (n + 1) * energy)).full_simplify()


def intersections(a, b):
    '''
    https://stackoverflow.com/a/40368603/6663588
    '''
    ranges = []
    i = j = 0
    while i < len(a) and j < len(b):
        a_left, a_right = a[i]
        b_left, b_right = b[j]

        if a_right < b_right:
            i += 1
        else:
            j += 1

        if a_right >= b_left and b_right >= a_left:
            end_pts = sorted([a_left, a_right, b_left, b_right])
            middle = [end_pts[1], end_pts[2]]
            ranges.append(middle)

    ri = 0
    while ri < len(ranges) - 1:
        if ranges[ri][1] == ranges[ri + 1][0]:
            ranges[ri:ri + 2] = [[ranges[ri][0], ranges[ri + 1][1]]]

        ri += 1

    return ranges


def get_matrix_variable(mat):
    for i in range(mat.nrows()):
        for j in range(mat.ncols()):
            vars = mat[i][j].variables()
            if len(vars):
                return vars[0]
    raise ValueError(f'No variable in the matrix. {mat}')


def substitute_by_real(mat, value):
    x = get_matrix_variable(mat)
    return mat.substitute({x: value}).apply_map(Real300)


def evaluate_det(mat, value):
    return substitute_by_real(mat, value).det()


root_folder = Path(__file__).parent
image_folder = root_folder / 'images'
dump_foler = root_folder / 'dumps'
os.makedirs(image_folder, exist_ok=True)
os.makedirs(dump_foler, exist_ok=True)

def save_and_show_graph(graph, filename):
    save(graph, str(image_folder / filename))
    os.system(f'display "{image_folder / filename}" &')


def newton_method(input_matrix,
                  derivative_matrix,
                  init_guess,
                  eps,
                  epsx,
                  maxiter=100):

    def evaluate_derivative_det(real_matrix, real_matrix_det, value):
        return real_matrix_det * (real_matrix.inverse() * substitute_by_real(
            derivative_matrix, value)).trace()

    current_guess = init_guess
    iterates = [current_guess]
    real_matrix = substitute_by_real(input_matrix, current_guess)
    fc = real_matrix.det()
    first_evaluation = fc
    for _ in range(maxiter - 1):
        # if abs(fc) < eps:
        #     break
        if fc == 0:
            break
        diff = fc / evaluate_derivative_det(real_matrix, fc, current_guess)
        if abs(diff) == infinity:
            raise ValueError('get infinity')
        current_guess -= diff
        iterates.append(current_guess)
        if len(iterates) >= 5 and max(iterates[-5:]) - min(
                iterates[-5:]) < epsx:
            break
        real_matrix = substitute_by_real(input_matrix, current_guess)
        fc = real_matrix.det()
    return current_guess, first_evaluation, iterates


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


def update_ranges_by_probing(hankel, current_ranges, accepted_range_size):
    depth = hankel.nrows()
    begin('weight_hankel')
    weight_hankel = matrix(depth, depth, [
        hankel[i][j] / hankel[0][j] / hankel[i][0]
        for i in range(depth)
        for j in range(depth)
    ])
    end('weight_hankel')

    interested_areas = {
        14: [[-0.12952, -0.12501], [-0.035683, -0.027594],
             [-0.024358, -0.0033258]],
    }
    if depth in interested_areas:
        for x_range in interested_areas[depth]:
            show_interested_area(depth, weight_hankel, x_range)

    begin('derivative_hankel')
    derivative_hankel = weight_hankel.apply_map(lambda f: f.derivative(energy))
    end('derivative_hankel')

    new_ranges = []
    for lower, upper in current_ranges:
        lower, upper = Real300(lower), Real300(upper)

        if lower == upper:
            if evaluate_det(hankel, lower) >= 0:
                # It may be a singular point, so need to evaluate on a original hankel.
                new_ranges.append([lower, lower])
            continue

        if upper - lower < accepted_range_size:
            new_ranges.append([lower, upper])
            continue

        denominator = 10
        step = (upper - lower) / denominator
        operators = [None, operator.lt, operator.ge]  # index: [_, 1, -1]
        low = None
        prev_x, current_x = None, lower
        has_zero_end = False
        sign = evaluate_det(weight_hankel, current_x)

        evaluated_value = evaluate_det(weight_hankel, current_x)
        while evaluated_value == 0 and current_x < lower + step:
            prev_x, current_x = current_x, current_x + step / denominator
            evaluated_value = evaluate_det(weight_hankel, current_x)
            # assume evaluated_value will not be zero again after shifting
            if evaluated_value < 0:
                new_ranges.append([lower, prev_x])
        if evaluated_value > 0:
            low = lower
        if evaluated_value == 0:
            raise ValueError('evaluated_value is still zero after shifting')
        sign = 1 if evaluated_value >= 0 else -1

        for i in range(1, denominator + 1):
            prev_x, current_x = current_x, lower + step * i
            evaluated_value = evaluate_det(weight_hankel, current_x)
            if evaluated_value == 0:
                if i == denominator:
                    has_zero_end = True
                current_x -= step / denominator
                evaluated_value = evaluate_det(weight_hankel, current_x)
            if operators[sign](evaluated_value, 0):
                # flip sign
                # The middle point may locate at the peak, so we try
                candidates = [(prev_x + current_x) / 2, prev_x, current_x]
                for guess in candidates:
                    begin('newton method')
                    root, first_evaluation, iterates = newton_method(
                        weight_hankel, derivative_hankel, guess, 10 ^ -50,
                        accepted_range_size)
                    end('newton method')
                    print(f'first_evaluation:{Real20(first_evaluation)}')
                    print(f'iterates: {len(iterates)}')
                    if len(iterates) > 20:
                        pprint([Real20(i) for i in iterates])
                    if root >= prev_x and root <= current_x:
                        break
                else:
                    print(
                        f'root out of range. root:{Real20(root)} range:{[Real20(prev_x), Real20(current_x)]}'
                    )
                    show_interested_area(depth, weight_hankel,
                                         [prev_x, current_x])
                    # pick middle temporary
                    root = (prev_x + current_x) / 2
                if sign == 1:
                    if low is None:
                        raise ValueError(
                            f'None low 1 i:{i} lower:{Real20(lower)} upper:{Real20(upper)} current_ranges:{[[Real20(low), Real20(high)] for low, high in current_ranges]}'
                        )
                    new_ranges.append([low, root])
                    low = None
                else:
                    low = root
                sign *= -1
        if low is not None:
            new_ranges.append([low, upper])
        elif has_zero_end:
            new_ranges.append([upper, upper])
    return new_ranges


function_graph = Graphics()
positive_range_graph = Graphics()
intersection_range_graphs = [Graphics()]


def update_ranges_by_solving_det(hankel, current_ranges, hue_value):
    global function_graph
    global positive_range_graph
    depth = hankel.nrows()
    begin('calculate det')
    f = hankel.det()
    end('calculate det')
    # print(f)

    begin('simplify')
    f = f.simplify_full()
    # f = expand(f)
    end('simplify')
    print(f)

    begin('solve inequation')
    positive_solutions = solve(f >= 0, energy)
    end('solve inequation')

    begin('find intersections')
    positive_ranges = []
    for solution in positive_solutions:
        if len(solution) > 2:
            print(f'more than two expression')
        if len(solution) == 1:
            if solution[0].operator() in [operator.le, operator.lt]:
                positive_ranges.append((xmin, solution[0].rhs()))
            elif solution[0].operator() in [operator.ge, operator.gt]:
                positive_ranges.append((solution[0].rhs(), xmax))
            elif solution[0].operator() == operator.eq:
                positive_ranges.append((solution[0].rhs(), solution[0].rhs()))
            else:
                print(f'Cannot recognize operator: {solution[0]}')
        elif len(solution) == 2:
            positive_ranges.append((solution[0].rhs(), solution[1].rhs()))
    # print(positive_ranges)

    # I don't know why there are complex numbers in the results.
    real_positive_ranges = []
    for positive_range in positive_ranges:
        try:
            real_positive_ranges.append([Real300(x) for x in positive_range])
        except TypeError as e:
            print(e)
    positive_ranges = real_positive_ranges
    # pprint([[Real20(low), Real20(high)] for low, high in positive_ranges])
    end('find intersections')

    function_graph += plot(f, ymax=5, ymin=-5, hue=hue_value, detect_poles=True)
    for positive_range in positive_ranges:
        dots = [(x, depth) for x in positive_range]
        positive_range_graph += line(dots, hue=hue_value)
        positive_range_graph += points(dots, color='red')
    return intersections(positive_ranges, current_ranges)


def save_current_ranges(depth, current_ranges):
    filename = f'current_ranges-l_{l}-depth_{depth}.sobj'
    save(current_ranges,  str(dump_foler / filename))


def load_current_ranges(depth):
    filename = f'current_ranges-l_{l}-depth_{depth}.sobj'
    try:
        return load(str(dump_foler / filename))
    except FileNotFoundError:
        return None


min_depth = 2
max_depth = 5
xmin = -0.7
xmax = 0.5
current_ranges = [[xmin, xmax]]
offset = 0
accepted_range_size = 10 ^ -5

for depth in range(min_depth, max_depth + 1):
    print(f'depth:{depth}')
    begin(f'depth:{depth}')

    hue_value = (depth - min_depth) / (max_depth - min_depth + 1)

    stored_current_ranges = load_current_ranges(depth)
    if stored_current_ranges is not None:
        current_ranges = stored_current_ranges
    else:
        begin('generate matrix')
        hankel = matrix.hankel(
            [expected_distance(i + offset, l) for i in range(depth)], [
                expected_distance(i + offset + depth - 1, l)
                for i in range(1, depth)
            ], SR)
        end('generate matrix')

        if depth >= 8:
            current_ranges = update_ranges_by_probing(hankel, current_ranges,
                                                      accepted_range_size)
        else:
            current_ranges = update_ranges_by_solving_det(
                hankel, current_ranges, hue_value)
        pprint([[Real20(low), Real20(high)] for low, high in current_ranges])
        save_current_ranges(depth, current_ranges)

    if depth % 10 == 1:
        intersection_range_graphs.append(Graphics())
    for current_range in current_ranges:
        dots = [(x, depth) for x in current_range]
        intersection_range_graphs[0] += line(dots, hue=hue_value)
        intersection_range_graphs[0] += points(dots, hue=hue_value)
        low, high = current_range
        if high - low <= 10 ^ -4:
            # Do not show the singular point in graph for high depths.
            continue
        for i, _ in enumerate(intersection_range_graphs[1:], 1):
            intersection_range_graphs[i] += line(dots, hue=hue_value)
            intersection_range_graphs[i] += points(dots, hue=hue_value)

    end(f'depth:{depth}')
    print()


filename_prefix = f'depth_{min_depth}_{max_depth}-l_{l}'
save_and_show_graph(function_graph, f'{filename_prefix}.png')
save_and_show_graph(positive_range_graph, f'{filename_prefix}-positive.png')
for i, graph in enumerate(intersection_range_graphs):
    save_and_show_graph(graph, f'{filename_prefix}-intersection_{i}.png')
