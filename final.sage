import operator
from functools import cache
from pprint import pprint

import sage.libs.ecl
from sage.all_cmdline import *

from benchmark import *
from graph import *
from symbolic_matrix import *
from precalculation import *
from precision import *

sage.libs.ecl.ecl_eval("(ext:set-limit 'ext:heap-size 0)")

l = 1
poly_int = QQ['energy']
(energy,) = poly_int._first_ngens(1)


@cache
def expected_distance(n, l):
    if n <= -2:
        return 0
    if n == 0:
        return 1
    if n == -1:
        return -2 * energy

    # `.reduce()` doesn't change anything, so it's no need to call it.
    return (
        - 4 * (2 * (n + 1) - 1) * expected_distance(n - 1, l) \
        - n * ((n + 1) * (n - 1) - 4 * l * (l + 1)) * expected_distance(n - 2, l)
        ) / (8 * (n + 1) * energy)


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


def update_ranges_by_probing(hankel, current_ranges, accepted_range_size):
    depth = hankel.nrows()

    # There is no difference after calling `.reduce()` on each element of
    # `weight_hankel` and `derivative_hankel`, so applying `.reduce()` is
    # useless
    begin('weight_hankel')
    weight_hankel = matrix(depth, depth,
                           [(hankel[i][j] / hankel[0][j] / hankel[i][0])
                            for i in range(depth)
                            for j in range(depth)])
    end('weight_hankel')

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
                        weight_hankel, derivative_hankel, guess,
                        accepted_range_size, (prev_x, current_x))
                    end('newton method')
                    print(f'first_evaluation:{Real20(first_evaluation)}')
                    print(f'iterates: {len(iterates)}')
                    if len(iterates) > 20:
                        pprint([Real20(i) for i in iterates])
                    if root is not None and root >= prev_x and root <= current_x:
                        break
                else:
                    print(
                        f'root out of range. root:{Real20(root)} range:{[Real20(prev_x), Real20(current_x)]}'
                    )
                    plot_determinant(weight_hankel, [prev_x, current_x])
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


def update_ranges_by_solving_det(hankel, current_ranges, search_range,
                                 hue_value):
    global function_graph
    global positive_range_graph
    xmin, xmax = search_range
    depth = hankel.nrows()
    begin('calculate det')
    # `algorithm='hessenberg'` doesn't improve efficiency in this case.
    f = hankel.det()
    end('calculate det')
    # print(f)

    begin('solve inequation')
    f = SR(f)
    positive_solutions = solve(f >= 0, SR.var(energy))
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


def find_possible_energy(depth_range, search_range):
    min_depth, max_depth = depth_range
    current_ranges = [list(search_range)]
    accepted_range_size = 10 ^ -5

    for depth in range(min_depth, max_depth + 1):
        print(f'depth:{depth}')
        begin(f'depth:{depth}')

        hue_value = (depth - min_depth) / (max_depth - min_depth + 1)

        stored_current_ranges = load_current_ranges(l, depth)
        if stored_current_ranges is not None:
            current_ranges = stored_current_ranges
        else:
            begin('generate functions')
            functions = [expected_distance(i, l) for i in range(depth * 2 - 1)]
            end('generate functions')
            # Calculating determinant is faster in the fraction field of
            # `QQ[]` than in `SR` when the variables are not substituted by
            # the numbers.
            begin('generate matrix')
            hankel = matrix.hankel(functions[:depth],
                                   functions[depth:depth * 2 - 1],
                                   poly_int.fraction_field())
            end('generate matrix')

            # It costs about 17 seconds for depth 10 to calculate the
            # determinant without replacing the variable with numbers.
            if depth >= 10:
                current_ranges = update_ranges_by_probing(
                    hankel, current_ranges, accepted_range_size)
            else:
                current_ranges = update_ranges_by_solving_det(
                    hankel, current_ranges, search_range, hue_value)
            pprint([[Real20(low), Real20(high)] for low, high in current_ranges
                   ])
            save_current_ranges(l, depth, current_ranges)

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


def plot_determinant(mat,
                     x_range,
                     title=None,
                     filename=None,
                     power_scale=Real300(0.001),
                     plot_points=20,
                     mode='native'):
    '''Plot the determinant value of the symbolic matrix.

    This method is useful for picking the enough precision. If the curve
    fluctuates like a noise, a higher precision should be applied.
    '''

    def valuate(x):
        det_value = evaluate_det(mat, Real300(x))
        if power_scale != 1:
            det_value = abs(det_value)**power_scale * sign(det_value)
        return det_value

    x_min, x_max = [Real300(x) for x in x_range]
    if title is None:
        title = f'l_{l}-depth_{mat.nrows()}-{[Real20(x) for x in x_range]}'
    if filename is None:
        filename = f'{title}.png'

    if mode == 'native':
        graph = plot(
            valuate,
            (x_min, x_max),
            plot_points=plot_points,
            title=title,
        )
    elif mode == 'interval':
        step = (x_max - x_min) / (plot_points - 1)
        x_list = [x_min + step * i for i in range(plot_points)]
        points = [(x, valuate(x)) for x in x_list]
        graph = list_plot(
            points,
            title=title,
            plotjoined=True,
        )
    else:
        raise ValueError(f'mode:{mode} is not supported.')
    save_and_show_graph(graph, filename)


if __name__ == '__main__':
    find_possible_energy((2, 5), (-0.7, 0.5))

    interested_areas = {
        # 14: [[-0.12952, -0.12501], [-0.035683, -0.027594],
        #      [-0.024358, -0.0033258]],
        # 60: [(-0.0070218, -1.4842e-7)],
        # 80: [(-0.0070218, -1.4842e-7)],
    }

    for depth, x_ranges in interested_areas.items():
        functions = [expected_distance(i, l) for i in range(depth * 2 - 1)]

        hankel = matrix.hankel(functions[:depth],
                               functions[depth:depth * 2 - 1],
                               poly_int.fraction_field())

        begin('weight_hankel')
        weight_hankel = matrix(depth, depth,
                               [(hankel[i][j] / hankel[0][j] / hankel[i][0])
                                for i in range(depth)
                                for j in range(depth)])
        end('weight_hankel')

        for x_range in x_ranges:
            display_x_range = [Real20(x) for x in x_range]
            title = f'weighted-l_{l}-depth_{depth}-{display_x_range}'
            plot_determinant(weight_hankel, x_range, title=title)
