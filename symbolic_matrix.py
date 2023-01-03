from precision import *


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
