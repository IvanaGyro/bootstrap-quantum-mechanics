from sage.matrix.matrix_generic_dense import Matrix_generic_dense
from sage.matrix.matrix_symbolic_dense import Matrix_symbolic_dense
from sage.rings.infinity import infinity
from sage.rings.real_mpfr import RealNumber
from sage.symbolic.expression import Expression

from precision import *

_zero = Real300(0)
_infinity = Real300(infinity)


def get_matrix_variable(mat: Matrix_symbolic_dense) -> Expression:
    for i in range(mat.nrows()):
        for j in range(mat.ncols()):
            vars = mat[i][j].variables()
            if len(vars):
                return vars[0]
    raise ValueError(f'No variable in the matrix. {mat}')


def substitute_by_real(mat: Matrix_symbolic_dense,
                       value: RealNumber) -> Matrix_generic_dense:
    x = get_matrix_variable(mat)
    return mat.substitute({x: value}).apply_map(Real300)


def evaluate_det(mat: Matrix_symbolic_dense, value: RealNumber) -> RealNumber:
    return substitute_by_real(mat, value).det()


def evaluate_derivative_det(mat: Matrix_symbolic_dense,
                            derivative_mat: Matrix_symbolic_dense,
                            value: RealNumber):
    # The determinant value is cached by Sagemath if the matrix isn't changed.
    return mat.det() * (mat.inverse() *
                        substitute_by_real(derivative_mat, value)).trace()


def newton_method(mat: Matrix_symbolic_dense,
                  derivative_mat: Matrix_symbolic_dense,
                  init_guess: RealNumber,
                  eps: RealNumber,
                  epsx: RealNumber,
                  maxiter: int = 100):
    current_guess = init_guess
    iterates = [current_guess]
    real_matrix = substitute_by_real(mat, current_guess)
    fc = real_matrix.det()
    first_evaluation = fc
    for _ in range(maxiter - 1):
        # if abs(fc) < eps:
        #     break
        # Comparisons between RealField numbers of the same precision are 6
        # times faster than comparisons between RealField numbers and Python
        # integers.
        if fc == _zero:
            break
        diff = fc / evaluate_derivative_det(real_matrix, derivative_mat,
                                            current_guess)
        if abs(diff) == _infinity:
            raise ValueError('get infinity')
        current_guess -= diff
        iterates.append(current_guess)
        if len(iterates) >= 5 and max(iterates[-5:]) - min(
                iterates[-5:]) < epsx:
            break
        real_matrix = substitute_by_real(mat, current_guess)
        fc = real_matrix.det()
    return current_guess, first_evaluation, iterates
