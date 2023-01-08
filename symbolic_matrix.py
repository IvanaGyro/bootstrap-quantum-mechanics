from sage.matrix.matrix_generic_dense import Matrix_generic_dense
from sage.rings.infinity import infinity
from sage.rings.real_mpfr import RealNumber

from precision import *

_zero = Real300(0)
_infinity = Real300(infinity)


def substitute_by_real(mat: Matrix_generic_dense,
                       value: RealNumber) -> Matrix_generic_dense:
    x = mat.base_ring().gen()
    return mat.substitute({x: value})


def evaluate_det(mat: Matrix_generic_dense, value: RealNumber) -> RealNumber:
    return substitute_by_real(mat, value).det()


def evaluate_derivative_det(mat: Matrix_generic_dense,
                            derivative_mat: Matrix_generic_dense,
                            value: RealNumber):
    # The determinant value is cached by Sagemath if the matrix isn't changed.
    return mat.det() * (mat.inverse() *
                        substitute_by_real(derivative_mat, value)).trace()


def newton_method(mat: Matrix_generic_dense,
                  derivative_mat: Matrix_generic_dense,
                  init_guess: RealNumber,
                  epsx: RealNumber,
                  x_range: tuple[RealNumber, RealNumber],
                  maxiter: int = 100):
    x_min, x_max = x_range
    current_guess = init_guess
    iterates = [current_guess]
    real_matrix = substitute_by_real(mat, current_guess)
    fc = real_matrix.det()
    first_evaluation = fc
    for _ in range(maxiter - 1):
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
        if len(iterates) >= 5:
            try_max = max(iterates[-5:])
            try_min = min(iterates[-5:])
            if try_max < x_min or try_min > x_max:
                current_guess = None
                break
            if try_max - try_min < epsx:
                break
        real_matrix = substitute_by_real(mat, current_guess)
        fc = real_matrix.det()
    return current_guess, first_evaluation, iterates
