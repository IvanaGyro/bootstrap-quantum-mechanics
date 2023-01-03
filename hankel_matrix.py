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
