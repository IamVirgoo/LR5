from matplotlib import pyplot as plt
import numpy as np

import devtools
from devtools import *
from sympy import *

x = Symbol('x')


def lagrange_basis(j, table):
    number = len(table[0])
    ch = 1
    zn = 1
    for i in range(number):
        if i == j:
            continue
        ch *= (x - table[0][i])
        zn *= (table[0][j] - table[0][i])
    return ch / zn


def lagrange_polynomial(table):
    number = len(table[0])
    alph = 0
    for i in range(number):
        alph += lagrange_basis(i, table) * table[1][i]
    return alph


def func(i, table):
    alph = lagrange_polynomial(table)
    x = Symbol('x')
    return alph.subs(x, i)


if __name__ == '__main__':
    grid = devtools.table

    x_1 = grid[0]
    y_1 = grid[1]

    x_plot = np.arange(grid[0][0], grid[0][-1], 0.1)
    y_plot = [func(i, grid) for i in x_plot]

    print("Многочлен Лагранжа: ", lagrange_polynomial(grid))

    plt.plot(x_plot, y_plot, label='Многочлен Лагранжа')
    plt.scatter(x_1, y_1, label='Точки графика')
    plt.show()
