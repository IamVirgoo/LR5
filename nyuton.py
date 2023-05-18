from matplotlib import pyplot as plt
from sympy import *
from copy import deepcopy

import devtools
from devtools import *
import numpy as np


def nyuton_func(N, arg):
    x = Symbol('x')
    return N.subs(x, arg)


def nyuton(table):
    x = Symbol('x')

    n = len(table[0])
    table_n = deepcopy(table)
    table_n = np.array(table_n)

    for i in range(1, n):
        table_n = np.append(table_n, np.array([0 for i in range(n)]))
    table_n = table_n.reshape(n + 1, n)

    for i in range(2, n + 1):
        for j in range(0, n - i + 1):
            table_n[i][j] = table_n[i - 1][j + 1] - table_n[i - 1][j]
            table_n[i][j] /= table_n[0][j + i - 1] - table_n[0][j]

    N = 0
    x_step = 1

    for i in range(1, len(table_n)):
        N += table_n[i][0] * x_step
        x_step *= (x - table_n[0][i - 1])

    x_plot = np.arange(table[0][0], table[0][-1], 0.1)
    y_plot = [nyuton_func(N, i) for i in x_plot]

    print(table_n)

    plt.plot(x_plot, y_plot, label='Многочлен Ньютона')
    plt.scatter(table[0], table[1], label='Точки графика')
    plt.show()


nyuton(devtools.table)
