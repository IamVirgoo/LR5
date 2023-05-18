from matplotlib import pyplot as plt
from sympy import *
import numpy as np
import devtools


def function(arg, a, b):
    x = Symbol('x')
    return ((ln(x)) ** (a / b) * sin(x)).subs(x, arg)


def N_func(N, arg):
    x = Symbol('x')
    return N.subs(x, arg)


def splain(n, m, c, d, a_, b_):
    m *= n

    h = (d - c) / (m)
    ind = [i for i in np.arange(c, d + h / 2, h)]
    table = [[0 for i in range(m + 1)] for i in range(2)]
    table[0] = ind
    for i in range(m + 1):
        table[1][i] = function(table[0][i], a_, b_)

    a = [table[1][i] for i in range(m)]
    b = [0 for i in range(m)]
    for i in range(m - 1):
        b[i + 1] = 2 * (table[1][i + 1] - table[1][i]) / h - b[i]
    C = [(b[i + 1] - b[i]) / (2 * h) for i in range(m - 1)]
    C.append((table[1][m] - table[1][m - 1]) / h ** 2 - b[m - 1] / h)

    for i in range(m):
        x = Symbol('x')
        N = a[i] + b[i] * (x - table[0][i]) + C[i] * (x - table[0][i]) ** 2
        x_plot = np.arange(c + i * h, c + i * h + h, 0.01)
        y_plot = [N_func(N, i) for i in x_plot]
        plt.plot(x_plot, y_plot)

    x_plot = np.arange(c, d, 0.01)
    y_plot = [function(i, a_, b_) for i in x_plot]

    plt.plot(x_plot, y_plot, label='График заданной функции')
    plt.scatter(ind, [function(i, a_, b_) for i in ind], label='Точки графика')
    plt.legend(fontsize=8)
    plt.show()


if __name__ == "__main__":
    data = devtools.data
    table = devtools.table
    splain(data["n"], data["m"], data["c"], data["d"], data["a"], data["b"])
