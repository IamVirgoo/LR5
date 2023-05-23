from matplotlib import pyplot as plt
import numpy as np
from sympy import *
import devtools


# Метод пятый.
def F(arg):
    x = Symbol('x')
    # return ((ln(x))**(13/4)*sin(3*x)).subs(x, arg)
    return np.sin(np.pi * arg)


def N_func(N, arg):
    x = Symbol('x')
    return N.subs(x, arg)


def gen_matrix(N):
    diag = 4
    eps = 1
    matrix = np.array([np.array([float(0) for _ in range(N)]) for _ in range(N)])
    for i in range(N):
        for j in range(N):
            if i == j:
                matrix[i][j] = diag
            if abs(i - j) == 1:
                matrix[i][j] = eps
    return matrix


def splain(n, m, c, d):
    m *= n

    h = (d - c) / m

    ind = [i for i in np.arange(c, d + h / 2, h)]

    table = np.array([np.array([0. for i in range(m + 1)]) for i in range(2)])
    table[0] = ind

    for i in range(m + 1):
        table[1][i] = F(table[0][i])

    print(table)

    a = [table[1][i] for i in range(m)]
    arg_matrix = gen_matrix(m - 1)
    func_matrix = np.array([3 / h ** 2 * (table[1][i + 2] - 2 * table[1][i + 1] + table[1][i]) for i in range(m - 1)])[..., None]
    C = [0] + list(np.linalg.solve(arg_matrix, func_matrix).reshape(1, m - 1)[0])
    D = [(C[i + 1] - C[i]) / (3 * h) for i in range(m - 1)] + [-C[-1] / (3 * h)]
    b = [(table[1][i + 1] - table[1][i]) / h - C[i] * h - D[i] * h ** 2 for i in range(m)]

    print(f"Коэфициент a: {a}")
    print(f"Коэфициент b: {b}")
    print(f"Коэфициент c: {C}")
    print(f"Коэфициент d: {D}")

    for i in range(m):
        x = Symbol('x')
        N = a[i] + b[i] * (x - table[0][i]) + C[i] * (x - table[0][i]) ** 2 + D[i] * (x - table[0][i]) ** 3
        x_plot = np.arange(c + i * h, c + i * h + h, 0.1)
        y_plot = [N_func(N, i) for i in x_plot]
        plt.plot(x_plot, y_plot)

    print(len(ind))

    x_plot = np.arange(c, d, 0.01)
    y_plot = [F(i) for i in x_plot]
    plt.plot(x_plot, y_plot, label='График заданной функции')
    plt.scatter(ind, [F(i) for i in ind], label='Точки графика')
    plt.legend(fontsize=8)
    plt.show()


if __name__ == "__main__":
    data = devtools.data_dev
    table = devtools.table_dev
    splain(data["n"], data["m"], data["c"], data["d"])
