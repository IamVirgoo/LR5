from matplotlib import pyplot as plt
import numpy as np
from sympy import *


def f(arg):
    x = Symbol('x')
    #return ((ln(x))**(13/4)*sin(3*x)).subs(x, arg)
    return np.sin(np.pi * arg)


def create_grid(c, d, m, n, f):
    # Найдем длину каждого подотрезка
    delta_x = (d - c) / n

    # Создадим список узлов сетки
    grid_nodes = [c + i * delta_x for i in range(n + 1)]

    # Разобьем этот список на n списков по m элементов в каждом
    grid = [[grid_nodes[j + i * n] for j in range(m)] for i in range(n)]

    # Вычислим значения функции f(x) в узлах сетки
    values = [[f(x) for x in row] for row in grid]

    # Построим интерполяционные квадратические сплайны
    splines = []
    for i in range(n - 1):
        for j in range(m - 2):
            x1, x2, x3 = grid[i][j], grid[i][j + 1], grid[i + 1][j + 1]
            y1, y2, y3 = values[i][j], values[i][j + 1], values[i + 1][j + 1]

            a = (y1 - 2 * y2 + y3) / ((x1 - x2) ** 2 - 2 * (x2 - x3) ** 2 + (x3 - x1) ** 2)
            b = (y3 - y2) / (x3 - x2) - a * (x3 + x2)
            c = y2 - a * x2 ** 2 - b * x2

            splines.append((a, b, c, x2, x3))

    # Построим график функции
    x = [c + i * (d - c) / 100 for i in range(101)]
    y = [f(xi) for xi in x]
    plt.plot(x, y, 'r-', label='f(x)')

    # Построим интерполяционный сплайн
    for spline in splines:
        a, b, c, x1, x2 = spline
        xs = [x for x in grid[0] if x >= x1 and x <= x2]
        ys = [a * x ** 2 + b * x + c for x in xs]
        plt.plot(xs, ys, 'b-', label='Spline')

    # Отобразим график
    plt.legend()
    plt.show()

    return grid, values, splines


grid, values, splines = create_grid(0, 4, 4, 3, f)

print('Сетка:')
for row in grid:
    print(row)

print('Значения функции в узлах сетки:')
for row in values:
    print(row)

print('Квадратичные сплайны:')
for spline in splines:
    print(f'({spline[0]:.2f})x^2 + ({spline[1]:.2f})x + ({spline[2]:.2f}) на [{spline[3]:.2f}, {spline[4]:.2f}]')
