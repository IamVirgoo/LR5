from matplotlib import pyplot as plt
import numpy as np

import devtools
from sympy import *

x = Symbol('x')


"""
Вычисляется количество точек в таблице (переменная number).
Инициализируются переменные ch и zn равными единице. Они будут использоваться для вычисления числителя и знаменателя базисной функции соответственно.
В цикле проходятся все точки кроме i-ой точки, где i = j.
Вычисляется числитель базисной функции Лагранжа ch как произведение (x - table[0][i]). Где table[0][i] - это значение x в i-ой точке из таблицы.
Вычисляется знаменатель базисной функции Лагранжа zn как произведение (table[0][j] - table[0][i]). Где table[0][j] - это значение x в j-ой точке из таблицы.
Возвращается значение базисной функции Лагранжа, вычисленное как ch / zn
"""


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


"""
Функция принимает в качестве входных данных таблицу точек данных table, 
где каждый столбец соответствует отдельной переменной, а каждая строка соответствует другой точке данных. 
Первая строка таблицы содержит значения x, а вторая строка содержит соответствующие значения y. 
Код перебирает каждое значение x в таблице и вычисляет соответствующее значение y для этого x, 
используя базисную функцию Лагранжа и заданные точки данных. 
Затем вычисленное значение y умножается на базисную функцию Лагранжа и добавляется к текущей сумме. 
Наконец, он возвращает сумму как значение интерполированного многочлена в данной точке.
"""


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
    grid = devtools.table_dev


    x_1 = grid[0]
    y_1 = grid[1]

    x_plot = np.arange(grid[0][0], grid[0][-1], 0.1)
    y_plot = [func(i, grid) for i in x_plot]

    for i in range(len(x_plot)):
        print(f"Значению X: {x_plot[i]} соответсвует значение Y: {y_plot[i]}")

    print("Многочлен Лагранжа: ", lagrange_polynomial(grid))

    plt.plot(x_plot, y_plot, label='Многочлен Лагранжа')
    plt.scatter(x_1, y_1, label='Точки графика')
    plt.show()
