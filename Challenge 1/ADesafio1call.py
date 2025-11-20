# Desafio 1 - Artigo 7
# Aluno: Leonardo Aparecido Ferreira Souza - RA: 2140543
# Guardando funções

# Importando as bibliotecas
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from ADesafio1 import pvirk, rkf45, rkf5
import numpy as np
from math import cos

# Definindo a função a ser utilizada
def dfun(y, x):
    return -50 * (y - cos(x))

# Definindo os parâmetros a serem aplicados
y0 = 1.0
x0 = 0.0
x_end = 1.0
n = 101

x = np.linspace(x0, x_end, n)

# Definindo as funções para plot dos gráficos
y1 = odeint(dfun, y0, x)
y2 = pvirk(dfun, y0, x, 'rk1')
y3 = pvirk(dfun, y0, x, 'rk2')
y4 = pvirk(dfun, y0, x, 'rk3')
y5 = pvirk(dfun, y0, x, 'rk4')
y6 = rkf45(dfun, y0, x)
y7 = rkf5(dfun, y0, x)

# Gráfico teste
plt.plot(x, y1, label = 'odeint')
plt.plot(x, y7, label = 'rkf5')
plt.legend()
plt.show()

# Gráfico comparação entre métodos
plt.plot(x, y1, label = 'odeint')
plt.plot(x, y2, label = 'rk1')
plt.plot(x, y3, label = 'rk2')
plt.plot(x, y4, label = 'rk3')
plt.plot(x, y5, label = 'rk4')
plt.plot(x, y6, label = 'rkf45')
plt.plot(x, y7, label = 'rkf5')
plt.legend()
plt.show()
