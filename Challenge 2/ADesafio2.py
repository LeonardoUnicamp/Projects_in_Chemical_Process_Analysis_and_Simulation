# Desafio 2 - Artigo 2
# Aluno: Leonardo Aparecido Ferreira Souza - RA: 2140543

from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
from math import exp, pi

# Dados iniciais:
A0 = 3.447e+59
Cpa = 1007
Cpb = 2500
Cpv = 1867
D0 = 0
E0 = 1
Ea = 364070
Em = 141.3
F = 2
Fa = 60
h = 360000
L = 10              # Escolher esse valor
P = 760
R = 8.314
t0 = 0
tf = 200
T0 = 30
Tf = 60
Ta = 30
Ts = 30
W = 1
X0 = 0.002
Xm = 0.250
Yq = 8.366e+6
alpha = 0.0782
phi = 5             # Escolher esse valor
Entalp = 2500900
mi = 0.0782         # Escolher esse valor
pb = 178.5
n = 1001

# Condições iniciais:
Y0 = [X0, E0, D0, T0]
t = np.linspace(t0, tf, n)

# Definindo as Equações:
def enzima(Y, t):
    X, E, D, T = Y
    k = A0 * exp(-Ea / (R * (T0 + 273.15)))
    Pvap = exp(18.3036 - 3816.44/((T + 273.15) - 46.13))
    H = 0.6246 * Pvap/(P - Pvap)
    Ha = 0.6246 * Pvap/(P - Pvap)
    B = pb * L * pi * (phi/2)**2
    A = L * phi * pi
    Tw = Ts - F * (T - Ts)
    dXdt = mi * X * (1 - X/Xm)
    dEdt = alpha * (E + D) * (1 - ((E + D) / Em)) - k * E
    dDdt = k * E
    dTdt = Fa * Cpa*(Ta - T) + Fa*(Entalp*(Ha - H) + Cpv*(Ha*Ta - H*T)) - h*A*(T-Tw) + B*Yq*dXdt
    return [dXdt, dEdt, dDdt, dTdt]

# Solução:
X = odeint(enzima, Y0, t)
crescim = X[:, 0]
atv = X[:, 1]
desnat = X[:, 2]
temp = X[:, 3]

# Gráfico 1:
plt.plot(t, crescim)
plt.xlabel('Tempo (h)')
plt.ylabel('Taxa de crescimento')
plt.show()

# Gráfico 2:
plt.plot(t, atv)
plt.xlabel('Tempo (h)')
plt.ylabel('Atividade da protease')
plt.show()

# Gráfico 3:
plt.plot(t, desnat)
plt.xlabel('Tempo (h)')
plt.ylabel('Desnaturação da enzima')
plt.show()

# Gráfico 4:
plt.plot(t, temp)
plt.xlabel('Tempo (h)')
plt.ylabel('Temperatura ($^\circ$C)')
plt.show()
