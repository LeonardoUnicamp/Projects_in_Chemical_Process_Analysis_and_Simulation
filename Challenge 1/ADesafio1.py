# Desafio 1 - Artigo 7
# Aluno: Leonardo Aparecido Ferreira Souza - RA: 2140543

def pvirk(fun, y0, x, method = 'rk4'):
    if method == 'rk1':
        y = rk1(fun, y0, x)
    elif method == 'rk2':
        y = rk2(fun, y0, x)
    elif method == 'rk3':
        y = rk3(fun, y0, x)
    elif method == 'rk4':
        y = rk4(fun, y0, x)
    else:
        print('Método inválido!')
        return None
    return y

# Definindo o método RK de 1º Ordem
def rk1(fun, y0, x):
    n = len(x)
    y = [0] * n
    y[0] = y0
    for i in range(0, n-1):
        h = x[i+1] - x[i]
        y[i+1] = y[i] + h * fun(y[i], x[i])
    return y

# Definindo o método RK de 2º Ordem
def rk2(fun, y0, x):
    n = len(x)
    y = [0] * n
    y[0] = y0
    for i in range(0, n-1):
        h = x[i+1] - x[i]
        k1 = fun(y[i], x[i])
        k2 = fun(y[i] + h * k1, x[i] + h)
        y[i+1] = y[i] + h * (k1 + k2) / 2.0
    return y

# Definindo o método RK de 3º Ordem
def rk3(fun, y0, x):
    n = len(x)
    y = [0] * n
    y[0] = y0
    for i in range(0, n-1):
        h = x[i+1] - x[i]
        k1 = fun(y[i], x[i])
        k2 = fun(y[i] + h / 2.0 * k1, x[i] + h / 2.0)
        k3 = fun(y[i] + 3.0 / 4.0 * h * k2, x[i] + 3.0 / 4.0 * h)
        y[i+1] = y[i] + h * (2.0 * k1 + 3.0 * k2 + 4.0 * k3) / 9.0
    return y

# Definindo o método RK de 4º Ordem
def rk4(fun, y0, x):
    n = len(x)
    y = [0] * n
    y[0] = y0
    for i in range(0, n-1):
        h = x[i+1] - x[i]
        k1 = fun(y[i], x[i])
        k2 = fun(y[i] + h / 2.0 * k1, x[i] + h / 2.0)
        k3 = fun(y[i] + h / 2.0 * k2, x[i] + h / 2.0)
        k4 = fun(y[i] + h * k3, x[i] + h)
        y[i+1] = y[i] + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return y

# Definindo o método RKF de 4-5º Ordem
def rkf45(fun, y0, x, tol = 1e-7):
    n = len(x)
    y = [0] * n
    y[0] = y0
    h0 = 1e-3
    for i in range(0, n-1):
        x_iter = x[i]
        y_iter = y[i]
        h = h0
        while x_iter < x[i+1]:
            if x_iter + h > x[i+1]:
                h = x[i+1] - x_iter
            ynew, Erro = rkf45_passo(fun, y_iter, x_iter, h)
            q = ((tol * h) / (2 * Erro)) ** (1 / 4)
            if Erro < tol:
                y_iter = ynew
                x_iter = x_iter + h
            h = q * h
        y[i+1] = y_iter
    return y

# Definindo o passo para o método de RKF de 4-5º Ordem
def rkf45_passo(fun, y, x, h):
    k1 = fun(y, x)
    k2 = fun(y + 1 / 4 * h * k1, x + 1 / 4 * h)
    k3 = fun(y + h / 32 * (3 * k1 + 9 * k2), x + 3 / 8 * h)
    k4 = fun(y + h / 2197 * (1932 * k1 - 7200 * k2 + 7296 * k3), x + 12 / 13 * h)
    k5 = fun(y + h * (439 / 216 * k1 - 8 * k2 + 3680 / 513 * k3 - 845 / 4104 * k4), x + h)
    k6 = fun(y + h * (- 8 / 27 * k1 + 2 * k2 - 3544 / 2565 * k3 + 1859 / 4104 * k4 - 11 / 40 * k5), x + h / 2)

    y_plus = y + h * (25 / 216 * k1 + 1408 / 2565 * k3 + 2197 / 4104 * k4 - 1 / 5 * k5)
    Erro = abs(1 / 360 * k1 - 128 / 4275 * k3 - 2197 / 75240 * k4 + 1 / 50 * k5 + 2 / 55 * k6)

    return y_plus, Erro

# Definindo o método RKF de 5º Ordem
def rkf5(fun, y0, x, tol = 1e-3):
    n = len(x)
    y = [0] * n
    y[0] = y0
    h = 1e-3
    for i in range(0, n - 1):
        x_iter = x[i]
        y_iter = y[i]
        while x_iter < x[i + 1]:
            if x_iter + h > x[i + 1]:
                h = x[i + 1] - x_iter
            ynew, Erro = rkf5_passo(fun, y_iter, x_iter, h)
            q = ((tol * h) / (2 * Erro)) ** (1 / 4)
            if Erro < tol:
                y_iter = ynew
                x_iter = x_iter + h
            h = q * h
        y[i + 1] = y_iter
    return y

# Definindo o passo para o método de RKF de 5º Ordem
def rkf5_passo(fun, y, x, h):
    k0 = fun(y, x)
    k1 = fun(y + 1 / 6 * h * k0, x + 1 / 6 * h)
    k2 = fun(y + h / 75 * (4 * k0 + 16 * k1), x + 4 / 15 * h)
    k3 = fun(y + h * (5 / 6 * k0 - 8 / 3 * k1 + 5 / 2 * k2), x + 2 / 3 * h)
    k4 = fun(y + h * (- 8 / 5 * k0 + 144 / 25 * k1 - 4 * k2 + 16 / 25 * k3), x + 4 / 5 * h)
    k5 = fun(y + h * (361 / 320 * k0 - 18 / 5 * k1 + 407 / 128 * k2 - 11 / 80 * k3 + 55 / 128 * k4), x + h)
    k6 = fun(y + h * (- 11 / 640 * k0 + 0 * k1 + 11 / 256 * k2 - 11 / 160 * k3 + 11 / 256 * k4 + 0 * k5), x + 0 * h)
    k7 = fun(y + h * (93 / 640 * k0 - 18 / 5 * k1 + 803 / 256 * k2 - 11 / 160 * k3 + 99 / 256 * k4 + 0 * k5
                      + k6), x + h)


    y_plus = y + h * (31 / 384 * k0 + 0 * k1 + 1125 / 2816 * k2 + 9 / 32 * k3 + 125 / 768 * k4 + 5 / 66 * k5)
    Erro = abs(5 / 66 * h * (k0 + k5 - k6 - k7))

    return y_plus, Erro
