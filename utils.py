from sympy import Symbol, diff

alpha = 3
beta = 2
gamma = 1

x = Symbol('x')
p = 1 + x ** gamma
g = 1 + x
u = x ** alpha * (1 - x) ** beta
f = - diff((p * diff(u, x))) + g * u


def coeff(n, h, p, g):
    a_coeff, b_coeff, c_coeff = ([0] * n for _ in range(3))
    for i in range(n):
        a_coeff[i] = - float(p.subs(x, i * h) / (h * h))
        b_coeff[i] = - float((p.subs(x, (i + 1) * h) + p.subs(x, i * h) + h * h * g.subs(x, i * h)) / (h * h))
        c_coeff[i] = - float(p.subs(x, (i + 1) * h) / (h * h))

    return a_coeff, b_coeff, c_coeff


def ab_mass(n, h, f, a_coeff, b_coeff, c_coeff):
    P_coeff, Q_coeff = ([0] * (n + 1) for _ in range(2))
    P_coeff[0] = 0
    Q_coeff[0] = 0
    for i in range(1, n):
        P_coeff[i + 1] = c_coeff[i] / (b_coeff[i] - P_coeff[i] * a_coeff[i])
        Q_coeff[i + 1] = float((a_coeff[i] * Q_coeff[i] - f.subs(x, i * h)) / (b_coeff[i] - a_coeff[i] * P_coeff[i]))

    return P_coeff, Q_coeff


def mass(n, P_coeff, Q_coeff):
    x_values = [0] * (n + 1)
    x_values[0] = 0
    x_values[n] = 0
    for i in range(n - 1, 0, -1):
        x_values[i] = P_coeff[i + 1] * x_values[i + 1] + Q_coeff[i + 1]
    return x_values


def jacobi(n, h, x_values):
    y = [0.0] * (n + 1)
    tmp = [0.0] * (n + 1)
    iterations = 0
    ai = [0.0] * (n + 1)
    gi = [0.0] * (n + 1)
    fi = [0.0] * (n + 1)

    for i in range(n + 1):
        ai[i] = p.subs(x, i * h) / h / h
        gi[i] = g.subs(x, i * h) / h / h
        fi[i] = f.subs(x, i * h)

    while True:
        iterations += 1
        for i in range(1, n):
            tmp[i] = (fi[i] + ai[i] * y[i - 1] + ai[i + 1] * y[i + 1]) / (ai[i + 1] + ai[i] + gi[i] * h * h)

        norm = abs(y[0] - tmp[0])
        for i in range(1, n):
            if abs(y[i] - tmp[i]) > norm:
                norm = abs(y[i] - tmp[i])
            y[i] = tmp[i]

        if norm <= 10e-6:
            break

    print("")
    print("         i*H       Jacobi[i]           y[i]           J-y")
    for i in range(n + 1):
        print(f"{i * h:12.8f}   {y[i]:12.8f}   {x_values[i]:12.8f}   {abs(y[i] - x_values[i]):12.8f}")
    print("iter =", iterations)


def zadel(n, h, x_values):
    y = [0.0] * (n + 1)
    tmp = [0.0] * (n + 1)
    iterations = 0
    ai = [0.0] * (n + 1)
    gi = [0.0] * (n + 1)
    fi = [0.0] * (n + 1)

    for i in range(n + 1):
        ai[i] = p.subs(x, i * h) / h / h
        gi[i] = g.subs(x, i * h) / h / h
        fi[i] = f.subs(x, i * h)

    while True:
        iterations += 1
        for i in range(1, n):
            tmp[i] = (fi[i] + ai[i] * tmp[i - 1] + ai[i + 1] * y[i + 1]) / (ai[i + 1] + ai[i] + gi[i] * h * h)

        norm = abs(y[0] - tmp[0])
        for i in range(1, n):
            if abs(y[i] - tmp[i]) > norm:
                norm = abs(y[i] - tmp[i])
            y[i] = tmp[i]

        if norm <= 10e-6:
            break

    print("")
    print("         i*H       Zadel[i]           y[i]           Z-y")
    for i in range(n + 1):
        print(f"{i * h:12.8f}   {y[i]:12.8f}   {x_values[i]:12.8f}   {abs(y[i] - x_values[i]):12.8f}")
    print("iter =", iterations)


def relax(n, h, deep=2):
    ai = [0.0] * (n + 1)
    gi = [0.0] * (n + 1)
    fi = [0.0] * (n + 1)

    for i in range(n + 1):
        ai[i] = p.subs(x, i * h)
        gi[i] = g.subs(x, i * h)
        fi[i] = f.subs(x, i * h) * h * h

    iteration_list = []

    def w_opt(w):
        if round(w, 2) == 1.0 or w == 0.0:
            return -1

        iterations = 0
        y = [0.0] * (n + 1)
        tmp = [0.0] * (n + 1)

        while True:
            iterations += 1
            for i in range(1, n):
                tmp[i] = (1 - w) * y[i] + w * (ai[i] * tmp[i - 1] + ai[i + 1] * y[i + 1] + fi[i]) / (
                        ai[i + 1] + ai[i] + gi[i] * h * h)

            norm = abs(y[0] - tmp[0])
            for i in range(1, n):
                if abs(y[i] - tmp[i]) > norm:
                    norm = abs(y[i] - tmp[i])
                y[i] = tmp[i]

            if norm <= 10e-6:
                return iterations

    k = 0
    best_interval = 0, 20
    while k < deep:
        for ki in range(int(best_interval[0]), int(best_interval[1])):
            w = ki * (10 ** -(k+1))
            iterations = w_opt(w)
            print(f"w = {round(w, 2)}, iter = {iterations}")
            if iterations > 0:
                iteration_list.append((iterations, w))
        min_iters = min(iteration_list)
        best_interval = [min_iters[1] - (10 ** -(k+1)), min_iters[1] + (10 ** -(k+1))]
        best_interval = [best_interval[0] * 10, best_interval[1] * 10 + 1]
        print(best_interval)
        k += 1


def main():
    n = 10
    h = 1 / n

    a_coeff, b_coeff, c_coeff = coeff(n, h, p, g)
    P_coeff, Q_coeff = ab_mass(n, h, f, a_coeff, b_coeff, c_coeff)
    x_values = mass(n, P_coeff, Q_coeff)

    print(f"{'x':>12}   {'y(x)':>12}   {'u(x)':>12}  {'?':>12}")
    for i in range(n + 1):
        print(f"{round(i * h, 8):12}   ", end="")
        print(f"{round(x_values[i], 8):12}   ", end="")
        print(f"{round(float(u.subs(x, i * h)), 8):12}   ", end="")
        print(f"{round(abs(x_values[i] - float(u.subs(x, i * h))), 8):12}")

    # jacobi(n, h, x_values)
    # zadel(n, h, x_values)
    relax(n, h)


if __name__ == "__main__":
    main()
