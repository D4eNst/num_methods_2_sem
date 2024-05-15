import math

N = 10
H = 1.0 / N
Alpha = 1
Beta = 3
Gamma = 3
alphaDoubles = [0] * (N + 1)
betaDoubles = [0] * (N + 1)
aCoeff = [0] * N
bCoeff = [0] * N
cCoeff = [0] * N
xDoubles = [0] * (N + 1)


def P(x):
    return 1 + x ** Gamma


def G(x):
    return 1 + x


def U(x):
    return x ** Alpha * (1 - x) ** Beta


def F(x):
    return 24 * x ** 5 - 45 * x ** 4 + 24 * x ** 3 + 9 * x ** 2 - 18 * x + 6 + (1 + x) * (x ** Alpha * (1 - x) ** Beta)


def Coeff():
    for i in range(N):
        aCoeff[i] = -P(i * H) / (H ** 2)
        bCoeff[i] = -(P((i + 1) * H) + P(i * H) + H ** 2 * G(i * H)) / (H ** 2)
        cCoeff[i] = -P((i + 1) * H) / (H ** 2)


def AB_MASS():
    alphaDoubles[0] = 0
    betaDoubles[0] = 0
    for i in range(1, N):
        alphaDoubles[i + 1] = cCoeff[i] / (bCoeff[i] - alphaDoubles[i] * aCoeff[i])
        betaDoubles[i + 1] = (aCoeff[i] * betaDoubles[i] - F(i * H)) / (bCoeff[i] - aCoeff[i] * alphaDoubles[i])


def MASS():
    xDoubles[0] = 0
    xDoubles[N] = 0
    for i in range(N - 1, 0, -1):
        xDoubles[i] = alphaDoubles[i + 1] * xDoubles[i + 1] + betaDoubles[i + 1]

    for i in range(N + 1):
        print(f"{round(i * H, 8):12}   {round(xDoubles[i], 8):12}   {round(U(i * H), 8):12}   {round(abs(xDoubles[i] - U(i * H)), 8):12}")


def main():
    Coeff()
    AB_MASS()
    MASS()


if __name__ == "__main__":
    main()
