import numpy as np
import math
import matplotlib.pyplot as plt
# GRS80
a = 6378137
e2 = 0.00669438002290

# Krasowskiego
a_k = 6378245
e2_k = 0.00669342

# punkty
A = [50 + 15/60, 20 + 45/60]
B = [50, 20 + 45/60]
C = [50 + 15/60, 21 + 15/60]
D = [50, 21 + 15/60]
Mean = [50.125, 21]
M = [50.12527044890504, 21.00065108816731]


def phi_lambda_to_xyz(A: list):
    phi, lam = [np.deg2rad(i) for i in A]
    N = a / (np.sqrt(1-e2*(np.sin(phi)**2)))

    x = N * np.cos(phi) * np.cos(lam)
    y = N * np.cos(phi) * np.sin(lam)
    z = N*(1-e2) * np.sin(phi)
    print('geo na xyz', x, y, z)
    return x, y, z


def hirvionen(xyz: list):
    x, y, z = [i for i in xyz]
    r = np.sqrt(x**2 + y**2)

    epsilon = np.deg2rad(0.00005 / 3600)
    # krok 2
    phi = np.arctan((z / r) * ((1 - e2_k) ** -1))
    while True:
        # krok 3
        N = a_k / (np.sqrt(1-e2_k*(np.sin(phi)**2)))
        h = r / np.cos(phi) - N

        # krok 4
        phi_new = np.arctan((z/r*((1-e2_k)*(N / (N + h)))**-1))

        if np.fabs(phi_new - phi) >= epsilon:
            phi = phi_new
        else:
            break

    lam = np.arctan(y/x)

    # kontrola
    N = a_k / (np.sqrt(1 - e2_k * (np.sin(phi) ** 2)))
    h = r / np.cos(phi) - N
    x = (N+h) * np.cos(phi) * np.cos(lam)
    y = (N+h) * np.cos(phi) * np.sin(lam)
    z = (N*(1 - e2_k) + h) * np.sin(phi)

    return np.degrees(phi), np.degrees(lam)


print(hirvionen(phi_lambda_to_xyz(A)))


