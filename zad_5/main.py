import numpy as np
import math
import matplotlib.pyplot as plt
from prettytable import PrettyTable
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
    # print('geo na xyz', x, y, z)
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


def transformacja(xyz: list):
    x, y, z = [i for i in xyz]
    x_0 = -33.4297
    y_0 = 146.5746
    z_0 = 76.2865
    # m = 1 + 0.8407728 * (10**-6)
    m = 0.8407728 * (10**-6)
    e_x = np.deg2rad(-0.35867/3600)  # alfa
    e_y = np.deg2rad(-0.05283/3600)  # beta
    e_z = np.deg2rad(0.84354/3600)  # gamma

    alfa = e_x
    beta = e_y
    gamma = e_z
    first = np.array([[m, gamma, -beta],
                  [-gamma, m, alfa],
                  [beta, -alfa, m]])
    second = np.array([[x],
                       [y],
                       [z]])
    third = np.array([[x_0],
                      [y_0],
                      [z_0]])
    xyz_k = second + first.dot(second) + third
    return xyz_k


# funkcja z trzeciego zadania
def real_degrees(degrees):
    deg = int(degrees)
    minutes = int((degrees - deg)*60)
    seconds = (degrees - deg - minutes/60) * 3600
    seconds = round(seconds, 5)
    deg = f"{deg:02d}"
    minutes = f"{minutes:02d}"
    int_sec = f"{int(seconds):02d}"
    float_seconds = str(round(seconds-int(seconds), 5))[1:]

    return f'{deg}°{minutes}\'{int_sec}{float_seconds}"'


points = ['A', 'B', 'C', 'D', 'Mean', 'M']
xyz_grs80 = PrettyTable(['GRS80', 'X', 'Y', 'Z'])
transformed = PrettyTable(['Krasowskiego', 'X', 'Y', 'Z'])
geo_krasowskiego = PrettyTable(['Geo Krasowskiego', 'Phi', 'Lambda'])
for index, val in enumerate([A, B, C, D, Mean, M]):
    # XYZ GRS80
    xyz_grs80.add_row([points[index], phi_lambda_to_xyz(val)[0], phi_lambda_to_xyz(val)[1], phi_lambda_to_xyz(val)[2]])

    # XYZ Krasowskiego
    xyz = phi_lambda_to_xyz(val)
    kraso = transformacja(xyz)
    transposed = list(np.transpose(kraso)[0])
    transformed.add_row([points[index], transposed[0], transposed[1], transposed[2]])

    # Geo Krasowskiego
    phi, lam = hirvionen(transposed)
    geo_krasowskiego.add_row([points[index], real_degrees(phi), real_degrees(lam)])


print(xyz_grs80)
print(transformed)
print(geo_krasowskiego)



