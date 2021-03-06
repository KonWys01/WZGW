import numpy as np
from prettytable import PrettyTable
from shapely.geometry import Polygon
import math

# Parametry GRS80
e2 = 0.00669438002290
a = 6378137
b2 = a**2*(1-e2)
e2_prim = (a**2 - b2)/b2


A = [50 + 15/60, 20 + 45/60]
B = [50, 20 + 45/60]
C = [50 + 15/60, 21 + 15/60]
D = [50, 21 + 15/60]
Mean = [50.125, 21]
M = [50.12527044890504, 21.00065108816731]


def GK(A, l0=19):  # domyslne 19 dla 2000
    A = [np.deg2rad(i) for i in A]
    phi, lam = A
    l0 = np.deg2rad(l0)
    delta_lambda = lam - l0
    t = np.tan(phi)
    n2 = e2_prim*(np.cos(phi)**2)
    N = a / (np.sqrt(1-e2*(np.sin(phi)**2)))

    A_0 = 1 - e2/4 - (3*(e2**2))/64 - (5*(e2**3))/256
    A_2 = 3/8 * (e2 + (e2**2)/4 + (15*(e2**3))/128)
    A_4 = 15/256*(e2**2 + 3*(e2**3)/4)
    A_6 = 35/3072*(e2**3)
    sigma = a*(A_0*phi - A_2*np.sin(2*phi) + A_4*np.sin(4*phi) - A_6*np.sin(6*phi))
    x = sigma + (delta_lambda**2)/2*N*np.sin(phi)*np.cos(phi)*\
        (1+(delta_lambda**2)/12*(np.cos(phi)**2)*(5-t**2+9*n2+4*n2**2) +
         (delta_lambda**4)/360*(np.cos(phi)**4)*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*(t**2)))

    y = delta_lambda * N * np.cos(phi) * \
        (1 + (delta_lambda ** 2)/6 * (np.cos(phi)**2)*(1 - t**2 + n2) +
         (delta_lambda ** 4)/120*(np.cos(phi)**4)*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))

    return x, y


def from_GK(x, y, l0=19):
    A_0 = 1 - e2 / 4 - (3 * (e2 ** 2)) / 64 - (5 * (e2 ** 3)) / 256
    A_2 = 3 / 8 * (e2 + (e2 ** 2) / 4 + (15 * (e2 ** 3)) / 128)
    A_4 = 15 / 256 * (e2 ** 2 + 3 * (e2 ** 3) / 4)
    A_6 = 35 / 3072 * (e2 ** 3)

    l0 = np.deg2rad(l0)
    phi = x / (a * A_0)

    while True:
        sigma = a * (A_0 * phi - A_2 * np.sin(2 * phi) + A_4 * np.sin(4 * phi) - A_6 * np.sin(6 * phi))
        phi_new = phi + (x - sigma) / (a * A_0)
        if np.fabs(phi_new - phi) < np.deg2rad(0.000001 / 3600):
            break
        else:
            phi = phi_new

    N = a / (np.sqrt(1 - e2 * (np.sin(phi) ** 2)))
    M = a * (1 - e2) / (np.sqrt((1 - e2 * np.sin(phi) ** 2) ** 3))
    t = np.tan(phi)
    n2 = e2_prim * (np.cos(phi) ** 2)

    old_phi = phi

    phi = old_phi - ((y ** 2) * t) * (
                1 - ((y ** 2) / (12 * N ** 2)) * (5 + 3 * t ** 2 + n2 - 0 * n2 * (t ** 2) - 4 * n2 ** 2) + (
                    (y ** 4) / (360 * N ** 4)) * (61 + 90 * t ** 2 + 45 * t ** 4)) / (2 * M * N)
    lam = l0 + (y / (N * np.cos(old_phi))) * (
                1 - (y ** 2 / (6 * N ** 2)) * (1 + 2 * t ** 2 + n2) + (y ** 4 / (120 * N ** 4)) * (
                    5 + 28 * t ** 2 + 24 * t ** 4 + 6 * n2 + 8 * n2 * t ** 2))

    phi, lam = [np.degrees(i) for i in [phi, lam]]
    return [phi, lam, x, y]

def GK_to_1992(A, l0=19):  # domyslne 19 dla 92
    A = [np.deg2rad(i) for i in A]
    phi, lam = A
    l0 = np.deg2rad(l0)
    delta_lambda = lam - l0
    t = np.tan(phi)
    n2 = e2_prim*(np.cos(phi)**2)
    N = a / (np.sqrt(1-e2*(np.sin(phi)**2)))

    A_0 = 1 - e2/4 - (3*(e2**2))/64 - (5*(e2**3))/256
    A_2 = 3/8 * (e2 + (e2**2)/4 + (15*(e2**3))/128)
    A_4 = 15/256*(e2**2 + 3*(e2**3)/4)
    A_6 = 35/3072*(e2**3)
    sigma = a*(A_0*phi - A_2*np.sin(2*phi) + A_4*np.sin(4*phi) - A_6*np.sin(6*phi))
    x = sigma + (delta_lambda**2)/2*N*np.sin(phi)*np.cos(phi)*\
        (1+(delta_lambda**2)/12*(np.cos(phi)**2)*(5-t**2+9*n2+4*n2**2) +
         (delta_lambda**4)/360*(np.cos(phi)**4)*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*(t**2)))

    y = delta_lambda * N * np.cos(phi) * \
        (1 + (delta_lambda ** 2)/6 * (np.cos(phi)**2)*(1 - t**2 + n2) +
         (delta_lambda ** 4)/120*(np.cos(phi)**4)*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))

    m = 0.9993
    x_92 = m * x - 5300000
    y_92 = m * y + 500000

    return x_92, y_92, x, y


def GK_to_2000(A):
    lam = A[1]
    # po??udnik z lambdy
    if lam < 16.5:
        l0 = 15
        strefa = 5
    elif lam < 19.5:
        l0 = 18
        strefa = 6
    elif lam < 22.5:
        l0 = 21
        strefa = 7
    else:
        l0 = 24
        strefa = 8

    A = [np.deg2rad(i) for i in A]
    phi, lam = A
    l0 = np.deg2rad(l0)
    delta_lambda = lam - l0
    t = np.tan(phi)
    n2 = e2_prim*(np.cos(phi)**2)
    N = a / (np.sqrt(1-e2*(np.sin(phi)**2)))

    A_0 = 1 - e2/4 - (3*(e2**2))/64 - (5*(e2**3))/256
    A_2 = 3/8 * (e2 + (e2**2)/4 + (15*(e2**3))/128)
    A_4 = 15/256*(e2**2 + 3*(e2**3)/4)
    A_6 = 35/3072*(e2**3)
    sigma = a*(A_0*phi - A_2*np.sin(2*phi) + A_4*np.sin(4*phi) - A_6*np.sin(6*phi))
    x = sigma + (delta_lambda**2)/2*N*np.sin(phi)*np.cos(phi)*\
        (1+(delta_lambda**2)/12*(np.cos(phi)**2)*(5-t**2+9*n2+4*n2**2) +
         (delta_lambda**4)/360*(np.cos(phi)**4)*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*(t**2)))

    y = delta_lambda * N * np.cos(phi) * \
        (1 + (delta_lambda ** 2)/6 * (np.cos(phi)**2)*(1 - t**2 + n2) +
         (delta_lambda ** 4)/120*(np.cos(phi)**4)*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))

    m = 0.999923
    x_2000 = m * x
    y_2000 = m * y + strefa * 1000000 + 500000

    return x_2000, y_2000, x, y


def from_92(x_92, y_92, l0=19):
    A_0 = 1 - e2 / 4 - (3 * (e2 ** 2)) / 64 - (5 * (e2 ** 3)) / 256
    A_2 = 3 / 8 * (e2 + (e2 ** 2) / 4 + (15 * (e2 ** 3)) / 128)
    A_4 = 15 / 256 * (e2 ** 2 + 3 * (e2 ** 3) / 4)
    A_6 = 35 / 3072 * (e2 ** 3)

    l0 = np.deg2rad(l0)

    m = 0.9993
    x = (x_92 + 5300000) / m
    y = (y_92 - 500000) / m

    phi = x / (a*A_0)

    while True:
        sigma = a*(A_0*phi - A_2*np.sin(2*phi) + A_4*np.sin(4*phi) - A_6*np.sin(6*phi))
        phi_new = phi + (x - sigma)/(a*A_0)
        if np.fabs(phi_new - phi) < np.deg2rad(0.000001 / 3600):
            break
        else:
            phi = phi_new

    N = a / (np.sqrt(1-e2*(np.sin(phi)**2)))
    M = a*(1-e2)/(np.sqrt((1-e2*np.sin(phi)**2)**3))
    t = np.tan(phi)
    n2 = e2_prim * (np.cos(phi) ** 2)

    old_phi = phi

    phi = old_phi - ((y**2)*t)*(1 - ((y**2)/(12*N**2))*(5 + 3*t**2 + n2 - 0*n2*(t**2) - 4*n2**2) + ((y**4)/(360*N**4))*(61 + 90*t**2 + 45*t**4))/(2*M*N)
    lam = l0 + (y/(N*np.cos(old_phi)))*(1 - (y**2/(6*N**2))*(1 + 2*t**2 + n2) + (y**4/(120*N**4))*(5 + 28*t**2 + 24*t**4 + 6*n2 + 8*n2*t**2))

    phi, lam = [np.degrees(i) for i in [phi, lam]]
    return [phi, lam, x, y]


def from_2000(x_2000, y_2000):
    A_0 = 1 - e2 / 4 - (3 * (e2 ** 2)) / 64 - (5 * (e2 ** 3)) / 256
    A_2 = 3 / 8 * (e2 + (e2 ** 2) / 4 + (15 * (e2 ** 3)) / 128)
    A_4 = 15 / 256 * (e2 ** 2 + 3 * (e2 ** 3) / 4)
    A_6 = 35 / 3072 * (e2 ** 3)

    if y_2000 < 6000000:
        strefa = 5
        l0 = 15
    elif y_2000 < 7000000:
        strefa = 6
        l0 = 18
    elif y_2000 < 8000000:
        strefa = 7
        l0 = 21
    else:
        strefa = 8
        l0 = 24

    l0 = np.deg2rad(l0)

    m = 0.999923
    x = x_2000 / m
    y = (y_2000 - strefa*1000000 - 500000) / m

    phi = x / (a * A_0)

    while True:
        sigma = a * (A_0 * phi - A_2 * np.sin(2 * phi) + A_4 * np.sin(4 * phi) - A_6 * np.sin(6 * phi))
        phi_new = phi + (x - sigma) / (a * A_0)
        if np.fabs(phi_new - phi) < np.deg2rad(0.000001 / 3600):
            break
        else:
            phi = phi_new

    N = a / (np.sqrt(1 - e2 * (np.sin(phi) ** 2)))
    M = a * (1 - e2) / (np.sqrt((1 - e2 * np.sin(phi) ** 2) ** 3))
    t = np.tan(phi)
    n2 = e2_prim * (np.cos(phi) ** 2)

    old_phi = phi

    phi = old_phi - ((y ** 2) * t) * (
                1 - ((y ** 2) / (12 * N ** 2)) * (5 + 3 * t ** 2 + n2 - 0 * n2 * (t ** 2) - 4 * n2 ** 2) + (
                    (y ** 4) / (360 * N ** 4)) * (61 + 90 * t ** 2 + 45 * t ** 4)) / (2 * M * N)
    lam = l0 + (y / (N * np.cos(old_phi))) * (
                1 - (y ** 2 / (6 * N ** 2)) * (1 + 2 * t ** 2 + n2) + (y ** 4 / (120 * N ** 4)) * (
                    5 + 28 * t ** 2 + 24 * t ** 4 + 6 * n2 + 8 * n2 * t ** 2))

    phi, lam = [np.degrees(i) for i in [phi, lam]]
    return [phi, lam, x, y]


def surface_area(A, B):  # z zadania 3
    A = [math.radians(i) for i in A]
    B = [math.radians(i) for i in B]
    e = np.sqrt(e2)
    Phi_A = np.sin(A[0])/(1 - e2*(np.sin(A[0])**2)) + np.log((1+e*np.sin(A[0]))/(1-e*np.sin(A[0])))/(2*e)
    Phi_B = np.sin(B[0])/(1 - e2*(np.sin(B[0])**2)) + np.log((1+e*np.sin(B[0]))/(1-e*np.sin(B[0])))/(2*e)

    b2 = (a * np.sqrt(1 - e2))**2
    area = b2*(B[1] - A[1])/2*(Phi_A - Phi_B)
    return round(area, 6)


def pola_powierchni():
    # elipsoidalne
    elipsoidalne = surface_area(A, D)

    # GK
    p = Polygon([GK(A), GK(B), GK(D), GK(C)])
    # print([GK(A), GK(B), GK(C), GK(D)])
    gauss = p.area

    # 92
    p = Polygon([GK_to_1992(A)[:2], GK_to_1992(B)[:2], GK_to_1992(D)[:2], GK_to_1992(C)[:2]])
    pole_92 = p.area

    # 2000
    p = Polygon([GK_to_2000(A)[:2], GK_to_2000(B)[:2], GK_to_2000(D)[:2], GK_to_2000(C)[:2]])
    pole_2000 = p.area

    return elipsoidalne, gauss, pole_92, pole_2000


def znieksztalcenia_92(A: list):  # phi lam x_2000 y_2000
    phi, lam, x, y = [i for i in A]
    # phi, lam = [np.degrees(i) for i in [phi, lam]]

    N = a / (np.sqrt(1 - e2 * math.pow(np.sin(float(phi)), 2)))
    M = (a * (1 - e2)) / (np.sqrt(math.pow(1 - e2 * math.pow(np.sin(float(phi)), 2), 3)))
    # print(phi, lam, N, M)

    Q = np.sqrt(M * N)
    m = 1 + (y**2)/(2*Q**2) + (y**2)/(24*Q**4)

    m_ukladu = m * 0.9993  # 0.9993 dla 92
    z_dlugosci = (m_ukladu - 1) * 1000

    m_pola = (m**2) * 0.9993**2
    z_pola = (m_pola - 1) * 10000
    return m_ukladu, z_dlugosci, m_pola, z_pola


def znieksztalcenia_2000(A: list):  # phi lam x_2000 y_2000
    phi, lam, x, y = [i for i in A]
    # phi, lam = [np.degrees(i) for i in [phi, lam]]

    N = a / (np.sqrt(1 - e2 * math.pow(np.sin(float(phi)), 2)))
    M = (a * (1 - e2)) / (np.sqrt(math.pow(1 - e2 * math.pow(np.sin(float(phi)), 2), 3)))
    # print(phi, lam, N, M)

    Q = np.sqrt(M * N)
    m = 1 + (y**2)/(2*Q**2) + (y**2)/(24*Q**4)

    m_ukladu = m * 0.999923  # 0.9993 dla 92
    z_dlugosci = (m_ukladu - 1) * 1000

    m_pola = (m**2) * 0.999923**2
    z_pola = (m_pola - 1) * 10000
    return m_ukladu, z_dlugosci, m_pola, z_pola


def znieksztalcenia_GK(A: list):  # phi lam x_2000 y_2000
    phi, lam, x, y = [i for i in A]
    # phi, lam = [np.degrees(i) for i in [phi, lam]]
    print('/////////', phi, lam, x, y)

    N = a / (np.sqrt(1 - e2 * math.pow(np.sin(float(phi)), 2)))
    M = (a * (1 - e2)) / (np.sqrt(math.pow(1 - e2 * math.pow(np.sin(float(phi)), 2), 3)))
    # print(phi, lam, N, M)

    Q = np.sqrt(M * N)
    m = 1 + (y**2)/(2*Q**2) + (y**2)/(24*Q**4)

    m_ukladu = m * 1  # 0.9993 dla 92
    z_dlugosci = (m_ukladu - 1) * 1000

    m_pola = (m**2) * 1**2
    z_pola = (m_pola - 1) * 10000
    return m_ukladu, z_dlugosci, m_pola, z_pola


wspolrzedne = PrettyTable(['Punkt', 'X_gk', 'Y_gk', 'X_1992', 'Y_1992', 'X_2000', 'Y_2000'])
skala_dl = PrettyTable(['Punkt', 'mGK', 'K_GK(1km)', 'm_1992', 'K_1992(1km)', 'm_2000', 'K_2000(1km)'])
skala_pol = PrettyTable(['Punkt', 'm^2 GK', 'K^2 GK(1ha)', 'm^2 1992', 'K^2 1992(1ha)', 'm^2 2000', 'K^2 2000(1ha)'])
names = ['A', 'B', 'C', 'D', 'Mean', 'M']
for index, i in enumerate([A, B, C, D, Mean, M]):
    wspolrzedne.add_row([names[index],
                         f"{GK(i)[0]:.3f}",
                         f"{GK(i)[1]:.3f}",
                         f"{GK_to_1992(i)[0]:.3f}",
                         f"{GK_to_1992(i)[1]:.3f}",
                         f"{GK_to_2000(i)[0]:.3f}",
                         f"{GK_to_2000(i)[0]:.3f}"])

    m_gk, k_gk, m2_gk, k2_gk = znieksztalcenia_GK(from_GK(GK(i)[0], GK(i)[1]))
    m_92, k_92, m2_92, k2_92 = znieksztalcenia_92(from_92(GK_to_1992(i)[0], GK_to_1992(i)[1]))
    m_2000, k_2000, m2_2000, k2_2000 = znieksztalcenia_2000(from_2000(GK_to_2000(i)[0], GK_to_2000(i)[1]))
    print(znieksztalcenia_92(from_92(GK_to_1992(i)[0], GK_to_1992(i)[1])))
    skala_dl.add_row([names[index], f"{m_gk:.6f}", f"{k_gk:.3f}", f"{m_92:.6f}", f"{k_92:.3f}", f"{m_2000:.6f}", f"{k_2000:.3f}"])
    skala_pol.add_row([names[index], f"{m2_gk:.6f}", f"{k2_gk:.3f}", f"{m2_92:.6f}", f"{k2_92:.3f}", f"{m2_2000:.6f}", f"{k2_2000:.3f}"])
print(wspolrzedne)

pola = PrettyTable(['Elipsoidalne', 'GK', '1992', '2000'])
pola.add_row([f"{round(i/1000000, 12)} km^2" for i in pola_powierchni()])
print(pola)

print(skala_dl)
print(skala_pol)

print(znieksztalcenia_92(from_92(266221.513, 624724.859)))
