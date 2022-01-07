import math

import numpy as np

s = 51123
nr = 10
e2 = 0.00669438002290
a = 6378137


# dane wierzchołków
A = [50 + 15/60, 20 + 45/60]
B = [50, 20 + 45/60]
C = [50 + 15/60, 21 + 15/60]
D = [50, 21 + 15/60]


def vincent(A, B):
    A = [np.deg2rad(i) for i in A]
    B = [np.deg2rad(i) for i in B]
    b = a * np.sqrt(1 - e2)
    f = 1 - (b/a)
    delta_lambda = B[1] - A[1]
    Ua = np.arctan((1 - f)*np.tan(A[0]))
    Ub = np.arctan((1 - f)*np.tan(B[0]))
    L = delta_lambda

    while True:
        sin_sigma = np.sqrt(
            (np.cos(Ub) * np.sin(L)) ** 2 + (np.cos(Ua) * np.sin(Ub) - np.sin(Ua) * np.cos(Ub) * np.cos(L)) ** 2)
        cos_sigma = np.sin(Ua) * np.sin(Ub) + np.cos(Ua) * np.cos(Ub) * np.cos(L)
        sigma = np.arctan(sin_sigma / cos_sigma)

        sin_a = (np.cos(Ua)*np.cos(Ub)*np.sin(L))/sin_sigma
        cos2_a = 1 - (sin_a**2)
        cos2_sigma_m = cos_sigma - (2*np.sin(Ua)*np.sin(Ub))/cos2_a
        C = (f / 16)*cos2_a*(4 + f*(4-3*cos2_a))
        L_iterated = delta_lambda + (1 - C)*f*sin_a*(sigma + C*sin_sigma*(cos2_sigma_m+C*cos_sigma*(1-2*(cos2_sigma_m**2))))

        # print(sin_sigma)
        # print(cos_sigma)
        # print(sigma)
        # print(sin_a)
        # print(cos2_a)
        # print(cos2_sigma_m)
        # print(C)
        # print(L_iterated)
        # break
        if np.fabs(math.radians(L_iterated - L)) < (0.000001 / 3600):
            break
        else:
            L = L_iterated

    u2 = (a**2 - b**2)*cos2_a/(b**2)
    A = 1 + (u2/16384) * (4096 + u2*(-768+u2*(320 - 175*u2)))
    B = (u2/1024) * (256 + u2*(-128 + u2*(74 - 47*u2)))

    delta_sigma = B*sin_sigma*(cos2_sigma_m +  0.25*B*(cos_sigma*(-1+2*(cos2_sigma_m**2))
                                                       - 1/6*B*cos2_sigma_m*(-3 + 4*(sin_sigma**2))*(-3+4*(cos2_sigma_m**2))))

    s_AB = b*A*(sigma - delta_sigma)
    A_ab = np.arctan((np.cos(Ub)*np.sin(L))/(np.cos(Ua)*np.sin(Ub) - np.sin(Ua)*np.cos(Ub)*np.cos(L)))
    A_ba = np.arctan((np.cos(Ua)*np.sin(L))/(-np.sin(Ua)*np.cos(Ub) + np.cos(Ua)*np.sin(Ub)*np.cos(L))) + np.pi

    print("A_ab, a_ba", A_ab, A_ba)
    return s_AB, math.degrees(A_ab), math.degrees(A_ba)

# def kivioji(fi, lam, A, s):
#     fi = np.deg2rad(fi)
#     A = np.deg2rad(A)
#     lam = np.deg2rad(lam)
#
#     ds = 1.2
#     n = s // ds
#
#     M = (A * (1-e2))/(np.sqrt((1-e2*(np.sin(fi)**2)))**3)
#     N = a / np.sqrt(1 - e2*(np.sin(fi)**2))
#     przyrost_fi = (ds * np.cos(A)) / M
#     przyrost_azymut = np.sin(A)*np.tan(lam)*ds / N
#
#     srednia_szerokosc_geodezyjna = fi + przyrost_fi
#     sredni_azymut = A + przyrost_azymut
#
#     przyrost_fi_2 = ds * np.cos(sredni_azymut) /


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


# print(f"punkt średniej szerokości phi={(A[0]+C[0])/2}, lambda={(A[1]+B[1])/2}")
print(f"Vincent {vincent(A, D)}")
print()
print(real_degrees(1.1233245434))
