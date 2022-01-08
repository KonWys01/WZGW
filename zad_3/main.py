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
Mean = [50.125, 21]


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

    """poprawka azymutu"""
    nominator_A_ab = np.cos(Ub)*np.sin(L)
    denominator_A_ab = np.cos(Ua)*np.sin(Ub) - np.sin(Ua)*np.cos(Ub)*np.cos(L)
    A_ab = np.arctan(nominator_A_ab/denominator_A_ab)
    if nominator_A_ab > 0 and denominator_A_ab > 0:
        pass
    elif nominator_A_ab > 0 and denominator_A_ab < 0:
        A_ab = A_ab + np.pi
    elif nominator_A_ab < 0 and denominator_A_ab < 0:
        A_ab = A_ab + np.pi
    elif nominator_A_ab < 0 and denominator_A_ab > 0:
        A_ab = A_ab + 2*np.pi

    nominator_A_ba = np.cos(Ua)*np.sin(L)
    denominator_A_ba = -np.sin(Ua)*np.cos(Ub) + np.cos(Ua)*np.sin(Ub)*np.cos(L)
    A_ba = np.arctan(nominator_A_ba/denominator_A_ba)
    if nominator_A_ba > 0 and denominator_A_ba > 0:
        pass
    elif nominator_A_ba > 0 and denominator_A_ba < 0:
        A_ba = A_ba + 2*np.pi
    elif nominator_A_ba < 0 and denominator_A_ba < 0:
        A_ba = A_ba + 2*np.pi
    elif nominator_A_ba < 0 and denominator_A_ba > 0:
        A_ba = A_ba + 3*np.pi
    return s_AB, math.degrees(A_ab), math.degrees(A_ba)


def kivioji():
    # link do algorytmu: https://docer.pl/doc/x0e0vv1

    ds = 1030
    s_AB, A_ab = vincent(A, D)[:2]
    s_AB = s_AB / 2
    n = int(s_AB / ds)

    Phi_A, Lambda_A = [math.radians(i) for i in A]
    A_ab = math.radians(A_ab)
    for i in range(n+1):
        M = a*(1-e2)/(np.sqrt((1-e2*np.sin(Phi_A)**2)**3))
        N = a / (np.sqrt(1-e2*(np.sin(Phi_A)**2)))

        Phi_przyrost = ds * np.cos(A_ab) / M
        Az_przyrost = ds * np.sin(A_ab) * np.tan(Phi_A) / N

        mid_phi = Phi_A + 1/2*Phi_przyrost
        mid_az = A_ab + 1/2*Az_przyrost

        M = a*(1-e2)/(np.sqrt((1-e2*np.sin(mid_phi)**2)**3))
        N = a / (np.sqrt(1 - e2 * (np.sin(mid_phi) ** 2)))

        Phi_przyrost = ds*np.cos(mid_az)/M
        Lambda_przyrost = ds*np.sin(mid_az)/(N*np.cos(mid_phi))
        Az_przyrost = np.sin(mid_az)*np.tan(mid_phi)*ds/N

        Phi_A = Phi_A + Phi_przyrost
        Lambda_A = Lambda_A + Lambda_przyrost
        A_ab = A_ab + Az_przyrost

    return math.degrees(Phi_A), math.degrees(Lambda_A), math.degrees(A_ab)


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


print(f"punkt średniej szerokości phi={real_degrees((A[0]+D[0])/2)}, lambda={real_degrees((A[1]+D[1])/2)}")
print(f"Azymut AD: {real_degrees(vincent(A, D)[1])} ----- Azymut DA: {real_degrees(vincent(A, D)[2])}")
print(f"srodkowy punkt-> Phi:{real_degrees(kivioji()[0])}    lambda:{real_degrees(kivioji()[1])}    Azymut:{real_degrees(kivioji()[2])}")
print(f"Odleglosc miedzy punktem średniej szerokości, a środkowym: {round(vincent(Mean, kivioji())[0], 3)}m")