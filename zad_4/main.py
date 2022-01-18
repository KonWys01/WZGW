import numpy as np

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
    # południk z lambdy
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

    phi = x_92 / (a*A_0)

    while True:
        sigma = a*(A_0*phi - A_2*np.sin(2*phi) + A_4*np.sin(4*phi) - A_6*np.sin(6*phi))
        phi_new = phi + (x_92 - sigma)/(a*A_0)
        if np.fabs(phi_new - phi) < np.deg2rad(0.000001 / 3600):
            break
        else:
            phi = phi_new

    N = a / (np.sqrt(1-e2*(np.sin(phi)**2)))
    M = a*(1-e2)/(np.sqrt((1-e2*np.sin(phi)**2)**3))
    t = np.tan(phi)
    n2 = e2_prim * (np.cos(phi) ** 2)

    old_phi = phi
    phi = old_phi - ((y_92**2)*t)*(1 - ((y_92**2)/(12*N))*(5 + 3*t**2 + n2 - 0*n2*(t**2) - 4*n2**2) + ((y_92**4)/(360*N**4))*(61 + 90*t**2 + 45*t**4))/(2*M*N)
    lam = l0 + (y_92/(N*np.cos(old_phi)))*(1 - (y**2/(6*N**2))*(1 + 2*t**2 + n2) + (y**4/(120*N**4))*(5 + 28*t**2 + 24*t**4 + 6*n2 + 8*n2*t**2))

    A = [np.degrees(i) for i in [phi, lam]]
    return A



print(GK_to_1992(A, 19))
print(GK_to_1992(A))
print(GK_to_2000(A))
