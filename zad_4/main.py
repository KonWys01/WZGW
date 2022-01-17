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


def GK(A, l0=19):
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
    print(lam, l0, delta_lambda)
    x = sigma + (delta_lambda**2)/2*N*np.sin(phi)*np.cos(phi)*\
        (1+(delta_lambda**2)/12*(np.cos(phi)**2)*(5-t**2+9*n2+4*n2**2) +
         (delta_lambda**4)/360*(np.cos(phi)**4)*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*(t**2)))

    y = delta_lambda * N * np.cos(phi) * \
        (1 + (delta_lambda ** 2)/6 * (np.cos(phi)**2)*(1 - t**2 + n2) +
         (delta_lambda ** 4)/120*(np.cos(phi)**4)*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))

    return x, y


print(GK(A))
