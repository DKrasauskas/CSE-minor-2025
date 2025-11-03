import numpy as np

def f(base, p1, p2, p3, p4, alpha = 40):
    return (
             np.exp(-alpha * (base.i * base.hx - p1[0]) ** 2 - alpha * (base.j * base.hy - p1[1]) ** 2)
            + np.exp(-alpha * (base.i * base.hx - p2[0]) ** 2 - alpha * (base.j * base.hy - p2[1]) ** 2)
            + np.exp(-alpha * (base.i * base.hx - p3[0]) ** 2 - alpha * (base.j * base.hy - p3[1]) ** 2)
            + np.exp(-alpha * (base.i * base.hx - p4[0]) ** 2 - alpha * (base.j * base.hy - p4[1]) ** 2)
    )

def f_wave(base, p1, p2, p3, p4, alpha = 40, omega = 4 * np.pi):
    return  np.sin(omega*base.t)*(
              np.exp(-alpha * (base.i * base.hx - p1[0]) ** 2 - alpha * (base.j * base.hy - p1[1]) ** 2)
            + np.exp(-alpha * (base.i * base.hx - p2[0]) ** 2 - alpha * (base.j * base.hy - p2[1]) ** 2)
            + np.exp(-alpha * (base.i * base.hx - p3[0]) ** 2 - alpha * (base.j * base.hy - p3[1]) ** 2)
            + np.exp(-alpha * (base.i * base.hx - p4[0]) ** 2 - alpha * (base.j * base.hy - p4[1]) ** 2)
    )

def k_vectorized(base, params = None):
    xarr = base.i * base.hx < base.Lx * 0.5
    yarr = base.j * base.hy < base.Ly * 0.5
    K = np.empty_like(base.i, dtype=np.float64)
    for i in range(len(K)):
        for j in range(len(K[i])):
            if xarr[i, j]:
                if yarr[i, j]:
                    K[i, j] = 0.1
                else:
                     K[i, j] = 0.4
            else:
                if yarr[i, j]:
                    K[i, j] = 1.0
                else:
                    K[i, j]=  0.7
    return K

def k(base, i, j,params = None):

    if i * base.hx < base.Lx * 0.5:
        if j * base.hy < base.Ly * 0.5:
            return 0.1
        else:
            return  0.4
    else:
        if j * base.hy < base.Ly * 0.5:
            return 1.0
        else:
            return 0.7
    return K