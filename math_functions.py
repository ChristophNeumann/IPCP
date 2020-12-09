import numpy as np
import math

def rounding(x, isInt):
    """Returns the componentwise rounding x_check of some point x (list or numpy array),
       where x_check[i] = x[i], if isInt is false, and x_check[i] = round(x[i]), otherwise """
    x = np.array(x)
    x_check = np.array(x)
    x_check[isInt] = np.round(x[isInt])
    return x_check


def floor_g(a, g):
    result = a
    if g != 0:
        result = math.floor(a) - (math.floor(a) % g)
    return result


def gcd_vec(a):
    a = np.array(a)
    a_nz = a[np.nonzero(a)]
    if (any(a_nz == 1)):
        result = 1
    else:
        result = int(a_nz[0])
        for i in a[1:]:
            result = math.gcd(result, int(i))
            if result == 1:
                break
    return result


