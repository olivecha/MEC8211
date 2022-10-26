import numpy as np


def mms_exp(r):
    C = (10/(1 - np.exp(1)))*(1 - np.exp(r  / 0.5)) + (10*(r-0.5)) / (0.5*(1.0 - np.exp(1)))
    return C


def mms_source_exp(r):
    T1= (10 * 1e-10 * np.exp(r/0.5))/(0.5**2 * (1-np.exp(1)))
    T2 = - 1e-10 * (- (10 * np.exp(r/0.5))/(0.5*(1-np.exp(1))) + 10/(0.5 * (1 - np.exp(1))))
    return T1 + T2


def mms_cos(r):
    C = 5*(1 - np.cos((np.pi*r)/0.5))
    return C


def mms_source_cos(r):
    r = np.max([r, np.ones_like(r) * 1e-60], axis=0)
    T1 = - np.pi * 10 * (1e-10) * np.sin(np.pi * r / 0.5)/(2 * 0.5 *r)
    T2 = - np.pi**2 * 10 * (1e-10) * np.cos(np.pi * r / 0.5)/(2*0.5**2)
    return T1 + T2


def mms_poly(r):
    C = 8.0*10*r**3
    return C


def mms_source_poly(r):
    R = -72.0 * 10 * 1e-10 * r
    return R
