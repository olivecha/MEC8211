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


def mms_source_time(r, t):
    De = 1e-10
    Cs = 10.0
    R = 0.5
    r = np.max([r, np.ones_like(r)*1e-20], axis=0)
    t = np.max([t, np.ones_like(t)*1e-20], axis=0)
    T1 = - De*(np.pi**2*Cs*np.cos(np.pi*r/R)/(2*R**2) + 4*r*np.log(t) - 2*(R - r)*np.log(t))
    T2 = - De*(np.pi*Cs*np.sin(np.pi*r/R)/(2*R) + r**2*np.log(t) + 2*r*(-R + r)*np.log(t))/r
    T3 = r**2*(-R + r)/t
    return T1 + T2 + T3


def mms_time(r, t):
    Cs = 10.0
    R = 0.5
    return Cs*(1 - np.cos(np.pi*r/R))/2 + r**2*(-R + r)*np.log(t+1)
