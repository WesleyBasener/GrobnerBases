import numpy as np

def lex_compare(a:np.array, b:np.array):
    assert len(a) == len(b)
    for i in range(len(a)):
        if a[i] > b[i]:
            return 1
        if b[i] > a[i]:
            return -1
    return 0 

def gr_lex_compare(a:np.array, b:np.array):
    assert len(a) == len(b)
    if sum(a) > sum(b):
        return 1
    if sum(a) < sum(b):
        return -1
    return lex_compare(a,b)

def grev_lex_compare(a:np.array, b:np.array):
    assert len(a) == len(b)
    if sum(a) > sum(b):
        return 1
    if sum(a) < sum(b):
        return -1
    return -1*lex_compare(a,b)