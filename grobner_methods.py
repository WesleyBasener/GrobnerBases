import numpy as np
import copy

from polynomial import Polynomial

def cancel_leading_terms(f1:Polynomial, f2:Polynomial):
    M = np.array([max(f1.exp[0,i], f2.exp[0,i]) for i in range(len(f1.exp[0]))])

    M_over_L1 = Polynomial(
        coef=np.array([1/f1.coef[0]]),
        exp=copy.copy(np.array([M-f1.exp[0]])),
        ordering=f1.ordering
    )
    M_over_L2 = Polynomial(
        coef=np.array([1/f2.coef[0]]),
        exp=copy.copy(np.array([M-f2.exp[0]])),
        ordering=f2.ordering
    )

    return M_over_L1*f1 - M_over_L2*f2

def buchberger_algorithm(G:list[Polynomial]):
    for i in range(len(G)):
        for j in range(i+1, len(G)):
            S_ij = cancel_leading_terms(G[i], G[j])
            
            _, r = S_ij/G
            if r != 0:
                G.append(r)
                return buchberger_algorithm(G)
    return G

def reduce_grobner_basis(G:list[Polynomial]):
    for i in range(len(G)):
        _, G[i] = G[i] / [G[j] for j in range(len(G)) if j != i and G[j] != 0]
    
    G = [g for g in G if g != 0]
    for g in G:
        g.coef = 1/g.coef[0]*g.coef
    
    return G