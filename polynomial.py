import numpy as np
import math
import copy

class Polynomial:
    
    def __init__(self, coef:np.array, exp:np.array, ordering="lex"):
        self.coef = coef
        self.exp = exp
        self.ordering = ordering
        self.order()
        self.simplify()
        
    def order(self):
        if len(self.coef) > 1:
            if self.ordering == "lex":
                comp = lex_compare
            
            for i in range(len(self.coef)):
                max_idx = i
                for j in range(i+1, len(self.coef)):
                    if comp(self.exp[max_idx], self.exp[j]) < 0:
                        max_idx = j
                
                self.coef[i], self.coef[max_idx] = self.coef[max_idx], self.coef[i]
                self.exp[i], self.exp[max_idx] = copy.copy(self.exp[max_idx]), copy.copy(self.exp[i])
    
    def simplify(self):
        num_coef = len(self.coef)
        i = 0
        while i < num_coef:
            j = i+1
            while j < num_coef:
                if (self.exp[i] == self.exp[j]).all():
                    self.coef[i] += self.coef[j]
                    self.exp = np.delete(self.exp, (j), axis=0)
                    self.coef = np.delete(self.coef, (j), axis=0)
                    num_coef -= 1
                else:
                    j += 1
            i += 1

        num_coef = len(self.coef)
        i = 0
        while i < num_coef:
            if self.coef[i] == 0:
                self.exp = np.delete(self.exp, (i), axis=0)
                self.coef = np.delete(self.coef, (i), axis=0)
                num_coef -= 1
            else:
                i += 1

        
    def remove_leading(self):
        self.exp = np.delete(self.exp, (0), axis=0)
        self.coef = np.delete(self.coef, (0), axis=0)

def add(p_1:Polynomial, p_2:Polynomial):
    if p_1 is None:
        return p_2
    if p_2 is None:
        return p_1
    
    assert p_1.ordering == p_2.ordering

    coef = np.append(copy.copy(p_1.coef), copy.copy(p_2.coef))
    exp = np.append(copy.copy(p_1.exp), copy.copy(p_2.exp), axis=0)

    return Polynomial(coef, exp, p_1.ordering)

def subtract(p_1:Polynomial, p_2:Polynomial):
    if p_1 is None:
        return p_2
    if p_2 is None:
        return p_1
    
    assert p_1.ordering == p_2.ordering

    coef = np.append(copy.copy(p_1.coef), -1*copy.copy(p_2.coef))
    exp = np.append(copy.copy(p_1.exp), copy.copy(p_2.exp), axis=0)

    return Polynomial(coef, exp, p_1.ordering)

def multiply(p_1:Polynomial, p_2:Polynomial):
    if p_1 is None:
        return p_2
    if p_2 is None:
        return p_1
    
    assert p_1.ordering == p_2.ordering
    
    coef = np.zeros(len(p_1.coef) * len(p_2.coef))
    exp = np.zeros((p_1.exp.shape[0] * p_2.exp.shape[0], p_1.exp.shape[1]))

    for i in range(len(p_1.coef)):
        for j in range(len(p_2.coef)):
            idx = len(p_1.coef)*i + j
            coef[idx] = p_1.coef[i] * p_2.coef[j]
            exp[idx,:] = p_1.exp[i] + p_2.exp[j]

    return Polynomial(coef, exp, p_1.ordering)

def divide(f:Polynomial, G:list[Polynomial]):
    Q = [None for g in G]
    
    return divide_rec(copy.deepcopy(f), copy.deepcopy(G), Q, None)

def divide_rec(f:Polynomial, G:list[Polynomial], Q:list[Polynomial], r:Polynomial):
    if len(f.coef) <= 0:
        return Q, r
    for i, g in enumerate(G):
        if (f.exp[0] >= g.exp[0]).all():
            a = Polynomial(
                coef = np.array([f.coef[0]/g.coef[0]]),
                exp = np.array([copy.copy(f.exp[0] - g.exp[0])]),
                ordering=f.ordering
                )
            Q[i] = add(Q[i], a)
            f = subtract(f, multiply(a, g))
            return divide_rec(f, G, Q, r)
    lt = Polynomial(
                coef = np.array([f.coef[0]]),
                exp = np.array([f.exp[0]]),
                ordering=f.ordering
                )
    r = add(r, lt)
    f.remove_leading()
    return divide_rec(f, G, Q, r)

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

    return subtract(multiply(M_over_L1, f1), multiply(M_over_L2, f2))

def buchberger_algorithm(G:list[Polynomial]):
    for i in range(len(G)):
        for j in range(i+1, len(G)):
            S_ij = cancel_leading_terms(G[i], G[j])
            
            _, r = divide(S_ij, G)
            if r is not None:
                G.append(r)
                return buchberger_algorithm(G)
    return G

def reduce_grobner_basis(G:list[Polynomial]):
    for i in range(len(G)):
        _, G[i] = divide(G[i], [G[j] for j in range(len(G)) if j != i and G[j] is not None])
    
    G = [g for g in G if g is not None]
    for g in G:
        g.coef = 1/g.coef[0]*g.coef
    
    return G

def lex_compare(a:np.array, b:np.array):
    assert len(a) == len(b)
    for i in range(len(a)):
        if a[i] > b[i]:
            return 1
        if b[i] > a[i]:
            return -1
    return 0 