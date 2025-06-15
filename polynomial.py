import numpy as np
import math
import copy

class Polynomial:
    
    def __init__(self, coef:np.array, exp:np.array, ordering="lex"):
        self.coef = coef
        self.exp = exp
        self.ordering = ordering

        if self.ordering == "lex":
            self.comp = lex_compare
            
        elif self.ordering == "grlex":
            self.comp = gr_lex_compare

        elif self.ordering == "grevlex":
            self.comp = grev_lex_compare

        else:
            raise Exception(f"Ordering \"{ordering}\" does not exist.")

        self.order()
        self.simplify()
        if len(self.coef) == 0:
            self.coef = np.array([0])
            self.exp = np.zeros((1, exp.shape[1]))
        
    def order(self):
        if len(self.coef) > 1:
            for i in range(len(self.coef)):
                max_idx = i
                for j in range(i+1, len(self.coef)):
                    if self.comp(self.exp[max_idx], self.exp[j]) < 0:
                        max_idx = j
                
                self.coef[i], self.coef[max_idx] = self.coef[max_idx], self.coef[i]
                self.exp[i], self.exp[max_idx] = copy.copy(self.exp[max_idx]), copy.copy(self.exp[i])
    
    def simplify(self):
        if len(self.coef) == 1 and self.coef[0] == 0:
            self.exp = np.zeros((1,len(self.exp[0])))
        else:
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
        if len(self) == 1:
            self.coef = np.zeros(1)
            self.exp = np.zeros((1,len(self.exp[0])))
        else:
            self.exp = np.delete(self.exp, (0), axis=0)
            self.coef = np.delete(self.coef, (0), axis=0)

    def get_leading(self):
        return self.coef[0], self.exp[0]

    def __len__(self):
        return len(self.coef)

    def __eq__(self, other):
        if isinstance(other, (int, float, complex)):
            return len(self.coef) == 1 and self.coef[0] == other
                
        return (self.coef == other.coef).all() and (self.exp == other.exp).all()
    
    def __lt__(self, other):
        for i in len(self.exp):
            if len(other.exp) < i+1:
                return True
            if self.comp(self.exp[i], other.exp[i]) == 1:
                return True
            if self.comp(self.exp[i], other.exp[i]) == -1:
                return False
        return False

    def __gt__(self, other):
        return not (self > other)
    
    def __repr__(self):
        return np.concat((self.coef.reshape(1,len(self.coef)), self.exp)).T
    
    def __str__(self):
        ret = ""
        for i in range(len(self.coef)):
            term = ""
            if self.exp[i][0] > 0:
                if self.exp[i][0] != 1:
                    term = term + f"x^{int(self.exp[i][0])}"
                else:
                    term = term + "x"
            if len(self.exp[i]) > 1 and self.exp[i][1] > 0:
                if self.exp[i][1] != 1:
                    term = term + f"y^{int(self.exp[i][1])}"
                else:
                    term = term + "y"
            if len(self.exp[i]) > 2 and self.exp[i][2] > 0:
                if self.exp[i][2] != 1:
                    term = term + f"z^{int(self.exp[i][2])}"
                else:
                    term = term + "y"
            if self.coef[i] - int(self.coef[i]) == 0:
                coef = int(self.coef[i])
            else:
                coef = self.coef[i]
            if coef == 1 and term != "":
                    if i != 0:
                        ret = ret + "+" + term
                    else:
                        ret = ret + term
            elif coef == -1 and term != "":
                    ret = ret + "-" + term
            else:
                if i == 0:
                    ret = ret + f"{coef}" + term 
                if self.coef[i] > 0 and i != 0:
                    ret = ret + f"+{coef}" + term
                if self.coef[i] < 0 and i != 0:
                    ret = ret + f"{coef}" + term
        return ret
    
    def __add__(self, other):
        assert self.ordering == other.ordering

        coef = np.append(copy.copy(self.coef), copy.copy(other.coef))
        exp = np.append(copy.copy(self.exp), copy.copy(other.exp), axis=0)

        return Polynomial(coef, exp, self.ordering)

    def __sub__(self, other):
        assert self.ordering == other.ordering

        coef = np.append(copy.copy(self.coef), -1*copy.copy(other.coef))
        exp = np.append(copy.copy(self.exp), copy.copy(other.exp), axis=0)

        return Polynomial(coef, exp, self.ordering)
    
    def __mul__(self, other):
        assert self.ordering == other.ordering
    
        coef = np.zeros(len(self.coef) * len(other.coef))
        exp = np.zeros((self.exp.shape[0] * other.exp.shape[0], self.exp.shape[1]))

        for i in range(len(self.coef)):
            for j in range(len(other.coef)):
                idx = len(self.coef)*i + j
                coef[idx] = self.coef[i] * other.coef[j]
                exp[idx,:] = self.exp[i] + other.exp[j]

        return Polynomial(coef, exp, self.ordering)
    
    def __truediv__(self, G):
        if type(G) is not list:
            G = [G]
        
        for g in G:
            assert g != 0
            
        Q = [Polynomial(coef=np.array([0]), exp=np.zeros((1,len(g.exp[0])))) for g in G]
        r = Polynomial(coef=np.array([0]), exp=np.zeros((1,len(self.exp[0]))))

        return divide_rec(copy.deepcopy(self), copy.deepcopy(G), Q, r)
    
def divide_rec(f:Polynomial, G:list[Polynomial], Q:list[Polynomial], r:Polynomial):
    if f == 0:
        return Q, r
    for i, g in enumerate(G):
        if (f.exp[0] >= g.exp[0]).all():
            a = Polynomial(
                coef = np.array([f.coef[0]/g.coef[0]]),
                exp = np.array([copy.copy(f.exp[0] - g.exp[0])]),
                ordering=f.ordering
                )
            Q[i] = Q[i] + a
            f = f - a*g
            return divide_rec(f, G, Q, r)
    lt = Polynomial(
                coef = np.array([f.coef[0]]),
                exp = np.array([f.exp[0]]),
                ordering=f.ordering
                )
    r = r + lt
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