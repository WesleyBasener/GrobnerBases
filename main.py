import numpy as np
from polynomial import Polynomial, divide, buchberger_algorithm, reduce_grobner_basis
#Polynomial(coef=np.array([]), exp=np.array([]))

if __name__ == "__main__":
    
    f = Polynomial(coef=np.array([3,1]), exp=np.array([[2,4], [3,3]]))
    g = Polynomial(coef=np.array([1]), exp=np.array([[1,4]]))
    
    Q, r = divide(f, [g])
    
    f = Polynomial(coef=np.array([1,1,-1,1]), exp=np.array([[2,0], [1,0], [0,2], [0,1]]))
    g1 = Polynomial(coef=np.array([1,1]), exp=np.array([[1,1], [0,0]]))
    g2 = Polynomial(coef=np.array([1,1]), exp=np.array([[1,0], [0,1]]))
    
    Q, r = divide(f, [g1, g2])
    
    f1 = Polynomial(coef=np.array([1,-1,1]), exp=np.array([[3,1], [1,2], [0,0]]))
    f2 = Polynomial(coef=np.array([1,-1,-1]), exp=np.array([[2,2],[0,3],[0,0]]))

    G = buchberger_algorithm([f1,f2])
    G1 = reduce_grobner_basis(G)

    f1 = Polynomial(coef=np.array([1,-1,1]), exp=np.array([[1,3], [2,1], [0,0]]))
    f2 = Polynomial(coef=np.array([1,-1,-1]), exp=np.array([[2,2],[3,0],[0,0]]))

    G = buchberger_algorithm([f1,f2])
    G2 = reduce_grobner_basis(G)
    pass

    h1 = Polynomial(coef=np.array([1,1,1]), exp=np.array([[2,0],[1,5],[0,4]]))
    h2 = Polynomial(coef=np.array([1,-1,1,-1]), exp=np.array([]))