#Montgomery curves arithmetics
from sage.all import EllipticCurve
from utilities.fast_roots import *

def diffadd(P,Q,PQ):

    Sum = [0,0]
    t0 = P[0] + P[1]
    t1 = P[0] - P[1]
    t2 = Q[0] + Q[1]
    t3 = Q[0] - Q[1]
    t0 *= t3
    t1 *= t2
    Sum[0] = PQ[1]*(t0 + t1)**2
    Sum[1] = PQ[0]*(t0 - t1)**2
    return Sum

def double(P,A24):

    A = [4*A24[0] - 2*A24[1],A24[1]]
    Db = [0,0]
    
    t0 = (P[0]-P[1])**2
    t1 = (P[0]+P[1])**2
    Db[1] = 4*A[1]*t0
    Db[0] = Db[1]*t1
    t1 = t1 - t0
    t0 = (2*A[1]+A[0])*t1
    Db[1] += t0
    Db[1] *= t1

    return Db

def double_add(P,Q,PQ,A24):
    
    Sum = [0,0]
    Db = [0,0]
    
    t0 = P[0] + P[1]
    t1 = P[0] - P[1]
    t2 = Q[0] + Q[1]
    t3 = Q[0] - Q[1]
    Db[0] = t0**2
    Db[1] = t1**2
    t0 = t0 * t3
    t1 = t1 * t2
    t2 = Db[0] - Db[1]
    Db[1] *= A24[1]
    Db[0] *= Db[1]
    Sum[0] = A24[0] * t2
    Sum[1] = t0 - t1
    Db[1] += Sum[0]
    Sum[0] = t0 + t1
    Db[1] *= t2
    Sum[0] = Sum[0]**2 * PQ[1]
    Sum[1] = Sum[1]**2 * PQ[0]

    return Db,Sum


def scalar(P,A24,n):

    P0 = [P[0],P[1]]
    Q = double(P,A24)
    for i in range(n.bit_length() - 1):
        if n & (1 << (n.bit_length()-i-2)):
            Q,P = double_add(Q,P,P0,A24)
        
        else:
            P,Q = double_add(P,Q,P0,A24)

    return P


def Ladder_3pt(P,Q,PQ,A24,n):
    
    P0 = [P[0],P[1]]
    Q0 = [Q[0],Q[1]]
    PQ0 = [PQ[0],PQ[1]]

    prevbit = 0
    for i in range(n.bit_length()):
        bit =  n & (1 << i)
        swap = bit ^ prevbit
        prevbit = bit
        
        if swap == 1:
            P0 = [P[0],P[1]]
            P = [PQ[0],PQ[1]]
            PQ = [P0[0],P0[1]]

        Q,PQ = double_add(Q,PQ,P,A24)

    return PQ


def isog_point(P,G):

    phiP = [1,1]
    t0 = P[0] + P[1]
    t1 = P[0] - P[1]
    
    for i in range(len(G)):
        t2 = (G[i][0] - G[i][1]) * t0
        t3 = (G[i][0] + G[i][1]) * t1
        phiP[0] *= (t2 + t3)
        phiP[1] *= (t2 - t3)

    phiP[0] = phiP[0]**2 * P[0]
    phiP[1] = phiP[1]**2 * P[1]
    
    return phiP


def isog_curve(A24,G):

    ell = len(G)
    Aimage = [A24[0]**(2*ell+1),A24[1]**(2*ell+1)]

    Aimage2 = [1,1]
    for i in range(ell):
        Aimage2[0] *= (2*G[i][0])
        Aimage2[1] *= (G[i][1]-G[i][0])

    Aimage[0] *= Aimage2[0]**8
    Aimage[1] *= Aimage2[1]**8
    
    return Aimage


def isog_curve_3(G):

    Aimage = [0,0]
    
    X3 = G[0][0]
    Z3 = G[0][1]
    K1 = X3 - Z3
    K2 = X3 + Z3
    R1 = K1**2
    R2 = K2**2
    R3 = R1 + R2
    R4 = (K1 + K2)**2 - R3
    R3 = R4 + R2
    R4 = R4 + R1
    R5 = 2*(R1 + R4) + R2
    Aimage[0] = R5 * R3
    R5 = (2*(R2 + R3) + R1) * R4
    Aimage[1] = R5 - Aimage[0]
    Aimage[0] = Aimage[0] + Aimage[1]
    
    return Aimage



def get_A(P,Q,PQ):

    xP = P[0]/P[1]
    xQ = Q[0]/Q[1]
    xR = PQ[0]/PQ[1]
    A = (xR*xP+xR*xQ+xP*xQ-1)**2/(4*xP*xQ*xR)
    A += -(xP+xQ+xR)

    return [A+2,4]


def MontoEllip(A24,S):

    A = [4*A24[0] - 2*A24[1],A24[1]]
    E = EllipticCurve([0,A[0]/A[1],0,1,0])

    Points = []
    for i in range(len(S)):
        xP = S[i][0]/S[i][1]
        #yP = (xP**3+xP**2*A[0]/A[1]+xP).sqrt()
        yP = sqrt_Fp2(xP**3+xP**2*A[0]/A[1]+xP)
        P = E(xP,yP)
        Points.append(P)

    return Points


def j_inv_Mon(A24):

    A = (4*A24[0] - 2*A24[1])/A24[1]
    A2 = A^2
    
    j_A = 256*(A2 - 3)**3/(A2 - 4)

    return j_A


def xTriple(P,A24):
    P0 = [P[0],P[1]]
    P2 = double(P,A24)
    P = diffadd(P,P2,P0)

    return P

def xTriplee(P,A24,e):
    for i in range(e):
        P = xTriple(P,A24)

    return P


#The code in https://crypto.stackexchange.com/questions/58375/how-to-get-an-optimal-strategy-in-computing-isogenies-for-sidh-sike provided by Luca De Feo
def optimal_strategy_3(n, p, q):
    S = { 1: [] }
    C = { 1: 0 }
    for i in range(2, n+2):
        b, cost = min(((b, C[i-b] + C[b] + b*p + (i-b)*q) for b in range(1,i)), key=lambda t: t[1])
        S[i] = [b] + S[i-b] + S[b]
        C[i] = cost
    return S[n+1]


def isog_3_strategy(b,A24,points,P,strategy):

    S = [[b,P]]
    i = 0

    while len(S) != 0:

        hP = S.pop()
        h = hP[0]
        P = hP[1]

        if h == 1:
            kernel = [P]
            A24 = isog_curve_3(kernel)
            S2 = []
            while len(S) != 0:
                hP = S.pop(0)
                P = isog_point(hP[1],kernel)
                S2.append([hP[0]-1,P])

            S = S2.copy()
            for j in range(len(points)):
                points[j] = isog_point(points[j],kernel)

        elif (strategy[i] < h) and (0 < strategy[i]):
            S.append([h,P])
            P = xTriplee(P,A24,strategy[i])
            S.append([h-strategy[i],P])
            i += 1

    return [A24,points]


def get_PQb_and_shift(A,A1,Pa,Qa,Pa1,Qa1,tau,Fp2,a,b,p):
    checker = 1
    while checker == 1:
        xP = Fp2.random_element()
        Pb = [xP,1]
        Pb = scalar(Pb,A,(p+1)//(3**b))
        P = scalar(Pb,A,3**(b-1))
        
        xP = Fp2.random_element()
        Qb = [xP,1]
        Qb = scalar(Qb,A,(p+1)//(3**b))
        Q = scalar(Qb,A,3**(b-1))
        if (P[1] != 0) and (scalar(P,A,3)[1] == 0) and (Q[1] != 0) and (scalar(Q,A,3)[1] == 0) and (P[0]*Q[1] != P[1]*Q[0]):
            checker = 0

    #Generate the auxiliary points that we need for computing the isogeny.
    E,S = MontoEllip(A,[Pa,Qa,Pb,Qb])
    E1,S1 = MontoEllip(A1,[Pa1,Qa1])
    if (S[0].weil_pairing(S[1],2**(a+2)))**tau != S1[0].weil_pairing(S1[1],2**(a+2)):
        S[1] = (-1)*S[1]

    T1_0 = (2**a) * S[0]
    T1_1 = (2**a) * S1[0]
    #T1 = (T1_0,T1_1): a point for shift

    Pa_shift_ = S[0] + T1_0
    Qa_shift_ = S[1] + T1_0

    Pa1_shift_ = S1[0] + T1_1
    Qa1_shift_ = S1[1] + T1_1

    Pa_shift = [Pa_shift_[0],Pa_shift_[2]]
    Qa_shift = [Qa_shift_[0],Qa_shift_[2]]
    Pa1_shift = [Pa1_shift_[0],Pa1_shift_[2]]
    Qa1_shift = [Qa1_shift_[0],Qa1_shift_[2]]
    
    PQbell = S[2] + S[3]

    Pb_shift_ = S[2] + 3*T1_0
    Qb_shift_ = S[3] + 3*T1_0
    PQb_shift_ = PQbell + 3*T1_0
    
    PQb = [PQbell[0],PQbell[2]]

    Pb_shift = [Pb_shift_[0],Pb_shift_[2]]
    Qb_shift = [Qb_shift_[0],Qb_shift_[2]]
    PQb_shift = [PQb_shift_[0],PQb_shift_[2]]

    PQa = [A,Pa,Qa]
    PQa1 = [A1,Pa1,Qa1]
    PQa_shift = [Pa_shift,Qa_shift,Pa1_shift,Qa1_shift]
    
    PQball = [Pb,Qb,PQb]
    PQb_shift = [Pb_shift,Qb_shift,PQb_shift]

    return PQa,PQa1,PQa_shift,PQball,PQb_shift
