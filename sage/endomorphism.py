from sage.groups.generic import order_from_multiple, discrete_log
from sage.rings.finite_rings.integer_mod_ring import Zmod, ZZ
from sage.matrix.all import Matrix
from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
from montgomery import *

def half_endo(summands):
    def _eval(P):
        E = P.curve()
        assert P in E
        F = E.base_field()
        if (halves := P.division_points(2)):
            Q = halves[0]
        else:
            Q = E.change_ring(F.extension(2,'v'))(P)
        R = sum(endo._eval(Q) for endo in summands)
        return E(R)
    return _eval

def dlp(P, Q, R, p, Fp2):
    n = P.order()
    assert(n != 0)
    assert P.order() == Q.order()

    """
    TODO
    pfact = (p^2-1).factor(limit = 10000000)
    plist = [l for (l, e) in pfact]

    print("p list", plist)
    print("")

    print("p:")
    print(p)

    print("")
    print("P order: ", P.order())
    print("")
    print("R order: ", R.order())
    print("")
    o = order_from_multiple(R, p^2-1, plist)
    print("o: ", o)
    print("")
    print(o * R)
    print("--")
    print((o-1) * R)
    print("")
    """
    o = R.order()

    assert o.divides(P.order())
    e = Fp2(P.weil_pairing(Q, n))
    a = Fp2(R.weil_pairing(Q, n))
    a = discrete_log(a, e, n)
    b = Fp2(P.weil_pairing(R, n))
    b = discrete_log(b, e, n)
    assert a*P + b*Q == R
    return a, b

def matrix_of_isogeny(phi, order, P, Q, p, Fp2, E):
    imP, imQ = map(phi, (P, Q))
    vecP = dlp(P, Q, imP, p, Fp2)
    vecQ = dlp(P, Q, imQ, p, Fp2)
    mat = Matrix(Zmod(order), [vecP, vecQ]).transpose()

    # debugging:
    
    endo_1 = E.scalar_multiplication(1)
    endo_i = E.automorphisms()[-2]
    endo_j = E.frobenius_isogeny()
    endo_k = endo_i * endo_j

    print("")
    print("??????????-----------")
    print("P: ", P)
    print("")
    print("Q: ", Q)
    print("")
    print(endo_i(P))
    print("")
    print(imP)

    assert imP == ZZ(mat[0][0])*P + ZZ(mat[1][0])*Q
    assert imQ == ZZ(mat[0][1])*P + ZZ(mat[1][1])*Q
    return mat

def get_endomorphism_application_on_torsion_group(E, Fp2, m, P, Q):
    """
    E is elliptic curve with j-invariant 1728.
    P, Q are generators of E[m].
    """
    endo_1 = E.scalar_multiplication(1)
    # TODO: E.automorphisms()[-1] returns (x, y) |-> (-x, -sqrt(-1) * y), not (x, y) |-> (-x, sqrt(-1) * y).
    endo_i = E.automorphisms()[-2]
    endo_j = E.frobenius_isogeny()
    endo_k = endo_i * endo_j

    gen1 = endo_1#._eval
    gen2 = endo_i#._eval
    gen3 = half_endo([endo_i, endo_j])
    gen4 = half_endo([endo_1, endo_k])

    p = E.base_ring().characteristic()
 
    # TODO: remove p (use Fp2 characteristic) and E
    mat2 = matrix_of_isogeny(gen2, m, P, Q, p, Fp2, E)
    mat3 = matrix_of_isogeny(gen3, m, P, Q, p, Fp2, E)
    mat4 = matrix_of_isogeny(gen4, m, P, Q, p, Fp2, E)

    return mat2, mat3, mat4

def compute_and_print_endomorphism_matrices_montgomery(E, Fp2, Pxz, Qxz, l, e):
    A24 = [Fp2(1),2]
    # Full points (not Montgomery):
    P, Q = MontoEllip(A24, [Pxz, Qxz])

    x, y = P.xy()
    P = E(x, -y)

    x, y = Q.xy()
    Q = E(x, -y)

    print("")
    print("order:")
    print("")
    print(P.order())
    print("")
    print(Q.order())
    print("")
    print("")

    print("")
    print("======================================")
    print("")
    print("P:")
    # print(P[0]/P[2])
    print("")
    print(P[1]/P[2])
    print("")
    print("Q:")
    # print(Q[0]/Q[2])
    print("")
    print(Q[1]/Q[2])
    print("")
    print("")


    compute_and_print_endomorphism_matrices(E, Fp2, P, Q, l, e)

def compute_and_print_endomorphism_matrices(E, Fp2, P, Q, l, e):
    name = str(l) + "**" + str(e)
    PmQ = P - Q
    torsion = P.order()
 
    assert(P.order() == Q.order())
    assert(torsion == l**e)

    mat2, mat3, mat4 = get_endomorphism_application_on_torsion_group(E, Fp2, torsion, P, Q)

    print("Matrix that corresponds to the second generator (endomorphism) of maximal order on E[%s]:" % name)
    print("mat2:")
    print(mat2)
    print("")

    print("Matrix that corresponds to the third generator (endomorphism) of maximal order on E[%s]:" % name)
    print("mat3:")
    print(mat3)
    print("")

    print("Matrix that corresponds to the fourth generator (endomorphism) of maximal order on E[%s]:" % name)
    print("mat4:")
    print(mat4)
    print("")

    s = """{
        "PxRe": "%s",
        "PxIm": "%s",
        "PzRe": "1",
        "PzIm": "0",
        "QxRe": "%s",
        "QxIm": "%s",
        "QzRe": "1",
        "QzIm": "0",
        "PmQxRe": "%s",
        "PmQxIm": "%s",
        "PmQzRe": "1",
        "PmQzIm": "0",
        "power": %s,
        "action_gen_2": ["%s", "%s", "%s", "%s"],
        "action_gen_3": ["%s", "%s", "%s", "%s"],
        "action_gen_4": ["%s", "%s", "%s", "%s"]
    }""" % (P.x()[0], P.x()[1], Q.x()[0], Q.x()[1],
        PmQ.x()[0], PmQ.x()[1], e, mat2[0][0], mat2[0][1], mat2[1][0], mat2[1][1],
        mat3[0][0], mat3[0][1], mat3[1][0], mat3[1][1],
        mat4[0][0], mat4[0][1], mat4[1][0], mat4[1][1])
    

    print("")
    print("s:")
    print("")
    print(s)
    print("")
