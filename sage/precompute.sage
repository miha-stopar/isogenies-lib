from montgomery import *
from endomorphism import *

def get_lit_si_gamal_params(lam):
    a = lam * 3
    if lam == 128:
        b = 162
        c = 56
        f = 30
    elif lam == 192:
        b = 243
        c = 83
        f = 118
    elif lam == 256:
        b = 324
        c = 111
        f = 436
    else:
        raise ValueError("lam should be 128 or 192 or 256.")

    # For public key generation we need integer n: we compute an isogeny
    # of degree l_a**(2*a) - n**2
    n = 3

    return 2, 3, 5, a, b, c, f, n

def generate_lit_si_gamal_pubparam(a, b, c, sqrtminus, Fp, Fp2, p):

    #The curve of j-invariant 1728
    A = [Fp2(1), 2]
    
    checker = 1
    while checker == 1:
        xP = Fp.random_element()
        P = [xP,1]
        P = scalar(P,A,(p+1)//((2**(a+2))*(3**b)*(5**c)))
        Pa = scalar(P,A,((2**(a+1))*(3**b)*(5**c)))
        Pb = scalar(P,A,((2**(a+2))*(3**(b-1))*(5**c)))
        Pc = scalar(P,A,((2**(a+2))*(3**b)*(5**(c-1))))
        if (Pa[1] != 0) and (Pb[1] != 0) and (Pc[1] != 0):
            checker = 0

    Pa = scalar(P, A, (3**b)*(5**c))
    Pbc = scalar(P, A, 2**(a+2))
    Qbc = [Pbc[0]*(-1), Pbc[1]]

    #[2**a]Pa = (1,*).
    Pa_ = scalar(Pa,A,2**a)

    if Pa_[0] != Pa_[1]:
        # print("+++++++++++")
        # print("")
        Pa[0] *= -1

    checker = 1
    while checker == 1:
        xP = Fp2.random_element()
        Qa = [xP,1]
        Qa = scalar(Qa,A,(p+1)//(2**(a+2)))
        Q = scalar(Qa,A,2**(a+1))
        if (Q[0]**2 == (-1)*Q[1]**2) and (Q[1] != 0):
            checker = 0

    Pb = scalar(Pbc, A, 5**c)
    Qb = scalar(Qbc, A, 5**c)

    Pc = scalar(Pbc, A, 3**b)
    Qc = scalar(Qbc, A, 3**b)

    return [Pa, Qa, Pb, Qb, Pc, Qc]

def precompute_lit_si_gamal_128():
    l_a, l_b, l_c, a, b, c, f, n = get_lit_si_gamal_params(128)

    p = (l_a**(a+2))*(l_b**b)*(l_c**c)*f - 1
    print("p:", p)

    Fp = GF(p)
    Fp2 = GF(p**2, name='sqrtminus', modulus = x**2+1)
    Fp2.inject_variables()

    A24 = [Fp2(1),2]
    AA = [4*A24[0] - 2*A24[1],A24[1]]
    # E = EllipticCurve([0, AA[0]/AA[1], 0, 1, 0])

    E = EllipticCurve(Fp2, [0, AA[0]/AA[1], 0, 1, 0])

    # strategy = optimal_strategy_3(b-1,12,8)
    # strategy_for_2dim = optimal_strategy_3(a-2,42,27)

    Pa, Qa, Pb, Qb, Pc, Qc = generate_lit_si_gamal_pubparam(a,b,c,sqrtminus,Fp,Fp2,p)

    # debugging
    """
    Pa = E(239936457053032974426498257089211308595947018800040055914670686390096178988935714368840043826359069061322058177871644057978435193697588598099267263608202635120006132567477404434753796213923487763022733484560617362373039516722337813908, 69960074188124660477965832111637444065604854405019199880757198764430735276866295644931112989091942112861058248498493060491087212119238605937217836757108892693771899969455010596121623171459225191221548831691691692923164176864610525622,1279169858307977916207913236987243269187098539324921249002290563431130699610196110165078316781878400663908105270104897281859213851511782154644230708898212331680403569731960611457995104381900185374998458400815718631034622684917017830196)
    # Pa = [239936457053032974426498257089211308595947018800040055914670686390096178988935714368840043826359069061322058177871644057978435193697588598099267263608202635120006132567477404434753796213923487763022733484560617362373039516722337813908, 1279169858307977916207913236987243269187098539324921249002290563431130699610196110165078316781878400663908105270104897281859213851511782154644230708898212331680403569731960611457995104381900185374998458400815718631034622684917017830196]
    Pa = [Pa[0], Pa[2]]


    Qa = E(sqrtminus*847191550412458376631227213302434154830430062842477619842818281073196329664369080179329258264423813210797526868244525727277157388794058716859319021391575884757943197277555108126807489368191974265817856298769502377523647362609461514286 + 858241212916385629383425522748633394889281159392111931012180800772284129754830305889335753837708313482859677045195270698451185085701509417464479997120233388948624488384389573902889091844679816123745084436323288797606004593860858903778, sqrtminus*668358515306903345536198919830696421876561094159296253599574083808187924005390889533792349310824777171135546043340193175642241403645845871791337886653377611692811639337080462634281600067620304736581938721854101097000871146702045345216 + 178247857566682946740243525391472010226735200696911440360523207834625193942375059469792563231575773387507291835424133549592463081350772112051088456092500484550982378106596337476079915920417043677101127896782247797145579535007519275531, sqrtminus*1141223164172540687167040967324107115467916967310970583655567125475271011847180539104181347031614870272222788989918681388428931724124928228197563918094050174973578530782988115489755079687061573419615647267561259918783189397236828119212 + 706609241234531418177198631468747138110032898940734000468509058201833052278944971120921254803218993559996178074555748532039738203529252134013571084726872858804142853477640595192620487510969690613528203690223707660410795199994253411194)
    # Qa = [sqrtminus*847191550412458376631227213302434154830430062842477619842818281073196329664369080179329258264423813210797526868244525727277157388794058716859319021391575884757943197277555108126807489368191974265817856298769502377523647362609461514286 + 858241212916385629383425522748633394889281159392111931012180800772284129754830305889335753837708313482859677045195270698451185085701509417464479997120233388948624488384389573902889091844679816123745084436323288797606004593860858903778, sqrtminus*1141223164172540687167040967324107115467916967310970583655567125475271011847180539104181347031614870272222788989918681388428931724124928228197563918094050174973578530782988115489755079687061573419615647267561259918783189397236828119212 + 706609241234531418177198631468747138110032898940734000468509058201833052278944971120921254803218993559996178074555748532039738203529252134013571084726872858804142853477640595192620487510969690613528203690223707660410795199994253411194]
    Qa = [Qa[0], Qa[2]]
    """
    # end of debugging 

    print("")
    print("l_:", l_a)
    print("l_:", l_b)
    print("l_:", l_c)
    print("")
    print("l_:", l_a**(a+2))

    print("")
    print("a: ", a)
    print("b: ", b)
    print("c: ", c)
    print("")

    # compute_and_print_endomorphism_matrices_montgomery(E, Fp2, Pa, Qa, l_a, a+2)
    # compute_and_print_endomorphism_matrices_montgomery(E, Fp2, Pb, Qb, l_b, b)
    compute_and_print_endomorphism_matrices_montgomery(E, Fp2, Pc, Qc, l_c, c)

def precompute_sqisign():
    import re

    from sage.misc.banner import require_version
    if not require_version(10, 0, print_message=True):
        exit('')

    for l in open('sqisign_parameters.txt'):
        for k in ('lvl', 'p', 'B'):
            m = re.search(rf'^\s*{k}\s*=\s*([x0-9a-f]+)', l)
            if m:
                v = ZZ(m.groups()[0], 0)
                globals()[k] = v

    L = {l for l,_ in (p**2 - 1).factor(limit=B+5) if l <= B}
    assert 2 in L
    L.remove(2)
    f = (p+1).valuation(2)
    if (p-1).valuation(2) > f:
        raise NotImplementedError('2-power torsion is on twist')
    # exp3 = (p-1).valuation(3)
    exp3 = (p+1).valuation(3)
    if (p-1).valuation(3) > exp3:
        raise NotImplementedError('3-power torsion is on twist')
    Lpls = {l for l in L if (p+1).valuation(l) >= (p-1).valuation(l)}
    Lmin = L - Lpls
    Lpls, Lmin = map(sorted, (Lpls, Lmin))
    Epls = [(p+1).valuation(l) for l in Lpls]
    Emin = [(p-1).valuation(l) for l in Lmin]
    Tpls = prod(l**e for l,e in zip(Lpls,Epls))
    Tmin = prod(l**e for l,e in zip(Lmin,Emin))

    Dcom = (Tpls*Tmin).prime_to_m_part(2)
    Dchall = 2**((p+1).valuation(2))
    # prod(l**(p+1).valuation(l) for l in (2))

    T = Tpls * Tmin

    if p % 4 != 3:
        raise NotImplementedError('requires p ≡ 3 (mod 4)')

    pfact = (p^2-1).factor(limit = 10000000)
    plist = [l for (l,e) in pfact]

    Fp2.<i> = GF((p,2), modulus=[1,0,1])
    Fp4 = Fp2.extension(2,'u')
    E = EllipticCurve(Fp4, [1,0])
    assert E.j_invariant() == 1728
    assert E.is_supersingular()
    assert E.change_ring(Fp2).frobenius() == -p
    assert E.order() == (p^2-1)^2

    from sage.groups.generic import order_from_multiple
    x = Fp4.gen()
    while True:
        x += 1
        try:
            P = E.lift_x(x)
        except ValueError:
            continue
        o = order_from_multiple(P, p^2-1, plist)
        if (T<<f).divides(o):
            P *= o // (T<<f)
            P.set_order(T<<f)
            break
    x = Fp4.gen()
    while True:
        x += 1
        try:
            Q = E.lift_x(x)
        except ValueError:
            continue
        o = order_from_multiple(Q, p^2-1, plist)
        if not (T<<f).divides(o):
            continue
        Q *= o // (T<<f)
        # Q.set_order(T<<f)
        if order_from_multiple(P.weil_pairing(Q, T<<f), T<<f, operation='*') == T<<f:
            break

    basis_two_order = 1<<f
    basis_three_order = 3**exp3

    print("basis three order: ", basis_three_order)
    print("")
    print("")
    print("")

    P_two, Q_two = ZZ(P.order()/basis_two_order)*P, ZZ(Q.order()/basis_two_order)*Q

    P_three, Q_three = ZZ(P.order()/basis_three_order)*P, ZZ(Q.order()/basis_three_order)*Q

    print("exp3: ", exp3)
    
    print("P:")
    print(P)

    print("Q:")
    print(Q)

    print("")
    print("P three order:", P_three.order())
    print("")
    print("Q three order:", Q_three.order())
    print("")

    # compute_and_print_endomorphism_matrices(E, Fp2, P_two, Q_two, 2, f)
    compute_and_print_endomorphism_matrices(E, Fp2, P_three, Q_three, 3, exp3)


# precompute_lit_si_gamal_128()
precompute_sqisign()