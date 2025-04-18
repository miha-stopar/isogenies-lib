#! /usr/bin/env sage

proof.arithmetic(False)

all_params = [
    {
        "banner": '"toy"',
        "p": 1981610078855556669567692070387599255633681794129762129193163375155414287691918076633799828157027504892442372430192785973110025991848722744133585514692446519295999,
        "b": 212,
        "d1": 415709681189302716705144249507062872921,
        "d2": 732304420484977534689364942307555175375
    },
    {
        "banner": '128-bit security',
        "p": 47024755958377559820179144403517951365146739428682715587273315623743750596061899792398466078144401013185715861778663433311982543504793487377320505349521764013889738945943941940719288568009419068532558629959602395774408311593047625247205227520083266453579670098107898431221700871924390803878116783902285625847004174086777944844067095845116835903577335827518413153233468762445855819038719999,
        "b": 632,
        "d1": 164689827851839165514063091879356521575625930886412402829946019942764725290681,
        "d2": 1587938149066897359671186423809537261105533753596954710157509996890614699089375
    },
    {
        "banner": '192-bit security',
        "p": 3989795910464850234544092016437599262156405882930270748007942277640374329764909450701017382298119417718075228827492438972399189797430748073952975120211362411225539619165193856075955921868203757050356678140740906219251807603970371623231282116806721678274928745107855520069582715077651032774292776341154446758673963216798416628169733661537007586004889965092043219442232326309212914425730233888098063333187024715534353168291196289326542828281483050941456739621249443006377886553788444285609603976701383081322473844061620853271098247180307945494907820660461209043068817353444264450905367117823999,
        "b": 992,
        "d1": 39765229602356933973867546686064878297883961049797387362167989199717196797566944088914230300103697118143423572390375,
        "d2": 460749310898670172007489084722214513166222146216865908558161563778706645819854540437584937770395888833782995576259881
    },
    {
        "banner": '256-bit security',
        "p": 11733423912110592098943129659050844606312554792190272416139222134125260126256293542355894090597151104849782385250047852082095791782963264272422477754752401134729288466033921644065210001521912515472747792624732203648860768810715096574444930741598687096149331706319733114751349674043478808912170856316437502354928382889858902223487847577640593372756405284562741981855046196126797425256011759705930093547767416564454816320263152972258477167858609259892440669154323283277394752920085739054342835100168847323666376793204002796841113778608062619081219422060112959495683623244822159182084177209507036172874169153132065572209685923391796593565672929583719519758210937482419800925740083545418151266675591940500266972810416236580195709873069266211654682012168405578100031819604139619086515972862805157115534015009687340591299724908610158318996682276165670366162825669852127013302207734132139526996112663259337639407303847758149986487766915329594770573087851486904319999,
        "b": 1892,
        "d1": 25426338098014906735793319747762781169403427411554660597931515963102157443824002958298252819578513436218302905815941213389386934746236807610173562623396761,
        "d2": 14305606430026385444469452808214919079481036242301813774147071977803879227277453295742269489644251831348508579743882075254168335354403770189014901261564375
    }
]

def print_int(name, x, N):
    x = int(x)
    print('    const %s: [u64; N] = [' % name)
    mask = 2**64 - 1
    for i in range(0, N):
        if i % 3 == 0:
            if i != 0:
                print('')
            print('        ', end='')
        else:
            print(' ', end='')
        print('0x%016X,' % (x & mask), end='')
        x >>= 64
    print('')
    print('    ];')

def print_intbytes(name, n):
    nbitlen = len(Integer(n).bits())
    nlen = (nbitlen + 7) >> 3
    print('    const %s_BITLEN: usize = %d;' % (name, nbitlen))
    print('    const %s: [u8; %d] = [' % (name, nlen), end='')
    for i in range(0, nlen):
        if i % 12 == 0:
            print('')
            print('        ', end='')
        else:
            print(' ', end='')
        print('0x%02X,' % (n & 255), end='')
        n >>= 8
    print('')
    print('    ];')
    return nbitlen

# Check that all q-th roots of 1 in GF(p^2) have distinct least significant
# 63 bits in Montgomery representation (not counting conjugates).
def check_dlp(p, R, q):
    K2 = GF(p^2, name = 'u', modulus = [1, 0, 1])
    u = K2.gen()
    while True:
        g = K2.random_element()**((p**2 - 1)//q)
        if g != K2(1):
            break
    assert g**q == K2(1)
    tt = []
    x = K2(R)
    for i in range(0, (q + 1) >> 1):
        x0 = x.polynomial()[0]
        tt.append(int(x0) % 2**63)
        x *= g
    tt.sort()
    for i in range(0, len(tt) - 1):
        if tt[i] == tt[i + 1]:
            raise Exception("collision!")

# Check that all w-th roots of 1 in GF(p^2) (with w = 2^lw) have distinct
# least significant 63 bits in Montgomery representation (not counting
# conjugates).
# FIXME remove (unused)
def check_dlp_w(p, R, lw):
    K2 = GF(p^2, name = 'u', modulus = [1, 0, 1])
    u = K2.gen()
    w = 1 << lw
    while True:
        g = K2.random_element()**((p**2 - 1)//w)
        if g**(w//2) != K2(1):
            break
    assert g**w == K2(1)
    tt = []
    x = K2(R)
    for i in range(0, (w >> 1) + 1):
        x0 = x.polynomial()[0]
        tt.append(int(x0) % 2**63)
        x *= g
    tt.sort()
    for i in range(0, len(tt) - 1):
        if tt[i] == tt[i + 1]:
            raise Exception("collision!")

# Compute indices of powers of g to save for DLP with order n.
def simu_dlp_n_inner(dd, base, lg):
    dd.append(base)
    if lg == 1:
        return
    lg0 = lg >> 1
    lg1 = lg - lg0
    simu_dlp_n_inner(dd, base + lg1, lg0)
    simu_dlp_n_inner(dd, base + lg0, lg1)

def simu_dlp_n(lg):
    dd = []
    simu_dlp_n_inner(dd, 0, lg)
    dd.sort()
    dd2 = []
    for i in range(0, len(dd)):
        v = dd[i]
        j = len(dd2)
        if j == 0 or v != dd2[j - 1]:
            dd2.append(v)
    return dd2

def mkfp(params):
    print('')
    print('/// FESTA %s' % params['banner'])
    p = params['p']
    assert is_prime(p)
    assert (p & 3) == 3
    bitlen = len(p.bits())
    print('pub mod p%d {' % bitlen)
    print('    const BITLEN: usize = %d;' % bitlen)
    N = (bitlen + 63) >> 6
    assert N >= 2
    print_int('MODULUS', p, N)
    print_int('HALF_MODULUS', (p + 1) >> 1, N)
    K = Zmod(p)
    print_int('R_VAL', K(2**(64*N)), N)
    print_int('MINUS_R_VAL', -K(2**(64*N)), N)
    print_int('DR_VAL', K(2*(2**(64*N))), N)
    print_int('TR_VAL', K(3*(2**(64*N))), N)
    print_int('QR_VAL', K(4*(2**(64*N))), N)
    print_int('R2_VAL', K(2**(2*64*N)), N)
    print('    const P0I: u64 = %d;' % int(-1/Zmod(2**64)(p)))
    n1 = floor((2*bitlen - 34) / 31)
    n2 = 2*bitlen - 31*n1 - 2
    print_int('TFIXDIV_VAL', K(2**(33*n1 + 64 - n2 + 2*64*N)), N)
    print_int('TDEC_VAL', K(2**(64*(2*N - 1))), N)
    e = (p + 1) >> 2
    el = 0
    while (e & 31) == 0:
        el += 1
        e >>= 5
    eh = []
    while e != 0:
        eh.append(e & 31)
        e >>= 5
    ehlen = len(eh)
    print('    const SQRT_EH: [u8; %d] = [' % ehlen)
    for i in range(0, ehlen):
        if i % 16 == 0:
            if i != 0:
                print('')
            print('        ', end='')
        else:
            print(' ', end='')
        print('%2d,' % eh[i], end='')
    print('')
    print('    ];')
    print('    const SQRT_EL: usize = %d;' % el)
    p1 = p >> (bitlen - 32)
    print('    const P1: u64 = %d;' % p1)
    print('    const P1DIV_M: u64 = %d;' % (((2**32 - p1)*2**64) // p1))
    b = params['b']
    n = 1 << b
    d1 = params['d1']
    d2 = params['d2']
    assert((p + 1) % (n*d1*d2) == 0)
    print_intbytes('P1_N', n)
    print('    const P1_B: usize = P1_N_BITLEN - 1;')
    mb = b;
    mp = (p + 1) >> mb;
    while (mp & 1) == 0:
        mp >>= 1
        mb += 1
    print('    const P1_MB: usize = %d;' % mb)
    dlpk = simu_dlp_n(b)
    print('    const P1_N_DLPK: [usize; %d] = [' % len(dlpk))
    s = ''
    for v in dlpk:
        vs = '%d,' % v
        if len(s) + len(vs) + 1 > 70:
            print('        %s' % s)
            s = vs
        else:
            if len(s) != 0:
                s += ' '
            s += vs
    if len(s) != 0:
        print('        %s' % s)
    print('    ];')
    print_intbytes('P1_DIV_N', (p + 1) // n)
    print_intbytes('P1_D1', d1)
    print_intbytes('P1_DIV_D1', (p + 1) // d1)
    print_intbytes('P1_DIV_D1MB', mp // d1)
    print_intbytes('P1_D2', d2)
    print_intbytes('P1_DIV_D2', (p + 1) // d2)
    print_intbytes('P1_DIV_D2MB', mp // d2)
    d1f = d1.factor()
    d2f = d2.factor()
    print('    pub struct DFactor { pub q: u32, pub e: u32, pub qe: u64, pub qe0i: u64, qeR: u64, qeR2: u64, iqnn: u64, }')
    print('    const P1_DF: [DFactor; %d] = [' % (len(d1f) + len(d2f)));
    qnn = 1
    for i in range(0, len(d1f) + len(d2f)):
        if i < len(d1f):
            (q, e) = d1f[i]
            dd = d1
        else:
            (q, e) = d2f[i - len(d1f)]
            dd = d2
        if i == len(d1f):
            qnn = 1
        if q >= 2**32:
            raise Exception("q does not fit")
        qe = q**e
        if qe >= 2**63:
            raise Exception("q^e does not fit")
        # FIXME: enable it again (no need to run it each time, but it
        # should be checked at least once).
        # check_dlp(p, 2**(64*N), q)
        qe0i = int(-1/Zmod(2**64)(qe))
        qeR = (2**64) % qe
        qeR2 = (2**128) % qe
        iqnn = int((2**64)/Zmod(qe)(qnn))
        qnn *= qe
        print('        DFactor { q: %d, e: %d, qe: %d, qe0i: %d, qeR: %d, qeR2: %d, iqnn: %d, },' % (q, e, qe, qe0i, qeR, qeR2, iqnn))
    print('    ];')
    print('    const P1_DF2_START: usize = %d;' % len(d1f))
    # FIXME remove -- check_dlp_w(p, 2**(64*N), 10)
    K2 = GF(p^2, name = 'u', modulus = [1, 0, 1])
    u = K2.gen()
    nqr = 1
    while (nqr + u).is_square():
        nqr += 1
    print_int('NQR_RE_VAL', K(nqr*(2^(64*N))), N)
    print('')
    print('    crate::fpcore::define_fp_core!{}')
    print('')
    print('    #[cfg(test)]')
    print('    mod tests {')
    print('        crate::fpcore::define_fp_tests!{}')
    print('    }')
    print('}')

print('#![allow(non_snake_case)]')
print('#![allow(non_upper_case_globals)]')
for params in all_params:
    mkfp(params)

#print('/// FESTA "toy"')
#mkfp('', 1981610078855556669567692070387599255633681794129762129193163375155414287691918076633799828157027504892442372430192785973110025991848722744133585514692446519295999)
#print('/// FESTA 128-bit security')
#mkfp('', 34210858419072694252629449264236339819295038177882714279734807826024242779933767310484671308986825809289385747576331492065712983367251481506001257912971613815319572151252069276055733192521959054419063889264829669600156166216063520472565273824830236099118629038414206133311310523899208042810982314166085222399)
#print('/// FESTA 192-bit security')
#mkfp('', 52945065226886545768973412262479390252463512451138474478119480690859004166406736708967110190366993284118191586243411219856259204415138115875661467055185050685622389213167412708846596810421184696547775398893428317977541560707895797517701316631165769571652785712586932342751385675212660962695405018705285220729253344914818373148367547522105287467492192998434906136755787164801476942611816386264353431353922852337823621595529215999)
#print('/// FESTA 256-bit security')
#mkfp('', 71823257994386405102088344965024640946107191389256345350941319407791925784292379309508680212408876253488412320311703370668361876636273656086033150557244003102885858860867344942244472109132224932018147905353524234816387094036415097947580778900149924967476078631896836733419781908980517890439542776241083892771303231793224716283889758090485987373332470987921600744474393725634714736411953274168726495994098127299525993295124325295134829426216019285626834619215315964935866306425983226139299056779564914298522931880577451399701026413843721990218586202370602462293266118367261919320997887999)
#print('/// FESTA+ "toy"')
#mkfp('', 228238118367136035497409838269287214037738932235300025052747824536514814766656865587136311544616926166666794031975075019855961626403180816056638892933119999)
#print('/// FESTA+ 128-bit security')
#mkfp('', 15378815106438587572954061140041242770230069876075065610547931824302491404930307930492046312364521934740624284588456914599446932194469252675249456417439714788703637184654695393084739570987969645708838200958647767565406081610743819769079387718550552575999)
#print('/// FESTA+ 192-bit security')
#mkfp('', 403455541884594839321056491219231444164194064021623906066845159016372755355985245510271878297324306628087746874923430348381235799056427544149082442065044220186442861662587748821242225155682683408822275643308864056737233644844659158899965732826255731119491517969601440425991818083046738035857980744025000074340577429191420271320559017492087830858024072558637650806072409975894115716719509885461332177715199999)
#print('/// FESTA+ 256-bit security')
#mkfp('', 102248908371217583707785449640411678452483985914322329152230555558097829205592222336160658682107783705194570832762295910595714412431275388952504597896882637438992513454117503537216447460568252657355051879653013914282315103082605535744565054448675328150013293349012925424598345867930177658162044970370127017173408498852735146501013033593406767113971823510877233997607987060063974230977375864670796993624522576920433901583098505559032693364867525229436563281641541010513810401668388997069093715012867180325439401855804054378157755596799)
#print('// For tests only')
#print('#[cfg(test)]')
#mkfp('ptests', 115792089237316195423570985008687907853269984665640564039457584007908834671663, with_factors=False)
