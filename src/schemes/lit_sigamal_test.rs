
#[cfg(test)]
mod tests {
    use crate::{ec_lit, util::Big};

    #[test]
    fn lit_sigamal_128_test() {
        /*
        lam = 128
        a = 3 * lam
        b = 162
        c = 56
        f = 30

        p = (2**(a+2))*(3**b)*(5**c)*f - 1
        p = 1290217975993796939363993419446162388979006021159541007293712082644700121088673466685157498316158528176855539315411759315356741765308895915108991692829754882889263058278152142847999999999999999999999999999999999999999999999999999999999
        */
        let p = "1290217975993796939363993419446162388979006021159541007293712082644700121088673466685157498316158528176855539315411759315356741765308895915108991692829754882889263058278152142847999999999999999999999999999999999999999999999999999999999".big();
        let A = ec_lit::Fq::ZERO;
        let curve = ec_lit::Curve::new(&A);
        let lit_sigamal = ec_lit::LITSiGamal::new(curve, p, 3);
        lit_sigamal.generate_pub_key();
    }
}
