
#[cfg(test)]
mod tests {
    use crate::{ec_lit, util::{generate_random_range, Big}};
    use num_traits::Pow;

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
        let (a, b, c, f) = ec_lit::get_params(128);
        let A = ec_lit::Fq::ZERO;
        let curve = ec_lit::Curve::new(&A);
        let l_a = 2;
        let l_b = 3;
        let l_c = 5;
        let lit_sigamal = ec_lit::LITSiGamal::new(curve, p, l_a, l_b, l_c, a, b, c, f, 3);
        let (pub_key, alice_secret)= lit_sigamal.generate_pub_key();


        let mu = l_c * generate_random_range(0.big(), l_c.big() * l_c.big().pow(pub_key.power_c - 1) - 1) +
                    generate_random_range(1.big(), 4.big()); // TODO: check if this can be 4
        // mu = randint(0,5**(c-1)-1) * 5 + randint(1,4)
        let cipher = lit_sigamal.encrypt(&pub_key, mu);
        lit_sigamal.decrypt(&pub_key, &cipher, alice_secret);
    }
}
