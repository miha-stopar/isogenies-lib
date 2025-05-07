
#[cfg(test)]
mod tests {
    use crate::{ec_lit, util::{generate_random_range, Big}};
    use num_traits::Pow;

    #[test]
    fn lit_sigamal_test() {
        /*
        lam = 128
        a = 3 * lam
        b = 162
        c = 56
        f = 30

        p = (2**(a+2))*(3**b)*(5**c)*f - 1
        p = 1290217975993796939363993419446162388979006021159541007293712082644700121088673466685157498316158528176855539315411759315356741765308895915108991692829754882889263058278152142847999999999999999999999999999999999999999999999999999999999
        */
        let l_c = 5;
        let lit_sigamal = ec_lit::LITSiGamal::new(128);
        let (pub_key, alice_secret)= lit_sigamal.generate_pub_key();
        // let (pub_key, alice_secret)= lit_sigamal.generate_pub_key_dbg();

        let mu = l_c * generate_random_range(0.big(), l_c.big() * l_c.big().pow(pub_key.power_c - 1) - 1) +
                    generate_random_range(1.big(), 4.big()); // TODO: check if this can be 4
        let mu = 234234.big(); // TODO, dbg

        println!("");
        println!("mu: {:?}", mu);
        println!("");
        let cipher = lit_sigamal.encrypt(&pub_key, mu.clone());
        let mu_decrypted = lit_sigamal.decrypt(&cipher, alice_secret);

        let check = l_c.big().pow(pub_key.power_c) - mu.clone();

        if mu_decrypted == mu || check == mu_decrypted {
            println!("Decryption successful!");
        } else {
            println!("Decryption failed.");
        } 
        println!("");
        println!("mu: {:?}", mu_decrypted);
        println!("");
    }

}
