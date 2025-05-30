
#[cfg(test)]
mod tests {
    use crate::{ec_lit128, ec_lit192, ec_lit256, util::{generate_random_range, Big}};
    use num_traits::Pow;
    use rug::Integer;

    #[test]
    fn lit_sigamal_128_test() {
        let l_c = 5;
        let lit_sigamal = ec_lit128::LITSiGamal::new(128);
        let (pub_key, alice_secret)= lit_sigamal.generate_pub_key();
        let mu: Integer = l_c.big() * 5818.big() + 1;
        let cipher = lit_sigamal.encrypt(&pub_key, mu.clone());
        let mu_decrypted = lit_sigamal.decrypt(&cipher, alice_secret);

        let check = l_c.big().pow(pub_key.power_c) - mu.clone();

        if mu_decrypted == mu || check == mu_decrypted {
            println!("Decryption successful!");
        } else {
            println!("Decryption failed.");
        } 
    }

    #[test]
    fn lit_sigamal_192_test() {
        let l_c = 5;
        let lit_sigamal = ec_lit192::LITSiGamal::new(192);
        let (pub_key, alice_secret)= lit_sigamal.generate_pub_key();
        let mu: Integer = l_c.big() * 5818.big() + 1;
        let cipher = lit_sigamal.encrypt(&pub_key, mu.clone());
        let mu_decrypted = lit_sigamal.decrypt(&cipher, alice_secret);

        let check = l_c.big().pow(pub_key.power_c) - mu.clone();

        if mu_decrypted == mu || check == mu_decrypted {
            println!("Decryption successful!");
        } else {
            println!("Decryption failed.");
        } 
    }
    
    #[test]
    fn lit_sigamal_256_test() {
        let l_c = 5;
        let lit_sigamal = ec_lit256::LITSiGamal::new(256);
        let (pub_key, alice_secret)= lit_sigamal.generate_pub_key();
        let mu: Integer = l_c.big() * 5818.big() + 1;
        let cipher = lit_sigamal.encrypt(&pub_key, mu.clone());
        let mu_decrypted = lit_sigamal.decrypt(&cipher, alice_secret);

        let check = l_c.big().pow(pub_key.power_c) - mu.clone();

        if mu_decrypted == mu || check == mu_decrypted {
            println!("Decryption successful!");
        } else {
            println!("Decryption failed.");
        } 
    }

}
