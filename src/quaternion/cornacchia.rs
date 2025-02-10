use rug::{integer::IsPrime, Integer};

use crate::util::{sqrt_mod_p, valuation, Big};

/// Compute the product of complex numbers a and b.
fn complex_mul(re_a: Integer, im_a: Integer, re_b: Integer, im_b: Integer) -> (Integer, Integer) {
    let mut re = re_a.clone() * re_b.clone();
    re = re - im_a.clone() * im_b.clone();
    let im = re_a * im_b + im_a * re_b;

    (re, im)
}

/// Compute a * b^exp where a and b are complex numbers.
fn complex_mul_by_complex_power(
    re_a: Integer,
    im_a: Integer,
    re_b: Integer,
    im_b: Integer,
    exp: u64,
) -> (Integer, Integer) {
    let mut re = Integer::from(1);
    let mut im = Integer::from(0);
    for i in 0..64 {
        (re, im) = complex_mul(re.clone(), im.clone(), re.clone(), im.clone());
        if (exp >> (63 - i)) & 1 == 1 {
            (re, im) = complex_mul(re.clone(), im.clone(), re_b.clone(), im_b.clone());
        }
    }
    complex_mul(re_a.clone(), im_a.clone(), re.clone(), im.clone())
}

/// Compute `x` and `y` such that x^2 + n*y^2 = p where p is prime.
/// It assumes that there is a sqrt of -1 mod p.
pub fn cornacchia_prime(n: Integer, p: Integer) -> (Integer, Integer) {
    let zero = Integer::from(0);
    let one = Integer::from(1);
    let two = Integer::from(2);

    if p == two {
        if n == one {
            return (one.clone(), one);
        } else {
            return (zero.clone(), zero);
        }
    }

    // test coprimality (should always be ok)
    if p.clone().gcd(&n) == one {
        // get sqrt of -n mod p
        let r2_res = sqrt_mod_p(-n.clone(), p.clone());
        let mut r2 = match r2_res {
            Ok(r2) => r2,
            Err(err) => panic!("{}", err),
        };

        let mut prod = p.clone();
        let mut r1 = p.clone();

        let mut r0 = zero.clone();
        while prod >= p {
            (_, r0) = r2.div_rem(r1.clone());
            prod = r0.clone() * r0.clone();
            r2 = r1.clone();
            r1 = r0.clone();
        }

        // test if result is solution
        let mut a = p.clone() - prod.clone();
        (a, r2) = a.div_rem(n.clone());

        let x: Integer;
        let y = a.sqrt();
        if r2 == zero.clone() {
            x = r0.clone();
            a = y.clone() * y.clone();
            a = a * n.clone();
            prod = prod + a;

            if prod == p {
                return (x.abs(), y.abs());
            }
        }
    } else {
        panic!("cornacchia_prime: p and n not coprime")
    }

    (zero.clone(), zero)
}

/// Given (`re_c`, `im_c`) a solution for x^2 + y^2 = n, compute (re, im) - a solution
/// for x^2 + y^2 = `p`, and return the real and imaginary part of (`re_c` + i * `im_c`) * (re + i * im)^val.
/// Note that it holds (`re_c` + i * `im_c`) * (re + i * im)^val = n * `p`^`val`.
/// Arguments:
/// `re_c, im_c`: It holds `re_c`^2 + `im_c`^2 = n. We do not know n here, it is of form
/// p_1^e_1 * ... * p_k^e_k and we would like to compute (`cornacchia_extended` function)
/// x^2 + y^2 = p_1^e_1 * ... * p_k^e_k * p_{k+1}^e_{k+1} * ... * p_m^e_m = N.
/// `p`: A prime factor of N (p_{k+1} in the notation above).
/// `val`: A valuation of prime factor `p` (e_{k+1} in the notation above).
fn cornacchia_extended_prime_loop(
    re_c: Integer,
    im_c: Integer,
    p: Integer,
    val: u64,
) -> (Integer, Integer) {
    let zero = Integer::from(0);
    let one = Integer::from(1);

    let (re, im) = cornacchia_prime(one, p.clone());
    if re != zero && im != zero {
        let (r, i) = complex_mul_by_complex_power(re_c, im_c, re, im, val);
        return (r, i);
    }

    (zero.clone(), zero)
}

/// Find x and y such that x^2 + y^2 = n.
/// It allows to solve x^2 + y^2 = n for some composite numbers.
/// It uses a prime factor decomposition of n via trial division for primes in the list,
/// computes solutions for n's prime factors and then uses complex multiplication.
/// Since (x+iz)(x-iy) = x^2 + y^2, a solution xp,yp for p and xq,yq for q give a solution for pq
/// by computing (xp+iyp)*(xq+iyq).
/// Arguments:
/// `n`: the parameter defining the equation. To get an output if one exists, only 1 of its prime factors can exceed the largest prime in `prime_list`.
/// `prime_list`: List of consecutive primes starting from 2.
/// `how_many_primes`: How many primes from `prime_list` to be used.
/// `bad_primes_prod`: Assumed to be a product of small primes which are 3 mod 4. Used only to accelerate failure in case its gcd with n is not 1. If 0, it is ignored.
pub fn cornacchia_extended(
    n: Integer,
    prime_list: Vec<Integer>,
    how_many_primes: usize,
    bad_primes_prod: Integer,
) -> (Integer, Integer) {
    let zero = Integer::from(0);
    let one = Integer::from(1);

    // If a prime which is 3 mod 4 divides n, extended Cornacchia can't solve the equation
    if bad_primes_prod != zero {
        let q = n.clone().gcd(&bad_primes_prod);
        if q != one {
            return (zero.clone(), zero);
        }
    }

    let mut nodd = n.clone();
    let mut valuations = vec![0; how_many_primes];
    for i in 0..how_many_primes {
        let p = prime_list[i].clone();
        if p.clone().modulo(&Integer::from(4)) == one || i == 0 {
            let v: u64;
            (v, nodd) = valuation(nodd.clone(), p.clone());
            valuations[i] = v;
        }
    }

    if nodd.clone() % 4.big() == 1.big() {
        // we hope the 'unfactored' part is a prime 1 mod 4
        let is_prob_prime = nodd.is_probably_prime(30);
        if is_prob_prime == IsPrime::Probably || is_prob_prime == IsPrime::Yes || nodd == one {
            let mut x: Integer;
            let mut y: Integer;
            if nodd == one {
                // the unfactored part is 1
                x = one;
                y = zero.clone();
            } else {
                // the 'unfactored' part is prime, can use Cornacchia
                (x, y) = cornacchia_prime(one, nodd);
            }

            if !(x == zero && y == zero) {
                // no need to continue if failure here
                for i in 0..how_many_primes {
                    if valuations[i] != 0 {
                        (x, y) = cornacchia_extended_prime_loop(
                            x,
                            y,
                            prime_list[i].clone(),
                            valuations[i],
                        );
                    }
                }
                return (x, y);
            }
        } else {
            return (zero.clone(), zero);
        }
    } else {
        return (zero.clone(), zero);
    }

    (zero.clone(), zero)
}

#[cfg(test)]
mod tests {
    use num_prime::detail::SMALL_PRIMES;
    use rug::Integer;

    use super::*;

    #[test]
    fn cornacchia_prime_test() {
        let test = |n: Integer, p: Integer| {
            let (x, y) = cornacchia_prime(n.clone(), p.clone());
            assert!(x.clone() * x.clone() + n.clone() * y.clone() * y.clone() == p.clone());
        };

        test(Integer::from(1), Integer::from(5));
        test(Integer::from(1), Integer::from(41));
        test(Integer::from(3), Integer::from(7));
        test(Integer::from(1), Integer::from(29));
        test(Integer::from(1), Integer::from(1381));
    }

    #[test]
    fn complex_mul_by_complex_power_test() {
        let test = |re_a: Integer,
                    im_a: Integer,
                    re_b: Integer,
                    im_b: Integer,
                    exp: u64,
                    re_cmp: Integer,
                    im_cmp: Integer| {
            let (x, y) = complex_mul_by_complex_power(re_a, im_a, re_b, im_b, exp);
            assert!(x == re_cmp);
            assert!(y == im_cmp);
        };

        test(
            Integer::from(3),
            -Integer::from(4),
            Integer::from(1),
            Integer::from(2),
            0,
            Integer::from(3),
            -Integer::from(4),
        );
        test(
            Integer::from(3),
            -Integer::from(4),
            Integer::from(1),
            Integer::from(2),
            1,
            Integer::from(11),
            Integer::from(2),
        );
        test(
            Integer::from(3),
            -Integer::from(4),
            Integer::from(1),
            Integer::from(2),
            2,
            Integer::from(7),
            Integer::from(24),
        );
    }

    #[test]
    fn cornacchia_extended_prime_loop_test() {
        let test = |re_c: Integer,
                    im_c: Integer,
                    p: Integer,
                    val: u64,
                    re_cmp: Integer,
                    im_cmp: Integer| {
            let (x, y) = cornacchia_extended_prime_loop(re_c, im_c, p, val);
            assert!(x == re_cmp);
            assert!(y == im_cmp);
        };

        test(
            Integer::from(1),
            Integer::from(1),
            Integer::from(5),
            2,
            -Integer::from(1),
            Integer::from(7),
        );
        test(
            Integer::from(1),
            Integer::from(0),
            Integer::from(2),
            2,
            Integer::from(0),
            Integer::from(2),
        );
        test(
            -Integer::from(8),
            Integer::from(6),
            Integer::from(41),
            1,
            -Integer::from(64),
            -Integer::from(2),
        );
        test(
            Integer::from(31),
            Integer::from(298),
            Integer::from(17),
            1,
            -Integer::from(174),
            Integer::from(1223),
        );

        let re_c = -Integer::from_str_radix("286292335164776146355279307454016", 10).unwrap();
        let im_c = Integer::from_str_radix("29390034047986289563880081110912", 10).unwrap();
        let re_cmp =
            Integer::from_str_radix("17365063775578586201598356324038842078208", 10).unwrap();
        let im_cmp =
            -Integer::from_str_radix("120145123962752358003585692097427531473856", 10).unwrap();
        test(re_c, im_c, Integer::from(37), 11, re_cmp, im_cmp);
    }

    #[test]
    fn cornacchia_extended_test() {
        let mut primes = vec![];
        for p in SMALL_PRIMES {
            primes.push(Integer::from(p));
        }
        let bed_primes_prod = Integer::from(3 * 7 * 11 * 19);

        let test = |n: Integer, how_many_primes: usize| {
            let (x, y) = cornacchia_extended(
                n.clone(),
                primes.clone(),
                how_many_primes,
                bed_primes_prod.clone(),
            );
            assert!(x.clone() * x.clone() + y.clone() * y.clone() == n.clone());
        };

        let mut n = Integer::from_str_radix("158828126405", 10).unwrap();
        test(n, 26);

        n = Integer::from_str_radix("40014512245953376820", 10).unwrap();
        test(n, 26);

        n = Integer::from_str_radix("13588601138373462550928562093344019331858765165447160913867695800135179682916335386618687107041342641315567209049089762325285514045525398790223782031779675017956905089308664998960177516311434400605404622520285574320810358353904423423012810358448642090586663211232583668667454734852141666802903725775590093049043999776164077647933766930542685009426257611359117058369584270139317425716467676553177434794742587363282934876629123721404464944683384304420820986239594188922650094326930372850811146481429730669474269713207255284173968500838746261882707217336240321810667791730219009998376346822365457541240875528870506928821817510280272725782371744788162331496464288694956312256388500637865037597849249048688448700135956112736510275629520042877377629061438837268497413468138373466165224807808346472654909601360038557113093422495558353125000", 10).unwrap();
        test(n, 100);
    }
}
