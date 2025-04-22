use std::mem;

use gmp_mpfr_sys::gmp::{self, mpz_limbs_read, mpz_size};
use num::BigUint;
use num_bigint::{BigInt, RandBigInt, Sign};
use num_prime::RandPrime;
use num_traits::{One, Pow, Zero};
use rug::{integer::Order, ops::DivRounding, Integer};
use gmp::mpz_t;

use crate::{
    error::{NotCoprime, NotQuadraticResidueError},
    linalg::matrix::Matrix,
};
use anyhow::{bail, Result};

#[cfg(target_arch = "x86")]
pub fn core_cycles() -> u64 {
    use core::arch::x86::{_mm_lfence, _rdtsc};
    unsafe {
        _mm_lfence();
        _rdtsc()
    }
}

#[cfg(target_arch = "x86_64")]
pub fn core_cycles() -> u64 {
    use core::arch::x86_64::{_mm_lfence, _rdtsc};
    unsafe {
        _mm_lfence();
        _rdtsc()
    }
}

#[cfg(target_arch = "aarch64")]
pub fn core_cycles() -> u64 {
    use core::arch::asm;
    let mut x: u64;
    unsafe {
        asm!("dsb sy", "mrs {}, pmccntr_el0", out(reg) x);
    }
    x
}

use rand_core::{CryptoRng, RngCore};
use sha2::{Digest, Sha512};

// Fake RNG for benchmarks only. NOT ACTUALLY SECURE! DO NOT USE!
pub struct DRNG {
    buf: [u8; 64],
    ptr: usize,
}

impl DRNG {
    pub fn new() -> Self {
        // Self::from_seed(&core_cycles().to_le_bytes())
        Self::from_seed(&u64::MAX.to_le_bytes())
    }

    pub fn from_seed(seed: &[u8]) -> Self {
        let mut d = Self {
            buf: [0u8; 64],
            ptr: 0,
        };
        let mut sh = Sha512::new();
        sh.update(seed);
        d.buf[..].copy_from_slice(&sh.finalize());
        d
    }

    pub fn reseed(&mut self, seed: &[u8]) {
        let mut sh = Sha512::new();
        sh.update(&self.buf[32..]);
        sh.update(seed);
        self.buf[..].copy_from_slice(&sh.finalize());
        self.ptr = 0;
    }
}

impl RngCore for DRNG {
    fn next_u32(&mut self) -> u32 {
        let mut buf = [0u8; 4];
        self.fill_bytes(&mut buf);
        u32::from_le_bytes(buf)
    }

    fn next_u64(&mut self) -> u64 {
        let mut buf = [0u8; 8];
        self.fill_bytes(&mut buf);
        u64::from_le_bytes(buf)
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        let len = dest.len();
        let mut off = 0;
        while off < len {
            let mut clen = 32 - self.ptr;
            if clen > (len - off) {
                clen = len - off;
            }
            dest[off..off + clen].copy_from_slice(&self.buf[self.ptr..self.ptr + clen]);
            self.ptr += clen;
            off += clen;
            if self.ptr == 32 {
                let mut sh = Sha512::new();
                sh.update(&self.buf);
                self.buf[..].copy_from_slice(&sh.finalize());
                self.ptr = 0;
            }
        }
    }

    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.fill_bytes(dest);
        Ok(())
    }
}

impl CryptoRng for DRNG {}

/// Trait that implements functionality to get a BigInt from
/// commonly used types.
pub trait Big {
    /// Return a BigInt for the type.
    fn big(self) -> Integer;
}

impl Big for i32 {
    #[inline]
    fn big(self) -> Integer {
        Integer::from(self)
    }
}

impl Big for u32 {
    #[inline]
    fn big(self) -> Integer {
        Integer::from(self)
    }
}

impl Big for u64 {
    #[inline]
    fn big(self) -> Integer {
        Integer::from(self)
    }
}

impl Big for &str {
    #[inline]
    fn big(self) -> Integer {
        self.parse().expect("p should be a number")
    }
}

pub fn bytes_from_str(a: &str) -> Vec<u8> {
    let f = a.parse::<BigInt>().unwrap();
    let bs = f.to_bytes_le();
    bs.1
}

pub fn bits_from_big(s: Integer) -> Vec<u8> {
    let mut n_bits = s.to_digits::<u8>(Order::LsfLe)
        .iter()
        .flat_map(|byte| (0..8).map(move |i| (byte >> i) & 1))
        .collect::<Vec<u8>>();
    let actual_length = s.significant_bits() as usize;
    n_bits.truncate(actual_length);

    n_bits
}

pub fn binary_from_num(a: u32) -> Vec<u32> {
    let s = format!("{a:b}");
    let char_vec: Vec<char> = s.chars().collect();
    char_vec.iter().map(|x| x.to_digit(10).unwrap()).collect()
}

/// Return `a` * `u` + `b` * `v`.
pub fn linear_combination(
    u: Integer,
    a: &Matrix<Integer>,
    v: Integer,
    b: &Matrix<Integer>,
) -> Matrix<Integer> {
    a * u + b * v
}

/// Return the greatest common divisor `g` of `a` and `b`,
/// and return `t`, `v` such that g = t * a + v * b.
pub fn extended_gcd(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    let mut s = (BigInt::zero(), BigInt::one());
    let mut t = (BigInt::one(), BigInt::zero());
    let mut r = (a.clone(), b.clone());

    while !r.0.is_zero() {
        let q = r.1.clone() / r.0.clone();
        let f = |mut r: (BigInt, BigInt)| {
            mem::swap(&mut r.0, &mut r.1);
            r.0 = r.0 - q.clone() * r.1.clone();
            r
        };
        r = f(r);
        s = f(s);
        t = f(t);
    }

    if r.1 >= BigInt::zero() {
        (r.1, t.1, s.1)
    } else {
        (
            BigInt::zero() - r.1,
            BigInt::zero() - t.1,
            BigInt::zero() - s.1,
        )
    }
}

/// Return the modular inverse of `a` modulo `p`.
pub fn mod_inv(a: Integer, p: Integer) -> Result<Integer> {
    let (d, u, _) = a.clone().extended_gcd(p.clone(), Integer::new());
    if d != 1.big() {
        bail!(NotCoprime());
    }
    if u < 0.big() {
        Ok(u + p)
    } else {
        Ok(u)
    }
}

/// Round a/b to closest integer.
pub fn rounded_div(a: Integer, b: Integer) -> Integer {
    let abs_b = b.clone().abs();
    // q is of same sign as a*b (and 0 if a is 0)
    let mut sign_q = a.clone() * b.clone();
    let (mut q, mut r) = a.div_rem(b);
    r = r.abs();
    r = r.clone() + r;
    if r > abs_b {
        r = 0.big();
        if sign_q < r {
            sign_q = -1.big();
        } else {
            sign_q = 1.big();
        }
        q = q + sign_q;
    }

    q
}

/// Return the p-adic valuation `val` of `a`.
/// It holds `a` = `p`^`val` ** `ac`.
pub fn valuation(a: Integer, p: Integer) -> (u64, Integer) {
    let mut ac = a.clone();
    let mut r = Integer::from(0);
    let mut q = ac.clone();
    let mut val = 0;
    while r == Integer::from(0) {
        val += 1;
        ac = q.clone();
        (q, r) = ac.clone().div_rem(p.clone());
    }
    val -= 1;

    (val, ac)
}

/// Return the square root of `a` modulo a prime `p`.
/// Two cases are handled separately:
/// - p mod 4 == 3
/// - p mod 8 == 5
/// Otherwise, if p mod 8 == 1), the Shanks-Tonelli algorithm is applied.
pub fn sqrt_mod_p(a: Integer, p: Integer) -> Result<Integer> {
    let sqrt;

    let a_mod = a % p.clone();
    if a_mod.jacobi(&p) != 1 {
        bail!(NotQuadraticResidueError());
    }
    let pm1 = p.clone() - 1;
    let mut tmp: Integer;
    let one = Integer::from(1);
    let two = Integer::from(2);
    let three = Integer::from(3);
    let five = Integer::from(5);

    if p.clone() % 4 == three {
        tmp = p.clone() + one;
        tmp = tmp >> 2;
        sqrt = a_mod.pow_mod(&tmp, &p).unwrap();
    } else if p.clone() % 8 == five {
        tmp = p.clone() - one.clone();
        tmp = tmp >> 2;
        tmp = a_mod.clone().pow_mod(&tmp, &p).unwrap(); // a^{(p-1)/4}
        if tmp == one {
            tmp = p.clone() + three;
            tmp = tmp >> 3;
            sqrt = a_mod.pow_mod(&tmp, &p).unwrap(); // a^{(p+3)/8} mod p
        } else {
            tmp = p.clone() - five;
            tmp = tmp >> 3; // (p - 5) // 8
            let a4: Integer = a_mod.clone() * 4;
            tmp = a4.pow_mod(&tmp, &p).unwrap();

            let a2 = a_mod * 2;
            tmp = a2 * tmp;
            sqrt = tmp % p;
        }
    } else {
        // p % 8 == 1 -> Shanks-Tonelli
        let mut q: Integer = p.clone() - 1;

        let mut e = 0;

        while !q.get_bit(e) {
            e += 1;
        }
        let two_to_e = 2.pow(e);
        q = q.div_floor(&Integer::from(two_to_e));

        // p - 1 = q * 2^e

        // Find generator - non-quadratic residue
        let mut qnr = two.clone();
        while qnr.legendre(&p) != -1 {
            qnr = qnr + one.clone();
        }
        let mut z = qnr.pow_mod(&q, &p).unwrap(); // c = z^q

        // Initialize
        let mut y = a_mod.clone().pow_mod(&q, &p).unwrap(); // t = n^q

        tmp = q + one.clone();
        tmp = tmp >> 1; // tmp = (q + 1) / 2

        let mut x = a_mod.pow_mod(&tmp, &p).unwrap(); // x = a^(q + 1)/2 mod p
        let mut exp = one << e - 2;

        for _ in 0..e {
            let b = y.clone().pow_mod(&exp, &p).unwrap();

            if b == pm1 {
                x = x * z.clone();
                x = x % p.clone();

                y = y * z.clone();
                y = y * z.clone();
                y = y % p.clone();
            }

            z = z.pow_mod(&two, &p).unwrap();
            exp = exp >> 1;
        }

        sqrt = x;
    }

    Ok(sqrt)
}

/// Return the random number between `a` and `b`.
/// The inputs are `BigInt`.
/// TODO: use generate_random_range without this method.
fn generate_random_range_big_int(a: BigInt, b: BigInt) -> BigInt {
    let mut rng = rand::thread_rng();
    rng.gen_bigint_range(&a, &b)
}

/// Return the length of bits of `a`.
pub fn big_bits_len(a: Integer) -> usize {
    unsafe { gmp::mpz_sizeinbase(a.as_raw(), 2) }
}

/// Return the random number between `a` and `b`.
pub fn generate_random_range(a: Integer, b: Integer) -> Integer {
    let mut a_digits = a.to_digits::<u8>(Order::MsfLe);
    a_digits.reverse();
    let a1 = BigInt::from_bytes_le(Sign::Plus, &a_digits);

    let mut b_digits = b.to_digits::<u8>(Order::MsfLe);
    b_digits.reverse();
    let b1 = BigInt::from_bytes_le(Sign::Plus, &b_digits);

    let p = generate_random_range_big_int(a1, b1);
    let digits = p.to_u64_digits().1;
    Integer::from_digits(&digits, Order::LsfLe)
}

/// Return the bytes of big integer.
pub fn big_to_bytes(a: Integer) -> Vec<u8> {
    let mut digits = a.to_digits::<u8>(Order::MsfLe);
    digits.reverse();

    // TODO: constant-time
    if digits.len() == 0 {
        digits.push(0);
    }
    
    digits
}

/// Return the random prime of bit length `bitsize`.
/// It filters out the primes which are not 3 mod 4 when `requires3mod4` = true.
pub fn generate_random_prime(bitsize: usize, requires3mod4: bool) -> Integer {
    let mut rng = rand::thread_rng();
    let mut p: BigUint;

    loop {
        p = rng.gen_prime_exact(bitsize, None);
        // let r = p.to_str_radix(16);
        // let b = Integer::from_str_radix(&r, 16).unwrap();

        let digits = p.to_u64_digits();
        let b = Integer::from_digits(&digits, Order::LsfLe);

        if !requires3mod4 || b.clone().div_rem(Integer::from(4)).1 == Integer::from(3) {
            return b;
        }
    }
}

const SMALL_PRIMES_COUNT: usize = 1006;
/// First 1006 primes as 32 bit unsigned integers.
pub static SMALL_PRIMES: [u32; SMALL_PRIMES_COUNT] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
    431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
    661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
    809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929,
    937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039,
    1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153,
    1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279,
    1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409,
    1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499,
    1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613,
    1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741,
    1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873,
    1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999,
    2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113,
    2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251,
    2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371,
    2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477,
    2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647,
    2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731,
    2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857,
    2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001,
    3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163,
    3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299,
    3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407,
    3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539,
    3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659,
    3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793,
    3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919,
    3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051,
    4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201,
    4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327,
    4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463,
    4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603,
    4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733,
    4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903,
    4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009,
    5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153,
    5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303,
    5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441,
    5443, 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569,
    5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701,
    5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839, 5843,
    5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987,
    6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131,
    6133, 6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269,
    6271, 6277, 6287, 6299, 6301, 6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373,
    6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553,
    6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691,
    6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829,
    6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949, 6959, 6961, 6967,
    6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109,
    7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247,
    7253, 7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451,
    7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559,
    7561, 7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687,
    7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841,
    7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919, 7927, 7933, 7937, 7949, 7951, 7963,
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sqrt_mod() {
        let p = rug::Integer::from(11); // p = 3 \pmod 4
        for a in 1..10 {
            if rug::Integer::from(a).jacobi(&p) != 1 {
                // not quadratic residue"
                continue;
            }
            let s = sqrt_mod_p(rug::Integer::from(a), p.clone()).unwrap();
            assert!(s.pow_mod(&rug::Integer::from(2), &p).unwrap() == rug::Integer::from(a));
        }

        let p = rug::Integer::from(13); // p = 5 \pmod 8
        for a in 1..12 {
            if rug::Integer::from(a).jacobi(&p) != 1 {
                // not quadratic residue"
                continue;
            }
            let s = sqrt_mod_p(rug::Integer::from(a), p.clone()).unwrap();
            assert!(s.pow_mod(&rug::Integer::from(2), &p).unwrap() == rug::Integer::from(a));
        }

        // Shanks-Tonelli:
        let p = rug::Integer::from(17);
        for a in 1..17 {
            if rug::Integer::from(a).jacobi(&p) != 1 {
                // not quadratic residue"
                continue;
            }
            let s = sqrt_mod_p(rug::Integer::from(a), p.clone()).unwrap();
            assert!(s.pow_mod(&rug::Integer::from(2), &p).unwrap() == rug::Integer::from(a));
        }

        let p = rug::Integer::from_str_radix(
            "61168504797541074430810428585729838043884986145671831191169392989232455728209",
            10,
        )
        .unwrap();
        let a = rug::Integer::from_str_radix(
            "44229201317184224679354873652471145530196604133128431668556891662306113797688",
            10,
        )
        .unwrap();
        let s = sqrt_mod_p(a.clone(), p.clone()).unwrap();
        assert!(s.pow_mod(&rug::Integer::from(2), &p).unwrap() == a);
    }
}
