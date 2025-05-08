# README

## Schemes

The repository provides primitives for quadratic fields and quaternion algebras, including the [KLPT algorithm](https://eprint.iacr.org/2014/505.pdf) and the [Cornacchia algorithm](https://en.wikipedia.org/wiki/Cornacchia%27s_algorithm).

It also includes implementations of two isogeny-based schemes:
 * [KLaPoTi: An asymptotically efficient isogeny group action from 2-dimensional isogenies](https://eprint.iacr.org/2024/1844)
 * LIT-SiGamal.

## Underlying libraries

This work uses:
 * [two-isogenies](https://github.com/ThetaIsogenies/two-isogenies) for computing (2, 2)-isogenies, implemented by [Jack](https://github.com/giacomopope).
 * Elliptic curve operations based on [finite field macros](https://github.com/pornin/crrl/blob/main/src/backend/w64/gfgen.rs) provided by [Thomas Pornin](https://github.com/pornin).
 * Sage scripts developed by [Tomoki](https://tomoriya.work).

## Computing the isogeny chain

This work uses the library
[two-isogenies](https://github.com/ThetaIsogenies/two-isogenies) for computing (2, 2)-isogenies.
To use its `product_isogeny` function, a strategy must first be precomputed using the `strategy.py` script.
Computing the optimal strategy for the isogeny chain involves the following two steps:

### Benchmark the operations (multiplication, square, invert)
 
These costs are computed using a benchmark `benches/msi.rs`, just run (make sure this benchmark is set in the Cargo.toml):
```
cargo bench
```

Example output:
```
New Multiplication (1745 bit)
                        time:   [2.9330 µs 2.9375 µs 2.9424 µs]
                        change: [-2.4460% -1.9078% -1.4550%] (p = 0.00 < 0.05)
                        Performance has improved.
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe

Square (1745 bit)       time:   [2.3517 µs 2.3552 µs 2.3589 µs]
                        change: [-1.8670% -1.6079% -1.3693%] (p = 0.00 < 0.05)
                        Performance has improved.
Found 13 outliers among 100 measurements (13.00%)
  6 (6.00%) low mild
  4 (4.00%) high mild
  3 (3.00%) high severe

Invert (1745 bit)       time:   [84.972 µs 85.162 µs 85.382 µs]
                        change: [-2.8793% -2.6136% -2.3365%] (p = 0.00 < 0.05)

```

As costs we take the middle values of the `time` list.

### Run the strategy script

Open `strategy.py` and set the `data` list (length of chain, cost of multiplication, cost of squaring, cost of inversion):

```
data = [1733, 2937, 2355, 85162]
strat = optimised_strategy(*data)
print(strat)
```

Run `strategy.py`.

Use the computed list in the code for the computation of the isogeny chain.

## Add a new finite field

Different schemes use different finite fields. To add a new field, use the `sage/finite_fields/gen_fp.sage` script, developed by Thomas Pornin. You will obtain an output such as:

```
pub mod FpLit192 {
    const N: usize = 19;
    const BITLEN: usize = 1163;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x3379544C3A34DAA7, 0x21B7935908E42B2B, 0x906B4DE5EB6D771F, 0x2CBE8C1918D4648A, 0xC406D98E3B5F726B, 0xAA8AF0C84CDDE099, 0x95C3EC55817B9EED, 0x9B48764AB5CC8E22, 0x70DD6317E5C58C3F, 0x00000000000006B8
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x99BCAA261D1A6D54, 0x90DBC9AC84721595, 0x4835A6F2F5B6BB8F, 0x965F460C8C6A3245, 0xE2036CC71DAFB935, 0xD5457864266EF04C, 0x4AE1F62AC0BDCF76, 0xCDA43B255AE64711, 0x386EB18BF2E2C61F, 0x000000000000035C
    ];
...
...
pub mod FpLit192Ext {
    use super::Fp1163::Fp;
    const NQR_RE: Fp = Fp::new([
        0x22988628721B7E83, 0x53D0AC2596D0BB38, 0xA8985368C042AFA7, 0xF0D632336068F736, 0xE84067B84DD4496E, 0x41CF80B4A799CAB0, 0x867E5E7842839AEE, 0x7D932E615A9BBB9D, 0xB74DCDC3057CD06F, 0x17B69EF2262D80A5, 0xC6C447EFBB1B3350, 0x47E013B0A12578B9, 0xAAC4F769D556E43B, 0x5B14A7F174EA04A8, 0x9CE7CE8C561E6141, 0xA67B8D57AE2D1EA5, 0x6E245B0887DF133C, 0x229385E491780B15, 0x00000000000004BA
    ]);
...
```

Rename the finite field and paste the code into `fields.rs`. Finally, paste something like below to `lib.rs`:

```
pub mod ec_lit192 {
    pub type Fp = crate::fields::FpLit192::Fp;
    pub type Fq = crate::fields::FpLit192Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::ec::ec_helpers::define_ec_helpers! {}
    
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::lit_sigamal::define_litsigamal! {}
}
```

To generate the torsion points, use the `sage/precompute.sage` script and paste the output to `precomputed.json`.