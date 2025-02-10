# README

This work uses:
 * [two-isogenies](https://github.com/ThetaIsogenies/two-isogenies) for the computation of (2, 2)-isogenies.
 * [KLaPoTi: An asymptotically efficient isogeny group action from 2-dimensional isogenies](https://eprint.iacr.org/2024/1844) which provides the Rust implementation for the quaternion algebra,
the [KLPT algorithm](https://eprint.iacr.org/2014/505.pdf),
and the [Cornacchia algorithm](https://en.wikipedia.org/wiki/Cornacchia%27s_algorithm).
 * Sage scripts from [LIT-SiGamal: An efficient isogeny-based PKE based on a LIT diagram](https://eprint.iacr.org/2024/521).
 * 3-isogenies code by [Giacomo Pope](https://github.com/giacomopope). Thanks Jack!
 * (2, 2)-isogenies and 3-isogenies code uses the [macros for finite fields](https://github.com/pornin/crrl/blob/main/src/backend/w64/gfgen.rs) provided by [Thomas Pornin](https://github.com/pornin).

## Computing the isogeny chain

This work uses the library
[two-isogenies](https://github.com/ThetaIsogenies/two-isogenies) for the computation of (2, 2)-isogenies.
To use its function `product_isogeny` a strategy has to be precomputed using the script `strategy.py`:
To compute the optimal strategy for the computation of the isogeny chain, the following two steps are needed.

## Benchmark the operations (multiplication, square, invert)

 
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