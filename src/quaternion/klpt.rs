use anyhow::{bail, Result};
use num_traits::Pow;
use rand::Rng;
use rug::Integer;

use crate::{
    error::{KLPTNotSuccessful, RepresentIntegerNotFound},
    linalg::{kernel::right_kernel_mod_prime, matrix::Matrix},
    quaternion::{
        cornacchia::cornacchia_extended,
        dim2::{enumerate_cvp_filter, ConditionParams},
        param::{
            KLPT_GAMMA_EXPONENT_CENTER_SHIFT, KLPT_GAMMA_EXPONENT_INTERVAL_SIZE,
            KLPT_KEYGEN_LENGTH, KLPT_KEYGEN_NUM_GAMMA_TRIAL, KLPT_REPRES_NUM_GAMMA_TRIAL,
            KLPT_SIGNING_NUMBER_STRONG_APPROX,
        },
        quaternion_algebra::{QuatAlg, QuatAlgEl},
        quaternion_order::standard_maximal_extremal_order,
    },
    util::{big_bits_len, generate_random_range, mod_inv, sqrt_mod_p, Big},
};

use super::{
    quaternion_ideal::QuaternionIdeal,
    quaternion_order::{standard_extremal_from, QuaternionOrder},
};

/// Find a primitive quaternion element gamma of norm `n_gamma` (divided by some integer to make the element primitive)
/// inside the standard maximal extremal order.
pub fn represent_integer(
    n_gamma: Integer,
    qa: QuatAlg,
    order: QuaternionOrder,
    bad_primes_prod: Integer,
) -> Result<(QuatAlgEl, Integer), anyhow::Error> {
    let prime_list = vec![5.big()];
    let how_many_primes = 1;

    // adjusting the norm of gamma (multiplied by 4 to find a solution in the full maximal order)
    let adjusted_n_gamma = n_gamma.clone() * 4.big();
    let (mut sq_bound, _) = adjusted_n_gamma.clone().div_rem(-qa.p.clone());

    // TODO
    if sq_bound < 1000.big() {
        sq_bound = 1000.big();
    }

    let bound = sq_bound.clone().sqrt();

    let mut x = 0.big();
    let mut y = 0.big();
    let mut z = 0.big();
    let mut t = 0.big();

    let mut counter = 0;

    while counter < KLPT_REPRES_NUM_GAMMA_TRIAL {
        counter += 1;

        // we start by sampling the third coordinate

        z = generate_random_range(1.big(), bound.clone());

        // then, we sample the second coordinate
        // computing the second bound in temp as sqrt( (adjusted_n_gamma - coeffs[2]²)/p )
        let mut cornacchia_target = z.clone() * z.clone();
        let mut temp = sq_bound.clone() - cornacchia_target.clone();
        temp = temp.sqrt();
        if temp == 0.big() {
            continue;
        }

        // sampling the fourth coordinate
        t = generate_random_range(1.big(), temp.clone());

        // compute cornacchia_target = adjusted_n_gamma - p * (z² + t²)
        temp = t.clone() * t.clone();
        cornacchia_target = cornacchia_target.clone() + temp;
        cornacchia_target = cornacchia_target * (-qa.p.clone());
        cornacchia_target = adjusted_n_gamma.clone() - cornacchia_target;

        (x, y) = cornacchia_extended(
            cornacchia_target.clone(),
            prime_list.clone(),
            how_many_primes,
            bad_primes_prod.clone(),
        );

        if !(x == 0.big() && y == 0.big()) {
            // check that we can divide by two at least once
            // we must have x = t mod  2 and y = z mod 2
            if x.clone() % 2 == t.clone() % 2 && y.clone() % 2 == z.clone() % 2 {
                break;
            }
        }
    }

    if x == 0.big() && y == 0.big() && z == 0.big() && t == 0.big() {
        bail!(RepresentIntegerNotFound())
    }

    let elem = QuatAlgEl::new(x, y, z, t, 1.big(), qa.clone());
    let order_basis = order;
    let (primitive, content) = elem.factor_in_order(order_basis.lattice.clone());
    let coeffs_new = order_basis.lattice.basis * primitive;

    let gamma = QuatAlgEl::new_from_matrix(coeffs_new, order_basis.lattice.denom, qa);
    let n_gamma_new = adjusted_n_gamma / (content.clone() * content);

    Ok((gamma, n_gamma_new))
}

/// Find two integers C[0], C[1] such that `gamma` * j * (C[0] + i C[1]]) * `delta` is in the eichler order ZZ + `ideal`.
/// If `gamma` norm is divisible by the norm of `ideal` (indicated by `is_divisible`),
/// then the result will be contained in J.
/// Failure should not be possible. It assumes that the ideal norm is prime.
pub fn solve_combi_eichler(
    ideal: QuaternionIdeal,
    gamma: QuatAlgEl,
    delta: QuatAlgEl,
    qa: QuatAlg,
    is_divisible: bool,
) -> (Integer, Integer) {
    let i = QuatAlgEl::new(0.big(), 1.big(), 0.big(), 0.big(), 1.big(), qa.clone());
    let j = QuatAlgEl::new(0.big(), 0.big(), 1.big(), 0.big(), 1.big(), qa.clone());
    let gamma_j = gamma.clone() * j;
    let mut gamma_ji_delta = gamma_j.clone() * i * delta.clone();
    let mut gamma_j_delta = gamma_j * delta.clone();

    gamma_ji_delta = gamma_ji_delta.normalize();
    gamma_j_delta = gamma_j_delta.normalize();

    let mut bas2 = QuatAlgEl::new_from_matrix(
        ideal.lattice.basis.get_col(2),
        ideal.lattice.denom.clone(),
        qa.clone(),
    );
    let mut bas3 = QuatAlgEl::new_from_matrix(
        ideal.lattice.basis.get_col(3),
        ideal.lattice.denom.clone(),
        qa.clone(),
    );

    gamma_j_delta = gamma_j_delta * ideal.lattice.denom.clone() * gamma_ji_delta.denom.clone();
    gamma_ji_delta = gamma_ji_delta * ideal.lattice.denom.clone() * gamma_j_delta.denom.clone();

    bas2 = bas2 * gamma_ji_delta.denom.clone() * gamma_j_delta.denom.clone();
    bas3 = bas3 * gamma_ji_delta.denom.clone() * gamma_j_delta.denom.clone();

    let ideal_norm = ideal.norm.clone();

    if is_divisible {
        let mut m: Matrix<Integer> = Matrix::zeros(4, 4);
        m.set_col(0, gamma_j_delta.coeffs());
        m.set_col(1, gamma_ji_delta.coeffs());
        m.set_col(2, bas2.coeffs());
        m.set_col(3, bas3.coeffs());
        let k = right_kernel_mod_prime(m, ideal_norm);

        (k[0].clone(), k[1].clone())
    } else {
        let gamma_delta = ideal.lattice.denom.clone()
            * gamma_ji_delta.denom.clone()
            * gamma_j_delta.denom.clone();
        let mut gamma_delta_mat = Matrix::zeros(4, 1);
        gamma_delta_mat[(0, 0)] = gamma_delta;

        let mut m: Matrix<Integer> = Matrix::zeros(4, 5);
        m.set_col(0, gamma_j_delta.coeffs());
        m.set_col(1, gamma_ji_delta.coeffs());
        m.set_col(2, gamma_delta_mat);
        m.set_col(3, bas2.coeffs());
        m.set_col(4, bas3.coeffs());
        let k = right_kernel_mod_prime(m, ideal_norm);

        (k[0].clone(), k[1].clone())
    }
}

/// Decide if a given vector is suitable for the signing KLPT algorithm.
/// The algorithm tests whether `gamma` * `mu` is primitive in the end.
///
/// Arguments:
/// `vec`: the vector to be tested.
fn condition_keygen_klpt(vec: Matrix<Integer>, params: ConditionParams) -> (bool, QuatAlgEl) {
    let mut mu = QuatAlgEl::new(
        0.big(),
        0.big(),
        vec[(0, 0)].clone(),
        vec[(1, 0)].clone(),
        1.big(),
        params.qa.clone(),
    );
    let norm = mu.reduced_norm();
    let norm = norm.numer();
    let mut rhs = params.target_norm.clone() - norm;
    let mut tmp = params.n.clone() * params.n.clone();
    (rhs, tmp) = rhs.div_rem(tmp);
    assert!(tmp == 0.big());

    let prime_list_1mod4 = vec![
        5.big(),
        13.big(),
        17.big(),
        37.big(),
        41.big(),
        53.big(),
        61.big(),
        73.big(),
        89.big(),
        97.big(),
    ];
    let bad_primes_prod = "11638895555051853627".big();

    let (a, b) = cornacchia_extended(
        rhs.clone(),
        prime_list_1mod4.clone(),
        prime_list_1mod4.len(),
        bad_primes_prod,
    );

    if !(a == 0.big() && b == 0.big()) {
        let quat_temp = QuatAlgEl::new(
            a * params.n.clone(),
            b * params.n.clone(),
            0.big(),
            0.big(),
            1.big(),
            params.qa.clone(),
        );
        let temp = params.target_norm.clone() - norm;
        let norm = quat_temp.reduced_norm();
        let norm = norm.numer().clone();
        assert!(norm == temp);

        mu = mu + quat_temp;
        mu = mu / 2.big();

        let mut coord = Matrix::zeros(4, 1);
        let found = params.order.order.lattice.contains(&mut coord, mu.clone());
        if found {
            // NOTE: strong approximation is written with representation x + y * i + z * j + t * j * i,
            // while here we have x + y * i + z * j + t * i * j.
            // We need to change the sign of `mu.t`:
            mu.t = -mu.t;

            return (true, mu);
        }
    }

    (false, mu)
}

/// Compute the strong approximation.
/// This algorithm finds a quaternion element mu of norm equal to `n_mu`/4 such that
/// `mu` = lambda* (c0 + i * c1)*j  mod ( n * order).
/// The value of lambda has been computed in a way to make this computation possible
/// The value of mu is computed by the function condition (which depends on some additional params) and
/// must satisfy a set of constraints defined by condition.
/// The number of trials is max_tries.
/// Assumes n is either a prime or a product of big prime,
/// this means the probability to encounter non-invertible elements is considered negligible.
///
/// Arguments:
/// n_mu: 4 times target norm of mu.
/// order: special extremal order.
/// c0, c1: coefficients to be lifted.
/// norm: norm of the quaternion element corresponding to c0, c1.
/// lambda: an integer used in the computation (target vector = lambda * c0, lambda * c1 + ...).
/// n: modulo for the lift (ideal norm).
/// max_tries int, bound on the number of tries before we abort
/// condition a filter function returning whether a vector should be output or not
/// params extra parameters passed to `condition`. May be NULL.
/// qa: the quaternion algebra.
fn strong_approximation(
    n_mu: Integer,
    c0: Integer,
    c1: Integer,
    norm: Integer,
    lambda: Integer,
    n: Integer,
    max_tries: u32,
    params: ConditionParams,
) -> Result<QuatAlgEl> {
    let p_lambda_2 = -params.qa.p.clone() * lambda.clone() * 2.big();
    let coeff_c0 = (p_lambda_2.clone() * c0.clone()) % n.clone();
    let mut coeff_c1 = (p_lambda_2 * c1.clone() * params.order.q.clone()) % n.clone();
    coeff_c1 = mod_inv(coeff_c1, n.clone()).unwrap();
    let (mut cst_term, tmp) =
        (n_mu.clone() - lambda.clone() * lambda.clone() * norm).div_rem(n.clone());
    cst_term = cst_term % n.clone();

    assert!(tmp == 0.big());
    cst_term = cst_term % n.clone();
    if cst_term < 0.big() {
        cst_term += n.clone();
    }

    // computation of the lattice
    let mut lat = Matrix::zeros(2, 2);

    // first_column = n * ( 1, - coeff_c0 * coeff_c1 mod n)
    // we start by inverting coeff_c1
    lat[(0, 0)] = n.clone();
    lat[(1, 0)] = (-coeff_c1.clone() * coeff_c0) % n.clone();
    if lat[(1, 0)] < 0.big() {
        lat[(1, 0)] += n.clone();
    }
    lat[(1, 0)] = lat[(1, 0)].clone() * n.clone();

    // second_colum = ( 0, n² )
    lat[(0, 1)] = 0.big();
    lat[(1, 1)] = n.clone() * n.clone();

    // computation of the target vector = (lambda * C[0], lambda * C[1] + n * cst_term * coeff_c1 )
    let mut target = Matrix::zeros(2, 1);
    target[(0, 0)] = lambda.clone() * c0;
    target[(1, 0)] = lambda * c1 + ((cst_term * coeff_c1) % n.clone()) * n.clone();

    // computation of the norm log_bound
    let (bound, _) = n_mu.div_rem(params.qa.p.clone());

    let mut dist_bound = big_bits_len(bound);
    let n_bits = big_bits_len(n);

    // TODO
    if dist_bound == 0 {
        dist_bound = 3 * n_bits + 10;
    }

    // to avoid that dist_bound is too big
    if dist_bound > 3 * n_bits + 10 {
        dist_bound = 3 * n_bits + 10;
    }

    enumerate_cvp_filter(
        lat,
        target,
        params.order.q.clone(),
        params,
        Box::new(condition_keygen_klpt),
        dist_bound as u32,
        max_tries,
    )
}

/// SQISign keygen KLPT.
pub fn keygen_klpt(ideal: QuaternionIdeal, qa: QuatAlg) -> Result<QuatAlgEl> {
    // initializing the center of the sample interval as log(p) - log( n(lideal) ) + KLPT_gamma_exponent_center_shift

    let p_bits = big_bits_len(qa.p.clone());
    let i_bits = big_bits_len(ideal.norm.clone());
    let center = p_bits as u32 - i_bits as u32 + KLPT_GAMMA_EXPONENT_CENTER_SHIFT;

    let mut cnt = 0;
    while cnt < KLPT_KEYGEN_NUM_GAMMA_TRIAL {
        cnt += 1;
        // we start by choosing the exponent of gamma at random inside some interval (depending on the choice of parameters this interval might actually be of size 1)
        // sampling the exponent
        let k = rand::thread_rng().gen_range(
            center - KLPT_GAMMA_EXPONENT_INTERVAL_SIZE..=center + KLPT_GAMMA_EXPONENT_INTERVAL_SIZE,
        );

        let order = standard_maximal_extremal_order().order;
        let (ok, el) = klpt(ideal.clone(), qa.clone(), order, k, KLPT_KEYGEN_LENGTH);
        if ok {
            return Ok(el);
        }
    }

    bail!(KLPTNotSuccessful())
}

/// Return generator of the equivalent ideal J ∼ `ideal` of norm 2^keygen_length.
pub fn klpt(
    ideal: QuaternionIdeal,
    qa: QuatAlg,
    order: QuaternionOrder,
    k: u32,
    keygen_length: u32,
) -> (bool, QuatAlgEl) {
    let n_gamma = 2.big().pow(k) * ideal.norm.clone();

    let (gamma, _) = match represent_integer(n_gamma.clone(), qa.clone(), order.clone(), 1.big()) {
        Ok((gamma, gamma_norm)) => (gamma, gamma_norm),
        Err(_) => return (false, QuatAlgEl::zero(qa.clone())),
    };

    // computing the norm of the remaining part
    // norm_rem = 2^keygen_length * ideal.norm
    let norm_rem = 2.big().pow(keygen_length) * ideal.norm.clone();
    let (mut n_mu, _) = norm_rem.div_rem(n_gamma);

    // computing the linear combination mod n(ideal)
    let (c0, c1) = solve_combi_eichler(
        ideal.clone(),
        gamma.clone(),
        QuatAlgEl::one(qa.clone()),
        qa.clone(),
        true,
    );

    let el = QuatAlgEl::new(
        0.big(),
        0.big(),
        c0.clone(),
        c1.clone(),
        1.big(),
        qa.clone(),
    );
    let norm = el.reduced_norm();
    let norm = norm.numer();

    let mut tmp = match mod_inv(norm.clone(), ideal.norm.clone()) {
        Ok(tmp) => tmp,
        Err(_) => return (false, QuatAlgEl::zero(qa.clone())),
    };
    tmp = n_mu.clone() * tmp;

    let mut lambda = match sqrt_mod_p(tmp, ideal.norm.clone()) {
        Ok(lambda) => lambda,
        Err(_) => return (false, QuatAlgEl::zero(qa.clone())),
    };

    // multiply lambda by 2 and the final norm by 4
    n_mu = n_mu * 4.big();
    lambda = lambda * 2.big();

    let params = ConditionParams::new(
        n_mu.clone(),
        qa.clone(),
        standard_extremal_from(order),
        ideal.norm.clone(),
    );

    match strong_approximation(
        n_mu.clone(),
        c0,
        c1,
        norm.clone(),
        lambda,
        ideal.norm.clone(),
        KLPT_SIGNING_NUMBER_STRONG_APPROX,
        params,
    ) {
        Ok(mu) => {
            let n = mu.reduced_norm();
            let n = n.numer() * 4.big();
            assert!(n == n_mu);

            return (true, gamma * mu);
        }
        Err(_) => return (false, QuatAlgEl::zero(qa.clone())),
    };
}

#[cfg(test)]
mod tests {
    use crate::{
        quaternion::{lattice::Lattice, quaternion_ideal::QuaternionIdeal},
        util::{generate_random_prime, Big},
    };

    use super::*;

    #[test]
    fn solve_combi_eichler_does_not_divide_ideal_norm_test() {
        let p =
            "207382631875713544995127616290604630368466823479170154862332108637983515463219".big();
        let qa = QuatAlg::new(-p);
        let order = standard_maximal_extremal_order().order;

        let data = vec![
            "907954721436953965751684344669123774054".big(),
            0.big(),
            "185132649634685922350761388101184156624".big(),
            "83772611778951282667046924681891206665".big(),
            0.big(),
            "907954721436953965751684344669123774054".big(),
            "824182109658002683084637419987232567389".big(),
            "185132649634685922350761388101184156624".big(),
            0.big(),
            0.big(),
            1.big(),
            0.big(),
            0.big(),
            0.big(),
            0.big(),
            1.big(),
        ];
        let basis = Matrix::create_from_data(data, 4, 4);
        let lattice = Lattice::new(basis, 2.big());
        let norm = "453977360718476982875842172334561887027".big();
        let ideal = QuaternionIdeal::new(lattice, norm, Some(order.clone()));

        let shift_gamma = QuatAlgEl::new(
            "5580021485201764633877987597567810067670013693047407665985".big(),
            "2320547272845645173275609279354266959468621448870911538608".big(),
            "4198880840257667550".big(),
            -"23795982715478488575".big(),
            2.big(),
            qa.clone(),
        );
        let delta = QuatAlgEl::new(
            "3602374210419591670642144183339784121997032109771738234825".big(),
            "307266378194183534037377607187105689418803506859320686982".big(),
            "7008935480242594896".big(),
            -"25452515529334439587".big(),
            2.big(),
            qa.clone(),
        );

        let (c0, c1) = solve_combi_eichler(ideal, shift_gamma, delta, qa, false);

        assert!(c0 == "300237822977117692966734613589639386254".big());
        assert!(c1 == "323781918246312966715393479975786929730".big());
    }

    #[test]
    fn solve_combi_eichler_divide_ideal_norm_test() {
        let p =
            "138883194693555113851393415309153156354060925518152854193510459328759979281071".big();
        let qa = QuatAlg::new(-p);
        let order = standard_maximal_extremal_order().order;

        let data = vec![
            "960472930926720839036858276645084775238".big(),
            0.big(),
            "427027823721248968964167152616035659878".big(),
            "34502334730772480448609784473656523425".big(),
            0.big(),
            "960472930926720839036858276645084775238".big(),
            "925970596195948358588248492171428251813".big(),
            "427027823721248968964167152616035659878".big(),
            0.big(),
            0.big(),
            1.big(),
            0.big(),
            0.big(),
            0.big(),
            0.big(),
            1.big(),
        ];
        let basis = Matrix::create_from_data(data, 4, 4);
        let lattice = Lattice::new(basis, 2.big());
        let norm = "480236465463360419518429138322542387619".big();
        let ideal = QuaternionIdeal::new(lattice, norm, Some(order.clone()));

        let gamma = QuatAlgEl::new(
            "392059506136331759216329497687781599770403611042830646".big(),
            "20875836824227216931066173115268039512153844545157315".big(),
            "1912719849983461".big(),
            -"1094483383615038".big(),
            2.big(),
            qa.clone(),
        );
        let delta = QuatAlgEl::new(1.big(), 0.big(), 0.big(), 0.big(), 1.big(), qa.clone());

        let (c0, c1) = solve_combi_eichler(ideal, gamma, delta, qa, false);

        assert!(c0 == "35231624699750769236970466201175770119".big());
        assert!(c1 == "172626538854885491913317077688770512365".big());
    }

    #[test]
    fn solve_combi_eichler_test() {
        let exp = 128;
        let p = generate_random_prime(2 * exp, true);

        let qa = QuatAlg::new(-p);
        let m = generate_random_prime(exp, true);
        let mut temp = generate_random_prime(exp + 100, true);
        temp = m.clone() * temp;

        let order = standard_maximal_extremal_order().order;
        let (gamma1, _) = represent_integer(temp, qa.clone(), order.clone(), 1.big()).unwrap();

        temp = 2.big().pow(exp as u32 * 3);
        let (shift_gamma, temp) =
            represent_integer(temp, qa.clone(), order.clone(), 1.big()).unwrap();

        let (delta, _) = represent_integer(temp, qa.clone(), order.clone(), 1.big()).unwrap();

        let ideal = QuaternionIdeal::new_left_ideal(gamma1, m, order.clone(), qa.clone());

        solve_combi_eichler(ideal, shift_gamma, delta, qa, false);
    }

    #[test]
    fn keygen_klpt_test() {
        let exp = 128;
        let p = generate_random_prime(2 * exp, true);
        let qa = QuatAlg::new(-p.clone());
        let n_tau = generate_random_prime(exp / 2, true);
        let mut tmp = generate_random_prime(3 * exp, true);
        tmp = tmp * n_tau.clone();
        let order = standard_maximal_extremal_order().order;
        let (gamma_tau, _) = represent_integer(tmp, qa.clone(), order.clone(), 1.big()).unwrap();

        let ideal = QuaternionIdeal::new_left_ideal(
            gamma_tau.clone(),
            n_tau.clone(),
            order.clone(),
            qa.clone(),
        );

        let _ = keygen_klpt(ideal, qa.clone());
    }
}
