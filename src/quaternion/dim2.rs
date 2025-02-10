use anyhow::{bail, Result};
use num_traits::Pow;
use rug::Integer;

use crate::{
    error::ShortDim2VecNotFound,
    linalg::matrix::Matrix,
    quaternion::quaternion_algebra::QuatAlg,
    util::{rounded_div, Big},
};

use super::{
    quaternion_algebra::QuatAlgEl,
    quaternion_order::{standard_maximal_extremal_order, ExtremalMaximalOrder},
};

/// Parameters needed to decide whether the short vector found by KLPT is ok (for example,
/// whether it is in the order).
#[derive(Clone, Debug)]
pub struct ConditionParams {
    /// KLPT target norm
    pub target_norm: Integer,
    /// Quaternion algebra
    pub qa: QuatAlg,
    /// Quaternion order
    pub order: ExtremalMaximalOrder,
    /// KLPT input ideal norm
    pub n: Integer,
}

impl ConditionParams {
    /// Create new `ConditionParams`.
    pub fn new(target_norm: Integer, qa: QuatAlg, order: ExtremalMaximalOrder, n: Integer) -> Self {
        ConditionParams {
            target_norm,
            qa,
            order,
            n,
        }
    }

    /// Create `ConditionParams` using quaterion algebra `qa` and prime `p`.
    /// Used only for testing.
    pub fn new_with_p(qa: QuatAlg, p: Integer) -> Self {
        let order = standard_maximal_extremal_order();
        let target_norm = 0.big();
        ConditionParams {
            target_norm,
            qa,
            order,
            n: p,
        }
    }
}

/// Find a lattice vector close to a target vector satisfying extra conditions.
/// Given a target vector `target` = (x₀,y₀), enumerate vectors `v` = (x,y) in lattice `lat` with distance
/// (x₀ - x)² + `qf` (y₀ - y)² < 2^`dist_bound`. On each vector `condition(res, target - v, params)` is called:
/// if it returns a non-zero value, processing stops and the same value is returned;
/// if it returns 0 processing continues until `max_tries` vectors have been tried.
///
/// The implementation will first reduce the basis by finding a short vector with algorithm 1.3.14 (Gauss)
/// from Henri Cohen's "A Course in Computational Algebraic Number Theory" (Springer Verlag, in series "Graduate texts in Mathematics") from 1993.
/// Then a second one is added to this basis after reduction by projection on the first one.
/// A close vector is found using 2 projection as in https://cims.nyu.edu/~regev/teaching/lattices_fall_2004/ln/cvp.pdf (15.5.2023,16h15CEST) using the already reduced matrix.
/// Finally, short vectors are enumerated below an enumeration bound set to
/// (2 ^ `dist_bound` - (distance of target to close vector)).
/// Each short vector is added to the close vector, and then the distance to the target is measured.
/// If it is below 2^dist_bound, it is tested for condition.
/// If coundition returns 1, the algorithm terminates and returns 1, otherwise it returns 0
/// when max_tries vectors have been tried.
///
/// The enumeration bound is a heuristic to avoid starting the search with vectors which will be
/// too far from the target, since enumeration starts at the largest vectors.
/// It can in some cases lead to failures even though solutions do exist, an in other cases be
/// insufficient to get a valid result in reasonable time.
///
/// Enumeration uses algorithm 2.7.5 (Fincke-Pohst) from
/// Henri Cohen's "A Course in Computational Algebraic Number Theory" (Springer Verlag, in series "Graduate texts in Mathematics") from 1993
/// Slightly adapted to work without rational numbers and their square roots.
/// Therefore needing a test to make sure the bounds are respected, which is integrated in
/// test_cvp_condition.
pub fn enumerate_cvp_filter(
    lattice_basis: Matrix<Integer>,
    target: Matrix<Integer>,
    qf: Integer,
    params: ConditionParams,
    condition: Box<dyn Fn(Matrix<Integer>, ConditionParams) -> (bool, QuatAlgEl)>,
    dist_bound: u32,
    max_tries: u32,
) -> Result<QuatAlgEl> {
    let reduced_basis = short_basis(lattice_basis.clone(), qf.clone());
    let (target_minus_closest, _) =
        closest_vector(reduced_basis.clone(), target.clone(), qf.clone());
    let norm_bound = 2.big().pow(dist_bound);

    enumerate_short_vec(
        target_minus_closest,
        reduced_basis,
        qf,
        params,
        condition,
        norm_bound,
        max_tries,
    )
}

fn norm(a: Matrix<Integer>, qf: Integer) -> Integer {
    a[(0, 0)].clone() * a[(0, 0)].clone() + qf * a[(1, 0)].clone() * a[(1, 0)].clone()
}

fn lattice_bilinear(a: Matrix<Integer>, b: Matrix<Integer>, qf: Integer) -> Integer {
    a[(0, 0)].clone() * b[(0, 0)].clone() + qf * a[(1, 0)].clone() * b[(1, 0)].clone()
}

/// Compute exact solution for shortest vector in dimension 2, than take a second, orthogonal vector
/// (algorithm 3.1.14 Cohen)
fn short_basis(lattice_basis: Matrix<Integer>, qf: Integer) -> Matrix<Integer> {
    let mut a = lattice_basis.get_col(0);
    let mut b = lattice_basis.get_col(1);

    let mut norm_a = norm(a.clone(), qf.clone());
    let mut norm_b = norm(b.clone(), qf.clone());

    // exchange if needed
    if norm_a < norm_b {
        let tmp = a.clone();
        a = b;
        b = tmp;
        let tmp = norm_a.clone();
        norm_a = norm_b;
        norm_b = tmp;
    }

    let mut r;
    let mut norm_t;

    loop {
        let n = lattice_bilinear(a.clone(), b.clone(), qf.clone());
        r = rounded_div(n.clone(), norm_b.clone());

        norm_t =
            norm_a.clone() - 2 * n.clone() * r.clone() + r.clone() * r.clone() * norm_b.clone();

        if norm_b > norm_t {
            norm_a = norm_b;
            norm_b = norm_t;

            let mut t: Matrix<Integer> = Matrix::zeros(2, 1);
            t[(0, 0)] = a[(0, 0)].clone() - r.clone() * b[(0, 0)].clone();
            t[(1, 0)] = a[(1, 0)].clone() - r.clone() * b[(1, 0)].clone();
            a = b;
            b = t;
        } else {
            break;
        }
    }
    // output : now b is short: need to get 2nd short vector: idea: take shortest among t and a
    if norm_t < norm_a {
        a[(0, 0)] = a[(0, 0)].clone() - r.clone() * b[(0, 0)].clone();
        a[(1, 0)] = a[(1, 0)].clone() - r.clone() * b[(1, 0)].clone();
    }

    let mut reduced = Matrix::zeros(2, 2);
    reduced[(0, 0)] = b[(0, 0)].clone();
    reduced[(1, 0)] = b[(1, 0)].clone();
    reduced[(0, 1)] = a[(0, 0)].clone();
    reduced[(1, 1)] = a[(1, 0)].clone();

    reduced
}

/// Compute the rounded value of <a*,t>/<a*,a*>, where a* is a orthogonalised with respect to b.
// fn get_coefficient_with_orthogonalisation(const ibz_t *a0,const ibz_t *a1,const ibz_t *b0,const ibz_t *b1,const ibz_t *t0,const ibz_t *t1,const ibz_t *norm_q){
fn get_coefficient_with_orthogonalisation(
    reduced_basis: Matrix<Integer>,
    target: Matrix<Integer>,
    qf: Integer,
) -> Integer {
    let a = reduced_basis.get_col(1);
    let b = reduced_basis.get_col(0);
    let norm_b = norm(b.clone(), qf.clone());
    let mut bilinear = lattice_bilinear(a.clone(), b.clone(), qf.clone());
    let a_star = a * norm_b.clone() - b * bilinear;
    let norm_a_star = norm(a_star.clone(), qf.clone());
    bilinear = lattice_bilinear(a_star, target, qf) * norm_b;

    rounded_div(bilinear, norm_a_star)
}

/// Test version of the condition argument in the cvp enumeration algorithms.
/// Sets elem[0] and elem[2] to vec[0], elem[1] and elem[3] to vec[1] if vec[0] + vec[1] mod p is 2.
/// This defines a quadratic form in dimension 2 where coord1 is the first and coord2 the second coordinate of the vector on which it is evaluated.
pub fn test_cvp_condition(vec: Matrix<Integer>, params: ConditionParams) -> (bool, QuatAlgEl) {
    let mut sum = (vec[(0, 0)].clone() + vec[(1, 0)].clone()) % params.n.clone();
    if sum < 0.big() {
        sum += params.n;
    }
    let el = QuatAlgEl::new(
        vec[(0, 0)].clone(),
        vec[(1, 0)].clone(),
        vec[(0, 0)].clone(),
        vec[(1, 0)].clone(),
        2.big(),
        params.qa.clone(),
    );

    if sum == 2.big() {
        return (true, el);
    }

    (false, el)
}

fn qf_value_bound_generation(
    num_a: Integer,
    denom_a: Integer,
    num_b: Integer,
    denom_b: Integer,
) -> Integer {
    if denom_a <= 0.big() || denom_b == 0.big() {
        return 0.big();
    }
    let mut sqrt_num = num_a.sqrt() + 1.big();
    let sqrt_denom = denom_a.sqrt();
    let common_denom = denom_b.clone() * sqrt_denom.clone();
    sqrt_num = sqrt_num * denom_b + sqrt_denom * num_b;
    let (d, _) = sqrt_num.div_rem(common_denom);

    d + 1.big()
}

/// Return short vectors.
/// Uses algorithm 2.7.5 (Fincke-Pohst) from Henri Cohen's "A Course in Computational Algebraic Number Theory" (Springer Verlag, in series "Graduate texts in Mathematics") from 1993
/// Slightly adapted to work without rational numbers and their square roots.
/// Therefore a test to make sure the bounds are respected is needed.
fn enumerate_short_vec(
    target_minus_closest: Matrix<Integer>,
    reduced_basis: Matrix<Integer>,
    qf: Integer,
    params: ConditionParams,
    condition: Box<dyn Fn(Matrix<Integer>, ConditionParams) -> (bool, QuatAlgEl)>,
    norm_bound: Integer,
    max_tries: u32,
) -> Result<QuatAlgEl> {
    // Subtract norm of distance from bound to not enumerate too large vectors at beginning,
    // this is an heuristic which can fail, so in case it is really bad, we do not do it.
    let mut norm_bound_for_enumeration =
        norm_bound.clone() - norm(target_minus_closest.clone(), qf.clone());
    if norm_bound_for_enumeration <= 0.big() {
        norm_bound_for_enumeration = norm_bound.clone();
    }

    // Set qf_a, qf_b, qf_c such that for x,y, a vector represented in lattice basis,
    // ax^2 + bxy + cy^2 is the norm defined by qf
    let qf_b = 2.big()
        * lattice_bilinear(
            reduced_basis.get_col(0),
            reduced_basis.get_col(1),
            qf.clone(),
        );
    let qf_a = norm(reduced_basis.get_col(0), qf.clone());
    let qf_c = norm(reduced_basis.get_col(1), qf.clone());

    // discriminant computation
    let mut prod = qf_b.clone() * qf_b.clone();
    let disc = 4.big() * qf_a.clone() * qf_c.clone() - prod.clone();

    if disc == 0.big() {
        bail!(ShortDim2VecNotFound());
    } else {
        let two_a = 2.big() * qf_a.clone();
        let four_a2 = two_a.clone() * two_a.clone();
        let four_a2_norm_bound = four_a2.clone() * norm_bound_for_enumeration;
        let four_a3 = four_a2.clone() * qf_a;
        let mut four_a2_c_minus_b2 = prod; // equals qf_b^2
        prod = four_a2 * qf_c;
        four_a2_c_minus_b2 = prod - four_a2_c_minus_b2;

        let bound_y = qf_value_bound_generation(
            four_a2_norm_bound.clone(),
            four_a2_c_minus_b2.clone(),
            0.big(),
            1.big(),
        );

        let mut y = -bound_y.clone() - 1.big();
        let mut stop = false;
        let mut tries = 0;

        while !stop && y.clone() < bound_y && tries < max_tries {
            y = y + 1.big();
            let y2 = y.clone() * y.clone();
            prod = qf_b.clone() * y.clone();
            let mut var = -prod;
            prod = y2 * four_a2_c_minus_b2.clone();
            prod = prod + four_a2_norm_bound.clone(); //prod = 4a^2*norm_bound + (4ca^2 - b^2)y^2
            let bound_x = qf_value_bound_generation(
                prod.clone(),
                four_a3.clone(),
                var.clone(),
                two_a.clone(),
            );
            var = -var;
            let mut x =
                qf_value_bound_generation(prod.clone(), four_a3.clone(), var, two_a.clone());
            x = -x - 1.big();

            while !stop && x < bound_x && tries < max_tries {
                tries += 1;
                x = x + 1.big();

                // put x, y from lattice basis in canonical basis
                let mut proposal = Matrix::create_from_data(vec![x.clone(), y.clone()], 2, 1);
                proposal = reduced_basis.clone() * proposal;
                let sum = target_minus_closest.clone() - proposal;
                let norm = norm(sum.clone(), qf.clone());
                if norm <= norm_bound {
                    let (ok, el) = condition(sum.clone(), params.clone());
                    if ok {
                        return Ok(el);
                    }
                }

                stop = x == 0.big() && y == 0.big();
            }
        }
    }

    bail!(ShortDim2VecNotFound());
}

/// Find the closest vector to target in lattice, using a reduced basis given as argument basis.
/// Nearest plane algorithm as in https://cims.nyu.edu/~regev/teaching/lattices_fall_2004/ln/cvp.pdf, but without basis, just a projection.
fn closest_vector(
    reduced_basis: Matrix<Integer>,
    target: Matrix<Integer>,
    qf: Integer,
) -> (Matrix<Integer>, Matrix<Integer>) {
    let norm_a = norm(reduced_basis.get_col(0), qf.clone());

    // Use 2nd basis vector (the larger one), and orthogonalise it with respect to the first one:
    let coord1 =
        get_coefficient_with_orthogonalisation(reduced_basis.clone(), target.clone(), qf.clone());

    // Subtract projection from vector
    let mut target_minus_closest = target.clone();
    target_minus_closest = target_minus_closest - reduced_basis.get_col(1) * coord1.clone();

    // Use 1st basis vector (the smaller one)
    let mut coord0 = lattice_bilinear(
        target_minus_closest.clone(),
        reduced_basis.get_col(0),
        qf.clone(),
    );
    coord0 = rounded_div(coord0.clone(), norm_a);

    // Subtract projection from vector
    target_minus_closest = target_minus_closest - reduced_basis.get_col(0) * coord0.clone();
    let closest_coords_in_basis = Matrix::create_from_data(vec![coord0, coord1], 2, 1);

    (target_minus_closest, closest_coords_in_basis)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{linalg::matrix::Matrix, util::Big};

    #[test]
    fn short_basis_test() {
        let mut basis = Matrix::zeros(2, 2);
        basis[(0, 0)] = -3.big();
        basis[(0, 1)] = 2.big();
        basis[(1, 0)] = 1.big();
        basis[(1, 1)] = 2.big();
        let qf = 2.big();
        let reduced = short_basis(basis.clone(), qf);
        let mut cmp = Matrix::zeros(2, 2);
        cmp[(0, 0)] = -3.big();
        cmp[(0, 1)] = 2.big();
        cmp[(1, 0)] = 1.big();
        cmp[(1, 1)] = 2.big();
        assert!(reduced == cmp);

        basis[(0, 0)] = 3.big();
        basis[(0, 1)] = 5.big();
        basis[(1, 0)] = 7.big();
        basis[(1, 1)] = 9.big();
        let qf = 2.big();
        let reduced = short_basis(basis.clone(), qf);
        assert!(reduced == cmp);

        basis[(0, 0)] = 48.big();
        basis[(0, 1)] = 4.big();
        basis[(1, 0)] = 81.big();
        basis[(1, 1)] = 9.big();
        let qf = 9.big();
        let reduced = short_basis(basis.clone(), qf);
        cmp[(0, 0)] = 12.big();
        cmp[(0, 1)] = 4.big();
        cmp[(1, 0)] = 0.big();
        cmp[(1, 1)] = 9.big();
        assert!(reduced == cmp);

        basis[(0, 0)] = 364.big();
        basis[(0, 1)] = 0.big();
        basis[(1, 0)] = 1323546.big();
        basis[(1, 1)] = 266606.big();
        let qf = 12165.big();
        let reduced = short_basis(basis.clone(), qf);
        cmp[(0, 0)] = 92092.big();
        cmp[(0, 1)] = 10192.big();
        cmp[(1, 0)] = 2.big();
        cmp[(1, 1)] = 1054.big();
        assert!(reduced == cmp);
    }

    #[test]
    fn enumerate_cvp_filter_test() {
        let mut basis = Matrix::zeros(2, 2);
        basis[(0, 0)] = -3.big();
        basis[(0, 1)] = 2.big();
        basis[(1, 0)] = 1.big();
        basis[(1, 1)] = 2.big();

        let mut target = Matrix::zeros(2, 1);
        target[(0, 0)] = -3.big();
        target[(1, 0)] = 3.big();

        let p = 3.big();
        let mut qf = 2.big();
        let mut dist_bound = 10;
        let mut max_tries = 20000;

        let qa = QuatAlg::new(-p.clone());

        let params = ConditionParams::new_with_p(qa.clone(), p.clone());

        let cmp = QuatAlgEl::new(4.big(), 22.big(), 4.big(), 22.big(), 2.big(), qa.clone());
        let el = enumerate_cvp_filter(
            basis.clone(),
            target.clone(),
            qf,
            params.clone(),
            Box::new(test_cvp_condition),
            dist_bound,
            max_tries,
        );
        assert!(el.unwrap() == cmp);

        target[(0, 0)] = 1216765.big();
        target[(1, 0)] = 879777.big();

        basis[(0, 0)] = 364.big();
        basis[(0, 1)] = 0.big();
        basis[(1, 0)] = 1323546.big();
        basis[(1, 1)] = 266606.big();

        qf = 12165.big();
        dist_bound = 32;
        max_tries = 100;

        let cmp = QuatAlgEl::new(
            -18287.big(),
            -155.big(),
            -18287.big(),
            -155.big(),
            2.big(),
            qa.clone(),
        );
        let el = enumerate_cvp_filter(
            basis.clone(),
            target.clone(),
            qf,
            params.clone(),
            Box::new(test_cvp_condition),
            dist_bound,
            max_tries,
        );
        assert!(el.unwrap() == cmp);

        target[(0, 0)] = 1216765.big();
        target[(1, 0)] = -223409.big();

        basis[(0, 0)] = 364.big();
        basis[(0, 1)] = 2.big();
        basis[(1, 0)] = -132346.big();
        basis[(1, 1)] = 26606.big();

        qf = 12.big();
        dist_bound = 25;
        max_tries = 50;

        let cmp = QuatAlgEl::new(
            1449.big(),
            1343.big(),
            1449.big(),
            1343.big(),
            2.big(),
            qa.clone(),
        );
        let el = enumerate_cvp_filter(
            basis.clone(),
            target.clone(),
            qf,
            params.clone(),
            Box::new(test_cvp_condition),
            dist_bound,
            max_tries,
        );
        assert!(el.unwrap() == cmp);

        target[(0, 0)] = "24331603791763515691729329496192502951056898161849272968170597636916968011043891967931576902713229532952543949167870260".big();
        target[(1, 0)] = "320864169455876864410246920205472260746024223953596490297261861428291460634916982031632479909031611896628831462307645903".big();

        basis[(0, 0)] = "83823929950128301594915602697228677079260482500022964929273".big();
        basis[(0, 1)] = 0.big();
        basis[(1, 0)] = "26601400243014239200180507899382755112580169588401175904342542265512938583814884922817611265843666255091702318534544478".big();
        basis[(1, 1)] = "702645123228401649030937397460831326065409163240113467483219315690011711473572222698591407715178198721852976513892308529".big();

        qf = 1.big();
        dist_bound = 600;
        max_tries = 1000;

        let a = -"1879950577972003736005350319459510223264087780255173448508833053707834624991079524305041640".big();
        let b = "545356496340105514642782361550937077798117534298654683033730684848663459606736968920419716".big();
        let cmp = QuatAlgEl::new(a.clone(), b.clone(), a, b, 2.big(), qa);
        let el = enumerate_cvp_filter(
            basis.clone(),
            target.clone(),
            qf,
            params,
            Box::new(test_cvp_condition),
            dist_bound,
            max_tries,
        );
        assert!(el.unwrap() == cmp);
    }
}
