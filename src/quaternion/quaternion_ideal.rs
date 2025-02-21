use rug::Integer;

use crate::{
    quaternion::{quaternion_algebra::QuatAlgEl, quaternion_order::QuaternionOrder},
    util::Big,
};

use super::{lattice::Lattice, quaternion_algebra::QuatAlg};

#[derive(Clone, Debug, PartialEq)]
pub struct QuaternionIdeal {
    pub lattice: Lattice,
    pub norm: Integer,
    pub parent: Option<QuaternionOrder>,
}

impl QuaternionIdeal {
    pub fn new(
        lattice: Lattice,
        norm: Integer,
        parent: Option<QuaternionOrder>,
    ) -> QuaternionIdeal {
        Self {
            lattice,
            norm,
            parent,
        }
    }

    /// Create a (principal) left ideal in `order` generated by `x`.
    pub fn new_principal_left_ideal(x: QuatAlgEl, order: QuaternionOrder) -> QuaternionIdeal {
        let mat = x.matrix("right".to_owned());
        let basis = mat.basis * order.lattice.basis.clone();

        let mut lat = Lattice::new(basis, x.denom.clone() * order.lattice.denom.clone());
        lat = lat.reduce_denom();
        lat = lat.hnf();

        let norm = x.reduced_norm();
        assert!(norm.denom() == &1.big());
        let ideal = QuaternionIdeal::new(lat, norm.numer().clone(), Some(order));

        ideal
    }

    /// Create a left ideal in order generated by the primitive element `x` and the integer `n`.
    pub fn new_left_ideal_from_primitive(
        x: QuatAlgEl,
        n: Integer,
        order: QuaternionOrder,
    ) -> QuaternionIdeal {
        // Compute ideal generated by x
        let ideal = QuaternionIdeal::new_principal_left_ideal(x, order.clone());

        // Compute norm
        let norm = ideal.norm.gcd(&n);

        // Compute ideal generated by n (without reducing denominator)
        let basis = order.lattice.basis.clone() * n;
        let ideal_n_lattice = Lattice::new(basis, order.lattice.denom.clone());

        // Add lattices (reduces denominators)
        let lat = ideal.lattice.add(ideal_n_lattice);

        QuaternionIdeal::new(lat, norm, Some(order))
    }

    /// Make `x` primitive, remove content from `N`, then create ideal.
    ///
    /// Given `x` = n·y ∈ order with y primitive, given an integer `N`, create the ideal generated by y and N / gcd(n, N).
    /// Arguments:
    /// `x`: generating element, a primitive element (of order) obtained from it will be used for ideal generation.
    /// `N`: generating integer.
    /// `order`: maximal order whose left ideal is searched.
    /// `qa`: quaternion algebra.
    pub fn new_left_ideal(
        x: QuatAlgEl,
        n: Integer,
        order: QuaternionOrder,
        qa: QuatAlg,
    ) -> QuaternionIdeal {
        let (coord, imprim) = x.factor_in_order(order.lattice.clone());

        let coord_new = order.lattice.basis.clone() * coord;
        let prim = QuatAlgEl::new_from_matrix(coord_new, order.lattice.denom.clone(), qa);

        // imprim = gcd(imprimitive part of x, N)
        // n1 = N / imprim
        let gcd = imprim.gcd(&n);
        let n1 = n / gcd;

        QuaternionIdeal::new_left_ideal_from_primitive(prim, n1, order)
    }
}

#[cfg(test)]
mod tests {
    use crate::{linalg::matrix::Matrix, quaternion::quaternion_algebra::QuatAlg, util::Big};

    use super::*;

    #[test]
    fn new_principal_left_ideal_test() {
        let p = 367.big();
        let qa = QuatAlg::new(-p.clone());

        let basis = Matrix::zeros(4, 4);
        let mut lat = Lattice::new(basis.clone(), 6.big());
        lat.basis[(0, 0)] = 6.big();
        lat.basis[(1, 1)] = 6.big();
        lat.basis[(1, 2)] = 3.big();
        lat.basis[(2, 2)] = 3.big();
        lat.basis[(3, 3)] = 3.big();
        lat.basis[(0, 3)] = 3.big();
        let order = QuaternionOrder::new(lat);
        let gamma = QuatAlgEl::new(
            438.big(),
            400.big(),
            156.big(),
            -2.big(),
            2.big(),
            qa.clone(),
        );

        let ideal = QuaternionIdeal::new_principal_left_ideal(gamma, order);
        assert!(ideal.norm == 2321156.big());

        assert!(ideal.lattice.basis[(0, 0)] == 1160578.big());
        assert!(ideal.lattice.basis[(1, 0)] == 0.big());
        assert!(ideal.lattice.basis[(2, 0)] == 0.big());
        assert!(ideal.lattice.basis[(3, 0)] == 0.big());

        assert!(ideal.lattice.basis[(0, 1)] == 0.big());
        assert!(ideal.lattice.basis[(1, 1)] == 1160578.big());
        assert!(ideal.lattice.basis[(2, 1)] == 0.big());
        assert!(ideal.lattice.basis[(3, 1)] == 0.big());

        assert!(ideal.lattice.basis[(0, 2)] == 310126.big());
        assert!(ideal.lattice.basis[(1, 2)] == 182529.big());
        assert!(ideal.lattice.basis[(2, 2)] == 1.big());
        assert!(ideal.lattice.basis[(3, 2)] == 0.big());

        assert!(ideal.lattice.basis[(0, 3)] == 978049.big());
        assert!(ideal.lattice.basis[(1, 3)] == 310126.big());
        assert!(ideal.lattice.basis[(2, 3)] == 0.big());
        assert!(ideal.lattice.basis[(3, 3)] == 1.big());
    }

    #[test]
    fn new_left_ideal_from_primitive_test() {
        let p = 367.big();
        let qa = QuatAlg::new(-p.clone());

        let basis = Matrix::zeros(4, 4);
        let mut lat = Lattice::new(basis.clone(), 2.big());
        lat.basis[(0, 0)] = 2.big();
        lat.basis[(1, 1)] = 2.big();
        lat.basis[(1, 2)] = 1.big();
        lat.basis[(2, 2)] = 1.big();
        lat.basis[(3, 3)] = 1.big();
        lat.basis[(0, 3)] = 1.big();
        let order = QuaternionOrder::new(lat);

        let gamma = QuatAlgEl::new(
            219.big(),
            200.big(),
            78.big(),
            -1.big(),
            1.big(),
            qa.clone(),
        );
        let n = 31.big();

        let ideal = QuaternionIdeal::new_left_ideal_from_primitive(gamma, n.clone(), order);
        assert!(ideal.norm == n);
        assert!(ideal.lattice.denom == 2.big());

        assert!(ideal.lattice.basis[(0, 0)] == 62.big());
        assert!(ideal.lattice.basis[(1, 0)] == 0.big());
        assert!(ideal.lattice.basis[(2, 0)] == 0.big());
        assert!(ideal.lattice.basis[(3, 0)] == 0.big());

        assert!(ideal.lattice.basis[(0, 1)] == 0.big());
        assert!(ideal.lattice.basis[(1, 1)] == 62.big());
        assert!(ideal.lattice.basis[(2, 1)] == 0.big());
        assert!(ideal.lattice.basis[(3, 1)] == 0.big());

        assert!(ideal.lattice.basis[(0, 2)] == 2.big());
        assert!(ideal.lattice.basis[(1, 2)] == 1.big());
        assert!(ideal.lattice.basis[(2, 2)] == 1.big());
        assert!(ideal.lattice.basis[(3, 2)] == 0.big());

        assert!(ideal.lattice.basis[(0, 3)] == 61.big());
        assert!(ideal.lattice.basis[(1, 3)] == 2.big());
        assert!(ideal.lattice.basis[(2, 3)] == 0.big());
        assert!(ideal.lattice.basis[(3, 3)] == 1.big());
    }

    #[test]
    fn new_left_ideal_test() {
        let p = 103.big();
        let qa = QuatAlg::new(-p.clone());

        let basis = Matrix::zeros(4, 4);
        let mut lat = Lattice::new(basis.clone(), 2.big());

        lat.basis[(0, 0)] = 2.big();
        lat.basis[(1, 1)] = 2.big();
        lat.basis[(1, 2)] = 1.big();
        lat.basis[(2, 2)] = 1.big();
        lat.basis[(3, 3)] = 1.big();
        lat.basis[(0, 3)] = 1.big();

        let order = QuaternionOrder::new(lat.clone());

        let x = QuatAlgEl::new(6.big(), 10.big(), 14.big(), 22.big(), 1.big(), qa.clone());
        let n = 17.big();

        let ideal =
            QuaternionIdeal::new_left_ideal(x.clone(), n.clone(), order.clone(), qa.clone());
        let gen = x.get_primitive_in_order(order.lattice.clone(), qa);
        let ideal_cmp = QuaternionIdeal::new_left_ideal_from_primitive(gen, n, order);

        assert!(ideal == ideal_cmp);
    }
}
