use crate::{linalg::matrix::Matrix, util::Big};
use rug::Integer;

use super::lattice::Lattice;

#[derive(Clone, Debug, PartialEq)]
pub struct QuaternionOrder {
    pub lattice: Lattice,
}

impl QuaternionOrder {
    pub fn new(lattice: Lattice) -> Self {
        QuaternionOrder { lattice }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ExtremalMaximalOrder {
    pub order: QuaternionOrder,
    pub q: Integer, // absolute value of sqrt of i
}

impl ExtremalMaximalOrder {
    pub fn new(order: QuaternionOrder, q: Integer) -> Self {
        ExtremalMaximalOrder { order, q }
    }
}

/*
The standard maximal extremal order is a maximal order that is isomorphic
to the endomorphism ring of the curve of j-invariant 1728:
E: y^2 = x^3 + x, defined over F_{p^2}.
End(E) = <1, i, (i+j)/2, (1+k)/2> with i^2 = -1, j^2 = -p, k = ij.
Moreover, i = iota and j = pi, where
iota(x, y) = (-x, \sqrt{-1} * y),
pi(x, y) = (x^p, y^p)
*/
pub fn standard_maximal_extremal_order() -> ExtremalMaximalOrder {
    let basis = Matrix::zeros(4, 4);
    let mut lat = Lattice::new(basis.clone(), 2.big());
    lat.basis[(0, 0)] = 2.big();
    lat.basis[(0, 3)] = 1.big();
    lat.basis[(1, 1)] = 2.big();
    lat.basis[(1, 2)] = 1.big();
    lat.basis[(2, 2)] = 1.big();
    lat.basis[(3, 3)] = 1.big();
    let order = QuaternionOrder::new(lat);

    ExtremalMaximalOrder::new(order, 1.big())
}

pub fn standard_extremal_from(order: QuaternionOrder) -> ExtremalMaximalOrder {
    ExtremalMaximalOrder::new(order, 1.big())
}
