use crate::{linalg::matrix::Matrix, util::Big};
use num::rational::Ratio;
use rug::Integer;
use std::ops::{Add, Div, Index, Mul};

use super::lattice::Lattice;

/// Quaternion algebra with i^2=-1 and j^2 = p.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct QuatAlg {
    pub p: Integer,
}

impl QuatAlg {
    pub fn new(p: Integer) -> Self {
        QuatAlg { p }
    }
}

/// The elements are represented in basis (1, i, j, ij) of the quaternion algebra,
/// with i^2=-1 and j^2 = p.
#[derive(Clone, Debug, Default)]
pub struct QuatAlgEl {
    pub x: Integer,
    pub y: Integer,
    pub z: Integer,
    pub t: Integer,
    pub denom: Integer,
    pub algebra: QuatAlg,
}

impl QuatAlgEl {
    pub fn new(
        x: Integer,
        y: Integer,
        z: Integer,
        t: Integer,
        denom: Integer,
        algebra: QuatAlg,
    ) -> Self {
        QuatAlgEl {
            x: x,
            y: y,
            z: z,
            t: t,
            denom: denom,
            algebra: algebra,
        }
    }

    pub fn new_from_matrix(m: Matrix<Integer>, denom: Integer, algebra: QuatAlg) -> Self {
        QuatAlgEl {
            x: m[0].clone(),
            y: m[1].clone(),
            z: m[2].clone(),
            t: m[3].clone(),
            denom,
            algebra,
        }
    }

    pub fn zero(qa: QuatAlg) -> Self {
        QuatAlgEl::new(0.big(), 0.big(), 0.big(), 0.big(), 1.big(), qa.clone())
    }

    pub fn one(qa: QuatAlg) -> Self {
        QuatAlgEl::new(1.big(), 0.big(), 0.big(), 0.big(), 1.big(), qa.clone())
    }

    pub fn new_from_ratio_matrix(m: Matrix<Ratio<Integer>>, algebra: QuatAlg) -> Self {
        let lcm = Matrix::lcm(&m);
        let to_common_denom =
            |a: Ratio<Integer>| -> Integer { a.numer() * lcm.clone() / a.denom() };

        QuatAlgEl {
            x: to_common_denom(m[0].clone()),
            y: to_common_denom(m[1].clone()),
            z: to_common_denom(m[2].clone()),
            t: to_common_denom(m[3].clone()),
            denom: lcm,
            algebra,
        }
    }

    pub fn coeffs(&self) -> Matrix<Integer> {
        let mut coeffs = Matrix::zeros(1, 4);
        coeffs[0] = self.x.clone();
        coeffs[1] = self.y.clone();
        coeffs[2] = self.z.clone();
        coeffs[3] = self.t.clone();

        coeffs
    }

    pub fn coeffs_from_ratio(&self) -> Matrix<Ratio<Integer>> {
        let mut coeffs = Matrix::zeros(1, 4);
        coeffs[0] = Ratio::new(self.x.clone(), self.denom.clone());
        coeffs[1] = Ratio::new(self.y.clone(), self.denom.clone());
        coeffs[2] = Ratio::new(self.z.clone(), self.denom.clone());
        coeffs[3] = Ratio::new(self.t.clone(), self.denom.clone());

        coeffs
    }

    pub fn reduced_norm(&self) -> Ratio<Integer> {
        let n = self.x.clone() * self.x.clone() + self.y.clone() * self.y.clone()
            - self.algebra.p.clone() * self.z.clone() * self.z.clone()
            - self.algebra.p.clone() * self.t.clone() * self.t.clone();

        Ratio::new(n, self.denom.clone() * self.denom.clone())
    }

    pub fn conjugate(&self) -> QuatAlgEl {
        QuatAlgEl::new(
            self.x.clone(),
            -self.y.clone(),
            -self.z.clone(),
            -self.t.clone(),
            self.denom.clone(),
            self.algebra.clone(),
        )
    }

    pub fn inverse(self) -> QuatAlgEl {
        let n = self.reduced_norm();
        let conj = self.conjugate();
        let mut i = conj * n.denom().clone();
        i.denom = i.denom * n.numer();

        i.normalize()
    }

    /// Compute the coefficients of the element in a given basis.
    pub fn express_in_basis(&self, basis: Lattice) -> Matrix<Integer> {
        let mut elem_in_basis = Matrix::zeros(4, 1);
        let contains = basis.contains(&mut elem_in_basis, self.clone());
        assert!(contains);

        elem_in_basis
    }

    /// Factor the quaternion algebra element into its primitive and imprimitive parts
    /// according to the order basis.
    /// First write the element in the order basis: self = `order_basis` * elem_in_order_basis.
    /// Then compute `primitive` such that elem_in_order_basis = `content` * primitive.
    /// Output (`primitive`, `content`).
    pub fn factor_in_order(&self, order_basis: Lattice) -> (Matrix<Integer>, Integer) {
        let elem_in_order_basis = self.express_in_basis(order_basis);
        let content = elem_in_order_basis.gcd();

        (elem_in_order_basis.div(content.clone()), content)
    }

    /// Factor the quaternion algebra element into its primitive and imprimitive parts
    /// according to the order basis. Return the primitive part expressed in the standard basis.
    pub fn get_primitive_in_order(&self, order_basis: Lattice, qa: QuatAlg) -> QuatAlgEl {
        let (coord, _) = self.factor_in_order(order_basis.clone());
        let coord_new = order_basis.basis.clone() * coord;

        QuatAlgEl::new_from_matrix(coord_new, order_basis.denom.clone(), qa)
    }

    pub fn normalize(&self) -> QuatAlgEl {
        let mut gcd = self.x.clone().gcd(&self.y);
        gcd = gcd.gcd(&self.z);
        gcd = gcd.gcd(&self.t);
        gcd = gcd.gcd(&self.denom);
        if gcd > 1.big() {
            QuatAlgEl::new(
                self.x.clone() / gcd.clone(),
                self.y.clone() / gcd.clone(),
                self.z.clone() / gcd.clone(),
                self.t.clone() / gcd.clone(),
                self.denom.clone() / gcd,
                self.algebra.clone(),
            )
        } else {
            self.clone()
        }
    }

    /// Take the two given elements and return two equivalent ones with the same denominator.
    fn put_to_same_denom(a: QuatAlgEl, b: QuatAlgEl) -> (QuatAlgEl, QuatAlgEl) {
        let gcd = a.denom.clone().gcd(&b.denom);
        let mut new_a = QuatAlgEl::default();
        let mut new_b = QuatAlgEl::default();

        new_a.denom = a.clone().denom / gcd.clone();
        new_b.denom = b.clone().denom / gcd;

        new_a = a.clone() * new_b.clone().denom;
        new_b = b.clone() * new_a.clone().denom;

        new_a.denom = b.denom * new_a.denom;
        new_b.denom = a.denom * new_b.denom;

        (new_a, new_b)
    }

    /// Return a matrix where columns present 1 * self, i * self, j * self, k * self for
    /// action right and where columns present self * 1, self * i, self * j, self * k for
    /// action left.
    ///
    /// For example, when an ideal by an element x is to be created,
    /// M = x.matrix("right") returns the matrix with x, i * x, j * x, k * x in its columns
    /// (the columns present the basis of the ideal (x)).
    pub fn matrix(&self, action: String) -> Lattice {
        let mut m = Matrix::zeros(4, 4);
        m[(0, 0)] = self.x.clone();
        m[(1, 0)] = self.y.clone();
        m[(2, 0)] = self.z.clone();
        m[(3, 0)] = self.t.clone();

        if action == "right" {
            m[(0, 1)] = -self.y.clone();
            m[(1, 1)] = self.x.clone();
            m[(2, 1)] = -self.t.clone();
            m[(3, 1)] = self.z.clone();

            m[(0, 2)] = self.algebra.p.clone() * self.z.clone();
            m[(1, 2)] = -self.algebra.p.clone() * self.t.clone();
            m[(2, 2)] = self.x.clone();
            m[(3, 2)] = -self.y.clone();

            m[(0, 3)] = self.algebra.p.clone() * self.t.clone();
            m[(1, 3)] = self.algebra.p.clone() * self.z.clone();
            m[(2, 3)] = self.y.clone();
            m[(3, 3)] = self.x.clone();
        } else if action == "left" {
            m[(0, 1)] = -self.y.clone();
            m[(1, 1)] = self.x.clone();
            m[(2, 1)] = self.t.clone();
            m[(3, 1)] = -self.z.clone();

            m[(0, 2)] = self.algebra.p.clone() * self.z.clone();
            m[(1, 2)] = self.algebra.p.clone() * self.t.clone();
            m[(2, 2)] = self.x.clone();
            m[(3, 2)] = self.y.clone();

            m[(0, 3)] = self.algebra.p.clone() * self.t.clone();
            m[(1, 3)] = -self.algebra.p.clone() * self.z.clone();
            m[(2, 3)] = -self.y.clone();
            m[(3, 3)] = self.x.clone();
        }

        Lattice::new(m, self.denom.clone())
    }
}

impl PartialEq for QuatAlgEl {
    fn eq(&self, other: &Self) -> bool {
        self.coeffs() * other.denom.clone() == other.coeffs() * self.denom.clone()
    }
}

impl Add for QuatAlgEl {
    type Output = Self;

    fn add(self, other: QuatAlgEl) -> QuatAlgEl {
        let (new_a, new_b) = QuatAlgEl::put_to_same_denom(self, other);

        QuatAlgEl {
            x: new_a.x + new_b.x,
            y: new_a.y + new_b.y,
            z: new_a.z + new_b.z,
            t: new_a.t + new_b.t,
            denom: new_a.denom,
            algebra: new_a.algebra,
        }
    }
}

impl Index<usize> for QuatAlgEl {
    type Output = Integer;

    fn index(&self, ind: usize) -> &Integer {
        match ind {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            3 => &self.t,
            _ => panic!("cannot access by this index"),
        }
    }
}

impl<T> Mul<T> for QuatAlgEl
where
    Integer: From<T>,
    T: Clone,
{
    type Output = Self;

    fn mul(self, other: T) -> QuatAlgEl {
        QuatAlgEl {
            x: self.x * Integer::from(other.clone()),
            y: self.y * Integer::from(other.clone()),
            z: self.z * Integer::from(other.clone()),
            t: self.t * Integer::from(other),
            denom: self.denom,
            algebra: self.algebra,
        }
    }
}

impl<T> Div<T> for QuatAlgEl
where
    Integer: From<T>,
    T: Clone,
{
    type Output = Self;

    fn div(self, other: T) -> QuatAlgEl {
        let other_big = Integer::from(other.clone());
        assert!(other_big != 0.big());

        let coeffs = vec![
            self.x.clone(),
            self.y.clone(),
            self.z.clone(),
            self.t.clone(),
        ];
        let non_zero = coeffs.iter().filter(|x| **x != 0.big());
        let gcd = non_zero
            .into_iter()
            .fold(other_big.clone(), |acc, x| acc.gcd(x));

        if gcd == other_big {
            QuatAlgEl {
                x: self.x / Integer::from(other.clone()),
                y: self.y / Integer::from(other.clone()),
                z: self.z / Integer::from(other.clone()),
                t: self.t / Integer::from(other),
                denom: self.denom,
                algebra: self.algebra,
            }
        } else {
            QuatAlgEl {
                x: self.x,
                y: self.y,
                z: self.z,
                t: self.t,
                denom: self.denom * Integer::from(other.clone()),
                algebra: self.algebra,
            }
        }
    }
}

impl Mul<QuatAlgEl> for QuatAlgEl {
    type Output = Self;

    fn mul(self, other: QuatAlgEl) -> QuatAlgEl {
        let x = self.algebra.p.clone()
            * (self.z.clone() * other.z.clone() + self.t.clone() * other.t.clone())
            + self.x.clone() * other.x.clone()
            - self.y.clone() * other.y.clone();
        let y = -self.algebra.p.clone()
            * (self.z.clone() * other.t.clone() - self.t.clone() * other.z.clone())
            + self.x.clone() * other.y.clone()
            + self.y.clone() * other.x.clone();
        let z = self.x.clone() * other.z.clone() + self.z.clone() * other.x.clone()
            - self.y.clone() * other.t.clone()
            + self.t.clone() * other.y.clone();
        let t = self.x.clone() * other.t.clone() + self.t.clone() * other.x.clone()
            - self.z.clone() * other.y.clone()
            + self.y.clone() * other.z.clone();

        QuatAlgEl {
            x: x,
            y: y,
            z: z,
            t: t,
            denom: self.denom.clone() * other.denom,
            algebra: self.algebra,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::util::Big;

    use super::*;

    #[test]
    fn ops() {
        let zero: Integer = Integer::from(0u32);
        let one: Integer = Integer::from(1u32);
        let two = Integer::from(2u32);
        let three = Integer::from(3u32);
        let four = Integer::from(4u32);
        let five = Integer::from(5u32);
        let six = Integer::from(6u32);
        let seven = Integer::from(7u32);
        let nine = Integer::from(9u32);
        let ten = Integer::from(10u32);

        let p =
            "73743043621499797449074820543863456997944695372324032511999999999999999999999".big();
        let qa = QuatAlg { p: -p };

        let a = QuatAlgEl::new(
            one.clone(),
            one.clone(),
            one.clone(),
            two.clone(),
            two.clone(),
            qa.clone(),
        );
        let b = QuatAlgEl::new(
            two.clone(),
            two.clone(),
            one.clone(),
            two.clone(),
            three.clone(),
            qa.clone(),
        );
        let new_a = QuatAlgEl::new(
            three.clone(),
            three.clone(),
            three.clone(),
            six.clone(),
            six.clone(),
            qa.clone(),
        );
        let new_b = QuatAlgEl::new(
            four.clone(),
            four.clone(),
            two.clone(),
            four.clone(),
            six.clone(),
            qa.clone(),
        );
        let a_plus_b = QuatAlgEl::new(
            seven.clone(),
            seven.clone(),
            five.clone(),
            ten.clone(),
            six.clone(),
            qa.clone(),
        );

        let (computed_new_a, computed_new_b) = QuatAlgEl::put_to_same_denom(a.clone(), b.clone());

        assert_eq!(computed_new_a, new_a);
        assert_eq!(computed_new_b, new_b);

        let computed_a_plus_b = a.clone() + b;

        assert_eq!(computed_a_plus_b, a_plus_b);

        let c = QuatAlgEl::new(
            nine.clone(),                 // 9
            nine.clone() + three.clone(), // 12
            nine.clone() + nine.clone(),  // 18
            nine.clone() + three.clone(), // 12
            six.clone(),                  // 6
            qa.clone(),
        );
        let norm_c = QuatAlgEl::new(
            three.clone(),
            four.clone(),
            six.clone(),
            four.clone(),
            two.clone(),
            qa.clone(),
        );
        let computed_norm_c = c.normalize();
        assert_eq!(computed_norm_c, norm_c);

        let d = QuatAlgEl::new(
            nine.clone(),                 // 9
            ten.clone(),                  // 10
            nine.clone() + nine.clone(),  // 18
            nine.clone() + three.clone(), // 12
            six.clone(),                  // 6
            qa.clone(),
        );
        let norm_d = QuatAlgEl::new(
            nine.clone(),                 // 9
            ten.clone(),                  // 10
            nine.clone() + nine.clone(),  // 18
            nine.clone() + three.clone(), // 12
            six.clone(),                  // 6
            qa.clone(),
        );

        let computed_norm_d = d.normalize();
        assert_eq!(computed_norm_d, norm_d);

        let e = QuatAlgEl::new(
            four.clone(),
            six.clone(),
            four.clone(),
            six.clone(),
            three.clone(),
            qa.clone(),
        );
        let e_div_by_2 = QuatAlgEl::new(
            two.clone(),
            three.clone(),
            two.clone(),
            three.clone(),
            three.clone(),
            qa.clone(),
        );
        let computed_e_div_by_2 = e / two.clone();

        assert_eq!(computed_e_div_by_2, e_div_by_2);

        let algebra_one = QuatAlgEl::new(
            one.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            one.clone(),
            qa.clone(),
        );

        let a_inv = a.clone().inverse();
        assert_eq!(a.clone() * a_inv, algebra_one);

        let m = a.matrix("left".to_owned());
        println!("{:?}", m);
    }
}
