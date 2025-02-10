use std::ops::Div;

// use crate::quaternion::quaternion_algebra::{QuatAlg, QuatAlgEl};
use crate::util::Big;
use crate::{linalg::matrix::Matrix, util::linear_combination};
use num_traits::{One, Zero};
use rug::Integer;

use super::quaternion_algebra::{QuatAlg, QuatAlgEl};

#[derive(Clone, Debug)]
pub struct Lattice {
    pub basis: Matrix<Integer>,
    pub denom: Integer,
}

impl Lattice {
    pub fn new(basis: Matrix<Integer>, denom: Integer) -> Lattice {
        Self {
            basis: basis,
            denom: denom,
        }
    }

    pub fn reduce_denom(&self) -> Lattice {
        let mut gcd = self.basis.gcd();
        gcd = gcd.gcd(&self.denom);
        let basis = self.basis.div(gcd.clone());
        let denom = self.denom.clone().div(gcd);

        Lattice::new(basis, denom)
    }

    pub fn hnf(&self) -> Lattice {
        let generators =
            Matrix::<Integer>::concatenate(Matrix::<Integer>::zeros(4, 4), self.basis.clone());
        let mut a = vec![];
        for l in 0..8 {
            a.push(generators.get_col(l));
        }

        let basis = Lattice::hnf_core(a);
        let lat = Lattice::new(basis, self.denom.clone());
        lat.reduce_denom()
    }

    pub fn is_hnf(&self) -> bool {
        let mut ok = true;
        let mut ind = 0;
        for i in 0..4 {
            for j in 0..i {
                ok &= self.basis[(i, j)].is_zero();
            }
            // find first non 0 element of line
            let mut found = false;
            for j in i..4 {
                if found {
                    // all values are positive, and first non-0 is the largest of that line
                    ok &= self.basis[(i, j)].is_positive() || self.basis[(i, j)].is_zero();
                    ok &= self.basis[(i, ind)] > self.basis[(i, j)];
                } else {
                    if !self.basis[(i, j)].is_zero() {
                        found = true;
                        ind = j;
                        // must be non-negative
                        ok &= self.basis[(i, j)].is_positive() || self.basis[(i, j)].is_zero();
                    }
                }
            }
        }
        // check that first non-zero elements index per column is strictly increasing
        let linestart = 0;
        let mut i = 0;
        for j in 0..4 {
            while i < 4 && self.basis[(i, j)].is_zero() {
                i = i + 1;
            }
            if i != 4 {
                ok &= linestart < i + 1;
            }
            i = 0;
        }

        ok
    }

    /// Algorithm 2.4.5 in Henri Cohen's "A Course in Computational Algebraic Number Theory"
    /// (Springer Verlag, in series "Graduate texts in Mathematics") from 1993.
    pub fn hnf_core(mat: Vec<Matrix<Integer>>) -> Matrix<Integer> {
        let mut a = mat.clone();

        let mut i = 4;
        let mut j = 7;
        let mut k = 7;

        let mut d;
        let mut u;
        let mut v;
        let mut r;
        let mut coeff_1;
        let mut coeff_2;

        while i > 0 {
            while j != 0 {
                j -= 1;
                if a[j][i - 1] != Integer::zero() {
                    // assumption that ibz_xgcd outputs u,v which are small in absolute value is needed here
                    (d, u, v) = a[k][i - 1]
                        .clone()
                        .extended_gcd(a[j][i - 1].clone(), Integer::new());

                    if u.is_zero() {
                        (v, _) = a[k][i - 1].clone().div_rem(a[j][i - 1].clone());
                        u = Integer::one();
                        v = u.clone() - v;
                    }

                    let c = linear_combination(u, &a[k], v, &a[j]);
                    (coeff_1, _) = a[k][i - 1].clone().div_rem(d.clone());
                    (coeff_2, _) = a[j][i - 1].clone().div_rem(d.clone());
                    coeff_2 = -coeff_2;
                    a[j] = linear_combination(coeff_1, &a[j], coeff_2, &a[k]);
                    a[k] = c;
                }
            }

            let mut b = a[k][i - 1].clone();
            let min_one: Integer = -Integer::one();
            if Integer::is_negative(&b) {
                a[k] = a[k].clone() * min_one.clone();
                b = b * min_one.clone();
            }
            if b.is_zero() {
                k += 1;
            } else {
                for j in k + 1..8 {
                    (d, r) = a[j][i - 1].clone().div_rem(b.clone());
                    if Integer::is_negative(&r) {
                        r = Integer::one();
                        d = d - r;
                    }
                    r = Integer::one();
                    d = d * min_one.clone();
                    a[j] = linear_combination(r, &a[j], d, &a[k]);
                }
            }

            if i > 1 {
                k -= 1;
                j = k;
            }
            i -= 1;
        }

        let mut basis = Matrix::eye(4);
        for j in 4..8 {
            for i in 0..4 {
                basis[(i, j - 4)] = a[j][i].clone();
            }
        }

        basis
    }

    pub fn dual_without_hnf(&self) -> Lattice {
        let t = self.basis.transpose();
        let (inv, det) = Matrix::inverse_with_det_as_denom(&t);
        // dual_denom = det/lat_denom
        let inv_basis = inv * self.denom.clone();

        Lattice::new(inv_basis, det)
    }

    /// Compute the sum of two lattices. Since all lattices contain zero, the addition is equal to
    /// the lattice generated by the union of the two lattices.
    pub fn add(&self, other: Lattice) -> Lattice {
        let m1 = other.basis * self.denom.clone();
        let m2 = self.basis.clone() * other.denom.clone();

        let mut a = vec![];
        for i in 0..4 {
            a.push(m1.get_col(i));
        }
        for i in 0..4 {
            a.push(m2.get_col(i));
        }

        let basis = Lattice::hnf_core(a);
        let lat = Lattice::new(basis, self.denom.clone() * other.denom);
        lat.reduce_denom()
    }

    /// Compute the intersection of two lattices. The intersection is computed by taking
    /// the duals of the two lattices, by computing their intersection,
    /// and taking the dual of the intersection.
    pub fn intersect(&self, other: Lattice) -> Lattice {
        let dual1 = self.dual_without_hnf();
        let dual2 = other.dual_without_hnf();
        let sum = dual1.add(dual2);
        let res = sum.dual_without_hnf();
        res.hnf()
    }

    pub fn mul(&self, other: Lattice, qa: QuatAlg) -> Lattice {
        let mut a1 = vec![];
        let mut a2 = vec![];
        let denom = self.denom.clone() * other.denom.clone();
        for i in 0..4 {
            let lat1_vec =
                QuatAlgEl::new_from_matrix(self.basis.get_col(i), self.denom.clone(), qa.clone());
            for j in 0..4 {
                let lat2_vec = QuatAlgEl::new_from_matrix(
                    other.basis.get_col(j),
                    other.denom.clone(),
                    qa.clone(),
                );
                let mut m = lat1_vec.clone() * lat2_vec.clone();
                if m.denom != denom.clone() {
                    let d = denom.clone() / m.clone().denom;
                    m = m * d;
                }
                if i < 2 {
                    a1.push(m.coeffs());
                } else {
                    a2.push(m.coeffs());
                }
            }
        }

        let mut basis1 = Lattice::hnf_core(a1);
        let basis2 = Lattice::hnf_core(a2);

        a1 = vec![];
        for l in 0..4 {
            a1.push(basis1.get_col(l));
        }
        for l in 0..4 {
            a1.push(basis2.get_col(l));
        }
        basis1 = Lattice::hnf_core(a1);

        Lattice::new(basis1, denom)
    }

    /// Solve self * coord = x. Self needs to be in a HNF format and of a full rank.
    /// Return true if x is in the lattice, return false otherwise.
    pub fn contains(&self, coord: &mut Matrix<Integer>, x: QuatAlgEl) -> bool {
        // test if rank 4 lattice under HNF
        assert!(self.is_hnf());
        for i in 0..4 {
            assert!(!self.basis[(i, i)].is_zero());
        }
        let mut work_x = Matrix::zeros(4, 1);
        let mut column = Matrix::zeros(4, 1);
        for i in 0..4 {
            // put on same denominator, 1st part
            work_x[i] = x[i].clone() * self.denom.clone();
        }
        let mut res = true;
        for i in 0..4 {
            if res {
                // put on same denominator, 2nd part
                let prod = x.denom.clone() * self.basis[(3 - i, 3 - i)].clone();
                let mut r: Integer;
                (coord[3 - i], r) = work_x[3 - i].clone().div_rem(prod);
                if r.is_zero() {
                    for j in 0..4 {
                        // put on same denominator here also
                        column[j] = self.basis[(j, 3 - i)].clone() * x.denom.clone();
                    }
                    // negate quotient
                    r = -coord[3 - i].clone();
                    work_x = linear_combination(1.big(), &work_x, r, &column);
                } else {
                    res = false;
                }
            }
        }

        for i in 0..4 {
            res &= work_x[i].is_zero();
        }

        res
    }

    pub fn index(&self, overlat: Lattice) -> Integer {
        // index = overlat.denom**4
        let mut index = overlat.denom.clone() * overlat.denom;
        index = index.clone() * index;
        // index = overlat.denom**4 * det(self.basis)
        for i in 0..4 {
            index = index.clone() * self.basis[(i, i)].clone();
        }
        // tmp = self.denom**4
        let mut tmp = self.denom.clone() * self.denom.clone();
        tmp = tmp.clone() * tmp;
        // tmp = self.denom**4 * det(overlat.basis)
        for i in 0..4 {
            tmp = tmp.clone() * overlat.basis[(i, i)].clone();
        }
        index = index / tmp.clone();
        assert!(tmp != Integer::zero());

        index.abs()
    }
}

impl PartialEq for Lattice {
    fn eq(&self, other: &Self) -> bool {
        self.basis.clone() * other.denom.clone().abs()
            == other.basis.clone() * self.denom.clone().abs()
    }
}

#[cfg(test)]
mod tests {
    use rug::Integer;

    use crate::{quaternion::{quaternion_algebra::QuatAlg, quaternion_order::standard_maximal_extremal_order}, util::Big};

    use super::*;

    #[test]
    fn hnf_core_1() {
        let generators = Matrix::zeros(4, 8);
        let mut a = vec![];
        for l in 0..8 {
            a.push(generators.get_col(l));
        }

        a[2][0] = 2.big();
        a[4][0] = 4.big();
        a[2][3] = 5.big();
        a[7][3] = 6.big();
        a[7][1] = 7.big();
        a[3][1] = 8.big();
        a[1][1] = 9.big();
        a[6][0] = 10.big();
        a[5][0] = 11.big();
        a[0][0] = 12.big();

        /*
        matrix:
        12 0 2 0 4 11 10 0
        0 9 0 8 0 0 0 7
        0 0 0 0 0 0 0 0
        0 0 5 0 0 0 0 6
        */
        let mut cmp: Matrix<Integer> = Matrix::zeros(4, 4);
        cmp[(0, 1)] = 1.big();
        cmp[(1, 2)] = 1.big();
        cmp[(3, 3)] = 1.big();

        let basis = Lattice::hnf_core(a);
        assert_eq!(basis, cmp);
    }

    #[test]
    fn hnf_core_2() {
        let generators = Matrix::zeros(4, 8);
        let mut a = vec![];
        for l in 0..8 {
            a.push(generators.get_col(l));
        }

        a[2][0] = 2.big();
        a[4][0] = 4.big();
        a[2][3] = 5.big();
        a[7][3] = 6.big();
        a[7][1] = 7.big();
        a[3][1] = 8.big();
        a[1][1] = 9.big();
        a[6][0] = 10.big();
        a[5][0] = 11.big();
        a[0][0] = 12.big();
        a[0][2] = 2.big();
        a[5][2] = 1.big();

        /*
        matrix:
        12 0 2 0 4 11 10 0
        0 9 0 8 0 0 0 7
        2 0 0 0 0 1 0 0
        0 0 5 0 0 0 0 6
        */
        let mut cmp: Matrix<Integer> = Matrix::zeros(4, 4);
        cmp[(0, 0)] = 2.big();
        cmp[(1, 1)] = 1.big();
        cmp[(0, 2)] = 1.big();
        cmp[(2, 2)] = 1.big();
        cmp[(3, 3)] = 1.big();

        let basis = Lattice::hnf_core(a);
        assert_eq!(basis, cmp);
    }

    #[test]
    fn hnf_core_3() {
        let generators = Matrix::zeros(4, 8);
        let mut a = vec![];
        for l in 0..8 {
            a.push(generators.get_col(l));
        }

        a[0][0] = 4.big();
        a[1][1] = 5.big();
        a[2][0] = 3.big();
        a[2][2] = 3.big();
        a[3][3] = 7.big();
        a[4][0] = 1.big();
        a[5][1] = -2.big();
        a[5][2] = 1.big();
        a[6][2] = 1.big();
        a[7][0] = -1.big();
        a[7][3] = -3.big();

        /*
        matrix:
        4 0 3 0 1 0 0 -1
        0 5 0 0 0 -2 0 0
        0 0 3 0 0 1 1 0
        0 0 0 7 0 0 0 -3
        */
        let mut cmp: Matrix<Integer> = Matrix::zeros(4, 4);
        cmp[(0, 0)] = 1.big();
        cmp[(1, 1)] = 1.big();
        cmp[(2, 2)] = 1.big();
        cmp[(3, 3)] = 1.big();

        let basis = Lattice::hnf_core(a);
        assert_eq!(basis, cmp);
    }

    #[test]
    fn reduce_denom() {
        let mut basis1 = Matrix::eye(4);
        let mut basis2 = Matrix::eye(4);
        let s = 15.big();
        for i in 0..4 {
            for j in 0..4 {
                basis1[(i, j)] = (i + j) * s.clone();
                basis2[(i, j)] = (i + j) * 1.big();
            }
        }
        let lat1 = Lattice::new(basis1, 4.big() * s);
        let lat2 = Lattice::new(basis2, 4.big());

        let reduced = lat1.reduce_denom();
        assert_eq!(reduced.basis, lat2.basis);
    }

    #[test]
    fn dual_without_hnf() {
        let basis1 = Matrix::zeros(4, 4);
        let basis2 = Matrix::zeros(4, 4);
        let mut lat = Lattice::new(basis1, 6.big());
        let mut cmp = Lattice::new(basis2, 1.big());

        lat.basis[(0, 0)] = 1.big();
        lat.basis[(0, 3)] = -1.big();
        lat.basis[(1, 1)] = -2.big();
        lat.basis[(2, 2)] = 1.big();
        lat.basis[(2, 1)] = 1.big();
        lat.basis[(3, 3)] = -3.big();

        cmp.basis[(0, 0)] = 6.big();
        cmp.basis[(1, 1)] = 3.big();
        cmp.basis[(2, 2)] = 6.big();
        cmp.basis[(3, 3)] = 2.big();

        lat = lat.hnf();

        let mut dual = lat.dual_without_hnf();
        dual = dual.hnf();
        cmp = cmp.hnf();
        assert_eq!(dual, cmp);
        assert_ne!(dual, lat);

        dual = dual.dual_without_hnf();
        dual = dual.hnf();
        assert_eq!(dual, lat);
    }

    #[test]
    fn add() {
        let basis1 = Matrix::zeros(4, 4);
        let basis2 = Matrix::zeros(4, 4);
        let basis3 = Matrix::zeros(4, 4);
        let mut lat1 = Lattice::new(basis1, 4.big());
        let mut lat2 = Lattice::new(basis2, 6.big());
        let mut cmp = Lattice::new(basis3, 12.big());

        lat1.basis[(0, 0)] = 44.big();
        lat1.basis[(0, 2)] = 3.big();
        lat1.basis[(0, 3)] = 32.big();
        lat2.basis[(0, 0)] = 1.big();
        cmp.basis[(0, 0)] = 2.big();
        cmp.basis[(0, 2)] = 1.big();
        lat1.basis[(1, 1)] = 5.big();
        lat2.basis[(1, 1)] = 2.big();
        cmp.basis[(1, 1)] = 1.big();
        lat1.basis[(2, 2)] = 3.big();
        lat2.basis[(2, 2)] = 1.big();
        cmp.basis[(2, 2)] = 1.big();
        lat1.basis[(3, 3)] = 1.big();
        lat2.basis[(3, 3)] = 3.big();
        cmp.basis[(3, 3)] = 3.big();

        let sum = lat1.add(lat2.clone());
        assert_eq!(cmp.basis, sum.basis);
        assert_eq!(cmp.denom, sum.denom);

        cmp = lat2.clone();
        cmp = cmp.hnf();
        lat2 = lat2.add(lat2.clone());

        assert_eq!(lat2.basis, cmp.basis);
        assert_eq!(lat2.denom, cmp.denom);
    }

    #[test]
    fn intersect() {
        let basis1 = Matrix::zeros(4, 4);
        let basis2 = Matrix::zeros(4, 4);
        let basis3 = Matrix::zeros(4, 4);
        let mut lat1 = Lattice::new(basis1, 4.big());
        let mut lat2 = Lattice::new(basis2, 6.big());
        let mut cmp = Lattice::new(basis3, 2.big());

        lat1.basis[(0, 0)] = 4.big();
        lat1.basis[(0, 2)] = 3.big();
        lat2.basis[(0, 0)] = 1.big();
        lat2.basis[(0, 3)] = -1.big();
        lat1.basis[(1, 1)] = 5.big();
        lat2.basis[(1, 1)] = -2.big();
        lat1.basis[(2, 2)] = 3.big();
        lat2.basis[(2, 2)] = 1.big();
        lat2.basis[(2, 1)] = 1.big();
        lat1.basis[(3, 3)] = 7.big();
        lat2.basis[(3, 3)] = -3.big();

        lat1 = lat1.hnf();
        lat2 = lat2.hnf();

        cmp.basis[(0, 0)] = 2.big();
        cmp.basis[(0, 2)] = 1.big();
        cmp.basis[(1, 1)] = 10.big();
        cmp.basis[(2, 2)] = 3.big();
        cmp.basis[(3, 3)] = 7.big();

        let inter = lat1.intersect(lat2.clone());
        assert_eq!(inter, cmp);

        lat2 = lat1.intersect(lat2);
        assert_eq!(lat2, cmp);

        cmp = lat1.clone();
        lat1 = lat1.intersect(lat1.clone());
        assert_eq!(lat1, cmp);
    }

    #[test]
    fn mul() {
        let qa = QuatAlg { p: -19.big() };
        let basis1 = Matrix::zeros(4, 4);
        let basis2 = Matrix::zeros(4, 4);
        let basis3 = Matrix::zeros(4, 4);
        let mut lat1 = Lattice::new(basis1.clone(), 4.big());
        let mut lat2 = Lattice::new(basis2.clone(), 6.big());
        let mut cmp = Lattice::new(basis3, 24.big());

        lat1.basis[(0, 0)] = 44.big();
        lat1.basis[(0, 2)] = 3.big();
        lat1.basis[(0, 3)] = 32.big();
        lat2.basis[(0, 0)] = 1.big();
        cmp.basis[(0, 0)] = 1.big();
        lat1.basis[(1, 1)] = 5.big();
        lat2.basis[(1, 1)] = 2.big();
        cmp.basis[(1, 1)] = 1.big();

        lat1.basis[(2, 2)] = 3.big();
        lat2.basis[(2, 2)] = 1.big();
        cmp.basis[(2, 2)] = 1.big();
        lat1.basis[(3, 3)] = 1.big();
        lat2.basis[(3, 3)] = 3.big();
        cmp.basis[(3, 3)] = 1.big();

        let mut prod = lat1.mul(lat2, qa.clone());
        assert_eq!(prod, cmp);

        lat1 = Lattice::new(basis1.clone(), 4.big());
        lat2 = Lattice::new(basis2, 6.big());
        lat1.basis[(0, 0)] = 4.big();
        lat1.basis[(0, 2)] = 3.big();
        lat2.basis[(0, 0)] = 1.big();
        lat2.basis[(0, 3)] = -1.big();
        lat1.basis[(1, 1)] = 5.big();
        lat2.basis[(1, 1)] = -2.big();
        lat1.basis[(2, 2)] = 3.big();
        lat2.basis[(2, 2)] = 1.big();
        lat2.basis[(2, 1)] = 1.big();
        lat1.basis[(3, 3)] = 7.big();
        lat2.basis[(3, 3)] = -3.big();

        prod = lat1.mul(lat2.clone(), qa.clone());
        assert_eq!(prod, cmp);

        cmp = Lattice::new(basis1, 36.big());
        cmp.basis[(0, 0)] = 1.big();
        cmp.basis[(1, 1)] = 1.big();
        cmp.basis[(2, 2)] = 1.big();
        cmp.basis[(3, 3)] = 1.big();

        lat2 = lat2.clone().mul(lat2, qa);
        assert_eq!(lat2, cmp);
    }

    #[test]
    fn contains() {
        let basis = Matrix::zeros(4, 4);
        let mut lat = Lattice::new(basis.clone(), 4.big());

        lat.basis[(0, 0)] = 4.big();
        lat.basis[(0, 2)] = 3.big();
        lat.basis[(1, 1)] = 5.big();
        lat.basis[(2, 2)] = 3.big();
        lat.basis[(3, 3)] = 7.big();

        let qa = QuatAlg { p: -19.big() };
        let mut coord = Matrix::zeros(4, 1);
        let x = QuatAlgEl::new(1.big(), -2.big(), 26.big(), 9.big(), 3.big(), qa.clone());

        assert!(!lat.contains(&mut coord, x.clone()));

        lat.basis[(0, 0)] = 1.big();
        lat.basis[(0, 3)] = -1.big();
        lat.basis[(1, 1)] = -2.big();
        lat.basis[(2, 2)] = 1.big();
        lat.basis[(2, 1)] = 1.big();
        lat.basis[(3, 3)] = -3.big();
        lat.denom = 6.big();

        let lat_hnf = lat.hnf();
        let mut cmp = Matrix::zeros(4, 1);
        cmp[0] = 2.big();
        cmp[1] = -2.big();
        cmp[2] = 52.big();
        cmp[3] = 6.big();

        assert!(lat_hnf.contains(&mut coord, x.clone()));

        let l = (x.clone() * lat_hnf.denom).normalize().coeffs().transpose();
        let r = lat_hnf.basis.clone() * coord.clone();
        assert!(l == r);

        assert!(coord == cmp);

        let order = standard_maximal_extremal_order().order;
        let x = QuatAlgEl::new(1.big(), 0.big(), 0.big(), 0.big(), 1.big(), qa.clone());
        assert!(order.lattice.contains(&mut coord, x.clone()));
    }

    #[test]
    fn index() {
        let basis1 = Matrix::zeros(4, 4);
        let basis2 = Matrix::eye(4);
        let mut sublat = Lattice::new(basis1.clone(), 2.big());
        let overlat = Lattice::new(basis2.clone(), 2.big());

        sublat.basis[(0, 0)] = 2.big();
        sublat.basis[(0, 2)] = 1.big();
        sublat.basis[(1, 1)] = 4.big();
        sublat.basis[(1, 2)] = 2.big();
        sublat.basis[(1, 3)] = 3.big();
        sublat.basis[(2, 2)] = 1.big();
        sublat.basis[(3, 3)] = 1.big();

        let index = sublat.index(overlat);
        assert!(index == 8.big());
    }

    #[test]
    fn hnf() {
        for _ in 0..10000 {
            let basis = Matrix::zeros(4, 4);
            let mut lat = Lattice::new(basis.clone(), 6.big());
            let mut cmp = Lattice::new(basis.clone(), 6.big());

            lat.basis[(0, 0)] = 1.big();
            lat.basis[(0, 3)] = -1.big();
            lat.basis[(1, 1)] = -2.big();
            lat.basis[(2, 2)] = 1.big();
            lat.basis[(2, 1)] = 1.big();
            lat.basis[(3, 3)] = -3.big();
            cmp.basis[(0, 0)] = 1.big();
            cmp.basis[(1, 1)] = 2.big();
            cmp.basis[(2, 2)] = 1.big();
            cmp.basis[(3, 3)] = 3.big();

            lat = lat.hnf();
            assert!(lat == cmp);
        }
    }
}
