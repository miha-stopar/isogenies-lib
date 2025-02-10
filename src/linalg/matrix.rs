use core::ops::{Add, Mul, Neg, Sub};
use num::{rational::Ratio, One, Zero};
use num_bigint::BigInt;
use num_integer::Integer;
use std::{
    fmt,
    iter::Sum,
    ops::{Index, IndexMut},
};

#[derive(PartialEq, Debug, Clone)]
pub struct Scalar<T>(pub T);

#[derive(PartialEq, Debug, Clone)]
pub struct Matrix<T: Clone + PartialEq + Zero + One> {
    pub data: Vec<T>,
    pub n_rows: usize,
    pub n_cols: usize,
}

impl<T: Clone + PartialEq + Zero + One> Matrix<T> {
    /// Create a matrix from vector of vectors.
    pub fn new(mut v: Vec<Vec<T>>) -> Matrix<T> {
        let n_rows = v.len();
        let n_cols = v[0].len();
        let mut data = vec![];
        for row in v.iter_mut() {
            while row.len() > 0 {
                data.push(row.remove(0));
            }
        }
        Matrix {
            data,
            n_rows,
            n_cols,
        }
    }

    /// Create a matrix from a vector.
    pub fn create_from_data(data: Vec<T>, n_rows: usize, n_cols: usize) -> Matrix<T> {
        if data.len() != n_rows * n_cols {
            panic!("not compatible dimension!");
        }
        Matrix {
            data,
            n_rows,
            n_cols,
        }
    }

    /// Obtain the element at row `i` and column `j`.
    pub fn get(&self, i: usize, j: usize) -> T {
        self.data[i * self.n_cols + j].clone()
    }

    /// Obtain the `i`'th row.
    pub fn get_row(&self, i: usize) -> Matrix<T> {
        let mut row = vec![];
        for e in &self.data[i * self.n_cols..(i + 1) * self.n_cols] {
            row.push(e.clone());
        }
        let n_cols = row.len();
        Matrix::create_from_data(row, 1, n_cols)
    }

    /// Set the `r`'th row.
    pub fn set_row(&mut self, r: usize, row: Matrix<T>) {
        for i in 0..self.n_cols {
            self.data[r * self.n_cols + i] = row[i].clone();
        }
    }

    /// Set the `r`'th row by applying function `f` on the elements of the row.
    pub fn set_row_with_func<F>(&mut self, r: usize, mut f: F)
    where
        F: FnMut(&T) -> T,
    {
        for i in 0..self.n_cols {
            self.data[r * self.n_cols + i] = f(&self.data[r * self.n_cols + i])
        }
    }

    /// Obtain the `i`'th column.
    pub fn get_col(&self, i: usize) -> Matrix<T> {
        let col = self
            .iter()
            .filter(|(_, _, col)| *col == i)
            .map(|(e, _, _)| e.clone())
            .collect::<Vec<T>>();
        let n_rows = col.len();
        Matrix::create_from_data(col, n_rows, 1)
    }

    /// Set the `r`'th column.
    pub fn set_col(&mut self, r: usize, col: Matrix<T>) {
        for i in 0..self.n_rows {
            self.data[r + self.n_cols * i] = col[i].clone();
        }
    }

    /// Obtain the submatrix from row `r1` to row `r2-1` and from column `c1` to column `c2-1`.
    pub fn submatrix(&self, r1: usize, c1: usize, r2: usize, c2: usize) -> Matrix<T> {
        let mut data = vec![];
        for i in r1..r2 {
            for j in c1..c2 {
                data.push(self[(i, j)].clone())
            }
        }
        Matrix {
            data,
            n_rows: r2 - r1,
            n_cols: c2 - c1,
        }
    }

    /// Create an immutable iterator over the matrix.
    pub fn iter<'a>(&'a self) -> MatrixIterator<'a, T> {
        MatrixIterator::new(0, 0, 0, &self.data, self.n_rows, self.n_cols)
    }

    /// Create an mutable iterator over the matrix.
    pub fn iter_mut<'a>(&'a mut self) -> MatrixIteratorMut<'a, T> {
        MatrixIteratorMut::new(0, 0, &mut self.data, self.n_rows, self.n_cols)
    }

    /// Compute gcd of the elements in the matrix.
    pub fn gcd(&self) -> T
    where
        T: Integer,
    {
        let mut gcd = self[(0, 0)].clone();
        for i in 0..self.n_rows {
            for j in 0..self.n_cols {
                gcd = self[(i, j)].gcd(&gcd);
            }
        }
        gcd
    }

    /// Divide the elements of the matrix by d.
    pub fn div(&self, d: T) -> Matrix<T>
    where
        T: Integer,
    {
        let mut div_m = self.clone();
        for i in 0..self.n_rows {
            for j in 0..self.n_cols {
                let (n, r) = self[(i, j)].div_rem(&d);
                div_m[(i, j)] = n;
                assert!(r == T::zero());
            }
        }
        div_m
    }

    /// Transpose a copy of the matrix.
    pub fn transpose(&self) -> Matrix<T> {
        let mut data = vec![];
        for j in 0..self.n_cols {
            for i in 0..self.n_rows {
                data.push(self.data[j + i * self.n_cols].clone());
            }
        }
        Matrix::create_from_data(data, self.n_cols, self.n_rows)
    }

    /// Swap two elements.
    pub fn swap(&mut self, ind1: (usize, usize), ind2: (usize, usize)) {
        let el1 = self[ind1].clone();
        let el2 = self[ind2].clone();
        self[ind1] = el2;
        self[ind2] = el1;
    }

    /// Compute the determinant of the matrix.
    pub fn determinant(&self) -> T
    where
        T: Sub<Output = T>,
    {
        assert!(self.n_rows == self.n_cols);
        match self.n_rows {
            0 => T::one(),
            1 => self[(0, 0)].clone(),
            2 => {
                let m11 = self[(0, 0)].clone();
                let m12 = self[(0, 1)].clone();
                let m21 = self[(1, 0)].clone();
                let m22 = self[(1, 1)].clone();

                m11 * m22 - m21 * m12
            }
            3 => {
                let m11 = self[(0, 0)].clone();
                let m12 = self[(0, 1)].clone();
                let m13 = self[(0, 2)].clone();

                let m21 = self[(1, 0)].clone();
                let m22 = self[(1, 1)].clone();
                let m23 = self[(1, 2)].clone();

                let m31 = self[(2, 0)].clone();
                let m32 = self[(2, 1)].clone();
                let m33 = self[(2, 2)].clone();

                let minor_m12_m23 = m22.clone() * m33.clone() - m32.clone() * m23.clone();
                let minor_m11_m23 = m21.clone() * m33 - m31.clone() * m23;
                let minor_m11_m22 = m21 * m32 - m31 * m22;

                m11 * minor_m12_m23 - m12 * minor_m11_m23 + m13 * minor_m11_m22
            }
            _ => todo!(),
        }
    }

    /// Zip two matrices and apply function `F`.
    pub fn zip_mut_with<F>(&mut self, other: Matrix<T>, mut f: F)
    where
        F: FnMut(&mut T, T),
    {
        for i in 0..self.n_cols {
            for j in 0..self.n_rows {
                f(&mut self[(i, j)], other[(i, j)].clone());
            }
        }
    }

    /// Apply function `f` on the matrix.
    pub fn map<F, K>(&self, mut f: F) -> Matrix<K>
    where
        F: FnMut(&T) -> K,
        K: Clone + PartialEq + Zero + One,
    {
        let new_data: Vec<K> = self.iter().map(|e| f(e.0)).collect();
        Matrix::create_from_data(new_data, self.n_rows, self.n_cols)
    }

    /// Create a matrix of all zeros.
    pub fn zeros(n_rows: usize, n_cols: usize) -> Matrix<T> {
        let mut data = vec![];
        let zero = T::zero();
        for _ in 0..n_cols {
            for _ in 0..n_rows {
                data.push(zero.clone());
            }
        }
        Matrix::create_from_data(data, n_rows, n_cols)
    }

    /// Create a matrix of all rational zeros.
    pub fn zeros_rational(n_rows: usize, n_cols: usize) -> Matrix<Ratio<T>>
    where
        T: Integer,
    {
        let mut data = vec![];
        let zero_r = Ratio::new(T::zero(), T::one());
        for _ in 0..n_cols {
            for _ in 0..n_rows {
                data.push(zero_r.clone());
            }
        }
        Matrix::create_from_data(data, n_rows, n_cols)
    }

    /// Create an identity matrix.
    pub fn eye(n: usize) -> Matrix<T> {
        let mut data = vec![];
        let zero = T::zero();
        let one = T::one();
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    data.push(one.clone());
                } else {
                    data.push(zero.clone());
                }
            }
        }
        Matrix::create_from_data(data, n, n)
    }

    /// Concatenate two matrices.
    pub fn concatenate(m1: Matrix<T>, m2: Matrix<T>) -> Matrix<T> {
        assert!(m1.n_rows == m2.n_rows);
        let mut data = vec![];
        for i in 0..m1.n_rows {
            for j in 0..m1.n_cols {
                data.push(m1[(i, j)].clone());
            }
            for j in 0..m2.n_cols {
                data.push(m2[(i, j)].clone());
            }
        }
        Matrix::create_from_data(data, m1.n_rows, m1.n_cols + m2.n_cols)
    }

    /// Put one matrix on the top of another.
    pub fn stack(m1: Matrix<T>, m2: Matrix<T>) -> Matrix<T> {
        assert!(m1.n_cols == m2.n_cols);
        let mut data1 = m1.data.clone();
        let mut data2 = m2.data.clone();
        data1.append(&mut data2);
        Matrix::create_from_data(data1, m1.n_rows + m2.n_rows, m1.n_cols)
    }

    /// Get least common denominator.
    pub fn lcm(m: &Matrix<Ratio<T>>) -> T
    where
        T: Integer,
    {
        m.iter()
            .fold(T::one(), |accum, item| accum.lcm(item.0.denom()))
    }

    /// Create matrix with denominators cleared.
    pub fn clear_denoms(m: &Matrix<Ratio<T>>) -> Matrix<T>
    where
        T: Integer,
    {
        let lcm = Self::lcm(m);
        let data = m
            .data
            .iter()
            .map(|x| (x * lcm.clone()).numer().clone())
            .collect();
        Matrix::create_from_data(data, m.n_rows, m.n_cols)
    }

    pub fn inverse(m: &Matrix<T>) -> Matrix<Ratio<T>>
    where
        T: Integer + Neg<Output = T>,
    {
        assert!(m.n_rows == m.n_cols);
        match m.n_rows {
            2 => {
                let m11 = m[(0, 0)].clone();
                let m12 = m[(0, 1)].clone();
                let m21 = m[(1, 0)].clone();
                let m22 = m[(1, 1)].clone();

                let determinant = m11.clone() * m22.clone() - m21.clone() * m12.clone();
                assert!(determinant != T::zero());

                let mut m_inv = Matrix::zeros_rational(m.n_rows, m.n_cols);

                m_inv[(0, 0)] = Ratio::new(m22, determinant.clone());
                m_inv[(0, 1)] = Ratio::new(-m12, determinant.clone());
                m_inv[(1, 0)] = Ratio::new(-m21, determinant.clone());
                m_inv[(1, 1)] = Ratio::new(m11, determinant.clone());

                m_inv
            }
            3 => {
                let m11 = m[(0, 0)].clone();
                let m12 = m[(0, 1)].clone();
                let m13 = m[(0, 2)].clone();

                let m21 = m[(1, 0)].clone();
                let m22 = m[(1, 1)].clone();
                let m23 = m[(1, 2)].clone();

                let m31 = m[(2, 0)].clone();
                let m32 = m[(2, 1)].clone();
                let m33 = m[(2, 2)].clone();

                let minor_m12_m23 = m22.clone() * m33.clone() - m32.clone() * m23.clone();
                let minor_m11_m23 = m21.clone() * m33.clone() - m31.clone() * m23.clone();
                let minor_m11_m22 = m21.clone() * m32.clone() - m31.clone() * m22.clone();

                let determinant = m11.clone() * minor_m12_m23.clone()
                    - m12.clone() * minor_m11_m23.clone()
                    + m13.clone() * minor_m11_m22.clone();

                assert!(determinant != T::zero());

                let mut m_inv = Matrix::zeros_rational(m.n_rows, m.n_cols);

                m_inv[(0, 0)] = Ratio::new(minor_m12_m23, determinant.clone());
                m_inv[(0, 1)] = Ratio::new(
                    m13.clone() * m32.clone() - m33.clone() * m12.clone(),
                    determinant.clone(),
                );
                m_inv[(0, 2)] = Ratio::new(
                    m12.clone() * m23.clone() - m22.clone() * m13.clone(),
                    determinant.clone(),
                );
                m_inv[(1, 0)] = Ratio::new(-minor_m11_m23, determinant.clone());
                m_inv[(1, 1)] = Ratio::new(
                    m11.clone() * m33 - m31.clone() * m13.clone(),
                    determinant.clone(),
                );
                m_inv[(1, 2)] =
                    Ratio::new(m13 * m21.clone() - m23 * m11.clone(), determinant.clone());
                m_inv[(2, 0)] = Ratio::new(minor_m11_m22, determinant.clone());
                m_inv[(2, 1)] =
                    Ratio::new(m12.clone() * m31 - m32 * m11.clone(), determinant.clone());
                m_inv[(2, 2)] = Ratio::new(m11 * m22 - m21 * m12, determinant.clone());

                m_inv
            }
            _ => todo!(),
        }
    }

    fn minor_2(a11: T, a12: T, a21: T, a22: T) -> T
    where
        T: Integer + Neg<Output = T>,
    {
        a11 * a22 - a12 * a21
    }

    /// Helper 4x4 inverse function.
    fn pmp(a1: T, a2: T, b1: T, b2: T, c1: T, c2: T) -> T
    where
        T: Integer + Neg<Output = T>,
    {
        a1 * a2 - b1 * b2 + c1 * c2
    }

    /// Helper 4x4 inverse function.
    fn mpm(a1: T, a2: T, b1: T, b2: T, c1: T, c2: T) -> T
    where
        T: Integer + Neg<Output = T>,
    {
        b1 * b2 - a1 * a2 - c1 * c2
    }

    pub fn inverse_with_det_as_denom(m: &Matrix<T>) -> (Matrix<T>, T)
    where
        T: Integer + Neg<Output = T> + Zero + std::fmt::Debug,
    {
        let mut s: [T; 6] = [
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
        ];
        let mut c: [T; 6] = [
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
        ];
        for i in 0..3 {
            s[i] = Self::minor_2(
                m[(0, 0)].clone(),
                m[(0, i + 1)].clone(),
                m[(1, 0)].clone(),
                m[(1, i + 1)].clone(),
            );
            c[i] = Self::minor_2(
                m[(2, 0)].clone(),
                m[(2, i + 1)].clone(),
                m[(3, 0)].clone(),
                m[(3, i + 1)].clone(),
            );
        }
        for i in 0..2 {
            s[3 + i] = Self::minor_2(
                m[(0, 1)].clone(),
                m[(0, 2 + i)].clone(),
                m[(1, 1)].clone(),
                m[(1, 2 + i)].clone(),
            );
            c[3 + i] = Self::minor_2(
                m[(2, 1)].clone(),
                m[(2, 2 + i)].clone(),
                m[(3, 1)].clone(),
                m[(3, 2 + i)].clone(),
            );
        }
        s[5] = Self::minor_2(
            m[(0, 2)].clone(),
            m[(0, 3)].clone(),
            m[(1, 2)].clone(),
            m[(1, 3)].clone(),
        );
        c[5] = Self::minor_2(
            m[(2, 2)].clone(),
            m[(2, 3)].clone(),
            m[(3, 2)].clone(),
            m[(3, 3)].clone(),
        );

        //compute det
        let mut work_det: T = T::zero();
        for i in 0..6 {
            let prod = s[i].clone() * c[5 - i].clone();
            if i != 1 && i != 4 {
                work_det = work_det + prod;
            } else {
                work_det = work_det - prod;
            }
        }
        let mut work = Matrix::eye(4);
        if work_det != T::zero() {
            // compute transposed adjugate
            for j in 0..4 {
                for k in 0..2 {
                    if (k + j + 1) % 2 == 1 {
                        work[(j, k)] = Self::pmp(
                            m[(1 - k, (j == 0) as usize)].clone(),
                            c[6 - j - (j == 0) as usize].clone(),
                            m[(1 - k, 2 - (j > 1) as usize)].clone(),
                            c[4 - j - (j == 1) as usize].clone(),
                            m[(1 - k, 3 - (j == 3) as usize)].clone(),
                            c[3 - j - (j == 1) as usize - (j == 2) as usize].clone(),
                        );
                    } else {
                        work[(j, k)] = Self::mpm(
                            m[(1 - k, (j == 0) as usize)].clone(),
                            c[6 - j - (j == 0) as usize].clone(),
                            m[(1 - k, 2 - (j > 1) as usize)].clone(),
                            c[4 - j - (j == 1) as usize].clone(),
                            m[(1 - k, 3 - (j == 3) as usize)].clone(),
                            c[3 - j - (j == 1) as usize - (j == 2) as usize].clone(),
                        );
                    }
                }
                for k in 2..4 {
                    if (k + j + 1) % 2 == 1 {
                        work[(j, k)] = Self::pmp(
                            m[(3 - (k == 3) as usize, (j == 0) as usize)].clone(),
                            s[6 - j - (j == 0) as usize].clone(),
                            m[(3 - (k == 3) as usize, 2 - (j > 1) as usize)].clone(),
                            s[4 - j - (j == 1) as usize].clone(),
                            m[(3 - (k == 3) as usize, 3 - (j == 3) as usize)].clone(),
                            s[3 - j - (j == 1) as usize - (j == 2) as usize].clone(),
                        );
                    } else {
                        work[(j, k)] = Self::mpm(
                            m[(3 - (k == 3) as usize, (j == 0) as usize)].clone(),
                            s[6 - j - (j == 0) as usize].clone(),
                            m[(3 - (k == 3) as usize, 2 - (j > 1) as usize)].clone(),
                            s[4 - j - (j == 1) as usize].clone(),
                            m[(3 - (k == 3) as usize, 3 - (j == 3) as usize)].clone(),
                            s[3 - j - (j == 1) as usize - (j == 2) as usize].clone(),
                        );
                    }
                }
            }
        }

        (work, work_det)
    }
}

pub struct MatrixIteratorMut<'a, T: Clone> {
    row_idx: usize,
    col_idx: usize,
    data: &'a mut [T],
    n_rows: usize,
    n_cols: usize,
}

impl<'a, T: Clone> MatrixIteratorMut<'a, T> {
    pub fn new(
        row_idx: usize,
        col_idx: usize,
        data: &'a mut [T],
        n_rows: usize,
        n_cols: usize,
    ) -> MatrixIteratorMut<'a, T> {
        MatrixIteratorMut {
            row_idx,
            col_idx,
            data,
            n_rows,
            n_cols,
        }
    }
}

impl<'a, T: Clone> std::iter::Iterator for MatrixIteratorMut<'a, T> {
    type Item = (&'a mut T, usize, usize);
    fn next(&mut self) -> Option<Self::Item> {
        if self.col_idx == self.n_cols {
            self.row_idx += 1;
            self.col_idx = 0;
        }
        return if self.row_idx == self.n_rows {
            None
        } else {
            let col_idx = self.col_idx;
            self.col_idx += 1;
            let data = std::mem::replace(&mut self.data, &mut []);
            if let Some((v, rest)) = data.split_first_mut() {
                self.data = rest;
                Some((v, self.row_idx, col_idx))
            } else {
                None
            }
        };
    }
}

pub struct MatrixIterator<'a, T: Clone> {
    row_idx: usize,
    col_idx: usize,
    idx: usize,
    data: &'a [T],
    n_rows: usize,
    n_cols: usize,
}

impl<'a, T: Clone> MatrixIterator<'a, T> {
    pub fn new(
        row_idx: usize,
        col_idx: usize,
        idx: usize,
        data: &'a [T],
        n_rows: usize,
        n_cols: usize,
    ) -> MatrixIterator<'a, T> {
        MatrixIterator {
            row_idx,
            col_idx,
            idx,
            data,
            n_rows,
            n_cols,
        }
    }
}

impl<'a, T: Clone> std::iter::Iterator for MatrixIterator<'a, T> {
    type Item = (&'a T, usize, usize);
    fn next(&mut self) -> Option<Self::Item> {
        if self.col_idx == self.n_cols {
            self.row_idx += 1;
            self.col_idx = 0;
        }
        return if self.row_idx == self.n_rows {
            None
        } else {
            let col_idx = self.col_idx;
            let idx = self.idx;
            self.col_idx += 1;
            self.idx += 1;
            Some((&self.data[idx], self.row_idx, col_idx))
        };
    }
}

macro_rules! matrix_add {
    ($LHS:ty, $RHS:ty, $ScalarType:tt ) => {
        impl<$ScalarType: Add<Output = $ScalarType> + Clone + PartialEq + Zero + One> Add<$RHS>
            for $LHS
        {
            type Output = Matrix<$ScalarType>;
            fn add(self, rhs: $RHS) -> Self::Output {
                let mut res: Vec<$ScalarType> = vec![];
                for (e, row, col) in self.iter() {
                    res.push(e.clone() + rhs.get(row, col));
                }
                Matrix::create_from_data(res, self.n_rows, self.n_cols)
            }
        }
    };
}
matrix_add!(&Matrix<T>, &Matrix<T>, T);
matrix_add!(Matrix<T>, Matrix<T>, T);
matrix_add!(&Matrix<T>, Matrix<T>, T);
matrix_add!(Matrix<T>, &Matrix<T>, T);

macro_rules! matrix_mult {
    ($LHS:ty, $RHS:ty, $ScalarType:tt ) => {
        impl<
                $ScalarType: Sum
                    + Mul<Output = $ScalarType>
                    + Add<Output = $ScalarType>
                    + Clone
                    + PartialEq
                    + Zero
                    + One,
            > Mul<$RHS> for $LHS
        {
            type Output = Matrix<$ScalarType>;
            fn mul(self, rhs: $RHS) -> Self::Output {
                if rhs.n_rows != self.n_cols {
                    panic!("dimensions do not match!");
                }
                let mut res: Vec<$ScalarType> = vec![];
                for i in 0..self.n_rows {
                    for j in 0..rhs.n_cols {
                        res.push(
                            self.data[i * self.n_cols..(i + 1) * self.n_cols]
                                .iter()
                                .enumerate()
                                .map(|(rhs_row_idx, v)| {
                                    v.clone() * rhs.data[j + rhs_row_idx * rhs.n_cols].clone()
                                })
                                .sum(),
                        );
                    }
                }
                Matrix::create_from_data(res, self.n_rows, rhs.n_cols)
            }
        }
    };
}
matrix_mult!(&Matrix<T>, &Matrix<T>, T);
matrix_mult!(Matrix<T>, Matrix<T>, T);
matrix_mult!(&Matrix<T>, Matrix<T>, T);
matrix_mult!(Matrix<T>, &Matrix<T>, T);

macro_rules! matrix_neg {
    ($RHS:ty, $ScalarType:tt ) => {
        impl<$ScalarType: Neg<Output = $ScalarType> + Clone + PartialEq + Zero + One> Neg for $RHS {
            type Output = Matrix<$ScalarType>;
            fn neg(self) -> Self::Output {
                let mut res: Vec<$ScalarType> = vec![];
                for (e, _, _) in self.iter() {
                    res.push(-e.clone());
                }
                Matrix::create_from_data(res, self.n_rows, self.n_cols)
            }
        }
    };
}
matrix_neg!(&Matrix<T>, T);
matrix_neg!(Matrix<T>, T);

macro_rules! scalar_mult {
    ($LHS:ty, $RHS:ty, $ScalarType:tt ) => {
        impl<$ScalarType: Mul<Output = $ScalarType> + Clone> Mul<$RHS> for $LHS {
            type Output = Scalar<$ScalarType>;
            fn mul(self, rhs: $RHS) -> Self::Output {
                Scalar(self.0.clone() * rhs.0.clone())
            }
        }
    };
}
scalar_mult!(Scalar<T>, Scalar<T>, T);
scalar_mult!(&Scalar<T>, Scalar<T>, T);
scalar_mult!(&Scalar<T>, &Scalar<T>, T);
scalar_mult!(Scalar<T>, &Scalar<T>, T);

macro_rules! scalar_matrix_mult {
    ($LHS:ty, $RHS:ty, $ScalarType:tt ) => {
        impl<$ScalarType: Mul<Output = $ScalarType> + Clone + PartialEq + Zero + One> Mul<$RHS>
            for $LHS
        {
            type Output = Matrix<$ScalarType>;
            fn mul(self, rhs: $RHS) -> Self::Output {
                let mut res: Vec<$ScalarType> = vec![];
                for v in rhs.data.iter() {
                    res.push(self.0.clone() * v.clone());
                }
                Matrix::create_from_data(res, rhs.n_rows, rhs.n_cols)
            }
        }
    };
}
scalar_matrix_mult!(Scalar<T>, Matrix<T>, T);
scalar_matrix_mult!(Scalar<T>, &Matrix<T>, T);
scalar_matrix_mult!(&Scalar<T>, Matrix<T>, T);
scalar_matrix_mult!(&Scalar<T>, &Matrix<T>, T);

macro_rules! scalar_type_matrix_mult {
    ($LHS:ty, $ScalarType:tt ) => {
        impl<$ScalarType: Mul<Output = $ScalarType> + Clone + PartialEq + Zero + One>
            Mul<$ScalarType> for $LHS
        {
            type Output = Matrix<$ScalarType>;
            fn mul(self, rhs: $ScalarType) -> Self::Output {
                let mut res: Vec<$ScalarType> = vec![];
                for v in self.data.iter() {
                    res.push(rhs.clone() * v.clone());
                }
                Matrix::create_from_data(res, self.n_rows, self.n_cols)
            }
        }
    };
}
scalar_type_matrix_mult!(Matrix<T>, T);
scalar_type_matrix_mult!(&Matrix<T>, T);

macro_rules! matrix_subtraction {
    ($LHS:ty, $RHS:ty, $ScalarType:tt ) => {
        impl<$ScalarType: Sub<Output = $ScalarType> + Clone + PartialEq + Zero + One> Sub<$RHS>
            for $LHS
        {
            type Output = Matrix<$ScalarType>;
            fn sub(self, rhs: $RHS) -> Self::Output {
                let mut res: Vec<$ScalarType> = vec![];
                for (e, row, col) in self.iter() {
                    res.push(e.clone() - rhs.get(row, col));
                }
                Matrix::create_from_data(res, self.n_rows, self.n_cols)
            }
        }
    };
}
matrix_subtraction!(&Matrix<T>, &Matrix<T>, T);
matrix_subtraction!(Matrix<T>, &Matrix<T>, T);
matrix_subtraction!(&Matrix<T>, Matrix<T>, T);
matrix_subtraction!(Matrix<T>, Matrix<T>, T);

impl<T> Index<usize> for Matrix<T>
where
    T: Clone + PartialEq + Zero + One,
{
    type Output = T;

    fn index(&self, ind: usize) -> &T {
        &self.data[ind]
    }
}

impl<T> Index<(usize, usize)> for Matrix<T>
where
    T: Clone + PartialEq + Zero + One,
{
    type Output = T;

    fn index(&self, ind: (usize, usize)) -> &T {
        &self[ind.0 * self.n_cols + ind.1]
    }
}

impl<T> IndexMut<usize> for Matrix<T>
where
    T: Clone + PartialEq + Zero + One,
{
    fn index_mut(self: &mut Matrix<T>, ind: usize) -> &mut T {
        &mut self.data[ind]
    }
}

impl<T> IndexMut<(usize, usize)> for Matrix<T>
where
    T: Clone + PartialEq + Zero + One,
{
    fn index_mut(self: &mut Matrix<T>, ind: (usize, usize)) -> &mut T {
        &mut self.data[ind.0 * self.n_cols + ind.1]
    }
}

impl<T: Clone + PartialEq + Zero + One> IntoIterator for Matrix<T> {
    type Item = Matrix<T>;
    type IntoIter = MatrixIntoIterator<T>;

    fn into_iter(self) -> Self::IntoIter {
        MatrixIntoIterator {
            index: 0,
            matrix: self,
        }
    }
}

impl<T: Clone + PartialEq + Zero + One + Integer> Into<Matrix<Ratio<T>>> for Matrix<T> {
    fn into(self) -> Matrix<Ratio<T>> {
        self.map(|e| Ratio::from(e.clone()))
    }
}

pub struct MatrixIntoIterator<T: Clone + PartialEq + Zero + One> {
    matrix: Matrix<T>,
    index: usize,
}

impl<T: Clone + PartialEq + Zero + One> Iterator for MatrixIntoIterator<T> {
    type Item = Matrix<T>;
    fn next(&mut self) -> Option<Matrix<T>> {
        if self.index < self.matrix.n_rows {
            self.index += 1;
            Some(self.matrix.get_row(self.index - 1))
        } else {
            None
        }
    }
}

impl fmt::Display for Matrix<BigInt> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for i in 0..self.n_rows {
            let mut row = String::new();

            for j in 0..self.n_cols {
                let idx = i * self.n_cols + j;
                let el = &self.data[idx];
                row.push_str(&el.to_string());
                row.push_str(" ");
            }
            if i != 0 {
                write!(f, "\n")?;
            }
            write!(f, "{}", row)?;
        }

        Ok(())
    }
}

impl fmt::Display for Matrix<Ratio<BigInt>> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for i in 0..self.n_rows {
            let mut row = String::new();

            for j in 0..self.n_cols {
                let idx = i * self.n_cols + j;
                let el = &self.data[idx];
                let n = el.numer();
                let d = el.denom();
                let mut s = n.to_string();
                if d != &BigInt::one() {
                    s += &("/".to_owned() + &d.to_string());
                }
                row.push_str(&s);
                row.push_str(" ");
            }
            if i != 0 {
                write!(f, "\n")?;
            }
            write!(f, "{}", row)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iter() {
        let m = Matrix::new(vec![vec![0, 1], vec![2, 3]]);
        assert_eq!(
            m.iter().collect::<Vec<(&u32, usize, usize)>>(),
            vec![(&0, 0, 0), (&1, 0, 1), (&2, 1, 0), (&3, 1, 1)]
        );
    }

    #[test]
    fn test_iter_mut() {
        let mut m = Matrix::new(vec![vec![0, 1], vec![2, 3]]);
        assert_eq!(
            m.iter_mut().collect::<Vec<(&mut u32, usize, usize)>>(),
            vec![
                (&mut 0, 0, 0),
                (&mut 1, 0, 1),
                (&mut 2, 1, 0),
                (&mut 3, 1, 1)
            ]
        );
    }

    #[test]
    fn test_get_row_1() {
        let m = Matrix::new(vec![vec![0, 1], vec![2, 3]]);
        assert_eq!(m.get_row(0), Matrix::create_from_data(vec![0, 1], 1, 2));
    }

    #[test]
    fn test_get_row_2() {
        let m = Matrix::new(vec![vec![0, 1], vec![2, 3]]);
        assert_eq!(m.get_row(1), Matrix::create_from_data(vec![2, 3], 1, 2));
    }

    #[test]
    fn test_get_col_1() {
        let m = Matrix::new(vec![vec![0, 1], vec![2, 3]]);
        assert_eq!(m.get_col(0), Matrix::create_from_data(vec![0, 2], 2, 1));
    }

    #[test]
    fn test_get_col_2() {
        let m = Matrix::new(vec![vec![0, 1], vec![2, 3]]);
        assert_eq!(m.get_col(1), Matrix::create_from_data(vec![1, 3], 2, 1));
    }

    #[test]
    fn test_trans() {
        let m = Matrix::new(vec![vec![0, 1], vec![2, 3]]);
        assert_eq!(m.transpose(), Matrix::new(vec![vec![0, 2], vec![1, 3]]));
    }
}
