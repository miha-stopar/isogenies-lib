use crate::linalg::matrix::Matrix;
use crate::util::{mod_inv, Big};
use rug::Integer;

pub fn right_kernel_mod_prime(m: Matrix<Integer>, p: Integer) -> Matrix<Integer> {
    let mut mat = m.map(|x| x % p.clone());
    mat = mat.map(|x| {
        if x < &0.big() {
            x + p.clone()
        } else {
            x + 0.big()
        }
    });

    let columns = m.n_cols;
    let rows = m.n_rows;
    let mut c = vec![0; rows];
    let mut d = vec![0; columns];

    let mut kernel_dim = 0;
    let mut k = 0;

    while k < columns {
        let mut j = 0;
        while j < rows && (mat[(j, k)] == 0.big() || c[j] != 0) {
            j = j + 1;
        }
        if j == rows {
            kernel_dim += 1;
            d[k] = 0;
            k = k + 1;
        } else {
            let mut prod = mod_inv(mat[(j, k)].clone(), p.clone()).unwrap();
            prod = -prod;
            if prod < 0.big() {
                prod += p.clone();
            }
            mat[(j, k)] = p.clone() - 1;

            for s in k + 1..columns {
                mat[(j, s)] = mat[(j, s)].clone() * prod.clone();
                mat[(j, s)] = mat[(j, s)].clone() % p.clone();
            }
            for i in 0..rows {
                if i != j {
                    let var = mat[(i, k)].clone();
                    mat[(i, k)] = 0.big();
                    for s in k + 1..columns {
                        prod = mat[(j, s)].clone() * var.clone();
                        mat[(i, s)] = mat[(i, s)].clone() + prod;
                        mat[(i, s)] = mat[(i, s)].clone() % p.clone();
                    }
                }
            }
            c[j] = k + 1;
            d[k] = j + 1;
            k = k + 1;
        }
    }

    let mut kernel: Matrix<Integer> = Matrix::zeros(columns, 1);
    if kernel_dim == 1 {
        for k in 0..columns {
            // should be true exactly for 1 k, since kernel_dim = 1
            if d[k] == 0 {
                for s in 0..columns {
                    if s == k {
                        kernel[(s, 0)] = 1.big();
                    }
                    if d[s] > 0 {
                        kernel[(s, 0)] = mat[(d[s] - 1, k)].clone() % p.clone();
                    }
                }
            }
        }
    }

    kernel
}
