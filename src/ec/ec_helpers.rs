#![allow(non_snake_case)]

macro_rules! define_ec_helpers {
    () => {
        use std::fs;
        use crate::linalg::matrix::Matrix;
        // use std::ops::Div;
        use crate::util::{big_to_bytes, bytes_from_str, Big};
        use rug::integer::Order;
        use num_traits::Pow;
        use rand::prelude::*;
        use rand_chacha::ChaCha20Rng;

        /// The m-torsion group (denoted E[m]) is a group of points P for which it holds mP = 0.
        /// E[m] is isomorphic to Z_m x Z_m, let's denote P and Q the generators of E[m].
        /// If we know how the endomorphism maps the generators P and Q,
        /// we can construct a 2x2 matrix that corresponds to the endomorphism on E[m].
        pub fn from_endomorphism_to_matrix(coord: Matrix<Integer>, imprim: Integer, torsion: Integer, mat2: Matrix<Integer>, mat3: Matrix<Integer>, mat4: Matrix<Integer>) -> Matrix<Integer> {
            // 2x2 matrix that corresponds to:
            // coord[0] * I + coord[1] * mat2 + coord[2] * mat3 + coord[3] * mat4
            let mut mat: Matrix<Integer> = Matrix::zeros(2, 2);

            for i in 0..2 {
                mat[(i, i)] += coord[0].clone();
                for j in 0..2 {
                    mat[(i, j)] += mat2[(i, j)].clone() * coord[1].clone()
                        + mat3[(i, j)].clone() * coord[2].clone()
                        + mat4[(i, j)].clone() * coord[3].clone();
                    mat[(i, j)] *= imprim.clone();
                }
            }

            // reduction mod torsion
            for i in 0..2 {
                for j in 0..2 {
                    mat[(i, j)] = mat[(i, j)].clone() % torsion.clone();
                }
            }

            mat
        }

        pub fn get_montgomery_A24(P: &PointX, Q: &PointX, PmQ: &PointX) -> (Fq, Fq) {
            let xP = &P.X / &P.Z;
            let xQ = &Q.X / &Q.Z;
            let xR = &PmQ.X / &PmQ.Z;
            let two_fp = Fp::from_u32(2);
            let four_fp = Fp::from_u32(4);
            let two = Fq::new(&two_fp, &Fp::ZERO);
            let four = Fq::new(&four_fp, &Fp::ZERO);
            let mut A = (&xR * &xP + &xR * &xQ + &xP * &xQ - Fq::ONE).square() / (four.clone() * &xP * &xQ * &xR);
            A += -(&xP + &xQ + &xR);

            (A + two, four)
        }

        pub fn endomorphism_matrix_to_bytes(mat: Matrix<Integer>, torsion: Integer) -> Vec<Vec<u8>> {
            let mat00 = big_to_bytes(mat[(0, 0)].clone());
            let mat10 = big_to_bytes(mat[(1, 0)].clone());
            let mat01 = big_to_bytes(mat[(0, 1)].clone());
            let mat11 = big_to_bytes(mat[(1, 1)].clone());

            let mut mat00mat01 = mat[(0, 0)].clone() - mat[(0, 1)].clone();
            mat00mat01 = mat00mat01.clone() % torsion.clone();
            let mat00mat01 = big_to_bytes(mat00mat01);
            let mut mat10mat11 = mat[(1, 0)].clone() - mat[(1, 1)].clone();
            mat10mat11 = mat10mat11.clone() % torsion.clone();
            let mat10mat11 = big_to_bytes(mat10mat11);

            vec![mat00, mat01, mat10, mat11, mat00mat01, mat10mat11]
        }

        pub fn endomorphism_matrix_to_digits(mat: Matrix<Integer>, torsion: Integer) -> Vec<Integer> {
            // TODO
            let mat00 = mat[(0, 0)].clone();
            let mat10 = mat[(1, 0)].clone();
            let mat01 = mat[(0, 1)].clone();
            let mat11 = mat[(1, 1)].clone();

            let mut mat00mat01 = mat[(0, 0)].clone() - mat[(0, 1)].clone();
            mat00mat01 = mat00mat01.clone() % torsion.clone();
            let mut mat10mat11 = mat[(1, 0)].clone() - mat[(1, 1)].clone();
            mat10mat11 = mat10mat11.clone() % torsion.clone();

            vec![mat00, mat01, mat10, mat11, mat00mat01, mat10mat11]
        }

        pub fn load_torsion_info(fileName: &str, scheme: &str, torsion: u8) -> (PointX, PointX, PointX, Matrix<Integer>, Matrix<Integer>, Matrix<Integer>, u64) {
            let file = fs::File::open(fileName)
                .expect("file should open read only");
            let json: serde_json::Value =
                serde_json::from_reader(file).expect("JSON was not well-formatted");

            let torsion_name = format!("torsion_{}", torsion);
            let j = json[scheme][torsion_name].clone();
            let power = j["power"].as_u64().unwrap();

            let Px = Fq::new(
                &Fp::decode_reduce(&bytes_from_str(&j["PxRe"].as_str().unwrap())),
                &Fp::decode_reduce(&bytes_from_str(&j["PxIm"].as_str().unwrap())),
            );

            let Pz = Fq::new(
                &Fp::decode_reduce(&bytes_from_str(&j["PzRe"].as_str().unwrap())),
                &Fp::decode_reduce(&bytes_from_str(&j["PzIm"].as_str().unwrap())),
            );

            let Qx = Fq::new(
                &Fp::decode_reduce(&bytes_from_str(&j["QxRe"].as_str().unwrap())),
                &Fp::decode_reduce(&bytes_from_str(&j["QxIm"].as_str().unwrap())),
            );

            let Qz = Fq::new(
                &Fp::decode_reduce(&bytes_from_str(&j["QzRe"].as_str().unwrap())),
                &Fp::decode_reduce(&bytes_from_str(&j["QzIm"].as_str().unwrap())),
            );

            let PmQx = Fq::new(
                &Fp::decode_reduce(&bytes_from_str(&j["PmQxRe"].as_str().unwrap())),
                &Fp::decode_reduce(&bytes_from_str(&j["PmQxIm"].as_str().unwrap())),
            );

            let PmQz = Fq::new(
                &Fp::decode_reduce(&bytes_from_str(&j["PmQzRe"].as_str().unwrap())),
                &Fp::decode_reduce(&bytes_from_str(&j["PmQzIm"].as_str().unwrap())),
            );

            let P = PointX::new_xz(&Px, &Pz);
            let Q = PointX::new_xz(&Qx, &Qz);
            let PmQ = PointX::new_xz(&PmQx, &PmQz);

            // 2x2 matrix that corresponds to the second generator of the order
            let mut mat2: Matrix<Integer> = Matrix::zeros(2, 2);
            mat2[(0, 0)] = j["action_gen_2"][0].as_str().unwrap().big();
            mat2[(0, 1)] = j["action_gen_2"][1].as_str().unwrap().big();
            mat2[(1, 0)] = j["action_gen_2"][2].as_str().unwrap().big();
            mat2[(1, 1)] = j["action_gen_2"][3].as_str().unwrap().big();

            // 2x2 matrix that corresponds to the third generator of the order
            let mut mat3: Matrix<Integer> = Matrix::zeros(2, 2);
            mat3[(0, 0)] = j["action_gen_3"][0].as_str().unwrap().big();
            mat3[(0, 1)] = j["action_gen_3"][1].as_str().unwrap().big();
            mat3[(1, 0)] = j["action_gen_3"][2].as_str().unwrap().big();
            mat3[(1, 1)] = j["action_gen_3"][3].as_str().unwrap().big();

            // 2x2 matrix that corresponds to the fourth generator of the order
            let mut mat4: Matrix<Integer> = Matrix::zeros(2, 2);
            mat4[(0, 0)] = j["action_gen_4"][0].as_str().unwrap().big();
            mat4[(0, 1)] = j["action_gen_4"][1].as_str().unwrap().big();
            mat4[(1, 0)] = j["action_gen_4"][2].as_str().unwrap().big();
            mat4[(1, 1)] = j["action_gen_4"][3].as_str().unwrap().big();

            (P, Q, PmQ, mat2, mat3, mat4, power)
        }

        /// The m-torsion group (denoted E[m]) is a group of points P for which it holds mP = 0.
        /// E[m] is isomorphic to Z_m x Z_m, let's denote P and Q the generators of E[m].
        /// If we know how the endomorphism maps the generators P and Q,
        /// we can construct a 2x2 matrix that corresponds to the endomorphism on E[m].
        /// This function applies the endomorphism gamma to the points P = k * P0, Q = l * Q0.
        /// It returns gamma(P), gamma(Q), gamma(P - Q).
        pub fn apply_endomorphism_on_torsion_group(curve: &Curve, coord: Matrix<Integer>, imprim: Integer, torsion: Integer, mat2: Matrix<Integer>, mat3: Matrix<Integer>, mat4: Matrix<Integer>, P: &Point, Q: &Point, PmQ: &Point) -> (Point, Point, Point) {
            /*
            M = |mat00 mat10|
                |mat01 mat11|
            
            v = |P| = |k * P0|
                |Q|   |l * Q0|
            */

            let mat = from_endomorphism_to_matrix(coord, imprim, torsion.clone(), mat2, mat3, mat4);
            // TODO
            let mat_bytes = endomorphism_matrix_to_bytes(mat.clone(), torsion.clone());
            let mat_bytes1 = endomorphism_matrix_to_digits(mat, torsion.clone());

            let mat00 = &mat_bytes[0];
            let mat01 = &mat_bytes[1];
            let mat10 = &mat_bytes[2];
            let mat11 = &mat_bytes[3];
            let mat00mat01 = &mat_bytes[4];
            let mat10mat11 = &mat_bytes[5];

            let amat00 = &mat_bytes1[0];
            let amat01 = &mat_bytes1[1];
            let amat10 = &mat_bytes1[2];
            let amat11 = &mat_bytes1[3];
            let amat00mat01 = &mat_bytes1[4];
            let amat10mat11 = &mat_bytes1[5];


            let P_new1 = curve.mul(P, &mat00, mat00.len() * 8);
            let P_new2 = curve.mul(Q, &mat10, mat10.len() * 8);
            // mat00 * P + mat10 * Q
            let P_new = curve.add(&P_new1, &P_new2);

            /*
            let mut k_digits = amat00.to_digits::<u64>(Order::MsfLe);
            k_digits.reverse();
            let mut l_digits = amat10.to_digits::<u64>(Order::MsfLe);
            l_digits.reverse();
            let f: usize = 36; // TODO
            let aP_new = curve.xdblmul_bounded(&P, &k_digits, &Q, &l_digits, &PmQ, f);
            */

            // assert!(P_new.equals(&aP_new) == 0xFFFFFFFF);



            let Q_new1 = curve.mul(&P, &mat01, mat01.len() * 8);
            let Q_new2 = curve.mul(&Q, &mat11, mat11.len() * 8);
            // mat01 * P + mat11 * Q
            let Q_new = curve.add(&Q_new1, &Q_new2);
            
            let PmQ_new1 = curve.mul(P, &mat00mat01, mat00mat01.len() * 8);
            let PmQ_new2 = curve.mul(Q, &mat10mat11, mat10mat11.len() * 8);
            // (mat00 * P + mat10 * Q) - (mat10 * P + mat11 * Q) = (mat00 - mat10) * P + (mat10 - mat11) * Q
            let PmQ_new = curve.add(&PmQ_new1, &PmQ_new2);

            // TODO: replace with xdblmul
            (P_new, Q_new, PmQ_new)
        }

        pub fn generate_random_fp(curve: &Curve, order: Integer, f: Integer) -> Point {
            // p + 1 = order * f
            let mut rng = ChaCha20Rng::from_entropy();
            let mut X: Fp;
            let mut P: Point;
            let f_bytes = big_to_bytes(f);

            loop {
                X = Fp::rand(&mut rng);
                let FqX = Fq::new(&X, &Fp::ZERO);
                let Pxz = PointX::new_xz(&FqX, &Fq::ONE);
                let (Px, _) = curve.complete_pointX(&Pxz); // TODO
                P = curve.mul(&Px, &f_bytes, f_bytes.len() * 8);

                let bytes = big_to_bytes(order.clone() - 1);
                let Q = curve.mul(&P, &bytes, bytes.len() * 8);

                if Q.isinfinity() == 0x00000000 {
                    return P
                }
            }
        }

        pub fn generate_random_fq(curve: &Curve, order: Integer, f: Integer) -> Point {
            // p + 1 = order * f
            let mut rng = ChaCha20Rng::from_entropy();
            let mut X: Fq;
            let mut P: Point;
            let f_bytes = big_to_bytes(f);

            loop {
                X = Fq::rand(&mut rng);
                let Pxz = PointX::new_xz(&X, &Fq::ONE);
                let (Px, _) = curve.complete_pointX(&Pxz); // TODO
                P = curve.mul(&Px, &f_bytes, f_bytes.len() * 8);

                let bytes = big_to_bytes(order.clone() - 1);
                let Q = curve.mul(&P, &bytes, bytes.len() * 8);

                if Q.isinfinity() == 0x00000000 {
                    return P
                }
            }
        }

        pub fn is_jac_equal(P: &Point, Q: &Point) -> bool {
            let mut t0 = Q.Z.square();
            let mut t2 = P.X * t0;
            let mut t1 = P.Z.square();
            let t3 = Q.X * t1;
            t2 = t2 - t3;
            t0 = t0 * Q.Z;
            t0 = P.Y * t0;
            t1 = t1 * P.Z;
            t1 = Q.Y * t1;
            t0 = t0 - t1;
            
            t0.iszero() == 0xFFFFFFFF && t2.iszero() == 0xFFFFFFFF
        }

        pub fn jac_neg(T: &mut Point, P: &Point) {
            T.X = P.X.clone();
            T.Y = P.Y.clone().neg();
            T.Z = P.Z.clone();
        }

        pub fn jac_dbl(curve: &Curve, Q: &mut Point, P: &Point) {
            if P.X.iszero() == 0xFFFFFFFF && P.Z.iszero() == 0xFFFFFFFF {
                *Q = Point::INFINITY;
                return;
            }

            let mut t0 = P.X.square();
            let mut t1 = t0 + t0;
            t0 = t0 + t1; // 3x^2
            t1 = P.Z.square(); // z^2
            let mut t2 = P.X * curve.A;
            t2 = t2 + t2; // 2Ax
            t2 = t1 + t2; // 2Ax+z^2
            t2 = t1 * t2; // z^2(2Ax+z^2)
            t2 = t0 + t2; // 3x^2+z^2(2Ax+z^2)
            Q.Z = P.Y * P.Z;
            Q.Z = Q.Z + Q.Z;
            t0 = Q.Z.square();
            t0 = t0 * curve.A;
            t1 = P.Y.square();
            t1 = t1 + t1; // 2y^2
            let mut t3 = P.X + P.X;
            t3 = t1 * t3; // 4xy^2
            Q.X = t2.square();
            Q.X = Q.X - t0; // alpha^2-4Ay^2z^2
            Q.X = Q.X - t3;
            Q.X = Q.X - t3; // alpha^2-4Ay^2z^2-8x^2y^2
            Q.Y = t3 - Q.X;
            Q.Y = Q.Y * t2;
            t1 = t1.square();
            Q.Y = Q.Y - t1;
            Q.Y = Q.Y - t1;
        }
        
        pub fn jac_add(curve: &Curve, R: &mut Point, P: &Point, Q: &Point) {
            if is_jac_equal(P, Q) {
                jac_dbl(curve, R, P);
                return;
            }

            let mut T = Point::INFINITY;
            jac_neg(&mut T, P);
            if is_jac_equal(&mut T, Q) {
                *R = Point::INFINITY;
                return;
            }
            
            if P.X.iszero() == 0xFFFFFFFF && P.Z.iszero() == 0xFFFFFFFF {
                *R = Q.clone();
                return;
            } else if Q.X.iszero() == 0xFFFFFFFF && Q.Z.iszero() == 0xFFFFFFFF {
                *R = P.clone();
                return;
            }

            let mut t0 = P.Z.square();
            let mut t1 = t0 * P.Z;
            let mut t2 = Q.Z.square();
            let mut t3 = t2 * Q.Z;
            t1 = t1 * Q.Y;
            t3 = t3 * P.Y;
            t1 = t1 - t3;
            t0 = t0 * Q.X;
            t2 = t2 * P.X;
            let t4 = t0 - t2;
            t0 = t0 + t2;
            let mut t5 = P.Z * Q.Z;
            R.Z = t4 * t5;
            t5 = t5.square();
            t5 = curve.A * t5;
            t0 = t0 + t5;
            let t6 = t4.square();
            t5 = t0 * t6;
            R.X = t1.square();
            R.X = R.X - t5;
            t3 = t3 * t4;
            t3 = t3 * t6;
            t2 = t2 * t6;
            R.Y = t2 - R.X;
            R.Y = R.Y * t1;
            R.Y = R.Y - t3;
        }

        pub fn jac_triple(curve: &Curve, Q: &mut Point, P: &Point) {
            let mut R = Point::INFINITY;

            jac_dbl(curve, &mut R, P);
            jac_add(curve, Q, &mut R, P);
        }

        fn mul(curve: &Curve, P: &Point, k: u32) -> Point {
            let pw = k.big().pow(1);
            let k_bytes = big_to_bytes(pw);
            let R = curve.mul(&P, &k_bytes, k_bytes.len() * 8);

            R
        }

        /// Pohlig-Hellman algorithm for points in subgroup of order 3^f.
        /// Return the discrete logarithm of R with respect to P.
        pub fn dlog_3(curve: &Curve, P: &Point, R: &Point, f: usize) -> Integer {
            let mut P_table = vec![P.clone(); f];
            let mut R_table = vec![R.clone(); f];
            for i in (0..f-1).rev() {
                P_table[i] = mul(curve, &P_table[i+1], 3);
                R_table[i] = mul(curve, &R_table[i+1], 3);
            }

            let mut dlog = Integer::from(0);
            let mut parts: Vec<usize> = vec![];
            let mut part;

            let P1 = &P_table[0];
            let P2 = curve.double(&P1);
        
            for i in 0..f {
                let mut R_adapted = R_table[i];
                for j in 0..i {
                    let p = parts[j].clone();
                    if p == Integer::from(1) {
                        R_adapted = curve.sub(&R_adapted, &P_table[i-j]);
                    } else if p == Integer::from(2) {
                        let PP = curve.double(&P_table[i-j]);
                        R_adapted = curve.sub(&R_adapted, &PP);
                    }
                }

                let p = 3.big().pow(i as u32);
                if R_adapted.equals(&P1) == 0xFFFFFFFF {
                    part = 1;
                    dlog += p;
                } else if R_adapted.equals(&P2) == 0xFFFFFFFF {
                    part = 2;
                    dlog += Integer::from(2) * p;
                } else {
                    part = 0;
                }

                parts.push(part);
            }
        
            dlog
        }

        /// Pohlig-Hellman algorithm for points in subgroup of order 5^f.
        /// Return the discrete logarithm of R with respect to P.
        pub fn dlog_5(curve: &Curve, P: &Point, R: &Point, f: usize) -> Integer {
            let mut P_table = vec![P.clone(); f];
            let mut R_table = vec![R.clone(); f];
            for i in (0..f-1).rev() {
                P_table[i] = mul(curve, &P_table[i+1], 5);
                R_table[i] = mul(curve, &R_table[i+1], 5);
            }

            let mut dlog = Integer::from(0);
            let mut parts: Vec<usize> = vec![];
            let mut part;

            let P1 = &P_table[0];
            let P2 = curve.double(&P1);
            let P3 = curve.add(&P1, &P2);
            let P4 = curve.double(&P2);
        
            for i in 0..f {
                let mut R_adapted = R_table[i];
                for j in 0..i {
                    let p = parts[j].clone();
                    if p == Integer::from(1) {
                        R_adapted = curve.sub(&R_adapted, &P_table[i-j]);
                    } else if p == Integer::from(2) {
                        let PP = curve.double(&P_table[i-j]);
                        R_adapted = curve.sub(&R_adapted, &PP);
                    } else if p == Integer::from(3) {
                        let mut PP = curve.double(&P_table[i-j]);
                        PP = curve.add(&P_table[i-j], &PP);
                        R_adapted = curve.sub(&R_adapted, &PP);
                    } else if p == Integer::from(4) {
                        let mut PP = curve.double(&P_table[i-j]);
                        PP = curve.double(&PP);
                        R_adapted = curve.sub(&R_adapted, &PP);
                    }
                }

                let p = 5.big().pow(i as u32);
                if R_adapted.equals(&P1) == 0xFFFFFFFF {
                    part = 1;
                    dlog += p;
                } else if R_adapted.equals(&P2) == 0xFFFFFFFF {
                    part = 2;
                    dlog += Integer::from(2) * p;
                } else if R_adapted.equals(&P3) == 0xFFFFFFFF {
                    part = 3;
                    dlog += Integer::from(3) * p;
                } else if R_adapted.equals(&P4) == 0xFFFFFFFF {
                    part = 4;
                    dlog += Integer::from(4) * p;
                } else {
                    part = 0;
                }

                parts.push(part);
            }
        
            dlog
        }
    };
} // End of macro: define_ec_helpers

pub(crate) use define_ec_helpers;
