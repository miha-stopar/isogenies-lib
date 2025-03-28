#![allow(non_snake_case)]

// NOTE:
// The majority of this code was written by Thomas Pornin, which is part of a
// larger work in progress implementing another project which is still ongoing.
// It is similar in form to the macros in the cryptographic library for rust:
// crrl. https://github.com/pornin/crrl

// A macro to create the following Elliptic Curve types:
// - Point:            Type of a point on an elliptic curve (X : Y : Z) with coordinates as elements of Fq
// - PointX:           Type of a point on the Kummer line of an elliptic curve (X : Z) with coordinates as elements of Fq
// - Curve:            Type of an elliptic curve in Montgomery form: E : y^2 = x(x^2 + Ax + 1)
// - CouplePoint       Type of a pair of `Point` P = (P1, P2) in E1 x E2
// - EllipticProduct   Type of a pair of `Curve` representing the product of two curves E1 x E2

// Macro expections:
//    Fq    type of field element
macro_rules! define_ec_core {
    () => {
        // use core::ops::Add;
        use core::ops::Neg;
        use rand_core::{CryptoRng, RngCore};
        use std::fmt;
        use crate::ec::mp::{mp_sub, select_ct, select_ct_arr, swap_ct, mp_shiftr};

        /// Curve point.
        /// Points do not know which curve they are on! The caller must ensure
        /// that only proper curve points are used on a given curve; the Rust
        /// type system does not enforce it.
        #[derive(Clone, Copy, Debug)]
        pub struct Point {
            pub(crate) X: Fq,
            pub(crate) Y: Fq,
            pub(crate) Z: Fq,
        }

        impl Point {
            /// The point-at-infinity (neutral element of the group law).
            pub const INFINITY: Self = Self {
                X: Fq::ZERO,
                Y: Fq::ONE,
                Z: Fq::ZERO,
            };

            /// Create a new point: WARNING no check is made
            pub fn new_xy(X: &Fq, Y: &Fq) -> Self {
                Self {
                    X: *X,
                    Y: *Y,
                    Z: Fq::ONE,
                }
            }

            pub fn new_xyz(X: &Fq, Y: &Fq, Z: &Fq) -> Self {
                Self {
                    X: *X,
                    Y: *Y,
                    Z: *Z,
                }
            }

            /// 0xFFFFFFFF for the point-at-infinity; 0x00000000 otherwise.
            pub fn isinfinity(self) -> u32 {
                self.Z.iszero()
            }

            /// Returns None for the point-at-infinity; otherwise, the (x,y)
            /// affine coordinates.
            pub fn to_xy_vartime(self) -> Option<(Fq, Fq)> {
                if self.Z.iszero() != 0 {
                    return None;
                }
                let t = self.Z.invert();
                Some((self.X * &t, self.Y * &t))
            }

            /// Get the (x,y) affine coordinates. For the point-at-infinity,
            /// this returns (0,0).
            pub fn to_xy(self) -> (Fq, Fq) {
                let t = self.Z.invert();
                (self.X * &t, self.Y * &t)
            }

            /// Returns the X and Z coordinates of the projective point
            pub fn to_xz(self) -> (Fq, Fq) {
                (self.X, self.Z)
            }

            /// Copy rhs into self if ctl == 0xFFFFFFFF.
            /// Do nothing is ctl == 0x00000000.
            /// ctl MUST be either 0xFFFFFFFF or 0x00000000.
            pub fn set_cond(&mut self, rhs: &Self, ctl: u32) {
                self.X.set_cond(&rhs.X, ctl);
                self.Y.set_cond(&rhs.Y, ctl);
                self.Z.set_cond(&rhs.Z, ctl);
            }

            /// Negate this point.
            pub fn set_neg(&mut self) {
                self.Y.set_neg();
            }

            /// Negate this point if ctl == 0xFFFFFFFF.
            /// Do nothing is ctl == 0x00000000.
            /// ctl MUST be either 0xFFFFFFFF or 0x00000000.
            pub fn set_condneg(&mut self, ctl: u32) {
                self.Y.set_condneg(ctl);
            }

            /// Return 0xFFFFFFFF if self and rhs represent the same point.
            /// Otherwise, return 0x00000000.
            pub fn equals(self, rhs: &Self) -> u32 {
                // P1 == P2 if and only if:
                //    P1 == inf AND P2 == inf
                //  OR:
                //    P1 != inf AND P2 != inf AND X1*Z2 = X2*Z1 AND Y1*Z2 = Y2*Z1
                let lz = self.Z.iszero();
                let rz = rhs.Z.iszero();
                let vx = (&self.X * &rhs.Z).equals(&(&rhs.X * &self.Z));
                let vy = (&self.Y * &rhs.Z).equals(&(&rhs.Y * &self.Z));
                (lz & rz) | (!lz & !rz & vx & vy)
            }
        }

        impl fmt::Display for Point {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(
                    f,
                    "Elliptic Curve Point: ({} : {} : {})",
                    self.X, self.Y, self.Z
                )
            }
        }

        impl Neg for Point {
            type Output = Point;

            #[inline(always)]
            fn neg(self) -> Point {
                let mut r = self;
                r.set_neg();
                r
            }
        }

        impl Neg for &Point {
            type Output = Point;

            #[inline(always)]
            fn neg(self) -> Point {
                let mut r = *self;
                r.set_neg();
                r
            }
        }

        /// Special X-only representation of a point (or a pair of points,
        /// since two Y coordinates may match a given X).
        #[derive(Clone, Copy, Debug)]
        pub struct PointX {
            pub(crate) X: Fq,
            pub(crate) Z: Fq,
        }

        impl PointX {
            /// The point-at-infinity (neutral element of the group law).
            #[allow(dead_code)]
            const INFINITY: Self = Self {
                X: Fq::ZERO,
                Z: Fq::ZERO,
            };

            #[allow(dead_code)]
            fn from_point(P: &Point) -> Self {
                Self { X: P.X, Z: P.Z }
            }

            #[allow(dead_code)]
            pub fn isinfinity(self) -> u32 {
                self.Z.iszero()
            }

            pub fn new_xz(X: &Fq, Z: &Fq) -> Self {
                Self { X: *X, Z: *Z }
            }

            #[allow(dead_code)]
            fn equals(self, rhs: &PointX) -> u32 {
                let inf1 = self.isinfinity();
                let inf2 = rhs.isinfinity();
                let e = (&self.X * &rhs.Z).equals(&(&rhs.X * &self.Z));
                (inf1 & inf2) | (!inf1 & !inf2 & e)
            }

            /// Returns the X and Z coordinates of the projective point
            pub fn to_xz(self) -> (Fq, Fq) {
                (self.X, self.Z)
            }
        }

        impl fmt::Display for PointX {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "Elliptic Curve Point: ({} : {})", self.X, self.Z)
            }
        }

        /// Curve y^2 = x^3 + A*x^2 + x, for a given constant A
        /// (special case of a Montgomery curve).
        #[derive(Clone, Copy, Debug)]
        pub struct Curve {
            pub A: Fq, // curve parameter
            pub A24: Fq,   // (A+2)/4
        }

        impl Curve {
            /// Create a new curve instance, with the provided constant.
            pub fn new(A: &Fq) -> Self {
                // We check that the curve is not singular, i.e. A^2 != 4.
                // FIXME: do we want to keep that check?
                assert!(A.equals(&Fq::TWO) == 0);
                assert!((A + &Fq::TWO).iszero() == 0);

                Self {
                    A: *A,
                    A24: (A + &Fq::TWO).half().half(),
                }
            }

            pub fn new_fromA24(A24: &Fq) -> Self {
                Self {
                    A: Fq::ONE, // TODO, but not used
                    A24: *A24,
                }
            }

            /// Set the point to the provided affine coordinate.
            /// IMPORTANT: this function does NOT check that the point is
            /// really part of the curve.
            pub fn set_xy_nocheck(self, P: &mut Point, x: &Fq, y: &Fq) {
                P.X = *x;
                P.Y = *y;
                P.Z = Fq::ONE;
            }

            /// Set the point P to the provided projective coordinates (use
            /// Z = 0 for the point-at-infinity).
            /// IMPORTANT: this function does NOT check that the point is
            /// really part of the curve.
            pub fn set_xyz_nocheck(self, P: &mut Point, x: &Fq, y: &Fq, z: &Fq) {
                P.X = *x;
                P.Y = *y;
                P.Z = *z;
            }

            /// Set the point P to the provided affine coordinates. This
            /// function checks that the point is on the curve. If the
            /// point is on the curve, then this function returns 0xFFFFFFFF;
            /// otherwise, P is set to the infinity point and this function
            /// returns 0x00000000.
            pub fn set_xy(self, P: &mut Point, x: &Fq, y: &Fq) -> u32 {
                let mut t = x + &self.A;
                t *= &x.square();
                t += x;
                let r = y.square().equals(&t);
                self.set_xy_nocheck(P, x, y);
                P.Z.set_cond(&Fq::ZERO, !r);
                r
            }

            /// Set the point P to the provided projective coordinates (use
            /// Z = 0 for the point-at-infinity). This function checks that
            /// the point is on the curve. If the point is on the curve, then
            /// this function returns 0xFFFFFFFF; otherwise, P is set to the
            /// infinity point and this function returns 0x00000000.
            pub fn set_xyz(self, P: &mut Point, x: &Fq, y: &Fq, z: &Fq) -> u32 {
                let mut t1 = (x + z).square();
                t1 += (&self.A - &Fq::TWO) * &(x * z);
                t1 *= x;
                let t2 = &y.square() * z;
                let r = t1.equals(&t2) | z.iszero();
                self.set_xyz_nocheck(P, x, y, z);
                P.Z.set_cond(&Fq::ZERO, !r);
                r
            }

            /// Get a Point instance set to the provided affine coordinates.
            /// IMPORTANT: this function does NOT check whether the point is on
            /// curve or not.
            pub fn point_xy_nocheck(self, x: &Fq, y: &Fq) -> Point {
                Point {
                    X: *x,
                    Y: *y,
                    Z: Fq::ONE,
                }
            }

            /// Get a Point instance set to the provided projective coordinates
            /// (use Z = 0 for the point-at-infinity).
            /// IMPORTANT: this function does NOT check whether the point is on
            /// curve or not.
            pub fn point_xyz_nocheck(self, x: &Fq, y: &Fq, z: &Fq) -> Point {
                Point {
                    X: *x,
                    Y: *y,
                    Z: *z,
                }
            }

            /// Get a Point instance from the provided affine coordinates.
            /// This returns None if the coordinates do not designate a point
            /// on the curve.
            /// CT: whether the point is on the curve or not may leak
            pub fn point_xy_vartime(self, x: &Fq, y: &Fq) -> Option<Point> {
                let mut P = Point::INFINITY;
                if self.set_xy(&mut P, x, y) == 0 {
                    None
                } else {
                    Some(P)
                }
            }

            /// Get a Point instance from the provided projective coordinates
            /// (use Z = 0 for the point-at-infinity). This returns None if the
            /// coordinates do not designate a point on the curve.
            /// CT: whether the point is on the curve or not may leak
            pub fn point_xyz_vartime(self, x: &Fq, y: &Fq, z: &Fq) -> Option<Point> {
                let mut P = Point::INFINITY;
                if self.set_xyz(&mut P, x, y, z) == 0 {
                    None
                } else {
                    Some(P)
                }
            }

            /// P3 <- P1 + P2
            pub fn add_into(self, P3: &mut Point, P1: &Point, P2: &Point) {
                // Complete routine, to handle all edge cases:
                //   if Z1 == 0:            # P1 == inf
                //       return P2
                //   if Z2 == 0:            # P2 == inf
                //       return P1
                //   L <- Y2*Z1 - Y1*Z2
                //   T <- X2*Z1 - X1*Z2
                //   if T == 0:             # x1 == x2
                //       if L == 0:         # ... and y1 == y2: doubling case
                //           L <- 3*X1^2 + 2*A*X1*Z1 + Z1^2
                //           T <- 2*Y1*Z1
                //       else:              # ... but y1 != y2, thus P2 = -P1
                //           return inf
                //   U <- Z1*Z2*L^2 - (X1*Z2 + X2*Z1 + A*Z1*Z2)*T^2
                //   X3 <- U*T
                //   Y3 <- L*(X1*Z2*T^2 - U) - Y1*Z2*T^3
                //   Z3 <- Z1*Z2*T^3
                //
                // Constant-time processing:
                //   Cases P1 == inf and P2 == inf are handled at the end.
                //   (L,T) are always computed for both normal and doubling cases.
                //   If P1 == -P2 then we can let T == 0 and L != 0, this will
                //   properly lead to Z3 == 0.
                //
                // Formulas from https://eprint.iacr.org/2015/1060 are faster
                // but do not cover the case when P1 - P2 is a point of order 2,
                // which can happen in all generality.
                //
                // : current formulas have cost 16M+5S; this can probably
                // be improved. Main issues to tackle:
                //   - Multiplications by A are expensive (since A can be any value)
                //   - There are three points of order 2; this makes finding
                //     complete formulas challenging.

                // T = X2*Z1 - X1*Z2
                // L = Y2*Z1 - Y1*Z2
                let x1z2 = &P1.X * &P2.Z;
                let x2z1 = &P2.X * &P1.Z;
                let mut T = &x2z1 - &x1z2;
                let y1z2 = &P1.Y * &P2.Z;
                let y2z1 = &P2.Y * &P1.Z;
                let mut L = &y2z1 - &y1z2;

                // Alternate (T,L) for doubling:
                //   Td = 2*Y1*Z1
                //   Ld = 3*X1^2 + 2*A*X1*Z1 + Z1^2
                let dbl = T.iszero() & L.iszero();
                let Td = (&P1.Y * &P1.Z).mul2();
                let x1x1 = P1.X.square();
                let z1z1 = P1.Z.square();
                let dx1z1 = &(&P1.X + &P1.Z).square() - &x1x1 - &z1z1;
                let Ld = &x1x1.mul3() + &z1z1 + &self.A * &dx1z1;
                T.set_cond(&Td, dbl);
                L.set_cond(&Ld, dbl);

                // U = L^2*Z1*Z2 - (X1*Z2 + X2*Z1 + A*Z1*Z2)*T^2
                let T2 = T.square();
                let T3 = &T * &T2;
                let z1z2 = &P1.Z * &P2.Z;
                let U = &(&L.square() * &z1z2) - &(&(&x1z2 + &x2z1 + &(&self.A * &z1z2)) * &T2);

                // X3 = U*T
                // Y3 = L*(X1*Z2*T^2 - U) - Y1*Z2*T^3
                // Z3 = Z1*Z2*T^3
                P3.X = &U * &T;
                P3.Y = &(&L * &(&(&x1z2 * &T2) - &U)) - &(&y1z2 * &T3);
                P3.Z = &z1z2 * &T3;

                // Corrective action in case one of the inputs was the
                // point-at-infinity.
                let inf1 = P1.Z.iszero();
                let inf2 = P2.Z.iszero();
                P3.set_cond(&P2, inf1);
                P3.set_cond(&P1, inf2);
            }

            /// P1 <- P1 + P2
            pub fn addto(self, P1: &mut Point, P2: &Point) {
                let mut P3 = Point::INFINITY;
                self.add_into(&mut P3, P1, P2);
                *P1 = P3;
            }

            /// Return P1 + P2 as a new point
            pub fn add(self, P1: &Point, P2: &Point) -> Point {
                let mut P3 = Point::INFINITY;
                self.add_into(&mut P3, P1, P2);
                P3
            }

            /// P3 <- P1 - P2
            pub fn sub_into(self, P3: &mut Point, P1: &Point, P2: &Point) {
                let mut nP2 = *P2;
                nP2.set_neg();
                self.add_into(P3, P1, &nP2);
            }

            /// P1 <- P1 - P2
            pub fn subfrom(self, P1: &mut Point, P2: &Point) {
                let mut nP2 = *P2;
                nP2.set_neg();
                self.addto(P1, &nP2);
            }

            /// Return P1 - P2 as a new point
            pub fn sub(self, P1: &Point, P2: &Point) -> Point {
                let mut nP2 = *P2;
                nP2.set_neg();
                self.add(P1, &nP2)
            }

            #[inline(always)]
            pub fn double_from_coords(self, X: &Fq, Y: &Fq, Z: &Fq) -> (Fq, Fq, Fq) {
                // Doubling formulas in cost 6M+6S
                // These formulas are complete.
                // Formulas from https://eprint.iacr.org/2015/1060 would be
                // more expensive, because multiplications by A are not cheap
                // in the general case.
                //
                // V <- X^2 - Z^2
                // M <- X^2 + Z^2
                // X' <- 2*Y*Z*V^2
                // Y' <- V*(M*(M + 2*A*X*Z) + 4*X^2*Z^2)
                // Z' <- 8*(Y*Z)^3
                let xx = X.square();
                let zz = Z.square();
                let dxz = &(X + Z).square() - &xx - &zz;
                let dyz = (Y * Z).mul2();
                let v = &xx - &zz;
                let m = &xx + &zz;
                let X2 = &dyz * &v.square();
                let Y2 = &v * (&(&m * &(&m + &(&self.A * &dxz))) + &dxz.square());
                let Z2 = &dyz * &dyz.square();

                (X2, Y2, Z2)
            }

            /// P3 <- 2*P1
            pub fn double_into(self, P3: &mut Point, P1: &Point) {
                let (X2, Y2, Z2) = self.double_from_coords(&P1.X, &P1.Y, &P1.Z);
                P3.X = X2;
                P3.Y = Y2;
                P3.Z = Z2;
            }

            // This is essentially reuse of the above, but allowing the point
            // itself to be mutated... Maybe it would be better to redesign the
            // below functions to stop the code duplication...
            pub fn double_self(self, P1: &mut Point) {
                let (X2, Y2, Z2) = self.double_from_coords(&P1.X, &P1.Y, &P1.Z);
                P1.X = X2;
                P1.Y = Y2;
                P1.Z = Z2;
            }

            /// Return 2*P as a new point
            pub fn double(self, P: &Point) -> Point {
                let mut P3 = Point::INFINITY;
                self.double_into(&mut P3, P);
                P3
            }

            /// Return [2^n]*P as a new point
            pub fn double_iter(self, P: &Point, n: usize) -> Point {
                let mut P3 = *P;
                for _ in 0..n {
                    self.double_self(&mut P3);
                }
                P3
            }

            /// Compute the x-only double of a given point and return
            /// the X-coords
            #[inline(always)]
            pub fn x_dbl_coords(self, X: &Fq, Z: &Fq) -> (Fq, Fq) {
                let mut V1 = (&*X + &*Z).square();
                let V2 = (&*X - &*Z).square();
                let X_new = &V1 * &V2;
                V1 -= &V2;
                let mut Z_new = V1;
                Z_new *= &self.A24;
                Z_new += &V2;
                Z_new *= &V1;

                (X_new, Z_new)
            }

            #[inline(always)]
            pub fn xadd(self, Xp: &Fq, Zp: &Fq, X0: &Fq, Z0: &Fq, X1: &mut Fq, Z1: &mut Fq) {
                let V1 = &(X0 - Z0) * &(&*X1 + &*Z1);
                let V2 = &(X0 + Z0) * &(&*X1 - &*Z1);
                *X1 = Zp * &(&V1 + &V2).square();
                *Z1 = Xp * &(&V1 - &V2).square();
            }

            #[inline(always)]
            pub fn xdbl(self, X: &mut Fq, Z: &mut Fq) {
                let mut V1 = (&*X + &*Z).square();
                let V2 = (&*X - &*Z).square();
                *X = &V1 * &V2;
                V1 -= &V2;
                *Z = V1;
                *Z *= &self.A24;
                *Z += &V2;
                *Z *= &V1;
            }

            /// Return [2^n]*P as a new point
            pub fn x_double_iter(self, P: &PointX, n: usize) -> PointX {
                let mut P3 = *P;
                for _ in 0..n {
                    self.xdbl(&mut P3.X, &mut P3.Z);
                }
                P3
            }

            #[inline(always)]
            pub fn xtriple(self, X: &mut Fq, Z: &mut Fq) {
                let X0 = *X;
                let Z0 = *Z;
                self.xdbl(X, Z);
                self.xadd(&X0, &Z0, &X0, &Z0, X, Z);
            }

            /// P3 <- n*P
            /// Integer n is encoded as unsigned little-endian, with length
            /// nbitlen bits. Bits beyond that length are ignored.
            pub fn mul_into(self, P3: &mut Point, P: &Point, n: &[u8], nbitlen: usize) {
                // Montgomery ladder: see https://eprint.iacr.org/2017/212

                #[inline(always)]
                fn xadd_aff(_curve: &Curve, Xp: &Fq, X0: &Fq, Z0: &Fq, X1: &mut Fq, Z1: &mut Fq) {
                    let V1 = &(X0 - Z0) * &(&*X1 + &*Z1);
                    let V2 = &(X0 + Z0) * &(&*X1 - &*Z1);
                    *X1 = (&V1 + &V2).square();
                    *Z1 = Xp * &(&V1 - &V2).square();
                }

                // We will need the complete 2*P at the end, to handle some
                // special cases of the formulas.
                let dP = self.double(P);
                let mut X0 = Fq::ONE;
                let mut Z0 = Fq::ZERO;
                let mut X1 = P.X;
                let mut Z1 = P.Z;
                let mut cc = 0u32;
                if nbitlen > 21 {
                    // If n is large enough then it is worthwhile to
                    // normalize the source point to affine.
                    // We do not care if P = inf, since that is handled at
                    // the end in the corrective steps.
                    let Xp = &P.X / &P.Z;
                    for i in (0..nbitlen).rev() {
                        let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                        Fq::condswap(&mut X0, &mut X1, ctl ^ cc);
                        Fq::condswap(&mut Z0, &mut Z1, ctl ^ cc);
                        xadd_aff(&self, &Xp, &X0, &Z0, &mut X1, &mut Z1);
                        self.xdbl(&mut X0, &mut Z0);
                        cc = ctl;
                    }
                } else {
                    for i in (0..nbitlen).rev() {
                        let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                        Fq::condswap(&mut X0, &mut X1, ctl ^ cc);
                        Fq::condswap(&mut Z0, &mut Z1, ctl ^ cc);
                        self.xadd(&P.X, &P.Z, &X0, &Z0, &mut X1, &mut Z1);
                        self.xdbl(&mut X0, &mut Z0);
                        cc = ctl;
                    }
                }
                Fq::condswap(&mut X0, &mut X1, cc);
                Fq::condswap(&mut Z0, &mut Z1, cc);

                // Special cases:
                //  - ladder fails if P = (0,0) (a point of order 2)
                //  - y is not reconstructed correctly if P has order 2,
                //    or if (n+1)*P = P, -P or infinity.
                let z0z = Z0.iszero();
                let z1z = Z1.iszero();
                let x1ex = (&X1 * &P.Z).equals(&(&P.X * &Z1));

                // (X0/Z0) is the X coordinate of P0 = n*P
                // (X1/Z1) is the X coordinate of P1 = (n + 1)*P
                // We recompute the Y coordinate of n*P (formulas from
                // Okeya and Sakurai).
                let xxzz = &(&P.X * &X0) + &(&P.Z * &Z0);
                let xpz0 = &P.X * &Z0;
                let x0zp = &X0 * &P.Z;
                let zz = &P.Z * &Z0;
                let zzdA = &self.A.mul2() * &zz;
                let u = &(&xxzz * &(&xpz0 + &x0zp + &zzdA)) - &(&zzdA * &zz);
                let v = &P.Y.mul2() * &zz * &Z1;
                P3.X = &X0 * &v;
                P3.Y = &(&u * &Z1) - &(&(&xpz0 - &x0zp).square() * &X1);
                P3.Z = &Z0 * &v;

                // Fix result for the special cases.
                //  P = inf                          -> inf
                //  P != inf, 2*P = inf              -> inf or P (depending on n_0)
                //  2*P != inf, P0 = inf             -> inf
                //  2*P != inf, P0 != inf, P1 = inf  -> -P
                //  2*P != inf, P0 != inf, P1 = -P   -> -2*P
                let order1 = P.Z.iszero();
                let order2 = !order1 & P.Y.iszero();
                let z0inf = !order1 & !order2 & z0z;
                let z1inf = !order1 & !order2 & !z0z & z1z;
                let p1mp = !order1 & !order2 & !z0z & !z1z & x1ex;

                let n_odd = ((n[0] as u32) & 1).wrapping_neg();
                P3.Z.set_cond(&Fq::ZERO, order1 | (order2 & !n_odd) | z0inf);
                P3.set_cond(&P, z1inf | (order2 & n_odd));
                P3.set_cond(&dP, p1mp);
                P3.set_condneg(z1inf | p1mp);
            }

            /// Return n*P as a new point.
            /// Integer n is encoded as unsigned little-endian, with length
            /// nbitlen bits. Bits beyond that length are ignored.
            pub fn mul(self, P: &Point, n: &[u8], nbitlen: usize) -> Point {
                let mut P3 = Point::INFINITY;
                self.mul_into(&mut P3, P, n, nbitlen);
                P3
            }

            /// P3 <- n*P
            /// CT: constant-time for the points, but the integer n may leak.
            pub fn mul_small_into(self, P3: &mut Point, P: &Point, n: u64) {
                match n {
                    0 => {
                        *P3 = Point::INFINITY;
                    }
                    1 => {
                        *P3 = *P;
                    }
                    2 => {
                        self.double_into(P3, P);
                    }
                    3 => {
                        self.add_into(P3, &self.double(P), P);
                    }
                    _ => {
                        self.mul_into(P3, P, &n.to_le_bytes(), (64 - n.leading_zeros()) as usize);
                    }
                }
            }

            /// Return n*P as a new point.
            /// CT: constant-time for the points, but the integer n may leak.
            pub fn mul_small(self, P: &Point, n: u64) -> Point {
                let mut P3 = Point::INFINITY;
                self.mul_small_into(&mut P3, P, n);
                P3
            }

            /// Set P to a random curve point.
            pub fn set_rand_point<T: CryptoRng + RngCore>(self, rng: &mut T, P: &mut Point) {
                // This function cannot actually return the point-at-infinity;
                // this is not a problem as long as the curve order is larger
                // than 2^128.
                P.Z = Fq::ONE;
                loop {
                    P.X.set_rand(rng);
                    P.Y = &(&(&(&P.X + &self.A) * &P.X) * &P.X) + &P.X;
                    if P.Y.legendre() >= 0 {
                        P.Y.set_sqrt();

                        // Randomly chooses the square root to use.
                        let mut tmp = [0u8; 1];
                        rng.fill_bytes(&mut tmp);
                        let ctl = 0u32.wrapping_sub((tmp[0] as u32) & 1);
                        P.Y.set_condneg(ctl);
                        return;
                    }
                }
            }

            /// Return a new random curve point.
            pub fn rand_point<T: CryptoRng + RngCore>(self, rng: &mut T) -> Point {
                let mut P = Point::INFINITY;
                self.set_rand_point(rng, &mut P);
                P
            }

            /// Complete an X-only point into a full point;
            /// (an error is returned if there is no matching Y coordinate).
            /// On error, P3 is set to the point-at-infinity.
            fn complete_pointX_into(self, P3: &mut Point, P: &PointX) -> u32 {
                let XZ = &P.X * &P.Z;
                let V = &(&P.X + &P.Z).square() + &(&(&self.A - &Fq::TWO) * &XZ);
                P3.X = XZ;
                P3.Y = &V * &XZ;
                let ok = P3.Y.set_sqrt();
                P3.Z = P.Z.square();

                // Set to inf on error.
                P3.Z.set_cond(&Fq::ZERO, !ok);

                ok
            }

            /// Complete an X-only point into a full point;
            /// On error, the output point is set to the point-at-infinity.
            pub fn complete_pointX(self, P: &PointX) -> (Point, u32) {
                let mut P3 = Point::INFINITY;
                let ok = self.complete_pointX_into(&mut P3, P);
                (P3, ok)
            }

            fn select_point(self, P: &Point, Q: &Point, ctl: u64) -> Point {
                let mut S = Point::INFINITY;
                let ctl32 = ctl as u32;
                
                S.X = Fq::select(&P.X, &Q.X, ctl32);
                S.Y = Fq::select(&P.Y, &Q.Y, ctl32);
                S.Z = Fq::select(&P.Z, &Q.Z, ctl32);

                S
            }

            fn swap_points(self, P: &mut Point, Q: &mut Point, ctl: u64) {
                let ctl32 = ctl as u32;
                Fq::condswap(&mut P.X, &mut Q.X, ctl32);
                Fq::condswap(&mut P.Y, &mut Q.Y, ctl32);
                Fq::condswap(&mut P.Z, &mut Q.Z, ctl32);
            }

            pub fn xdblmul_bounded(
                self,
                P: &Point,
                k: &[u64],
                Q: &Point,
                l: &[u64],
                PQ: &Point,
                f: usize, // TODO
            ) -> Point {
                // TODO: we use Point, but we could use PointX (and remove some unnecessary copying
                // of Y in select_point, swap_points...)
                let mut sigma = [0u64; 2];
                let mut pre_sigma = 0u64;
                let mut evens = 0;
                let mut mevens = 0;
                let mut bitk0 = k[0] & 1;
                let mut bitl0 = l[0] & 1;
                let maskk = 0u64.wrapping_sub(bitk0); // Parity masks
                let maskl = 0u64.wrapping_sub(bitl0);
            
                sigma[0] = bitk0 ^ 1; // 1 if k is even, 0 if k is odd.
                sigma[1] = bitl0 ^ 1;

                evens = sigma[0] + sigma[1];
                mevens = 0u64.wrapping_sub(evens & 1);

                // TODO:
                const NWORDS_ORDER: usize = 4;
            
                sigma[0] &= mevens;
                sigma[1] = (sigma[1] & mevens) | (1 & !mevens);
            
                // Convert even scalars to odd
                let mut one = [0u64; NWORDS_ORDER];
                one[0] = 1;
            
                let mut k_t = k.to_vec();
                let mut l_t = l.to_vec();

                mp_sub(&mut k_t, k, &one, NWORDS_ORDER);
                mp_sub(&mut l_t, l, &one, NWORDS_ORDER);

                let k_t_c = k_t.clone();
                let l_t_c = l_t.clone();
 
                select_ct_arr(&mut k_t, &k_t_c, k, maskk, NWORDS_ORDER);
                select_ct_arr(&mut l_t, &l_t_c, l, maskl, NWORDS_ORDER);

                const BITS: usize = 256;

                let mut r = [0u64; 2 * BITS];
                
                // Scalar recoding
                for i in 0..BITS {
                    let maskk = 0u64.wrapping_sub(sigma[0] ^ pre_sigma);
                    swap_ct(&mut k_t, &mut l_t, maskk, NWORDS_ORDER); 

                    let bs1_ip1 = if i == BITS - 1 { 0 } else { mp_shiftr(&mut k_t, 1, NWORDS_ORDER) };
                    let bs2_ip1 = if i == BITS - 1 { 0 } else { mp_shiftr(&mut l_t, 1, NWORDS_ORDER) };
            
                    let bs1_i = k_t[0] & 1;
                    let bs2_i = l_t[0] & 1;

                    r[2 * i] = bs1_i ^ bs1_ip1;
                    r[2 * i + 1] = bs2_i ^ bs2_ip1;
            
                    pre_sigma = sigma[0];
                    let maskk = 0u64.wrapping_sub(r[2 * i + 1]);

                    let temp = select_ct(sigma[0], sigma[1], maskk);
                    sigma[1] = select_ct(sigma[1], sigma[0], maskk);
                    sigma[0] = temp;
                }
                
                // Point initialization
                let maskk = 0u64.wrapping_sub(sigma[0]);

                // Do not use Point::INFINITY here, we need X = 1 for xadd operation below (otherwise
                // the addition of a point with INFINITY gives INFINITY).
                let mut R0 = Point::new_xyz(&Fq::ONE, &Fq::ONE, &Fq::ZERO);

                let mut R1 = self.select_point(P, Q, maskk);
                let mut R2 = self.select_point(Q, P, maskk);

                let mut diff1a = R1;
                let mut diff1b = R2;

                self.xadd(&PQ.X, &PQ.Z, &R1.X, &R1.Z, &mut R2.X, &mut R2.Z);

                let mut diff2a = R2;
                let mut diff2b = PQ.clone(); // TODO: why clone needed?

                let mut T0 = R0;
                let mut T1 = R1;
                let mut T2 = R2;

                // Main loop
                for i in (0..BITS).rev() {
                    // let apply = i <= f + 2 + (BITS - TORSION_PLUS_EVEN_POWER);
                    // TODO
                    let apply = i <= 2 + BITS;
            
                    let h = r[2 * i] + r[2 * i + 1];
                    let maskk = 0u64.wrapping_sub(h & 1);

                    if apply {
                        T0 = self.select_point(&R0, &R1, maskk);
                    }

                    let maskk = 0u64.wrapping_sub(h >> 1);
            
                    if apply {
                        T0 = self.select_point(&T0, &R2, maskk);
                        self.xdbl(&mut T0.X, &mut T0.Z);
                    }
            
                    let maskk = 0u64.wrapping_sub(r[2 * i + 1]);
                    if apply {
                        T1 = self.select_point(&R0, &R1, maskk);
                        T2 = self.select_point(&R1, &R2, maskk);
                    }
            
                    self.swap_points(&mut diff1a, &mut diff1b, maskk);

                    if apply {
                        self.xadd(&diff1a.X, &diff1a.Z, &T2.X, &T2.Z, &mut T1.X, &mut T1.Z);

                        let R2_old = R2.clone();
                        self.xadd(&diff2a.X, &diff2a.Z, &R0.X, &R0.Z, &mut R2.X, &mut R2.Z);
                        T2 = R2.clone(); // TODO
                        R2 = R2_old; // TODO
                    }
            
                    let maskk = 0u64.wrapping_sub(h & 1);
                    self.swap_points(&mut diff2a, &mut diff2b, maskk);
            
                    R0 = T0;
                    R1 = T1;
                    R2 = T2;
                }

                let mut S = self.select_point(&R0, &R1, mevens);
            
                let maskk = 0u64.wrapping_sub(bitk0 & bitl0);
                S = self.select_point(&S, &R2, maskk);

                let Sxz = PointX::new_xz(&S.X, &S.Z);
                let ok = self.complete_pointX_into(&mut S, &Sxz); // TODO: check ok

                S
            }

            // S <- 2*S
            fn line_double(self, pc: &mut PairingContext, S: &mut Point) {
                // Line slope is L/T. Special cases:
                //   S = inf                 -> current value is not changed
                //   S != inf and 2*S = inf  -> line is vertical
                // When 2*S = inf, we set the slope to (1,0). To make the
                // formulas also correct when S = inf, we enforce X != 0 in
                // that case.
                let sinf = S.isinfinity();
                S.X.set_cond(&Fq::ONE, sinf);
                let XX = S.X.square();
                let ZZ = S.Z.square();
                let dXZ = &(&S.X + &S.Z).square() - &XX - &ZZ;
                let mut L = &XX.mul3() + &(&self.A * &dXZ) + &ZZ;
                let YY = S.Y.square();
                let T = &(&S.Y + &S.Z).square() - &YY - &ZZ;
                let tz = T.iszero();
                L.set_cond(&Fq::ONE, tz);

                // Apply the line on Q+R2 and R2.
                let n1 = &(&L * &(&(&pc.xq * &S.Z) - &S.X)) - &(&T * &(&(&pc.yq * &S.Z) - &S.Y));
                let d3 = &(&L * &(&(&pc.xr * &S.Z) - &S.X)) - &(&T * &(&(&pc.yr * &S.Z) - &S.Y));

                // S' = 2*S
                // If S' = inf, we enforce X = 1.
                let V = (&XX - &ZZ).square();
                let dYY = YY.mul2();
                let TT = T.square();
                S.X = &V * &T;
                S.Y = &(&L * &(&(&dYY * &dXZ) - &V)) - &(&dYY * &TT);
                S.Z = &T * &TT;
                S.X.set_cond(&Fq::ONE, tz);

                // Apply the vertical line going through S'. If S' = inf,
                // this should be a no-operation; since we enforced X = 1
                // in that case, the formulas still work.
                let n2 = &(&pc.xr * &S.Z) - &S.X;
                let d4 = &(&pc.xq * &S.Z) - &S.X;

                // Update vn/vd.
                pc.vn.set_square();
                pc.vn *= &(&n1 * &n2);
                pc.vd.set_square();
                pc.vd *= &(&d3 * &d4);
            }

            fn weil_pairing_2exp(self, e: usize, P: &Point, Q: &Point) -> (Fq, u32) {
                // For order n = 2^e, we use R1 = P-Q and R2 = P+Q. Thus,
                // P+R1 = 2*P-Q, and Q+R2 = P+2*Q.
                // All Miller steps in f_P(A_Q) after the first iteration use
                // two line functions:
                //   - line from (2^k)*P to -(2^(k+1))*P, for k >= 1
                //   - line from -(2^(k+1))*P to +(2^(k+1))*P, for k >= 1
                // If one of these lines goes through R2 or Q+R2 (thus
                // triggering a zero), then this means that (2^k+1)*P or
                // (2^k-1)*P (for some integer k >= 1) is equal to either
                // +Q or -Q. P has order n = 2^e; 2^k+1 and 2^k-1 are odd,
                // thus invertible modulo n, which means that P and Q are
                // colinear and the Weil pairing should ultimately be 1. We
                // thus handle these cases properly by normalizing zeros to
                // a final result of 1. The same reasoning goes for f_Q(A_P).
                //
                // Remaining cases are for the computation of f1_P(A_Q) and
                // f1_Q(A_P), and the first iteration of the Miller loop.
                // We potentially get a zero when some point is on a given
                // line; here is the list of all cases to account for (for
                // f_P(A_Q) and f_Q(A_P)):
                //
                //    P = R2             Q = inf
                //    P = Q+R2           2*Q = inf
                //    R1 = R2            2*Q = inf
                //    R1 = Q+R2          3*Q = inf (impossible)
                //    -(P+R1) = R2       3*P = inf (impossible)
                //    -(P+R1) = Q+R2     Q = -3*P (colinear)
                //    P+R1 = R2          P = 2*Q (colinear)
                //    P+R1 = Q+R2        3*P = 3*Q (colinear)
                //
                //    Q = R1             P = 2*Q (colinear)
                //    Q = P+R1           2*P = 2*Q
                //    R2 = R1            2*Q = inf
                //    R2 = P+R1          P = -2*Q (colinear)
                //    -(Q+R2) = R1       2*P = -3*Q (colinear)
                //    -(Q+R2) = P+R1     Q = 3*P (colinear)
                //    Q+R2 = R1          2*P = 3*Q (colinear)
                //    Q+R2 = P+R1        P = 3*Q (colinear)
                //
                // The left column lists cases that trigger a zero in the
                // accumulator; the right column gives the corresponding
                // conditions on points P and Q, using R1 = P-Q and R2 = P+Q.
                // Those marked "impossible" cannot happen for points of
                // n-torsion (and the algorithm explicitly verifies that
                // both points are n-torsion). "Colinear" means that one of
                // the point is expressed as an integer multiple of the other,
                // leading to a Weil pairing of 1, and the normalization step
                // of zeros to 1 correctly handles these. The only remaining
                // cases that must be handled as corrective steps are the
                // following:
                //
                //    P = inf
                //    Q = inf
                //    2*P = inf
                //    2*Q = inf
                //    2*(P-Q) = inf
                //
                // When P = inf or Q = inf, w(P,Q) = 1.
                //
                // For the order-2 cases: let (U,V) be a basis of the
                // n-torsion group, so that any n-torsion point P can be
                // uniquely written as P = a*U + b*V, for two integers
                // (a,b) taken modulo n. The points of 2-torsion then
                // correspond to (0,0) (the point-at-infinity), (n/2,0),
                // (0,n/2) and (n/2,n/2). Then, if Q = a'*U + b'*V, we
                // have:
                //    w(P,Q) = g^(a*b'-a'*b)
                // for g = w(U,V) (which is a generator of the n-th roots of
                // unity in GF(p^2)). In particular, if 2*P = inf, then
                // both a and b are in {0,n/2}, and we have:
                //
                //    w(P,Q) = w((a/(n/2),b/(n/2)), (n/2)*Q)
                //
                // with (n/2)*Q being itself a point of 2-torsion. This yields
                // either 1 or -1:
                //
                //    If (n/2)*Q != inf and (n/2)*Q != P, then w(P,Q) = -1
                //    Otherwise, w(P,Q) = 1.
                //
                // We compute (n/2)*P and (n/2)*Q as part of Miller's algorithm,
                // so we can handle these cases with negligible overhead.
                //
                // For the final case of 2*(P-Q) = 0, we can remark that:
                //
                //    w(P,Q) = w(P-Q,Q)*w(Q,Q) = w(P-Q,Q)
                //
                // so that we can handle that case similarly to that of points
                // of order 2. The point P-Q is R1, we compute it explicitly;
                // testing whether it has order 2 is then simply checking whether
                // is Y coordinate is zero.

                let R1 = self.sub(P, Q);
                let R2 = self.add(P, Q);
                let PR1 = self.add(P, &R1);
                let QR2 = self.add(Q, &R2);

                // Get the slopes for the lines that double P and Q, respectively.
                let dA = self.A.mul2();
                let Lp = &P.X.square().mul3() + &(&(&(&dA * &P.X) + &P.Z) * &P.Z);
                let Tp = &P.Y.mul2() * &P.Z;
                let Lq = &Q.X.square().mul3() + &(&(&(&dA * &Q.X) + &Q.Z) * &Q.Z);
                let Tq = &Q.Y.mul2() * &Q.Z;

                // f1_P(A_Q)
                let n1 = &(&QR2.X * &PR1.Z) - &(&PR1.X * &QR2.Z);
                let d2 = &(&R2.X * &PR1.Z) - &(&PR1.X * &R2.Z);
                let L = &(&R1.Y * &P.Z) - &(&P.Y * &R1.Z);
                let T = &(&R1.X * &P.Z) - &(&P.X * &R1.Z);
                let n3 = &(&T * &(&(&R2.Y * &P.Z) - &(&P.Y * &R2.Z)))
                    - &(&L * &(&(&R2.X * &P.Z) - &(&P.X * &R2.Z)));
                let d4 = &(&T * &(&(&QR2.Y * &P.Z) - &(&P.Y * &QR2.Z)))
                    - &(&L * &(&(&QR2.X * &P.Z) - &(&P.X * &QR2.Z)));
                let f1pn = &n1 * &n3;
                let f1pd = &d2 * &d4;

                // f1_Q(A_P)
                let n1 = -&n1;
                let d2 = &(&R1.X * &QR2.Z) - &(&QR2.X * &R1.Z);
                let L = &(&R2.Y * &Q.Z) - &(&Q.Y * &R2.Z);
                let T = &(&R2.X * &Q.Z) - &(&Q.X * &R2.Z);
                let n3 = &(&T * &(&(&R1.Y * &Q.Z) - &(&Q.Y * &R1.Z)))
                    - &(&L * &(&(&R1.X * &Q.Z) - &(&Q.X * &R1.Z)));
                let d4 = &(&T * &(&(&PR1.Y * &Q.Z) - &(&Q.Y * &PR1.Z)))
                    - &(&L * &(&(&PR1.X * &Q.Z) - &(&Q.X * &PR1.Z)));
                let f1qn = &n1 * &n3;
                let f1qd = &d2 * &d4;

                // Get the affine coordinates for P, Q, R1, R2, P+R1 and Q+R2,
                // and also normalize lamb2p, lamb2q, f1p and f1q; this makes
                // everything else faster.
                let g1 = &P.Z * &Q.Z;
                let g2 = &g1 * &R1.Z;
                let g3 = &g2 * &R2.Z;
                let g4 = &g3 * &PR1.Z;
                let g5 = &g4 * &QR2.Z;
                let g6 = &g5 * &Tp;
                let g7 = &g6 * &Tq;
                let g8 = &g7 * &f1pd;
                let mut gg = (&g8 * &f1qd).invert();
                let if1qd = &gg * &g8;
                gg *= &f1qd;
                let if1pd = &gg * &g7;
                gg *= &f1pd;
                let iTq = &gg * &g6;
                gg *= &Tq;
                let iTp = &gg * &g5;
                gg *= &Tp;
                let iZqr2 = &gg * &g4;
                gg *= &QR2.Z;
                let iZpr1 = &gg * &g3;
                gg *= &PR1.Z;
                let iZr2 = &gg * &g2;
                gg *= &R2.Z;
                let iZr1 = &gg * &g1;
                gg *= &R1.Z;
                let iZq = &gg * &P.Z;
                let iZp = &gg * &Q.Z;

                let (xp, yp) = (&P.X * &iZp, &P.Y * &iZp);
                let (xq, yq) = (&Q.X * &iZq, &Q.Y * &iZq);
                let (xr1, yr1) = (&R1.X * &iZr1, &R1.Y * &iZr1);
                let (xr2, yr2) = (&R2.X * &iZr2, &R2.Y * &iZr2);
                let (xpr1, ypr1) = (&PR1.X * &iZpr1, &PR1.Y * &iZpr1);
                let (xqr2, yqr2) = (&QR2.X * &iZqr2, &QR2.Y * &iZqr2);
                let lamb2p = &Lp * &iTp;
                let lamb2q = &Lq * &iTq;
                let f1p = &f1pn * &if1pd;
                let f1q = &f1qn * &if1qd;

                // Set the pairing context for f_P(A_Q). In that context,
                // (xq,yq) is Q+R2, and (xr,yr) is R2.
                let mut pc = PairingContext {
                    xp: xp,
                    yp: yp,
                    xq: xqr2,
                    yq: yqr2,
                    xr: xr2,
                    yr: yr2,
                    f1: f1p,
                    lamb2: lamb2p,
                    xq_xp: &xqr2 - &xp,
                    yq_yp: &yqr2 - &yp,
                    xr_xp: &xr2 - &xp,
                    yr_yp: &yr2 - &yp,
                    xp_A: &xp + &self.A,
                    vn: f1p,
                    vd: Fq::ONE,
                };

                // Compute f_P(A_Q).
                // We keep the x coordinate of (n/2)*P in hnpX/hnpZ.
                let mut S = *P;
                for _ in 1..e {
                    self.line_double(&mut pc, &mut S);
                }
                let hnpX = S.X;
                let hnpZ = S.Z;
                self.line_double(&mut pc, &mut S);
                let okp = S.isinfinity();
                let vn1 = pc.vn;
                let vd1 = pc.vd;

                // Prepare the pairing context for f_Q(A_P).
                let mut pc = PairingContext {
                    xp: xq,
                    yp: yq,
                    xq: xpr1,
                    yq: ypr1,
                    xr: xr1,
                    yr: yr1,
                    f1: f1q,
                    lamb2: lamb2q,
                    xq_xp: &xpr1 - &xq,
                    yq_yp: &ypr1 - &yq,
                    xr_xp: &xr1 - &xq,
                    yr_yp: &yr1 - &yq,
                    xp_A: &xq + &self.A,
                    vn: f1q,
                    vd: Fq::ONE,
                };

                // Compute f_Q(A_P).
                // We keep the x coordinate of (n/2)*Q in hnqX/hnqZ.
                let mut S = *Q;
                for _ in 1..e {
                    self.line_double(&mut pc, &mut S);
                }
                let hnqX = S.X;
                let hnqZ = S.Z;
                self.line_double(&mut pc, &mut S);
                let okq = S.isinfinity();
                let vn2 = pc.vn;
                let vd2 = pc.vd;

                // w = wn/wd = f_P(A_Q) / f_Q(A_P)
                let wn = &vn1 * &vd2;
                let wd = &vd1 * &vn2;

                // Handling of points of order 1:
                //  - If P = inf or Q = inf, value is 1.
                //  - If P-Q = inf or P+Q = inf, value is 1.
                //  - If 2*P-Q = inf or P+2*Q = inf, value is 1.
                let p_z = P.Z.iszero();
                let q_z = Q.Z.iszero();
                let r1_z = R1.Z.iszero();
                let r2_z = R2.Z.iszero();
                let pr1_z = PR1.Z.iszero();
                let qr2_z = QR2.Z.iszero();
                let mut set1 = p_z | q_z | r1_z | r2_z | pr1_z | qr2_z;

                // Handling of points of order 2:
                //  - If P != inf and 2*P = inf, value is:
                //       -1 if (n/2)*Q != inf and (n/2)*Q != P
                //       1 otherwise
                //  - If Q != inf and 2*Q = inf, value is:
                //       -1 if (n/2)*P != inf and (n/2)*P != Q
                //       1 otherwise
                //  - If P-Q != inf and 2*(P-Q) = inf, value is:
                //       -1 if (n/2)*Q != inf and (n/2)*Q != P-Q
                //       1 otherwise
                let p_t = P.Y.iszero() & !p_z;
                let q_t = Q.Y.iszero() & !q_z;
                let r1_t = R1.Y.iszero() & !r1_z;
                let hnP_z = hnpZ.iszero();
                let hnQ_z = hnqZ.iszero();
                let hnQeqP =
                    (hnQ_z & p_z) | ((&hnqX * &P.Z).equals(&(&P.X * &hnqZ)) & p_t & !hnQ_z);
                let hnPeqQ =
                    (hnP_z & q_z) | ((&hnpX * &Q.Z).equals(&(&Q.X * &hnpZ)) & q_t & !hnP_z);
                let hnQeqR1 =
                    (hnQ_z & r1_z) | ((&hnqX * &R1.Z).equals(&(&R1.X * &hnqZ)) & r1_t & !hnQ_z);
                let mut neg1 = !set1
                    & ((p_t & !hnQ_z & !hnQeqP)
                        | (q_t & !hnP_z & !hnPeqQ)
                        | (r1_t & !hnQ_z & !hnQeqR1));
                set1 |= p_t | q_t | r1_t;

                // Handling of zeros:
                //  - If none of the cases above was encountered, a zero in
                //    either wn or wd implies colinearity, and the result is 1.
                set1 |= !set1 & (wn.iszero() | wd.iszero());

                // An error is reported if P or Q was not n-torsion. On error,
                // we force the result to 1.
                let ok = okp & okq;
                set1 |= !ok;
                neg1 &= ok;
                let mut w = &wn / &wd;
                w.set_cond(&Fq::ONE, set1);
                w.set_condneg(set1 & neg1);
                (w, ok)
            }
        }

        impl fmt::Display for Curve {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "Montgomery Curve with coefficient: {}, A24: {}", self.A, self.A24)
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CouplePoint {
            P1: Point,
            P2: Point,
        }

        impl CouplePoint {
            /// Return the pair of points at infinity: O1, O2 on E1 x E2
            pub const INFINITY: Self = Self {
                P1: Point::INFINITY,
                P2: Point::INFINITY,
            };

            /// Create a CouplePoint given a pair of points P1, P2 on E1 x E2
            pub fn new(P1: &Point, P2: &Point) -> Self {
                Self { P1: *P1, P2: *P2 }
            }

            /// Return the points P1, P2
            pub fn points(self) -> (Point, Point) {
                (self.P1, self.P2)
            }
        }

        /// Print debugging, not used within computations.
        impl fmt::Display for CouplePoint {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "Couple point with points:\n{}\n{}", self.P1, self.P2)
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CouplePointX {
            P1: PointX,
            P2: PointX,
        }

        impl CouplePointX {
            /// Create a CouplePointX given a pair of points P1, P2 on E1 x E2
            pub fn new(P1: &PointX, P2: &PointX) -> Self {
                Self { P1: *P1, P2: *P2 }
            }

            /// Return the points P1, P2
            pub fn points(self) -> (PointX, PointX) {
                (self.P1, self.P2)
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct EllipticProduct {
            E1: Curve,
            E2: Curve,
        }

        impl EllipticProduct {
            /// Create an EllipticProduct given a pair of elliptic curves of
            /// type Curve
            pub fn new(E1: &Curve, E2: &Curve) -> Self {
                Self { E1: *E1, E2: *E2 }
            }

            /// Return the pair of curves as a tuple
            pub fn curves(self) -> (Curve, Curve) {
                (self.E1, self.E2)
            }

            /// Addition of elements (P1, P2) and (Q1, Q2) on E1 x E2 is defined
            /// as (P1 + Q1, P2 + Q2). This function calls the add function for
            /// the pair of curves on the EllipticProduct
            pub fn add(self, C1: &CouplePoint, C2: &CouplePoint) -> CouplePoint {
                let mut C3 = CouplePoint::INFINITY;
                C3.P1 = self.E1.add(&C1.P1, &C2.P1);
                C3.P2 = self.E2.add(&C1.P2, &C2.P2);
                C3
            }

            /// Doubles the pair of points (P1, P2) on E1 x E2 as ([2]P1, [2]P2)
            pub fn double(self, C: &CouplePoint) -> CouplePoint {
                let mut C3 = *C;
                C3.P1 = self.E1.double(&C3.P1);
                C3.P2 = self.E2.double(&C3.P2);
                C3
            }

            pub fn x_double(self, C: &CouplePointX) -> CouplePointX {
                let mut P1 = C.P1;
                let mut P2 = C.P2;
                self.E1.xdbl(&mut P1.X, &mut P1.Z);
                self.E2.xdbl(&mut P2.X, &mut P2.Z);

                CouplePointX::new(&P1, &P2)
            }

            /// Repeatedly doubles the pair of points (P1, P2) on E1 x E2 to get
            /// ([2^n]P1, [2^n]P2)
            pub fn double_iter(self, C: &CouplePoint, n: usize) -> CouplePoint {
                let mut C3 = *C;
                C3.P1 = self.E1.double_iter(&C3.P1, n);
                C3.P2 = self.E2.double_iter(&C3.P2, n);
                C3
            }

            pub fn x_double_iter(self, C: &CouplePointX, n: usize) -> CouplePointX {
                let mut C3 = *C;
                C3.P1 = self.E1.x_double_iter(&C3.P1, n);
                C3.P2 = self.E2.x_double_iter(&C3.P2, n);
                C3
            }
        }

        // Implementation of 3-isogenies

        fn triple_e_point_iter_into(E: &Curve, P: &mut PointX, e: usize) {
            #[inline(always)]
            fn xTPL(E: &Curve, XP: &mut Fq, ZP: &mut Fq) {
                let mut X = *XP;
                let mut Z = *ZP;
                let XP0 = *XP;
                let ZP0 = *ZP;

                E.xdbl(&mut X, &mut Z);
                E.xadd(&XP0, &ZP0, &X, &Z, XP, ZP);
            }
            let mut X = P.X;
            let mut Z = P.Z;
            for _ in 0..e {
                xTPL(E, &mut X, &mut Z);
            }
            P.X = X;
            P.Z = Z;
        }

        /// Given a point P = (XP : ZP) of order 3, computes the
        /// 3-isogeny codomain with coefficient A represented as
        /// (A + 2C) / 4C (where A24_num = A + 2C, A24_denom = 4C)
        /// along with constants K1, K2 used for computing images
        #[inline]
        pub fn three_isogeny_codomain(P: &PointX) -> (Fq, Fq, Fq, Fq) {
            let K1 = &P.X - &P.Z;
            let K2 = &P.X + &P.Z;
            let R1 = K1.square();
            let R2 = K2.square();

            let mut R3 = R2 + R1;
            let mut R4 = K1 + K2;
            R4 = R4.square();
            R4 = R4 - R3;
            R3 = R4 + R2;
            R4 = R4 + R1;
            let mut R5 = R1 + R4;
            R5.set_mul2();
            R5 = R5 + R2;
            let mut A24_num = R5 * R3;
            R5 = R2 + R3;
            R5.set_mul2();
            R5 = R5 + R1;
            R5 = R5 * R4;
            let A24_denom = R5 - A24_num;
            A24_num += A24_denom; // TODO: simply R5?

            (A24_num, A24_denom, K1, K2)
        }

        /// Given constants (K1, K2) along with the point Q = (XQ : ZQ)
        /// compute the image of this point in place
        #[inline(always)]
        pub fn three_isogeny_image(K1: &Fq, K2: &Fq, Q: &mut PointX) {
            let mut t0 = &Q.X + &Q.Z;
            let mut t1 = &Q.X - &Q.Z;
            t0 *= K1;
            t1 *= K2;
            let mut t2 = &t0 + &t1;
            t0 = &t1 - &t0;
            t2.set_square();
            t0.set_square();
            Q.X *= &t2;
            Q.Z *= &t0;
        }

        /// 3^e isogeny chain using kernel
        /// Compute an isogeny between elliptic products, use an optimised
        /// strategy for all steps assuming doubling is always more expensive
        /// that images, which is not true for gluing.
        pub fn three_isogeny_chain<const N: usize>(
            E: &Curve,
            K: &PointX,
            eval_points: &[PointX; N],
            n: usize,
            strategy: &[usize],
        ) -> (Curve, [PointX; N]) {
            let mut kernel_pts = vec![*K];
            let mut image_points = *eval_points;

            let mut strat_idx = 0;
            let mut level: Vec<usize> = vec![0];
            let mut prev: usize;
            let mut kernel_len: usize;

            // For initalisation
            let mut S: PointX;
            let mut E_curr = *E;

            let mut A24_num;
            let mut A24_denom;
            let mut K1: Fq;
            let mut K2: Fq;

            for k in 0..n {
                prev = level.iter().sum();
                kernel_len = kernel_pts.len();

                // Recover the point from the list
                S = kernel_pts[kernel_len - 1];

                while prev != (n - 1 - k) {
                    // Add the next strategy to the level
                    level.push(strategy[strat_idx]);

                    // Triple the points according to the strategy
                    triple_e_point_iter_into(&E_curr, &mut S, strategy[strat_idx]);

                    // Add the point to the image points
                    kernel_pts.push(S);

                    // Update the strategy bookkeepping
                    prev += strategy[strat_idx];
                    strat_idx += 1;
                }

                // Clear out the used kernel point and update level
                kernel_pts.pop();
                level.pop();

                // Compute the codomain constants
                (A24_num, A24_denom, K1, K2) = three_isogeny_codomain(&S);
                E_curr = Curve::new_fromA24(&(A24_num/A24_denom));

                // Push all kernel points through the isogeny
                for ker in kernel_pts.iter_mut() {
                    three_isogeny_image(&K1, &K2, ker);
                }
                // Push the image points through the isogeny
                for imP in image_points.iter_mut() {
                    three_isogeny_image(&K1, &K2, imP);
                }
            }

            (E_curr, image_points)
        }

        /// Weil pairing
        
        // Context for a pairing computation. We are computing f_P(A_Q)
        // with A_Q = <Q+R2> - <R2>, and f_P being the rational function
        // of divisor n*<P+R1> - n*<R1>.
        //
        // Contents:
        //    xp, yp      affine coordinates of P
        //    xq, yq      affine coordinates of Q+R2
        //    xr, yr      affine coordinates of R2
        //    f1          f1_P(A_Q)
        //    lamb2       slope of tangent on P
        //    xq_xp       xq - xp
        //    yq_yp       yq - yp
        //    xr_xp       xr - xp
        //    yr_yp       yr - yp
        //    xp_A        xp + A
        //    vn, vd      current value (fraction vn/vd)
        //
        // The context assumes that the following special cases are handled
        // elsewhere by the caller:
        //    P = inf
        //    Q+R2 = inf
        //    R2 = inf
        //    2*P = inf
        struct PairingContext {
            xp: Fq,
            yp: Fq,
            xq: Fq,
            yq: Fq,
            xr: Fq,
            yr: Fq,
            f1: Fq,
            lamb2: Fq,
            xq_xp: Fq,
            yq_yp: Fq,
            xr_xp: Fq,
            yr_yp: Fq,
            xp_A: Fq,
            vn: Fq,
            vd: Fq,
        } 
        
        
    };
} // End of macro: define_ec_core

pub(crate) use define_ec_core;
