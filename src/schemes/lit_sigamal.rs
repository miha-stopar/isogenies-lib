macro_rules! define_litsigamal {
    () => {
        use crate::quaternion::quaternion_algebra::{QuatAlg, QuatAlgEl};
        use crate::quaternion::klpt::represent_integer;
        use crate::quaternion::quaternion_order::standard_maximal_extremal_order;
        use crate::util::{generate_random_range};
        use crate::ec_lit;
        use rug::integer::Order;

        pub fn get_params(lam: u32) -> (u32, u32, u32, u32) {
            let a = lam * 3;
            let b: u32;
            let c: u32;
            let f: u32;

            if lam == 128 {
                b = 162;
                c = 56;
                f = 30;
            } else if lam == 192 {
                b = 243;
                c = 83;
                f = 118;
            } else if lam == 256 {
                b = 324;
                c = 111;
                f = 436;
            } else {
                panic!("lam should be 128 or 192 or 256");
            }

            (a, b, c, f)
        }

        /// LITSiGamal struct
        #[derive(Clone, Debug)]
        pub struct LITSiGamal {
            curve: Curve,
            p: Integer,
            l_a: u32,
            l_b: u32,
            l_c: u32,
            a: u32,
            b: u32,
            c: u32,
            f: u32,
            scalar_without_b: Integer, // p + 1 = scalar_without_b * l_b**b
            n: u32,
        }

        impl LITSiGamal {
            pub fn new(curve: Curve, p: Integer, l_a: u32, l_b: u32, l_c: u32, a: u32, b: u32, c: u32, f: u32, n: u32) -> Self {
                let scalar_without_b = l_a.big().pow(a + 2) * l_c.big().pow(c) * f;

                Self {
                    curve, p, l_a, l_b, l_c, a, b, c, f, scalar_without_b, n
                }
            }

            pub fn generate_pub_key(&self) { 
                let fileName = "src/schemes/precomputed.json";
                let (PaX, QaX, PmQaX, mat2_2, mat2_3, mat2_4, power_a) = load_torsion_info(fileName, "lit-sigamal-128", 2);
                let (PbX, QbX, PmQbX, mat3_2, mat3_3, mat3_4, power_b) = load_torsion_info(fileName, "lit-sigamal-128", 3);
                let (PcX, QcX, PmQcX, mat5_2, mat5_3, mat5_4, power_c) = load_torsion_info(fileName, "lit-sigamal-128", 5);

                let power_a = power_a as u32;
                let power_b = power_b as u32;
                let power_c = power_c as u32;
    
                let (Pa, ok1) = self.curve.complete_pointX(&PaX); // TODO: check ok
                assert!(ok1 == 0xFFFFFFFF);

                let (mut Qa, ok2) = self.curve.complete_pointX(&QaX); // TODO: check ok

                let (Pb, ok1) = self.curve.complete_pointX(&PbX); // TODO: check ok
                assert!(ok1 == 0xFFFFFFFF);
                let (Qb, ok2) = self.curve.complete_pointX(&QbX); // TODO: check ok

                let (mut Pc, ok1) = self.curve.complete_pointX(&PcX); // TODO: check ok
                let (Qc, ok2) = self.curve.complete_pointX(&QcX); // TODO: check ok
                let (PmQc, ok3) = self.curve.complete_pointX(&PmQcX); // TODO: check ok
 
                let l_a = 2;
                let l_b = 3;
                let l_c = 5;
 
                let mut a1 = generate_random_range(0.big(), l_a.big().pow(power_a) - 1);
                let mut a2 = generate_random_range(0.big(), l_a.big().pow(power_a - 1)) * 2 + 1;
                // TODO: remove hardcoded a1, a2
                a1 = "67821213147103272121758535926151039005030912112029727951085328323311503610545609587044506281792717705600807906711952".big();
                a2 = "108275711446779976206444241046503643726530007828495696334502130746308494739140748101110547204481156360035409488865063".big();
                // Randomize Qa
                // Qa = a1 * Pa + a2 * Qa
                // TODO: use dblmul
                /*
                Note: don't do this here because in apply_ we need the original Qa
                let mut bytes = big_to_bytes(a1.clone());
                let a1Pa = self.curve.mul(&Pa, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(a2.clone());
                let a2Qa = self.curve.mul(&Qa, &bytes, bytes.len() * 8);
                Qa = self.curve.add(&a1Pa, &a2Qa);
                */

                let mut c1: Integer;
                let mut c2: Integer;

                let five = Integer::from(5);
                let zero = Integer::ZERO;
                loop {
                    c1 = generate_random_range(0.big(), l_c.big().pow(power_c) - 1);
                    c2 = generate_random_range(0.big(), l_c.big().pow(power_c) - 1);
                    // TODO: remove c1, c2
                    c1 = Integer::from(1);
                    c2 = Integer::from(1);
                    if c1.clone().modulo(&five) != zero && c2.clone().modulo(&five) != zero {
                        break;
                    }
                }
                // Pc = c1 * Pc + c2 * Qc
                // TODO: remove
                /*
                let mut bytes = big_to_bytes(c1.clone());
                let c1Pc = self.curve.mul(&Pc, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(c2.clone());
                let c2Qc = self.curve.mul(&Qc, &bytes, bytes.len() * 8);
                let Pc = self.curve.add(&c1Pc, &c2Qc);
                */

                // let mut k_digits = c1.to_digits::<u64>(Order::MsfLe);
                // k_digits.reverse();
                let k_digits: Vec<u64> = [1, 0, 0, 0].to_vec(); // TODO: remove, JUST DEBUGGING

                // let mut l_digits = c2.to_digits::<u64>(Order::MsfLe);
                // l_digits.reverse();
                let l_digits: Vec<u64> = [1, 0, 0, 0].to_vec(); // TODO: remove, JUST DEBUGGING

                let f: usize = 36; // TODO
                let Pc = self.curve.xdblmul_bounded(&Pc, &k_digits, &Qc, &l_digits, &PmQc, f);

                let Rx = PointX::new_xz(&Pc.X, &Pc.Z);

                // N = (l_a**(2*a) - n**2) * l_b**b
                let tau = l_a.big().pow(power_a * 2) - self.n * self.n;
                let N = tau.clone() * l_b.big().pow(power_b);
                let gamma = self.generate_gamma(N);

                // TODO: remove, just for debugging
                let qa = QuatAlg::new(-self.p.clone());
                // let gamma = QuatAlgEl::new("101064425814894250297266458592805333959492115142777713483064160339396739459362670629011189775195355480432667912928029565585343417641447712715314376809237973".big(), "93533154661012604105678420765706116710130443509048775286300216094985327466258160413386761474232055531385270204373695404255540646882449862261480847728263402".big(), "4753400567059888737315986212971719000".big(), "20561456972084008333059843341602645649".big(), 2.big(), qa.clone());
                // let gamma = QuatAlgEl::new("2011185045392199242481466685241103431236289913375605487988529329864805134553257030694326259418924424563378327667027885353428994775042308664556073155839862".big(), "10853501034135357035889336862977087605953185135487653652185963537341394638859475493620188339725879895509904749276277210794053785439623299471028262804950610".big(), "-10310443381405150524300858171174973139".big(), "-5988292950171597878820895074482233260".big(), 2.big(), qa.clone());
                // let gamma = QuatAlgEl::new("8".big(), "4".big(), "6".big(), "2".big(), 2.big(), qa.clone());

                let gamma = QuatAlgEl::new("8425099297649819352166379388758583438269420404259906871066684285709655835562922681747598532398853172000034459777928781756277965902285271705188053172406860".big(), "10672475422448534529312519343800163927437511879170590252915490695777092409588869040995562812407263134190749362059361293778411214787787911660929131931909741".big(), "7138531157822206556546664362792862977".big(), "6507363534081204197854356661058574173".big(), 1.big(), qa.clone());
                
                println!("");
                println!("");
                println!("");
                println!("gamma: {:?}", gamma);
                println!("");

                let order = standard_maximal_extremal_order().order;
                let (mut coord, imprim) = gamma.factor_in_order(order.lattice.clone());

                println!("coord: {:?}", coord);
                println!("");
                println!("imprim: {:?}", imprim);
                println!("");
                println!("");
 
                let torsion_a = l_a.big().pow(power_a);
                // let (Pa_gamma, _, _) = apply_endomorphism_on_torsion_group(&self.curve, coord.clone(), imprim.clone(), torsion_a, mat2_2, mat2_3, mat2_4, &Pa, &Qa);
                
                let (Pa_gamma, Qa_gamma, _) = apply_endomorphism_on_torsion_group(&self.curve, coord.clone(), imprim.clone(), torsion_a, mat2_2, mat2_3, mat2_4, &Pa, &Qa);

                println!("");
                println!("Pa_gamma");
                println!("{}", Pa_gamma.X / Pa_gamma.Z);
                println!("");

                /*
                println!("Qa_gamma");
                println!("{}", Qa_gamma.X / Qa_gamma.Z);
                println!("");
                */


                let bytes1 = big_to_bytes(a1);
                let a1Pa_gamma = self.curve.mul(&Pa_gamma, &bytes1, bytes1.len() * 8);
                let bytes2 = big_to_bytes(a2);
                let a2Qa_gamma = self.curve.mul(&Qa_gamma, &bytes2, bytes2.len() * 8);
                let Qa_rand_gamma = self.curve.add(&a1Pa_gamma, &a2Qa_gamma);


                let a1Pa = self.curve.mul(&Pa, &bytes1, bytes1.len() * 8);
                let a2Qa = self.curve.mul(&Qa, &bytes2, bytes2.len() * 8);
                let Qa_rand = self.curve.add(&a1Pa, &a2Qa);
                let Qa_rand_X = PointX::new_xz(&Qa_rand.X, &Qa_rand.Z);


                let torsion_b = l_b.big().pow(power_b);
                let (Pb_gamma, Qb_gamma, PmQb_gamma) = apply_endomorphism_on_torsion_group(&self.curve, coord.clone(), imprim.clone(), torsion_b.clone(), mat3_2, mat3_3, mat3_4, &Pb, &Qb);

                let torsion_c = l_c.big().pow(power_c);
                let (Pc_gamma, Qc_gamma, PmQc_gamma) = apply_endomorphism_on_torsion_group(&self.curve, coord, imprim, torsion_c, mat5_2, mat5_3, mat5_4, &Pc, &Qc);

                /*
                let mut bytes = big_to_bytes(c1);
                let c1Pc = self.curve.mul(&Pc_gamma, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(c2);
                let c2Qc = self.curve.mul(&Qc_gamma, &bytes, bytes.len() * 8);
                let R = self.curve.add(&c1Pc, &c2Qc);

                let Rx = PointX::new_xz(&R.X, &R.Z);
                */

                let dlog1 = ec_lit::dlog_3(&self.curve, &Pb_gamma, &Qb_gamma, 162);
                let dlog2 = ec_lit::dlog_3(&self.curve, &Qb_gamma, &Pb_gamma, 162);
 
                println!("dlog1: {}", dlog1);
                println!("");
                println!("dlog2: {}", dlog2);
                println!("");
                println!("");
                println!("");
                println!("");

                // TODO: set [2^a]Pa_gamma = (1,*)

                let mut bytes = big_to_bytes(dlog2); // TODO
                let dlog_Qb = self.curve.mul(&Qb, &bytes, bytes.len() * 8);
                let kernel1 = self.curve.sub(&Pb, &dlog_Qb); // TODO: PointX directly

                let kernel1x = PointX::new_xz(&kernel1.X, &kernel1.Z);

                println!("kernel1 X: {}", kernel1.X / kernel1.Z);
                println!("");
                println!("kernel Y: {}", kernel1.Y / kernel1.Z);
                println!("");

                // TODO 
                let n = 162;

                // Precomputed with strategy.py
                // let strategy: [usize; 35] = [21, 14, 5, 3, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1];
                // let strategy: [usize; 161] = [55, 38, 34, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 13, 8, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1];

                // From Sage LITSiGamal:
                let strategy: [usize; 161] = [65, 37, 23, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 9, 5, 4, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1];

                // endomorphism gamma is of order (l_a**(2*a) - n**2) * l_b**b
                // let's take the isogeny gamma1 with kernel: ker(gamma) \cap E[l_b**b]

                // TODO: eval_points are [Pa, Qa, R]

                // TODO: remove, just for testing, instead use the kernel of 
                // TODO: fix QaX above
                let eval_points = [PaX, Qa_rand_X, Rx]; 
                let (codomain, image_points) = ec_lit::three_isogeny_chain(&self.curve, &kernel1x, &eval_points, n, &strategy);
                // let (codomain, image_points) = ec_lit::three_isogeny_chain(&self.curve, &FF, &eval_points, n, &strategy);

                let (mut Pa_isog3X, mut Qa_isog3X, mut R_isog3X) = (image_points[0], image_points[1], image_points[2]);

                let (Pa, ok1) = codomain.complete_pointX(&Pa_isog3X);
                let (mut Qa, ok1) = codomain.complete_pointX(&Qa_isog3X);

                let (Pa_shift, Pa_shift1) = self.get_PQb_and_shift(&codomain, &self.curve, &Pa, &mut Qa, &Pa_gamma, &Qa_rand_gamma, torsion_b, tau);

                println!("++++++=================++++++++");
                println!("Pa {}", Pa.X / Pa.Z);
                println!("");

                println!("Qa {}", Qa_rand_X.X / Qa_rand_X.Z);
                println!("");
               
                println!("");
                println!("");
                println!("");
                println!("codomain: {}", codomain);
                println!("");
                println!("Pa after isog 3");
                println!("{}", image_points[0].X / image_points[0].Z);
                println!("");

                println!("Qa after isog 3");
                println!("{}", image_points[1].X / image_points[1].Z);
                println!("");

                println!("");
                println!("Pa1");
                println!("{}", Pa_gamma.X / Pa_gamma.Z);
                println!("");

                println!("Qa1");
                println!("{}", Qa_rand_gamma.X / Qa_rand_gamma.Z);
                println!("");

                println!("");
                println!("Pb_gamma");
                println!("{}", Pb_gamma.X / Pb_gamma.Z);
                println!("");

                println!("Qb_gamma");
                println!("{}", Qb_gamma.X / Qb_gamma.Z);
                println!("");

                println!("R");
                println!("{}", image_points[2].X / image_points[2].Z);
                println!("");

                /*
                println!("{}", Qa_new);
                println!("");
                */

                // TODO: compute the discrete logarithm s such that Qa_new = s * Pa_new
                // K := Qa_new - s * Pa_new
                // It holds: gamma(K) = 0
                // K is the kernel of l_b**b isogeny, let's denote this isogeny by gamma1

                // It holds:
                // gamma = gamma1 \circ phi
                // where gamma1: E -> E1 is the isogeny with kernel K and
                // phi: E1 -> E is the isogeny of degree l_a**(2*a) - n**2 
                
                /*
                E
                
                |
    gamma1 (degree l_b**b)            gamma: E -> E (degree l_a**(2*a) - n**2)
                |
                V

                E1                                                            E
                */


                let ell_product = EllipticProduct::new(&codomain, &self.curve);

                // apply endomorphism [n] - we assume n = 3
                assert!(self.n == 3);
                codomain.xtriple(&mut Pa_isog3X.X, &mut Pa_isog3X.Z);
                codomain.xtriple(&mut Qa_isog3X.X, &mut Qa_isog3X.Z);

                /*
                bytes = big_to_bytes(Integer::from(self.n));
                let (Pa_isog3, ok1) = self.curve.complete_pointX(&Pa_isog3X);
                let Pa_mul3 = self.curve.mul(&Pa_isog3, &bytes, bytes.len() * 8);
                println!("=========");
                println!("{}", Pa_mul3.X / Pa_mul3.Z);
                */

                let Pa_gammaX = PointX::new_xz(&Pa_gamma.X, &Pa_gamma.Z);
                let Qa_gammaX = PointX::new_xz(&Qa_gamma.X, &Qa_gamma.Z);

                println!("??????? 2222 ?????????");
                println!("");
                // println!("{}", Pa_mul3.X / Pa_mul3.Z);
                println!("{}", Pa_isog3X.X / Pa_isog3X.Z);
                println!("");
                // println!("{}", Qa_mul3.X / Qa_mul3.Z);
                // println!("{}", Qa_isog3X.X / Qa_isog3X.Z);
                println!("{}", Pa_gammaX.X / Pa_gammaX.Z);
                println!("");
                println!("");
                println!("");
                
                let P1P2 = CouplePointX::new(&Pa_isog3X, &Pa_gammaX);
                let Q1Q2 = CouplePointX::new(&Qa_isog3X, &Qa_gammaX);

                // TODO:
                let image_points = vec![
                    P1P2
                    // CouplePoint::new(&self.two_dim.P, &self.two_dim.Q),
                    // CouplePoint::new(&self.two_dim.omegaP, &self.two_dim.omegaQ),
                ];

                // Precomputed with strategy.py
                // TODO: derive 383 from self.a
                let strategy: [usize; 383] = [144, 89, 55, 34, 27, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 8, 6, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 55, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1];

                compute_isogeny(
                    &ell_product,
                    &P1P2,
                    &Q1Q2,
                    // &image_points,
                    self.a as usize,
                );

                /*
                let (product, points) = product_isogeny(
                    &ell_product,
                    &P1P2,
                    &Q1Q2,
                    &image_points,
                    self.a as usize,
                    &strategy,
                );

                println!("");
                println!("product E1: {}", product.E1);
                println!("");
                println!("product E2: {}", product.E2);
                println!("");
                */

                /*
                for i in 0..6 {
                    let (product, points) = product_isogeny(
                        &ell_product,
                        &P1P2,
                        &Q1Q2,
                        &image_points,
                        self.a as usize,
                        &strategy,
                    );
                }
                */

            }

            fn get_PQb_and_shift(&self, curve_1: &Curve, curve_2: &Curve, Pa: &Point, Qa: &mut Point,
                    Pa1: &Point, Qa1: &Point, torsion_b: Integer, tau: Integer) -> (Point, Point) {
                let Pb1 = generate_random_fq(curve_1, torsion_b.clone(), self.scalar_without_b.clone());
                let Qb1 = generate_random_fq(curve_1, torsion_b, self.scalar_without_b.clone());

                let t = self.l_a.big().pow(self.a + 2);
                let mut bytes = big_to_bytes(t);
                let Pa_check = curve_1.mul(&Pa, &bytes, bytes.len() * 8);

                println!("");
                println!("check: {}", Pa_check.X / Pa_check.Z);
                println!("");

                // let (foo, ok) = codomain.weil_pairing_2exp((self.a+2).try_into().unwrap(), &Pa_isog3, &Qa_isog3);
                let t = (self.a+2).try_into().unwrap();
                let (mut w1, ok) = curve_1.weil_pairing_2exp(t, Pa, Qa);
                bytes = big_to_bytes(tau);
                w1.set_pow_simple(&bytes);

                let (w2, ok) = curve_2.weil_pairing_2exp(t, Pa1, Qa1);
                println!("ok: {}", ok);
                println!("");

                if w1.equals(&w2) == 0xFFFFFFFF{
                    println!("========= ========= ======== ===========");
                    println!("========= ========= ======== ===========");
                    println!("========= ========= ======== ===========");
                    println!("========= ========= ======== ===========");
                } else {
                    Qa.Y.set_neg();
                    println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    println!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                }

                // TODO: this is again computed using X coordinates in compute_isogeny, perhaps
                // use this result
                let shift_1 = curve_1.double_iter(Pa, self.a as usize);
                let shift_2 = curve_2.double_iter(Pa1, self.a as usize);

                let Pa_shift = curve_1.add(Pa, &shift_1);
                let Pa1_shift = curve_2.add(Pa1, &shift_2);

                println!("??????? 1111 ?????????");
                println!("");
                // println!("{}", Pa_mul3.X / Pa_mul3.Z);
                println!("{}", Pa.X / Pa.Z);
                println!("");
                // println!("{}", Qa_mul3.X / Qa_mul3.Z);
                println!("{}", Qa.X / Qa.Z);
                println!("");
                println!("");
                println!("");

                (Pa_shift, Pa1_shift)
            }

            fn generate_gamma(&self, N: Integer) -> QuatAlgEl { 
                let mut bad_prod_primes = "140227657289781369".big();
                bad_prod_primes *= "8695006970070847579".big();
                bad_prod_primes *= "4359375434796427649".big();
                bad_prod_primes *= "221191130330393351".big();
                bad_prod_primes *= "1516192381681334191".big();
                bad_prod_primes *= "5474546011261709671".big();

                let qa = QuatAlg::new(-self.p.clone());
                let order = standard_maximal_extremal_order().order;

                let (gamma, _) =
                    represent_integer(N.clone(), qa.clone(), order.clone(), bad_prod_primes).unwrap();

                // TODO: remove
                assert!(*gamma.reduced_norm().numer() == N);

                gamma
            }
        }
    };
} // End of macro: define_litsigamal

pub(crate) use define_litsigamal;
