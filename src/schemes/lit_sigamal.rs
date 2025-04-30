macro_rules! define_litsigamal {
    () => {
        use crate::quaternion::quaternion_algebra::{QuatAlg, QuatAlgEl};
        use crate::quaternion::klpt::represent_integer;
        use crate::quaternion::quaternion_order::standard_maximal_extremal_order;
        use crate::util::{generate_random_range};
        use crate::ec_lit;
        use std::time::Instant;

        #[derive(Clone, Debug)]
        pub struct PubKeyPoints {
            pub Pa: PointX,
            pub Qa: PointX,
            pub Pb: PointX,
            pub Qb: PointX,
            pub PQb: PointX,
            pub R: PointX,
        }

        impl PubKeyPoints {
            pub fn new(
                Pa: PointX,
                Qa: PointX,
                Pb: PointX,
                Qb: PointX,
                PQb: PointX,
                R: PointX,
            ) -> Self {
                Self {
                    Pa,
                    Qa,
                    Pb,
                    Qb,
                    PQb,
                    R,
                }
            }
        }

        #[derive(Clone, Debug)]
        pub struct PubKey {
            pub points: PubKeyPoints,
            pub points1: PubKeyPoints,
            pub power_a: u32,
            pub power_b: u32,
            pub power_c: u32,
            pub tau: Integer,
        }

        impl PubKey {
            pub fn new(
                points: PubKeyPoints,
                points1: PubKeyPoints,
                power_a: u32,
                power_b: u32,
                power_c: u32,
                tau: Integer,
            ) -> Self {
                Self {
                    points,
                    points1,
                    power_a,
                    power_b,
                    power_c,
                    tau,
                }
            }
        }

        #[derive(Clone, Debug)]
        pub struct Cipher {
            pub points: CipherPoints,
            pub points1: CipherPoints,
            pub A24_num: Fq,
            pub A24_denom: Fq,
            pub A24_1_num: Fq,
            pub A24_1_denom: Fq,
        }

        impl Cipher {
            pub fn new(
                points: CipherPoints,
                points1: CipherPoints,
                A24_num: Fq,
                A24_denom: Fq,
                A24_1_num: Fq,
                A24_1_denom: Fq,
            ) -> Self {
                Self {
                    points,
                    points1,
                    A24_num,
                    A24_denom,
                    A24_1_num,
                    A24_1_denom,
                }
            }
        }

        #[derive(Clone, Debug)]
        pub struct CipherPoints {
            pub Pa: PointX,
            pub Qa: PointX,
            pub R: PointX,
        }

        impl CipherPoints {
            pub fn new(
                Pa: PointX,
                Qa: PointX,
                R: PointX,
            ) -> Self {
                Self {
                    Pa,
                    Qa,
                    R,
                }
            }
        }

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

            pub fn generate_pub_key(&self) -> (PubKey, (Integer, Integer, Integer)) {
                let start = Instant::now();

                let fileName = "src/schemes/precomputed.json";
                let (PaX, QaX, PmQaX, mat2_2, mat2_3, mat2_4, power_a) = load_torsion_info(fileName, "lit-sigamal-128", 2);
                let (PbX, QbX, PmQbX, mat3_2, mat3_3, mat3_4, power_b) = load_torsion_info(fileName, "lit-sigamal-128", 3);
                let (PcX, QcX, PmQcX, mat5_2, mat5_3, mat5_4, power_c) = load_torsion_info(fileName, "lit-sigamal-128", 5);

                let power_a = power_a as u32;
                let power_b = power_b as u32;
                let power_c = power_c as u32;

                let mut curve = self.curve.clone();
    
                // TODO: prepare completions
                let (Pa, ok1) = curve.complete_pointX(&PaX); // TODO: check ok
                assert!(ok1 == 0xFFFFFFFF);
                let (mut Qa, ok2) = curve.complete_pointX(&QaX); // TODO: check ok
                let (PmQa, ok2) = curve.complete_pointX(&PmQaX); // TODO: check ok

                let (Pb, ok1) = curve.complete_pointX(&PbX); // TODO: check ok
                assert!(ok1 == 0xFFFFFFFF);
                let (Qb, ok2) = curve.complete_pointX(&QbX); // TODO: check ok
                let (PmQb, ok3) = curve.complete_pointX(&PmQbX); // TODO: check ok

                let (mut Pc, ok1) = curve.complete_pointX(&PcX); // TODO: check ok
                let (Qc, ok2) = curve.complete_pointX(&QcX); // TODO: check ok
                let (PmQc, ok3) = curve.complete_pointX(&PmQcX); // TODO: check ok
 
                let l_a = 2;
                let l_b = 3;
                let l_c = 5;
 
                let mut a1 = generate_random_range(0.big(), l_a.big().pow(power_a) - 1);
                let mut a2 = generate_random_range(0.big(), l_a.big().pow(power_a - 1)) * l_a + 1;
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
                    // c1 = Integer::from(1);
                    // c2 = Integer::from(1);
                    if c1.clone().modulo(&five) != zero && c2.clone().modulo(&five) != zero {
                        break;
                    }
                }
                // Pc = c1 * Pc + c2 * Qc
                // TODO: remove
                let mut bytes = big_to_bytes(c1.clone());
                let c1Pc = self.curve.mul(&Pc, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(c2.clone());
                let c2Qc = self.curve.mul(&Qc, &bytes, bytes.len() * 8);
                // let Pc1 = self.curve.add(&c1Pc, &c2Qc);
                Pc = self.curve.add(&c1Pc, &c2Qc);

                /*
                let mut k_digits = c1.to_digits::<u64>(Order::MsfLe);
                k_digits.reverse();
                // TODO

                // let k_digits: Vec<u64> = [1, 0, 0, 0].to_vec(); // TODO: remove, JUST DEBUGGING

                let mut l_digits = c2.to_digits::<u64>(Order::MsfLe);
                l_digits.reverse();
                while l_digits.len() < 4 {
                    l_digits.push(0);
                }
                // let l_digits: Vec<u64> = [1, 0, 0, 0].to_vec(); // TODO: remove, JUST DEBUGGING

                let f: usize = 36; // TODO
                let Pc = self.curve.xdblmul_bounded(&Pc, &k_digits, &Qc, &l_digits, &PmQc, f);
                 */


                // assert!(Pc.equals(&Pc1) == 0xFFFFFFFF);

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
                
                let order = standard_maximal_extremal_order().order;
                let (mut coord, imprim) = gamma.factor_in_order(order.lattice.clone());

                println!("1: {:?}", start.elapsed());
                let second_part = Instant::now();

                let torsion_a = l_a.big().pow(power_a);
                // let (Pa_gamma, _, _) = apply_endomorphism_on_torsion_group(&self.curve, coord.clone(), imprim.clone(), torsion_a, mat2_2, mat2_3, mat2_4, &Pa, &Qa);
                
                let (Pa_gamma, Qa_gamma, R_gamma) = apply_endomorphism_on_torsion_group(&curve, coord.clone(), imprim.clone(), torsion_a, mat2_2, mat2_3, mat2_4, &Pa, &Qa, &PmQa);

                let mut Pa1_to_be_mapped = PointX::new_xz(&Pa_gamma.X, &Pa_gamma.Z);
                let mut Qa1_to_be_mapped = PointX::new_xz(&Qa_gamma.X, &Qa_gamma.Z);
                let mut Ra1_to_be_mapped = PointX::new_xz(&R_gamma.X, &R_gamma.Z);

                let bytes1 = big_to_bytes(a1);
                let a1Pa_gamma = curve.mul(&Pa_gamma, &bytes1, bytes1.len() * 8);
                let bytes2 = big_to_bytes(a2);
                let a2Qa_gamma = curve.mul(&Qa_gamma, &bytes2, bytes2.len() * 8);
                let Qa_rand_gamma = curve.add(&a1Pa_gamma, &a2Qa_gamma);


                let a1Pa = curve.mul(&Pa, &bytes1, bytes1.len() * 8);
                let a2Qa = curve.mul(&Qa, &bytes2, bytes2.len() * 8);
                let Qa_rand = curve.add(&a1Pa, &a2Qa);
                let Qa_rand_X = PointX::new_xz(&Qa_rand.X, &Qa_rand.Z);


                let torsion_b = l_b.big().pow(power_b);
                let (Pb_gamma, Qb_gamma, PmQb_gamma) = apply_endomorphism_on_torsion_group(&curve, coord.clone(), imprim.clone(), torsion_b.clone(), mat3_2, mat3_3, mat3_4, &Pb, &Qb, &PmQb);

                let torsion_c = l_c.big().pow(power_c);
                let (Pc_gamma, Qc_gamma, PmQc_gamma) = apply_endomorphism_on_torsion_group(&curve, coord, imprim, torsion_c, mat5_2, mat5_3, mat5_4, &Pc, &Qc, &PmQc);

                /*
                let mut bytes = big_to_bytes(c1);
                let c1Pc = self.curve.mul(&Pc_gamma, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(c2);
                let c2Qc = self.curve.mul(&Qc_gamma, &bytes, bytes.len() * 8);
                let R = self.curve.add(&c1Pc, &c2Qc);

                let Rx = PointX::new_xz(&R.X, &R.Z);
                */

                println!("2: {:?}", second_part.elapsed());
                let third_part = Instant::now();

                let dlog1 = ec_lit::dlog_3(&curve, &Pb_gamma, &Qb_gamma, 162);
                let dlog2 = ec_lit::dlog_3(&curve, &Qb_gamma, &Pb_gamma, 162);
 
                // TODO: set [2^a]Pa_gamma = (1,*)

                let mut bytes = big_to_bytes(dlog2); // TODO
                let dlog_Qb = curve.mul(&Qb, &bytes, bytes.len() * 8);
                let kernel1 = curve.sub(&Pb, &dlog_Qb); // TODO: PointX directly

                let kernel1x = PointX::new_xz(&kernel1.X, &kernel1.Z);

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
                let (mut codomain, mut image_points) = ec_lit::three_isogeny_chain(&curve, &kernel1x, eval_points.to_vec(), n, &strategy);
                let (mut Pa_isog3X, mut Qa_rand_isog3X, R_isog3X) = (image_points[0], image_points[1], image_points[2]);

                // TODO
                let mut Pa_to_be_mapped = Pa_isog3X.clone(); // to be mapped by 3-isogeny chain
                let mut Qa_to_be_mapped = Qa_rand_isog3X.clone();
                let mut Ra_to_be_mapped = R_isog3X.clone();

                let (Pa, ok1) = codomain.complete_pointX(&Pa_isog3X);
                let (mut Qa_rand, ok1) = codomain.complete_pointX(&Qa_rand_isog3X);

                let (Pa_shift, Qa_shift, Pa1_shift, Qa1_shift, Pb, Qb, PQb, Pb_shift, Qb_shift, PQb_shift) = self.get_PQb_and_shift(&codomain, &curve, &Pa, &mut Qa_rand, &Pa_gamma, &Qa_rand_gamma, torsion_b.clone(), tau.clone());

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

                let ell_product = EllipticProduct::new(&codomain, &curve);
 
                let Pa_gammaX = PointX::new_xz(&Pa_gamma.X, &Pa_gamma.Z);
                let Qa_rand_gammaX = PointX::new_xz(&Qa_rand_gamma.X, &Qa_rand_gamma.Z);

                let mut Pa_shiftX = PointX::new_xz(&Pa_shift.X, &Pa_shift.Z);
                let mut Qa_shiftX = PointX::new_xz(&Qa_shift.X, &Qa_shift.Z);

                let Pa1_shiftX = PointX::new_xz(&Pa1_shift.X, &Pa1_shift.Z);
                let Qa1_shiftX = PointX::new_xz(&Qa1_shift.X, &Qa1_shift.Z);

                let mut Pb_to_be_mapped = PointX::new_xz(&Pb.X, &Pb.Z);
                let mut Qb_to_be_mapped = PointX::new_xz(&Qb.X, &Qb.Z);
                let mut PQb_to_be_mapped = PointX::new_xz(&PQb.X, &PQb.Z);

                let Pb_shiftX = PointX::new_xz(&Pb_shift.X, &Pb_shift.Z);
                let Qb_shiftX = PointX::new_xz(&Qb_shift.X, &Qb_shift.Z);
                let PQb_shiftX = PointX::new_xz(&PQb_shift.X, &PQb_shift.Z);

                println!("3: {:?}", third_part.elapsed());
                let fourth_part = Instant::now();

                // Precomputed with strategy.py
                // TODO: derive 383 from self.a
                // let n = 383;
                // let strategy: [usize; 383] = [144, 89, 55, 34, 27, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 8, 6, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 55, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1];
                let strategy_2: [usize; 382] = [154, 93, 55, 33, 20, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 22, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 38, 22, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 65, 37, 21, 12, 7, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1];

                let points = compute_isogeny(
                    &ell_product,
                    &mut Pa_isog3X,
                    &mut Qa_rand_isog3X,
                    &Pa_gammaX,
                    &Qa_rand_gammaX,
                    &mut Pa_shiftX,
                    &mut Qa_shiftX,
                    &Pa1_shiftX,
                    &Qa1_shiftX,
                    &Pb_to_be_mapped,
                    &Qb_to_be_mapped,
                    &PQb_to_be_mapped,
                    &Pb_shiftX,
                    &Qb_shiftX,
                    &PQb_shiftX,
                    self.a as usize, // 384
                    &strategy_2,
                );

                let mut Pb1_to_be_mapped = points[0];
                let mut Qb1_to_be_mapped = points[1];
                let mut PQb1_to_be_mapped = points[2];

                let f_mul = l_b.big().pow(power_b - 1);
                let mut backtracking_check = PointX::INFINITY;

                println!("4: {:?}", fourth_part.elapsed());
                let fifth_part = Instant::now();

                for _ in 0..1 { // TODO: 0..6
                    let mut no_backtracking = false;
                    let mut s = 0.big();
                    let mut kernel = PointX::INFINITY;
                    while !no_backtracking {

                        s = generate_random_range(0.big(), l_b.big().pow(power_b - 1)) * l_b + 
                            generate_random_range(1.big(), 2.big()); // TODO: check if this is 1 or 2

                        kernel = codomain.ladder_3pt(&Pb_to_be_mapped, &Qb_to_be_mapped, &PQb_to_be_mapped, s.clone());
                        let mut check = kernel.clone();

                        codomain.xmul(&mut check, f_mul.clone());

                        if check.equals(&backtracking_check) != 0xFFFFFFFF {
                            no_backtracking = true;
                        }
                    }

                    let kernel1 = curve.ladder_3pt(&Pb1_to_be_mapped, &Qb1_to_be_mapped, &PQb1_to_be_mapped, s);

                    let eval_points = [Pa_to_be_mapped, Qa_to_be_mapped, Ra_to_be_mapped, Qb_to_be_mapped];
                    (codomain, image_points) = ec_lit::three_isogeny_chain(&codomain, &kernel, eval_points.to_vec(), n, &strategy);
                    (Pa_to_be_mapped, Qa_to_be_mapped, Ra_to_be_mapped, backtracking_check) = (image_points[0], image_points[1], image_points[2], image_points[3]);
                    codomain.xmul(&mut backtracking_check, f_mul.clone());

                    let eval_points = [Pa1_to_be_mapped, Qa1_to_be_mapped, Ra1_to_be_mapped]; 
                    (curve, image_points) = ec_lit::three_isogeny_chain(&curve, &kernel1, eval_points.to_vec(), n, &strategy);
                    (Pa1_to_be_mapped, Qa1_to_be_mapped, Ra1_to_be_mapped) = (image_points[0], image_points[1], image_points[2]);

                    let (Pa, ok1) = codomain.complete_pointX(&Pa_to_be_mapped);
                    let (mut Qa, ok1) = codomain.complete_pointX(&Qa_to_be_mapped);
                    let (Pa1, ok1) = codomain.complete_pointX(&Pa1_to_be_mapped);
                    let (Qa1, ok1) = codomain.complete_pointX(&Qa1_to_be_mapped);

                    let (Pa_shift, Qa_shift, Pa1_shift, Qa1_shift, Pb, Qb, PQb, Pb_shift, Qb_shift, PQb_shift) = self.get_PQb_and_shift(&codomain, &curve, &Pa, &mut Qa, &Pa1, &Qa1, torsion_b.clone(), tau.clone());

                    Pb_to_be_mapped = PointX::new_xz(&Pb.X, &Pb.Z);
                    Qb_to_be_mapped = PointX::new_xz(&Qb.X, &Qb.Z);
                    PQb_to_be_mapped = PointX::new_xz(&PQb.X, &PQb.Z);

                    // TODO
                    let mut Pa_shiftX = PointX::new_xz(&Pa_shift.X, &Pa_shift.Z);
                    let mut Qa_shiftX = PointX::new_xz(&Qa_shift.X, &Qa_shift.Z);

                    let Pa1_shiftX = PointX::new_xz(&Pa1_shift.X, &Pa1_shift.Z);
                    let Qa1_shiftX = PointX::new_xz(&Qa1_shift.X, &Qa1_shift.Z);

                    let Pb_shiftX = PointX::new_xz(&Pb_shift.X, &Pb_shift.Z);
                    let Qb_shiftX = PointX::new_xz(&Qb_shift.X, &Qb_shift.Z);
                    let PQb_shiftX = PointX::new_xz(&PQb_shift.X, &PQb_shift.Z);

                    let ell_product = EllipticProduct::new(&codomain, &curve);

                    let points = compute_isogeny(
                        &ell_product,
                        &mut Pa_to_be_mapped,
                        &mut Qa_to_be_mapped,
                        &Pa1_to_be_mapped,
                        &Qa1_to_be_mapped,
                        &mut Pa_shiftX,
                        &mut Qa_shiftX,
                        &Pa1_shiftX,
                        &Qa1_shiftX,
                        &Pb_to_be_mapped,
                        &Qb_to_be_mapped,
                        &PQb_to_be_mapped,
                        &Pb_shiftX,
                        &Qb_shiftX,
                        &PQb_shiftX,
                        self.a as usize,
                        &strategy_2,
                    );

                    Pb1_to_be_mapped = points[0];
                    Qb1_to_be_mapped = points[1];
                    PQb1_to_be_mapped = points[2];
                }

                println!("5: {:?}", fifth_part.elapsed());


                let alpha: Integer = generate_random_range(0.big(), l_c.big().pow(power_b - 1)) * l_c + 
                    generate_random_range(1.big(), 4.big()); // TODO: check if this is 1 or 4
                curve.xmul(&mut Ra1_to_be_mapped, alpha.clone());

                let alice_secret0: Integer = generate_random_range(0.big(), l_a.big().pow(power_a + 2) - 1) * l_a + 1;
                let alice_secret1: Integer = generate_random_range(0.big(), l_a.big().pow(power_a + 2) - 1) * l_a + 1;

                curve.xmul(&mut Pa1_to_be_mapped, alice_secret0.clone());
                curve.xmul(&mut Qa1_to_be_mapped, alice_secret1.clone());

                let pubkey_points = PubKeyPoints::new(
                    Pa_to_be_mapped,
                    Qa_to_be_mapped,
                    Pb_to_be_mapped,
                    Qb_to_be_mapped,
                    PQb_to_be_mapped,
                    Ra_to_be_mapped
                );

                let pubkey_points1 = PubKeyPoints::new(
                    Pa1_to_be_mapped,
                    Qa1_to_be_mapped,
                    Pb1_to_be_mapped,
                    Qb1_to_be_mapped,
                    PQb1_to_be_mapped,
                    Ra1_to_be_mapped
                );

                (PubKey::new(
                    pubkey_points,
                    pubkey_points1,
                    power_a,
                    power_b,
                    power_c,
                    tau,
                ), (alice_secret0, alice_secret1, alpha))
            }

            pub fn encrypt(&self, pub_key: &PubKey, mu: Integer) -> Cipher {
                let start = Instant::now();

                let (points, points1) = (&pub_key.points, &pub_key.points1);
                let (mut Pa, mut Qa, Pb, Qb, PQb, mut R) = (points.Pa, points.Qa, points.Pb, points.Qb, points.PQb, points.R);
                let (mut Pa1, mut Qa1, Pb1, Qb1, PQb1, mut R1) = (points1.Pa, points1.Qa, points1.Pb, points1.Qb, points1.PQb, points1.R);

                let A24 = get_montgomery_A24(&Pb, &Qb, &PQb);
                let A24_1 = get_montgomery_A24(&Pb1, &Qb1, &PQb1);

                let s = self.l_b * generate_random_range(0.big(), self.l_b.big().pow(pub_key.power_b) - 1) +
                            generate_random_range(1.big(), 2.big()); // TODO: check if this is 1 or 2

                let curve = Curve::new_fromA24(&A24.0, &A24.1);
                let curve1 = Curve::new_fromA24(&A24_1.0, &A24_1.1);
                let kernelx = curve.ladder_3pt(&Pb, &Qb, &PQb, s.clone());
                let kernel1x = curve1.ladder_3pt(&Pb1, &Qb1, &PQb1, s.clone());

                let n = 162;
                let strategy: [usize; 161] = [65, 37, 23, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 9, 5, 4, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1];

                let eval_points = [Pa, Qa, R];
                let (curve_new, image_points) = ec_lit::three_isogeny_chain(&curve, &kernelx, eval_points.to_vec(), n, &strategy);
                (Pa, Qa, R) = (image_points[0], image_points[1], image_points[2]);

                let eval_points = [Pa1, Qa1, R1];
                let (curve1_new, image_points) = ec_lit::three_isogeny_chain(&curve1, &kernel1x, eval_points.to_vec(), n, &strategy);
                (Pa1, Qa1, R1) = (image_points[0], image_points[1], image_points[2]);

                let beta = self.l_c * generate_random_range(0.big(), self.l_c.big().pow(pub_key.power_c - 1) - 1) +
                            generate_random_range(1.big(), 4.big()); // TODO: check if this can be 4

                curve.xmul(&mut R, beta.clone());
                curve1.xmul(&mut R1, beta * mu);

                let bob_secret_1: Integer = self.l_a * generate_random_range(0.big(), self.l_a.big().pow(pub_key.power_a + 2) - 1) + 1;
                let bob_secret_2: Integer = self.l_a * generate_random_range(0.big(), self.l_a.big().pow(pub_key.power_a + 2) - 1) + 1;

                curve.xmul(&mut Pa, bob_secret_1.clone());
                curve.xmul(&mut Qa, bob_secret_2.clone());

                curve1.xmul(&mut Pa1, bob_secret_1);
                curve1.xmul(&mut Qa1, bob_secret_2);

                let points = CipherPoints::new(
                    Pa,
                    Qa,
                    R,
                );

                let points1 = CipherPoints::new(
                    Pa1,
                    Qa1,
                    R1,
                );

                println!("encrypt: {:?}", start.elapsed());

                Cipher::new(
                    points,
                    points1,
                    curve_new.A24_num,
                    curve_new.A24_denom,
                    curve1_new.A24_num,
                    curve1_new.A24_denom,
                )
            }

            pub fn decrypt(&self, pub_key: &PubKey, cipher: &Cipher, alice_secret: (Integer, Integer, Integer)) {
                let start = Instant::now();

                let (points, points1) = (&cipher.points, &cipher.points1);
                let (mut Pa, mut Qa, mut R) = (points.Pa, points.Qa, points.R);
                let (mut Pa1, mut Qa1, mut R1) = (points1.Pa, points1.Qa, points1.R);

                let curve = Curve::new_fromA24(&cipher.A24_num, &cipher.A24_denom);
                let curve1 = Curve::new_fromA24(&cipher.A24_1_num, &cipher.A24_1_denom);

                let alpha = alice_secret.2;
                
                curve.xmul(&mut R, alpha);
                curve.xmul(&mut Pa, alice_secret.0);
                curve.xmul(&mut Qa, alice_secret.1);

                let (Pa_complete, _) = curve.complete_pointX(&Pa); // S[0]
                let (Qa_complete, _) = curve.complete_pointX(&Qa); // S[1]
                let (R_complete, _) = curve.complete_pointX(&R); // S[2]

                let (Pa1_complete, _) = curve1.complete_pointX(&Pa1); // S1[0]
                let (Qa1_complete, _) = curve1.complete_pointX(&Qa1);

                let t = self.l_a.big().pow(self.a);
                let mut bytes = big_to_bytes(t);
                let T1_0 = curve.mul(&Pa_complete, &bytes, bytes.len() * 8);
                let T1_1 = curve.mul(&Pa_complete, &bytes, bytes.len() * 8);

                // TODO: pairing check

                let Rs = curve.add(&R_complete, &T1_0);

                let Pa_shift = curve.add(&Pa_complete, &T1_0);
                let Qa_shift = curve.add(&Qa_complete, &T1_0);

                let Pa1_shift = curve.add(&Pa1_complete, &T1_1);
                let Qa1_shift = curve.add(&Qa1_complete, &T1_1);

                let Pa_shiftX = PointX::new_xz(&Pa_shift.X, &Pa_shift.Z);
                let Qa_shiftX = PointX::new_xz(&Qa_shift.X, &Qa_shift.Z);

                let Pa1_shiftX = PointX::new_xz(&Pa1_shift.X, &Pa1_shift.Z);
                let Qa1_shiftX = PointX::new_xz(&Qa1_shift.X, &Qa1_shift.Z);

                /*
                Rs = S[2] + T1_0
                R_shift = [Rs[0],Rs[2]]
                Pas = S[0] + T1_0
                Qas = S[1] + T1_0
                Pas1 = S1[0] + T1_1
                Qas1 = S1[1] + T1_1
                
                Pa_shift = [Pas[0],Pas[2]]
                Pa1_shift = [Pas1[0],Pas1[2]]
                Qa_shift = [Qas[0],Qas[2]]
                Qa1_shift = [Qas1[0],Qas1[2]]
                PQa_shift = [Pa_shift,Qa_shift,Pa1_shift,Qa1_shift]

                R_ = [R]
                R_shift_ = [R_shift]

                pointsR = compute_image_gamma0(PQa,PQa1,PQa_shift,R_,R_shift_,a,tau,sqrtminus,strategy_for_2dim)
                */

                let ell_product = EllipticProduct::new(&curve, &curve1);

                let strategy_2: [usize; 382] = [154, 93, 55, 33, 20, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 22, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 38, 22, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 65, 37, 21, 12, 7, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1];
                /*
                let points = compute_isogeny(
                        &ell_product,
                        &mut Pa,
                        &mut Qa,
                        &Pa1,
                        &Qa1,
                        &mut Pa_shiftX,
                        &mut Qa_shiftX,
                        &Pa1_shiftX,
                        &Qa1_shiftX,
                        &Pb_to_be_mapped,
                        &Qb_to_be_mapped,
                        &PQb_to_be_mapped,
                        &Pb_shiftX,
                        &Qb_shiftX,
                        &PQb_shiftX,
                        self.a as usize,
                        &strategy_2,
                    );
                */

                println!("decrypt: {:?}", start.elapsed());
            }

            fn get_PQb_and_shift(&self, curve_1: &Curve, curve_2: &Curve, Pa: &Point, Qa: &mut Point,
                    Pa1: &Point, Qa1: &Point, torsion_b: Integer, tau: Integer) -> (Point, Point, Point, Point, Point, Point, Point, Point, Point, Point) {
                let Pb1 = generate_random_fq(curve_1, torsion_b.clone(), self.scalar_without_b.clone());
                let Qb1 = generate_random_fq(curve_1, torsion_b, self.scalar_without_b.clone());

                // TODO: the following is only for debugging:
                /*
                let px = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("589347475779467028048732638707745169074099572448781051271349481559473727012677648788277674070708563815384409822118544148404313695739142076641733342195789664410029747217412725418639602788946160431051681659016463369120955499530350884793")),
                    &Fp::decode_reduce(&bytes_from_str("295541731903830973966435177821009510545806235315633443160350494708236538905366964828599192114390194343094738330899036091771606752115327511842840580779568519473668508004731498100416200716785989577694849244176022777068709128659699958651")),
                );

                let pz = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("89542995344710952948107022790100922686334678202826420161095874206526834236077000668814682919417119967989382371825802002025554791683974908178374722677745710441261906652561600152816296300222128057559398717792476238492791294646212107584")),
                    &Fp::decode_reduce(&bytes_from_str("1028199781928077642554616804903240208550974242321502258355144963022928866609527418284098186948125250284199941214438305372352745234972686884241582592073140746296427566777630409750957270938545946624116302693122101933302505969963727803020")),
                );
                let PX = PointX::new_xz(&px, &pz);
                let (mut Pb1, _) = curve_1.complete_pointX(&PX);
                Pb1.Y.set_neg();

                let qx = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("454140180413608233286394846003134365209314843053780796163906640030713711389055041986213108103635675206609008574032235475274293757590849226484000435247916029901267410923839011673395777802804904328624904321463415715607384207304067875730")),
                    &Fp::decode_reduce(&bytes_from_str("1243112448202039772158444090834398340074566723306283492337265843328542753962809102809078142471446057326889296295982219243530765736786289041817270190922322511604202629506602681892554629222996194581581219464530795781841550001797599410006")),
                );

                let qz = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("1052179633628854838348234790181091744705359231662359128521305087623977163897424533351031768785764629228680763061743717048694702723750786611598195928571348917667108362366464498799815518080430742521740861086356920849755947555393300743532")),
                    &Fp::decode_reduce(&bytes_from_str("1064275765062938759527616005874886327370713276008282087681071518299603510092884744948583012023499246802112511221842256640840835566271805426695601243375872188201751132299529941374828828373089760184532851360899361002657396410465738357472")),
                );

                let QX = PointX::new_xz(&qx, &qz);
                let (Qb1, _) = curve_1.complete_pointX(&QX);
                */
                // end of debugging

                let PQb1 = curve_1.sub(&Pb1, &Qb1);

                let t = self.l_a.big().pow(self.a + 2);
                let mut bytes = big_to_bytes(t);
                let Pa_check = curve_1.mul(&Pa, &bytes, bytes.len() * 8);

                // let (foo, ok) = codomain.weil_pairing_2exp((self.a+2).try_into().unwrap(), &Pa_isog3, &Qa_isog3);
                let t = (self.a+2).try_into().unwrap();
                let (mut w1, ok) = curve_1.weil_pairing_2exp(t, Pa, Qa);
                bytes = big_to_bytes(tau);
                w1.set_pow_simple(&bytes);

                let (w2, ok) = curve_2.weil_pairing_2exp(t, Pa1, Qa1);

                if w1.equals(&w2) == 0xFFFFFFFF {
                    
                } else {
                    Qa.Y.set_neg();
                }

                // TODO: this is again computed using X coordinates in compute_isogeny, perhaps
                // use this result
                let shift_1 = curve_1.double_iter(Pa, self.a as usize);
                let shift_2 = curve_2.double_iter(Pa1, self.a as usize);

                let Pa_shift = curve_1.add(Pa, &shift_1);
                let Qa_shift = curve_1.add(Qa, &shift_1);

                let Pa1_shift = curve_2.add(Pa1, &shift_2);
                let Qa1_shift = curve_2.add(Qa1, &shift_2);

                let shift_1_3 = curve_1.mul_small(&shift_1, 3);
                let Pb1_shift = curve_1.add(&Pb1, &shift_1_3);
                let Qb1_shift = curve_1.add(&Qb1, &shift_1_3);
                let PQb1_shift = curve_1.add(&PQb1, &shift_1_3);

                // TODO: output in some different form
                (Pa_shift, Qa_shift, Pa1_shift, Qa1_shift, Pb1, Qb1, PQb1, Pb1_shift, Qb1_shift, PQb1_shift)
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
