macro_rules! define_litsigamal {
    () => {
        use crate::quaternion::quaternion_algebra::{QuatAlg, QuatAlgEl};
        use crate::quaternion::klpt::represent_integer;
        use crate::quaternion::quaternion_order::standard_maximal_extremal_order;
        use crate::util::{generate_random_range};
        use crate::ec_lit;

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

                println!("++++++=================++++++++");
                println!("{}", PaX.X / PaX.Z);
                println!("");
                println!("{}", PaX);
                println!("");
    
                let (Pa, ok1) = self.curve.complete_pointX(&PaX); // TODO: check ok
                assert!(ok1 == 0xFFFFFFFF);

                let (mut Qa, ok2) = self.curve.complete_pointX(&QaX); // TODO: check ok

                let (Pb, ok1) = self.curve.complete_pointX(&PbX); // TODO: check ok
                assert!(ok1 == 0xFFFFFFFF);
                let (Qb, ok2) = self.curve.complete_pointX(&QbX); // TODO: check ok

                let (mut Pc, ok1) = self.curve.complete_pointX(&PcX); // TODO: check ok
                let (Qc, ok2) = self.curve.complete_pointX(&QcX); // TODO: check ok
 
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
                let mut bytes = big_to_bytes(a1);
                let a1Pa = self.curve.mul(&Pa, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(a2);
                let a2Qa = self.curve.mul(&Qa, &bytes, bytes.len() * 8);
                let Qa_rand = self.curve.add(&a1Pa, &a2Qa);
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
                // TODO: use dblmul
                /*
                let mut bytes = big_to_bytes(c1);
                let c1Pc = self.curve.mul(&Pc, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(c2);
                let c2Qc = self.curve.mul(&Qc, &bytes, bytes.len() * 8);
                Pc = self.curve.add(&c1Pc, &c2Qc);
                */

                // N = (l_a**(2*a) - n**2) * l_b**b
                let tau = l_a.big().pow(power_a * 2) - self.n * self.n;
                let N = tau * l_b.big().pow(power_b);
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
                let (Pa_gamma, Qa_gamma, PmQa_gamma) = apply_endomorphism_on_torsion_group(&self.curve, coord.clone(), imprim.clone(), torsion_a, mat2_2, mat2_3, mat2_4, &Pa, &Qa);

                //
                let mut bytes = big_to_bytes(a1);
                let a1Pa = self.curve.mul(&Pa_gamma, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(a2);
                let a2Qa = self.curve.mul(&Qa_gamma, &bytes, bytes.len() * 8);
                let Qa_rand_gamma = self.curve.add(&a1Pa, &a2Qa);

                let torsion_b = l_b.big().pow(power_b);
                let (Pb_gamma, Qb_gamma, PmQb_gamma) = apply_endomorphism_on_torsion_group(&self.curve, coord.clone(), imprim.clone(), torsion_b.clone(), mat3_2, mat3_3, mat3_4, &Pb, &Qb);

                let torsion_c = l_c.big().pow(power_c);
                let (Pc_gamma, Qc_gamma, PmQc_gamma) = apply_endomorphism_on_torsion_group(&self.curve, coord, imprim, torsion_c, mat5_2, mat5_3, mat5_4, &Pc, &Qc);

                let mut bytes = big_to_bytes(c1);
                let c1Pc = self.curve.mul(&Pc_gamma, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(c2);
                let c2Qc = self.curve.mul(&Qc_gamma, &bytes, bytes.len() * 8);
                let R = self.curve.add(&c1Pc, &c2Qc);

                let Rx = PointX::new_xz(&R.X, &R.Z);


                println!("Pb X: {}", Pb.X / Pb.Z);
                println!("");
                println!("Pb Y: {}", Pb.Y / Pb.Z);
                println!("");
                println!("");

                println!("Qb X: {}", Qb.X / Qb.Z);
                println!("");
                println!("Qb Y: {}", Qb.Y / Qb.Z);
                println!("");
                println!("");

                println!("");
                println!("");

                println!("Pa X: {}", Pa.X / Pa.Z);
                println!("");
                println!("Pa Y: {}", Pa.Y / Pa.Z);
                println!("");
                println!("");

                println!("Qa X: {}", Qa.X / Qa.Z);
                println!("");
                println!("Qa Y: {}", Qa.Y / Qa.Z);
                println!("");
                println!("");

                println!("Pc X: {}", Pc.X / Pc.Z);
                println!("");
                println!("Pc Y: {}", Pc.Y / Pc.Z);
                println!("");
                println!("");

                println!("Qc X: {}", Qc.X / Qc.Z);
                println!("");
                println!("Qc Y: {}", Qc.Y / Qc.Z);
                println!("");
                println!("");

                println!("========================");


                println!("Pb_gamma X: {}", Pb_gamma.X / Pb_gamma.Z);
                println!("");
                println!("Pb_gamma Y: {}", Pb_gamma.Y / Pb_gamma.Z);
                println!("");

                println!("Qb_gamma X: {}", Qb_gamma.X / Qb_gamma.Z);
                println!("");
                println!("Qb_gamma Y: {}", Qb_gamma.Y / Qb_gamma.Z);
                println!("");

                println!("");

                println!("Pa_gamma X: {}", Pa_gamma.X / Pa_gamma.Z);
                println!("");
                println!("Pa_gamma Y: {}", Pa_gamma.Y / Pa_gamma.Z);
                println!("");

                println!("Qa_rand_gamma X: {}", Qa_rand_gamma.X / Qa_rand_gamma.Z);
                println!("");
                println!("Qa_rand_gamma Y: {}", Qa_rand_gamma.Y / Qa_rand_gamma.Z);
                println!("");

                println!("");

                println!("R X: {}", R.X / R.Z);
                println!("");
                println!("R Y: {}", R.Y / R.Z);
                println!("");
 
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

                bytes = big_to_bytes(dlog2); // TODO
                let dlog_Qb = self.curve.mul(&Qb, &bytes, bytes.len() * 8);
                let kernel1 = self.curve.sub(&Pb, &dlog_Qb); // TODO: PointX directly

                let kernel1x = PointX::new_xz(&kernel1.X, &kernel1.Z);

                println!("kernel1 X: {}", kernel1.X / kernel1.Z);
                println!("");
                println!("kernel Y: {}", kernel1.Y / kernel1.Z);
                println!("");

                // TODO 
                let n = 36;

                // Precomputed with strategy.py
                let strategy: [usize; 35] = [21, 14, 5, 3, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1];

                // endomorphism gamma is of order (l_a**(2*a) - n**2) * l_b**b
                // let's take the isogeny gamma1 with kernel: ker(gamma) \cap E[l_b**b]

                // TODO: eval_points are [Pa, Qa, R]

                // TODO: remove, just for testing, instead use the kernel of 
                let eval_points = [PaX, QaX, Rx]; 
                let (codomain, image_points) = ec_lit::three_isogeny_chain(&self.curve, &kernel1x, &eval_points, n, &strategy);

                /*
                40462640216449207672168992988623511775336924066853855274546745222615865684900642222898108382455842618365117964778639196177347207085357604725664385427796579295592436532010871949667738164957851517164788313374798288498756871920103663953*sqrtminus
                + 970541812902597247581586118962672018735900274268797204095167395803112287493471491598698090128841200813831380874652397438640565161195302840313904161390870177508026650386592721124672238490254317774184133590356799364984120716582941706095

                1249722545006266441120006262083514632766815889237780158930324325671112785250730980429928525004234781246451849490611409020519930843445073584072794696443427504458825306819336321611899670903300187189923829452476229237796629934363670095154*sqrtminus
                + 220835060888732177714055873396130338211426525716735667126531717994944550825260327261201226324233211400965520569038349098555557429025991989295434705054844308790847237100964246142971652968315021176800459089181506969495264032362279571665
                */

                /*
                let Px = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("970541812902597247581586118962672018735900274268797204095167395803112287493471491598698090128841200813831380874652397438640565161195302840313904161390870177508026650386592721124672238490254317774184133590356799364984120716582941706095")),
                    &Fp::decode_reduce(&bytes_from_str("40462640216449207672168992988623511775336924066853855274546745222615865684900642222898108382455842618365117964778639196177347207085357604725664385427796579295592436532010871949667738164957851517164788313374798288498756871920103663953")),
                );
                let Pz = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("1249722545006266441120006262083514632766815889237780158930324325671112785250730980429928525004234781246451849490611409020519930843445073584072794696443427504458825306819336321611899670903300187189923829452476229237796629934363670095154")),
                    &Fp::decode_reduce(&bytes_from_str("220835060888732177714055873396130338211426525716735667126531717994944550825260327261201226324233211400965520569038349098555557429025991989295434705054844308790847237100964246142971652968315021176800459089181506969495264032362279571665")),
                );
               
                let FF = PointX::new_xz(&Px, &Pz);
                let (A24_plus, A24_minus, K1, K2) = three_isogeny_codomain(&FF);
                let A = &(&A24_plus + A24_minus).mul2() / &(&A24_plus - &A24_minus);

                let codomain = Curve::new(&A);

                println!("");
                println!("");
                println!("");
                println!("codomain: {}", codomain);
                println!("");
                println!("A24_plus: {}", A24_plus);
                println!("");

                println!("A24_minus: {}", A24_minus);
                println!("");
                */


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

                let Pbr = generate_random_fp(&self.curve, torsion_b.clone(), self.scalar_without_b.clone());
                let Qbr = generate_random_fq(&self.curve, torsion_b, self.scalar_without_b.clone());

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
