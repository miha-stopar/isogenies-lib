macro_rules! define_litsigamal {
    () => {
        use crate::quaternion::quaternion_algebra::{QuatAlg, QuatAlgEl};
        use crate::quaternion::klpt::represent_integer;
        use crate::quaternion::quaternion_order::standard_maximal_extremal_order;
        use crate::util::{generate_random_range};
        use crate::ec_lit;

        /// LITSiGamal struct
        #[derive(Clone, Debug)]
        pub struct LITSiGamal {
            curve: Curve,
            p: Integer,
            n: u32,
        }

        impl LITSiGamal {
            pub fn new(curve: Curve, p: Integer, n: u32) -> Self {
                Self {
                    curve, p, n
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
                let mut bytes = big_to_bytes(a1);
                let a1Pa = self.curve.mul(&Pa, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(a2);
                let a2Qa = self.curve.mul(&Qa, &bytes, bytes.len() * 8);
                Qa = self.curve.add(&a1Pa, &a2Qa);

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
                bytes = big_to_bytes(c1);
                let c1Pc = self.curve.mul(&Pc, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(c2);
                let c2Qc = self.curve.mul(&Qc, &bytes, bytes.len() * 8);
                Pc = self.curve.add(&c1Pc, &c2Qc);

                // N = (l_a**(2*a) - n**2) * l_b**b
                let tau = l_a.big().pow(power_a * 2) - self.n * self.n;
                let N = tau * l_b.big().pow(power_b);
                let gamma = self.generate_gamma(N);

                // TODO: remove, just for debugging
                let qa = QuatAlg::new(-self.p.clone());
                // let gamma = QuatAlgEl::new("101064425814894250297266458592805333959492115142777713483064160339396739459362670629011189775195355480432667912928029565585343417641447712715314376809237973".big(), "93533154661012604105678420765706116710130443509048775286300216094985327466258160413386761474232055531385270204373695404255540646882449862261480847728263402".big(), "4753400567059888737315986212971719000".big(), "20561456972084008333059843341602645649".big(), 2.big(), qa.clone());
                // let gamma = QuatAlgEl::new("2011185045392199242481466685241103431236289913375605487988529329864805134553257030694326259418924424563378327667027885353428994775042308664556073155839862".big(), "10853501034135357035889336862977087605953185135487653652185963537341394638859475493620188339725879895509904749276277210794053785439623299471028262804950610".big(), "-10310443381405150524300858171174973139".big(), "-5988292950171597878820895074482233260".big(), 2.big(), qa.clone());
                // let gamma = QuatAlgEl::new("8".big(), "4".big(), "6".big(), "2".big(), 2.big(), qa.clone());

                // let gamma = QuatAlgEl::new("2011185045392199242481466685241103431236289913375605487988529329864805134553257030694326259418924424563378327667027885353428994775042308664556073155839862".big(), "10853501034135357035889336862977087605953185135487653652185963537341394638859475493620188339725879895509904749276277210794053785439623299471028262804950610".big(), "-10310443381405150524300858171174973139".big(), "-5988292950171597878820895074482233260".big(), 2.big(), qa.clone());

                let gamma = QuatAlgEl::new("12438122981279082530580154370455515949426546571783775599841409735688020657743566670317349247383366128276797224590897687418693932406293124037751296924501724".big(), "33421844375087071646780030668302789396979173001517474713223550355563704785253646650376906179723303183092207272674743906841643325724028906856886927071952077".big(), "118287257238515240750857735462690960763".big(), "12842514859574082180489774512612831798".big(), 1.big(), qa.clone());
                // let gamma = QuatAlgEl::new("12438122981279082530580154370455515949426546571783775599841409735688020657743566670317349247383366128276797224590897687418693932406293124037751296924501724".big(), "0".big(), "0".big(), "0".big(), 1.big(), qa.clone());
                // let gamma = QuatAlgEl::new("12438122981279082530580154370455515949426546571783775599841409735688020657743566670317349247383366128276797224590897687418693932406293124037751296924501724".big(), "33421844375087071646780030668302789396979173001517474713223550355563704785253646650376906179723303183092207272674743906841643325724028906856886927071952077".big(), "0".big(), "0".big(), 1.big(), qa.clone());

                
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

                let torsion = l_a.big().pow(power_a);
                let (Pa_new, Qa_new, PmQa_new) = apply_endomorphism_on_torsion_group(&self.curve, coord.clone(), imprim.clone(), torsion, mat2_2, mat2_3, mat2_4, &Pa, &Qa);

                let torsion = l_b.big().pow(power_b);
                let (Pb_new, Qb_new, PmQb_new) = apply_endomorphism_on_torsion_group(&self.curve, coord, imprim, torsion, mat3_2, mat3_3, mat3_4, &Pb, &Qb);

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
                println!("");


                println!("Pb_new X: {}", Pb_new.X / Pb_new.Z);
                println!("");
                println!("Pb_new Y: {}", Pb_new.Y / Pb_new.Z);
                println!("");

                println!("Qb_new X: {}", Qb_new.X / Qb_new.Z);
                println!("");
                println!("Qb_new Y: {}", Qb_new.Y / Qb_new.Z);
                println!("");


                let dlog1 = ec_lit::dlog_3(&self.curve, &Pb_new, &Qb_new, 162);
                let dlog2 = ec_lit::dlog_3(&self.curve, &Qb_new, &Pb_new, 162);
                
                println!("dlog1: {}", dlog1);
                println!("");
                println!("dlog2: {}", dlog2);
                println!("");
                println!("");
                println!("");
                println!("");

                /*
                println!("====== mapped =======");
                println!("");
                println!("{}", Pa_new.X / Pa_new.Z);
                println!("");

                println!("");
                println!("{}", Pa_new.Y / Pa_new.Z);
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

                // TODO 
                let n = 36;

                // Precomputed with strategy.py
                let strategy: [usize; 35] = [21, 14, 5, 3, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1];

                // endomorphism gamma is of order (l_a**(2*a) - n**2) * l_b**b
                // let's take the isogeny gamma1 with kernel: ker(gamma) \cap E[l_b**b]

                // TODO: eval_points are [Pa, Qa, R]

                // TODO: remove, just for testing, instead use the kernel of 
                let eval_points = [PaX, QaX]; 
                let (codomain, image_points) = ec_lit::three_isogeny_chain(&self.curve, &PbX, &eval_points, n, &strategy);


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
