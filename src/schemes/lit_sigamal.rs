macro_rules! define_litsigamal {
    () => {
        use crate::quaternion::quaternion_algebra::{QuatAlg, QuatAlgEl};
        use crate::quaternion::klpt::represent_integer;
        use crate::quaternion::quaternion_order::standard_maximal_extremal_order;
        use crate::util::{generate_random_range};
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
        }

        impl PubKey {
            pub fn new(
                points: PubKeyPoints,
                points1: PubKeyPoints,
                power_a: u32,
                power_b: u32,
                power_c: u32,
            ) -> Self {
                Self {
                    points,
                    points1,
                    power_a,
                    power_b,
                    power_c,
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

        pub fn get_params(lam: u32) -> (Integer, u32, u32, u32, u32) {
            let p;
            let a = lam * 3;
            let b: u32;
            let c: u32;
            let f: u32;

            if lam == 128 {
                p = "1290217975993796939363993419446162388979006021159541007293712082644700121088673466685157498316158528176855539315411759315356741765308895915108991692829754882889263058278152142847999999999999999999999999999999999999999999999999999999999".big();
                b = 162;
                c = 56;
                f = 30;
            } else if lam == 192 {
                p = "105243372827141588334896706952131446393985610876498590866946033436085731820481224209718550756818498226506584579388540458536674694936506474952416076610303430661479448818071050554239639344046356810649842429481145316664146061991592319538503400256168957692064876729257164799999999999999999999999999999999999999999999999999999999999999999999999999999999999".big();
                b = 243;
                c = 83;
                f = 118;
            } else if lam == 256 {
                p = "40321823197322392735588260263395257580747909237527626646846604072674523904511201327512435027773599857137068322905929186031078022612322994889152431338934029723532461084420532453809231794766530701129753010640413009244929390629124476760178197654116840881679679002356836870380358444873109038766263029534615051169507870006410875301114902805114753409739451793407999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999".big();
                b = 324;
                c = 111;
                f = 436;
            } else {
                panic!("lam should be 128 or 192 or 256");
            }

            (p, a, b, c, f)
        }

        /// LITSiGamal struct
        #[derive(Clone, Debug)]
        pub struct LITSiGamal<'a> {
            param: u32,
            curve: Curve,
            p: Integer,
            l_a: u32,
            l_b: u32,
            l_c: u32,
            a: u32,
            b: u32,
            c: u32,
            scalar_without_b: Integer, // p + 1 = scalar_without_b * l_b**b
            n: u32,
            strategy: &'a [usize],
            strategy_2: &'a [usize],
        }

        impl<'a> LITSiGamal<'a> {
            pub fn new(param: u32) -> Self {
                let (p, a, b, c, f) = get_params(param);
                let n = 3;
                let A = Fq::ZERO;
                let curve = Curve::new(&A);
                let l_a = 2;
                let l_b = 3;
                let l_c = 5;

                let scalar_without_b = l_a.big().pow(a + 2) * l_c.big().pow(c) * f;
                
                let strategy: &[usize];
                let strategy_2: &[usize];
                if param == 128 {
                    // From Sage LITSiGamal, TODO: measure it!
                    strategy  = &[65, 37, 23, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 9, 5, 4, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1];

                    strategy_2 = &[154, 93, 55, 33, 20, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 22, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 38, 22, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 65, 37, 21, 12, 7, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1];
                } else if param == 192 {
                    // From Sage LITSiGamal, TODO: measure it!
                    strategy = &[92, 65, 37, 21, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 37, 21, 13, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1];

                    strategy_2 = &[229, 146, 86, 49, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 37, 21, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 60, 37, 21, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 23, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 9, 5, 4, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 93, 55, 33, 20, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 22, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 38, 22, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1];
                } else {
                    // From Sage LITSiGamal, TODO: measure it!
                    // strategy = &[124, 86, 49, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 37, 21, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 49, 28, 19, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 7, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1];
                    // strategy = &[110, 89, 55, 34, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 34, 21, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1];
                    strategy = &[111, 89, 55, 34, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 34, 22, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1];

                    // strategy_2 = &[309, 194, 113, 65, 37, 21, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 48, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 20, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 81, 48, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 20, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 33, 20, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 115, 81, 48, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 20, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 33, 20, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 49, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1];
                    strategy_2 = &[233, 177, 144, 89, 55, 34, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 55, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 55, 34, 33, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 12, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 89, 55, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1];
                }

                Self {
                    param, curve, p, l_a, l_b, l_c, a, b, c, scalar_without_b, n, strategy, strategy_2
                }
            }

            pub fn generate_pub_key(&self) -> (PubKey, (Integer, Integer, Integer)) {
                let start = Instant::now();

                let fileName = "src/schemes/precomputed.json";
                let s = &format!("lit-sigamal-{}", self.param);
                let (PaX, QaX, PmQaX, mat2_2, mat2_3, mat2_4, power_a) = load_torsion_info(fileName, s, 2);
                let (PbX, QbX, PmQbX, mat3_2, mat3_3, mat3_4, power_b) = load_torsion_info(fileName, s, 3);
                let (PcX, QcX, PmQcX, mat5_2, mat5_3, mat5_4, power_c) = load_torsion_info(fileName, s, 5);

                // Note: power_a = a + 2
                let power_a = power_a as u32;
                let power_b = power_b as u32;
                let power_c = power_c as u32;

                let mut curve = self.curve.clone();
    
                // TODO: prepare completions, also for Y!! Otherwise, the Y completion
                // in Sage and Rust might be different.
                let (mut Pa, ok) = curve.complete_pointX(&PaX);
                assert!(ok == 0xFFFFFFFF);
                let (mut Qa, ok) = curve.complete_pointX(&QaX);
                assert!(ok == 0xFFFFFFFF);
                let (PmQa, ok) = curve.complete_pointX(&PmQaX);
                assert!(ok == 0xFFFFFFFF);

                let (Pb, ok) = curve.complete_pointX(&PbX);
                assert!(ok == 0xFFFFFFFF);
                let (Qb, ok) = curve.complete_pointX(&QbX);
                assert!(ok == 0xFFFFFFFF);
                let (PmQb, ok) = curve.complete_pointX(&PmQbX);
                assert!(ok == 0xFFFFFFFF);

                let (Pc, ok) = curve.complete_pointX(&PcX);
                assert!(ok == 0xFFFFFFFF);
                let (Qc, ok) = curve.complete_pointX(&QcX);
                assert!(ok == 0xFFFFFFFF);
                let (PmQc, ok) = curve.complete_pointX(&PmQcX);
                assert!(ok == 0xFFFFFFFF);

                // Qa.Y.set_neg(); // TODO 1
                
                /*
                println!("");
                println!("====================");
                println!("");
                println!("Pa: {}", Pa.X / Pa.Z);
                println!("");
                println!("Pa: {}", Pa.Y / Pa.Z);
                println!("");

                println!("");
                println!("Qa: {}", Qa.X / Qa.Z);
                println!("");
                println!("Qa: {}", Qa.Y / Qa.Z);
                println!("");

                println!("");
                println!("Pb: {}", Pb.X / Pb.Z);
                println!("");
                println!("Pb: {}", Pb.Y / Pb.Z);
                println!("");

                println!("");
                println!("Qb: {}", Qb.X / Qb.Z);
                println!("");
                println!("Qb: {}", Qb.Y / Qb.Z);
                println!("");

                println!("");
                println!("Pc: {}", Pc.X / Pc.Z);
                println!("");
                println!("Pc: {}", Pc.Y / Pc.Z);
                println!("");

                println!("");
                println!("Qc: {}", Qc.X / Qc.Z);
                println!("");
                println!("Qc: {}", Qc.Y / Qc.Z);
                println!("");
                */
                 
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
                let Pc_rand = self.curve.add(&c1Pc, &c2Qc);
 
                let Rx = PointX::new_xz(&Pc_rand.X, &Pc_rand.Z);

                // N = (l_a**(2*a) - n**2) * l_b**b
                let tau = l_a.big().pow(self.a * 2) - self.n * self.n;
                let N = tau.clone() * l_b.big().pow(power_b);
                let mut gamma = self.generate_gamma(N.clone());

                // TODO: remove, just for debugging
                let qa = QuatAlg::new(-self.p.clone());
                // let gamma = QuatAlgEl::new("101064425814894250297266458592805333959492115142777713483064160339396739459362670629011189775195355480432667912928029565585343417641447712715314376809237973".big(), "93533154661012604105678420765706116710130443509048775286300216094985327466258160413386761474232055531385270204373695404255540646882449862261480847728263402".big(), "4753400567059888737315986212971719000".big(), "20561456972084008333059843341602645649".big(), 2.big(), qa.clone());
                // let gamma = QuatAlgEl::new("2011185045392199242481466685241103431236289913375605487988529329864805134553257030694326259418924424563378327667027885353428994775042308664556073155839862".big(), "10853501034135357035889336862977087605953185135487653652185963537341394638859475493620188339725879895509904749276277210794053785439623299471028262804950610".big(), "-10310443381405150524300858171174973139".big(), "-5988292950171597878820895074482233260".big(), 2.big(), qa.clone());
                // let gamma = QuatAlgEl::new("8".big(), "4".big(), "6".big(), "2".big(), 2.big(), qa.clone());

                // debugging:
                if self.param == 128 {
                    // gamma = QuatAlgEl::new("8425099297649819352166379388758583438269420404259906871066684285709655835562922681747598532398853172000034459777928781756277965902285271705188053172406860".big(), "10672475422448534529312519343800163927437511879170590252915490695777092409588869040995562812407263134190749362059361293778411214787787911660929131931909741".big(), "7138531157822206556546664362792862977".big(), "6507363534081204197854356661058574173".big(), 1.big(), qa.clone());
                } else if self.param == 192 {
                    // gamma = QuatAlgEl::new("662600986569554607086186139717380593600621700495654784370754813176584196881602237433761052273707611322728528748274174182542784592898311629798340420256581733000879541025286009940162273590609577777394911993777190222271995532700037707".big(), "1799354100419486712391339580362819447142359123820622055360650764778398997554824673132730837236860088729249818881893370453300960091771964990380640474206579364334639343768768257294433631204991561360624607148084574683793231827142040775".big(), "-125472024283789366450893441561831071654100276651462478219".big(), "268933293804616335725073763987973337907728658776001142".big(), 1.big(), qa.clone());
                } else {
                    // gamma = QuatAlgEl::new("138601276885266359348924423168250315165478982100976906860968440028802590202782345587413951105705556006524207959973750964087196404850926049290216371979961522797259549414688362603367850681343177750244976932214600226221154506142599743514660210665164334544479617325876989846049244977908157527833561789481301570878".big(), "203108580675992905414491090775344429818391560201936869394145342661635849317307986324805282886129261812822283285677764427289921365829467852342774407932637374237713787715658816530668155912380527928217313426997520452727693996824057379361319175112782539718488257668282926209879601449187438917780234308575133684498".big(), "-194209084683335020967895199072782968208342558998606564091220864772136318540".big(), "879697261687923968662903575211307786154010371542164913611531882979847514199".big(), 1.big(), qa.clone());
                }

                /*
                println!("");
                println!("");
                println!("gamma: ");
                println!("");
                println!("{:?}", gamma);
                println!("");
                println!("");
                */
                
                let order = standard_maximal_extremal_order().order;
                let (coord, imprim) = gamma.factor_in_order(order.lattice.clone());

                // println!("1 (before apply endomorphism): {:?}", start.elapsed());
                // let second_part = Instant::now();

                let torsion_a = l_a.big().pow(self.a + 2);
                // let (Pa_gamma, _, _) = apply_endomorphism_on_torsion_group(&self.curve, coord.clone(), imprim.clone(), torsion_a, mat2_2, mat2_3, mat2_4, &Pa, &Qa);
                
                let (mut Pa_gamma, Qa_gamma, R_gamma) = apply_endomorphism_on_torsion_group(&curve, coord.clone(), imprim.clone(), torsion_a, mat2_2, mat2_3, mat2_4, &Pa, &Qa, &PmQa);

                // TODO
                // let mut Qa1_to_be_mapped = PointX::new_xz(&Qa_gamma.X, &Qa_gamma.Z);
                // let mut Ra1_to_be_mapped = PointX::new_xz(&R_gamma.X, &R_gamma.Z);

                let bytes1 = big_to_bytes(a1);
                let a1Pa_gamma = curve.mul(&Pa_gamma, &bytes1, bytes1.len() * 8);
                let bytes2 = big_to_bytes(a2);
                let a2Qa_gamma = curve.mul(&Qa_gamma, &bytes2, bytes2.len() * 8);
                let mut Qa_rand_gamma = curve.add(&a1Pa_gamma, &a2Qa_gamma);

                let torsion_c = l_c.big().pow(power_c);
                let (Pc_gamma, Qc_gamma, PmQc_gamma) = apply_endomorphism_on_torsion_group(&curve, coord.clone(), imprim.clone(), torsion_c, mat5_2, mat5_3, mat5_4, &Pc, &Qc, &PmQc);

                let mut bytes = big_to_bytes(c1);
                let c1Pc = self.curve.mul(&Pc_gamma, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(c2);
                let c2Qc = self.curve.mul(&Qc_gamma, &bytes, bytes.len() * 8);
                let R = self.curve.add(&c1Pc, &c2Qc);
                let mut R_gamma = PointX::new_xz(&R.X, &R.Z);

                /*
                println!("");
                println!("+++++++++++++");
                println!("");
                println!("Pa1: {}", Pa_gamma.X / Pa_gamma.Z);
                println!("");
                println!("Pa1: {}", Pa_gamma.Y / Pa_gamma.Z);
                println!("");
                println!("Qa1: {}", Qa_rand_gamma.X / Qa_rand_gamma.Z);
                println!("");
                println!("Qa1: {}", Qa_rand_gamma.Y / Qa_rand_gamma.Z);
                println!("");
                println!("");
                println!("R1: {}", R_gamma.X / R_gamma.Z);
                println!("");
                println!("R1: {}", R.Y / R.Z);
                println!("");
                println!("");
                */

                // Set [2^a]Pa_gamma = (1,*)
                let t = self.l_a.big().pow(self.a); // TODO: define once
                let bytes = big_to_bytes(t);
                let Pa_gamma_4 = curve.mul(&Pa_gamma, &bytes, bytes.len() * 8);
                if Pa_gamma_4.X.equals(&Pa_gamma_4.Z) != 0xFFFFFFFF {
                    Pa_gamma.X.set_neg();
                    Qa_rand_gamma.X.set_neg();

                    // R.X.set_neg();
                    R_gamma.X.set_neg();
                    
                    let mut tmp = PointX::new_xz(&Pa_gamma.X, &Pa_gamma.Z);
                    (Pa_gamma, _) = curve.complete_pointX(&tmp);

                    tmp = PointX::new_xz(&Qa_rand_gamma.X, &Qa_rand_gamma.Z);
                    (Qa_rand_gamma, _) = curve.complete_pointX(&tmp);

                    // tmp = PointX::new_xz(&R_gamma.X, &R_gamma.Z);
                    // (R, _) = curve.complete_pointX(&tmp);

                    // println!("HHHHHMMMMMMMMMMMMMMMMMM");
                }
                let mut Pa1_to_be_mapped = PointX::new_xz(&Pa_gamma.X, &Pa_gamma.Z);

                /*
                println!("");
                println!("Pa1: {}", Pa_gamma.X / Pa_gamma.Z);
                println!("");
                println!("Qa1: {}", Qa_rand_gamma.X / Qa_rand_gamma.Z);
                println!("");
                println!("R1: {}", R_gamma.X / R_gamma.Z);
                println!("");
                println!("");
                */

                let mut Ra1_to_be_mapped = PointX::new_xz(&R_gamma.X, &R_gamma.Z);

                /*
                println!("Qa before rand");
                println!("Qa: {}", Qa.X / Qa.Z);
                println!("");
                println!("Qa: {}", Qa.Y / Qa.Z);
                println!("");
                */

                let a1Pa = curve.mul(&Pa, &bytes1, bytes1.len() * 8);
                let a2Qa = curve.mul(&Qa, &bytes2, bytes2.len() * 8);
                let Qa_rand = curve.add(&a1Pa, &a2Qa);
                let Qa_rand_X = PointX::new_xz(&Qa_rand.X, &Qa_rand.Z);

                /*
                println!("Qa after rand");
                println!("Qa: {}", Qa_rand.X / Qa_rand.Z);
                println!("");
                println!("Qa: {}", Qa_rand.Y / Qa_rand.Z);
                println!("");

                println!("============ 111 ==============");
                println!("Pa: {}", Pa.X / Pa.Z);
                println!("");
                println!("Pa: {}", Pa.Y / Pa.Z);
                println!("");
                println!("Qa: {}", Qa_rand.X / Qa_rand.Z);
                println!("");
                println!("Qa: {}", Qa_rand.Y / Qa_rand.Z);
                println!("");

                println!("");
                println!("Pb: {}", Pb.Y / Pb.Z);
                println!("");
                println!("Qb: {}", Qb.Y / Qb.Z);
                println!("");

                println!("");

                println!("");
                println!("");
                println!("Pc: {}", Pc_rand.X / Pc_rand.Z);
                println!("");
                println!("Pc: {}", Pc_rand.Y / Pc_rand.Z);
                println!("");
                println!("Qc: {}", Qc.Y / Qc.Z);
                println!("");
                */

                let mut Qa1_to_be_mapped = PointX::new_xz(&Qa_rand_gamma.X, &Qa_rand_gamma.Z);
            
                let torsion_b = l_b.big().pow(power_b);
                let torsion_b_minus = l_b.big().pow(power_b - 1);
                let (Pb_gamma, Qb_gamma, PmQb_gamma) = apply_endomorphism_on_torsion_group(&curve, coord.clone(), imprim.clone(), torsion_b.clone(), mat3_2, mat3_3, mat3_4, &Pb, &Qb, &PmQb);

                /* 
                println!("???????????????????");
                println!("");
                println!("Qb_gamma: {}", Qb_gamma.X / Qb_gamma.Z);
                println!("");
                println!("Qb_gamma: {}", Qb_gamma.Y / Qb_gamma.Z);
                println!("");

                println!("");
                println!("Pa_gamma: {}", Pa_gamma.X / Pa_gamma.Z);
                println!("");
                println!("Pa_gamma: {}", Pa_gamma.Y / Pa_gamma.Z);
                println!("");
                */

                // println!("2 (after apply endomorphism): {:?}", second_part.elapsed());
                // let third_part = Instant::now();

                // TODO: remove one dlog
                let dlog1 = dlog_3(&curve, &Pb_gamma, &Qb_gamma, self.b.try_into().unwrap());
                // let dlog2 = dlog_3(&curve, &Qb_gamma, &Pb_gamma, self.b.try_into().unwrap());
                // println!("3 (after dlog): {:?}", third_part.elapsed());
                // let fourth_part = Instant::now();

                /*
                let kernel1;
                if dlog1 != 0.big() {
                    println!("pubkey dlog 1");
                    let bytes = big_to_bytes(dlog1);
                    let dlog_Pb = curve.mul(&Pb, &bytes, bytes.len() * 8);
                    kernel1 = curve.sub(&Qb, &donsog_Pb);
                } else {
                    println!("pubkey dlog 2");
                    let bytes = big_to_bytes(dlog2);
                    let dlog_Qb = curve.mul(&Qb, &bytes, bytes.len() * 8);
                    kernel1 = curve.sub(&Pb, &dlog_Qb);
                }
                */
                let bytes = big_to_bytes(dlog1);
                let dlog_Pb = curve.mul(&Pb, &bytes, bytes.len() * 8);
                let kernel1 = curve.sub(&Qb, &dlog_Pb);

                let kernel1x = PointX::new_xz(&kernel1.X, &kernel1.Z);

                let n = self.strategy.len() + 1;

                // endomorphism gamma is of order (l_a**(2*a) - n**2) * l_b**b
                // let's take the isogeny gamma1 with kernel: ker(gamma) \cap E[l_b**b]

                /*
                println!("???????????????????");
                println!("");
                println!("Pa: {}", PaX.X / PaX.Z);
                println!("");
                println!("Qa rand: {}", Qa_rand_X.X / Qa_rand_X.Z);
                println!("");
                println!("R: {}", Rx.X / Rx.Z);
                println!("");
                println!("");
                */
 

                // TODO: eval_points are [Pa, Qa, R]
 
                // TODO: remove, just for testing, instead use the kernel of 
                // TODO: fix QaX above
                let eval_points = [PaX, Qa_rand_X, Rx]; 
                let (mut codomain, mut image_points) = three_isogeny_chain(&curve, &kernel1x, eval_points.to_vec(), n, self.strategy);
                let (mut Pa_isog3X, mut Qa_rand_isog3X, R_isog3X) = (image_points[0], image_points[1], image_points[2]);

                // TODO
                let mut Pa_to_be_mapped = Pa_isog3X.clone(); // to be mapped by 3-isogeny chain
                let mut Qa_to_be_mapped = Qa_rand_isog3X.clone();
                let mut Ra_to_be_mapped = R_isog3X.clone();

                let (mut Pa, _) = codomain.complete_pointX(&Pa_isog3X);
                let (mut Qa_rand, _) = codomain.complete_pointX(&Qa_rand_isog3X);

                // TODO 1: remove
                // Pa.Y.set_neg();

                /*
                println!("????????????/ after 3-isogeny");
                println!("");
                println!("Pa: {}", Pa.X / Pa.Z);
                println!("");
                println!("Pa: {}", Pa.Y / Pa.Z);
                println!("");
                println!("Qa: {}", Qa_rand.X / Qa_rand.Z);
                println!("");
                println!("Qa: {}", Qa_rand.Y / Qa_rand.Z);
                println!("");
                */

                let (Pa_shift, Qa_shift, Pa1_shift, Qa1_shift, Pb, Qb, PQb, Pb_shift, Qb_shift, PQb_shift) =
                    self.get_PQb_and_shift(&codomain, &curve, &Pa, &mut Qa_rand, &Pa_gamma, &Qa_rand_gamma, torsion_b_minus.clone(), tau.clone(), 1);

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

                // println!("4 (after 3-isogeny): {:?}", fourth_part.elapsed());
                // let fifth_part = Instant::now();

                let b_points = [Pb_to_be_mapped, Qb_to_be_mapped, PQb_to_be_mapped];
                let b_shift_points = [Pb_shiftX, Qb_shiftX, PQb_shiftX];

                /*
                println!("??????????????????? 1111");
                println!("");
                println!("Pa: {}", Pa_isog3X.X / Pa_isog3X.Z);
                println!("");
                println!("Qa rand: {}", Qa_rand_isog3X.X / Qa_rand_isog3X.Z);
                println!("");
                println!("R: {}", R_isog3X.X / R_isog3X.Z);
                println!("");
                println!("");

                println!("");
                println!("Pa1: {}", Pa_gamma.X / Pa_gamma.Z);
                println!("");
                println!("Qa1: {}", Qa_rand_gamma.X / Qa_rand_gamma.Z);
                println!("");
                println!("R1: {}", R_gamma.X / R_gamma.Z);
                println!("");
                println!("");


                println!("");
                println!("?????????????? 11111 2222222222");
                println!("");
                println!("Pb_to_be_mapped: {}", Pb_to_be_mapped.X / Pb_to_be_mapped.Z);
                println!("");
                println!("Qb_to_be_mapped: {}", Qb_to_be_mapped.X / Qb_to_be_mapped.Z);
                println!("");
                println!("PQb_to_be_mapped: {}", PQb_to_be_mapped.X / PQb_to_be_mapped.Z);
                println!("");

                println!("");
                println!("");
                println!("Pb_to_be_mapped shift: {}", Pb_shiftX.X / Pb_shiftX.Z);
                println!("");
                println!("Qb_to_be_mapped: {}", Qb_shiftX.X / Qb_shiftX.Z);
                println!("");
                println!("PQb_to_be_mapped: {}", PQb_shiftX.X / PQb_shiftX.Z);
                println!("");
                */

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
                    b_points.to_vec(),
                    b_shift_points.to_vec(),
                    self.a as usize,
                    self.strategy_2,
                );

                let mut Pb1_to_be_mapped = points[0];
                let mut Qb1_to_be_mapped = points[1];
                let mut PQb1_to_be_mapped = points[2];

                /*
                println!("");
                println!("?????????????? 2222222222");
                println!("");
                println!("Pb1_to_be_mapped: {}", Pb1_to_be_mapped.X / Pb1_to_be_mapped.Z);
                println!("");
                println!("Qb1_to_be_mapped: {}", Qb1_to_be_mapped.X / Qb1_to_be_mapped.Z);
                println!("");
                println!("PQb1_to_be_mapped: {}", PQb1_to_be_mapped.X / PQb1_to_be_mapped.Z);
                println!("");
                */

                let f_mul = l_b.big().pow(power_b - 1);
                let mut backtracking_check = PointX::INFINITY;

                // println!("5 (after (2,2)-isogeny): {:?}", fifth_part.elapsed());
                // let sixth_part = Instant::now();

                for _ in 0..6 {
                    let mut no_backtracking = false;
                    let mut s = 0.big();
                    let mut kernel = PointX::INFINITY;

                    while !no_backtracking {
                        s = generate_random_range(0.big(), l_b.big().pow(power_b - 1)) * l_b + 
                            generate_random_range(1.big(), 3.big());

                        kernel = codomain.ladder_3pt(&Pb_to_be_mapped, &Qb_to_be_mapped, &PQb_to_be_mapped, s.clone());

                        let mut check = kernel.clone();
                        codomain.xmul(&mut check, f_mul.clone());

                        if (check.X * backtracking_check.Z).equals(&(backtracking_check.X * check.Z)) != 0xFFFFFFFF || backtracking_check.equals(&PointX::INFINITY) == 0xFFFFFFFF {
                            no_backtracking = true;
                        }
                    }

                    let kernel1 = curve.ladder_3pt(&Pb1_to_be_mapped, &Qb1_to_be_mapped, &PQb1_to_be_mapped, s);

                    /*
                    println!("");
                    println!("kernel:");
                    println!("{}", kernel.X / kernel.Z);
                    println!("");
                    println!("");

                    println!("");
                    println!("kernel1:");
                    println!("{}", kernel1.X / kernel1.Z);
                    println!("");
                    println!("");

                    println!("");
                    println!("---------- before 2 x 3-isogeny ----------------");
                    println!("");
                    println!("Pa_to_be_mapped: {}", Pa_to_be_mapped.X / Pa_to_be_mapped.Z);
                    println!("");
                    println!("Qa_to_be_mapped: {}", Qa_to_be_mapped.X / Qa_to_be_mapped.Z);
                    println!("");
                    println!("PQa_to_be_mapped: {}", Ra_to_be_mapped.X / Ra_to_be_mapped.Z);
                    println!("");

                    println!("Qb_to_be_mapped: {}", Qb_to_be_mapped.X / Qb_to_be_mapped.Z);
                    println!("");
                    */

                    let eval_points = [Pa_to_be_mapped, Qa_to_be_mapped, Ra_to_be_mapped, Qb_to_be_mapped];
                    (codomain, image_points) = three_isogeny_chain(&codomain, &kernel, eval_points.to_vec(), n, &self.strategy);
                    (Pa_to_be_mapped, Qa_to_be_mapped, Ra_to_be_mapped, backtracking_check) = (image_points[0], image_points[1], image_points[2], image_points[3]);

                    codomain.xmul(&mut backtracking_check, f_mul.clone());

                    /*
                    println!("after first 3-isogeny");
                    println!("");
                    println!("Pa_to_be_mapped: {}", Pa_to_be_mapped.X / Pa_to_be_mapped.Z);
                    println!("");
                    println!("Qa_to_be_mapped: {}", Qa_to_be_mapped.X / Qa_to_be_mapped.Z);
                    println!("");
                    println!("PQa_to_be_mapped: {}", Ra_to_be_mapped.X / Ra_to_be_mapped.Z);
                    println!("");

                    println!("-------------");
                    println!("");
                    println!("Pa1_to_be_mapped: {}", Pa1_to_be_mapped.X / Pa1_to_be_mapped.Z);
                    println!("");
                    println!("Qa1_to_be_mapped: {}", Qa1_to_be_mapped.X / Qa1_to_be_mapped.Z);
                    println!("");
                    println!("PQa1_to_be_mapped: {}", Ra1_to_be_mapped.X / Ra1_to_be_mapped.Z);
                    println!("");
                    */

                    let eval_points = [Pa1_to_be_mapped, Qa1_to_be_mapped, Ra1_to_be_mapped]; 
                    (curve, image_points) = three_isogeny_chain(&curve, &kernel1, eval_points.to_vec(), n, &self.strategy);
                    (Pa1_to_be_mapped, Qa1_to_be_mapped, Ra1_to_be_mapped) = (image_points[0], image_points[1], image_points[2]);

                    /*
                    println!("after second 3-isogeny");
                    println!("");
                    println!("Pa1_to_be_mapped: {}", Pa1_to_be_mapped.X / Pa1_to_be_mapped.Z);
                    println!("");
                    println!("Qa1_to_be_mapped: {}", Qa1_to_be_mapped.X / Qa1_to_be_mapped.Z);
                    println!("");
                    println!("PQa1_to_be_mapped: {}", Ra1_to_be_mapped.X / Ra1_to_be_mapped.Z);
                    println!("");
                    */

                    let (Pa, _) = codomain.complete_pointX(&Pa_to_be_mapped);
                    let (mut Qa, _) = codomain.complete_pointX(&Qa_to_be_mapped);
                    let (Pa1, _) = curve.complete_pointX(&Pa1_to_be_mapped);
                    let (Qa1, _) = curve.complete_pointX(&Qa1_to_be_mapped);

                    let (Pa_shift, Qa_shift, Pa1_shift, Qa1_shift, Pb, Qb, PQb, Pb_shift, Qb_shift, PQb_shift) =
                        self.get_PQb_and_shift(&codomain, &curve, &Pa, &mut Qa, &Pa1, &Qa1, torsion_b_minus.clone(), tau.clone(), 2);

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

                    let b_points = [Pb_to_be_mapped, Qb_to_be_mapped, PQb_to_be_mapped];
                    let b_shift_points = [Pb_shiftX, Qb_shiftX, PQb_shiftX];
                
                    /*
                    println!("Pb_to_be_mapped: {}", Pb_to_be_mapped.X / Pb_to_be_mapped.Z);
                    println!("");
                    println!("Qb_to_be_mapped: {}", Qb_to_be_mapped.X / Qb_to_be_mapped.Z);
                    println!("");
                    println!("PQb_to_be_mapped: {}", PQb_to_be_mapped.X / PQb_to_be_mapped.Z);
                    println!("");

                    println!("");
                    println!("Pb1_to_be_mapped: {}", Pb1_to_be_mapped.X / Pb1_to_be_mapped.Z);
                    println!("");
                    println!("Qb1_to_be_mapped: {}", Qb1_to_be_mapped.X / Qb1_to_be_mapped.Z);
                    println!("");
                    println!("PQb1_to_be_mapped: {}", PQb1_to_be_mapped.X / PQb1_to_be_mapped.Z);
                    println!("");

                    println!("");
                    println!("Pb shift: {}", Pb_shiftX.X / Pb_shiftX.Z);
                    println!("");
                    println!("Qb shift: {}", Qb_shiftX.X / Qb_shiftX.Z);
                    println!("");
                    println!("PQb shift: {}", PQb_shiftX.X / PQb_shiftX.Z);
                    println!("");
                    println!("");

                    println!("");
                    println!("Pa shift: {}", Pa_shiftX.X / Pa_shiftX.Z);
                    println!("");
                    println!("Qa shift: {}", Qa_shiftX.X / Qa_shiftX.Z);
                    println!("");

                    println!("");
                    println!("Pa1 shift: {}", Pa1_shiftX.X / Pa1_shiftX.Z);
                    println!("");
                    println!("Qa1 shift: {}", Qa1_shiftX.X / Qa1_shiftX.Z);
                    println!("");
                    */

                    // TODO:
                    let mut Pa = Pa_to_be_mapped.clone();
                    let mut Qa = Qa_to_be_mapped.clone();

                    let points = compute_isogeny(
                        &ell_product,
                        &mut Pa,
                        &mut Qa,
                        &Pa1_to_be_mapped,
                        &Qa1_to_be_mapped,
                        &mut Pa_shiftX,
                        &mut Qa_shiftX,
                        &Pa1_shiftX,
                        &Qa1_shiftX,
                        b_points.to_vec(),
                        b_shift_points.to_vec(),
                        self.a as usize,
                        &self.strategy_2,
                    );

                    Pb1_to_be_mapped = points[0];
                    Qb1_to_be_mapped = points[1];
                    PQb1_to_be_mapped = points[2];
                }

                // println!("6 (after a series of (2,2)-isogenies): {:?}", sixth_part.elapsed());

                let alpha: Integer = generate_random_range(0.big(), l_c.big().pow(power_b - 1)) * l_c + 
                    generate_random_range(1.big(), 4.big()); // TODO: check if this is 1 or 4
                let alpha = 3.big(); // TODO, dbg
                curve.xmul(&mut Ra1_to_be_mapped, alpha.clone());

                let alice_secret0: Integer = generate_random_range(0.big(), l_a.big().pow(self.a + 2) - 1) * l_a + 1;
                let alice_secret1: Integer = generate_random_range(0.big(), l_a.big().pow(self.a + 2) - 1) * l_a + 1;

                // let alice_secret0 = (2*8 + 1).big(); // TODO: remove
                // let alice_secret1 = (2*7 + 1).big(); // TODO: remove

                curve.xmul(&mut Pa1_to_be_mapped, alice_secret0.clone());
                curve.xmul(&mut Qa1_to_be_mapped, alice_secret1.clone());

                /*
                println!("");
                println!("================= ????????????????????????????????????????????????");
                println!("");
                println!("Pa: {}", Pa_to_be_mapped.X / Pa_to_be_mapped.Z);
                println!("");

                println!("");
                println!("Qa: {}", Qa_to_be_mapped.X / Qa_to_be_mapped.Z);
                println!("");

                println!("");
                println!("Ra: {}", Ra_to_be_mapped.X / Ra_to_be_mapped.Z);
                println!("");

                println!("");
                println!("Pa1: {}", Pa1_to_be_mapped.X / Pa1_to_be_mapped.Z);
                println!("");

                println!("");
                println!("Qa1: {}", Qa1_to_be_mapped.X / Qa1_to_be_mapped.Z);
                println!("");

                println!("");
                println!("Ra1: {}", Ra1_to_be_mapped.X / Ra1_to_be_mapped.Z);
                println!("");

                println!("");
                println!("Pb_to_be_mapped: {}", Pb_to_be_mapped.X / Pb_to_be_mapped.Z);
                println!("");
                println!("Qb_to_be_mapped: {}", Qb_to_be_mapped.X / Qb_to_be_mapped.Z);
                println!("");
                println!("PQb_to_be_mapped: {}", PQb_to_be_mapped.X / PQb_to_be_mapped.Z);
                println!("");
                println!("");

                println!("");
                println!("========================== 2");
                println!("");
                println!("Pb1_to_be_mapped: {}", Pb1_to_be_mapped.X / Pb1_to_be_mapped.Z);
                println!("");
                println!("Qb1_to_be_mapped: {}", Qb1_to_be_mapped.X / Qb1_to_be_mapped.Z);
                println!("");
                println!("PQb1_to_be_mapped: {}", PQb1_to_be_mapped.X / PQb1_to_be_mapped.Z);
                println!("");
                */

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

                println!("pub key gen: {:?}", start.elapsed());

                (PubKey::new(
                    pubkey_points,
                    pubkey_points1,
                    power_a,
                    power_b,
                    power_c,
                ), (alice_secret0, alice_secret1, alpha))
            }

            pub fn generate_pub_key_dbg(&self) -> (PubKey, (Integer, Integer, Integer)) {

                let alice_secret0 = "55229142801563651167479381070741390294961084967249063319026542761781798602700093540824213202749429562504163971101963".big();
                let alice_secret1 = "296450948127633855881835738357114816850041462722352125163225898973449003875834390871812467152296118395134915003819243".big();
                let alpha = "1293729513794169818769802401092416653942".big();

                let mut X = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("263205817672033389536356441280802898652727220579613789487815894903424181914063001893156500498772503885454443139458492158271724536184331666045060682346746060824489129668970477949013515429597156928623540886902557963548504452514751412487")),
                    &Fp::decode_reduce(&bytes_from_str("1023272960850707793383137101084827385481369622026301268103519253667025397122133784990154670496541674440950389851464961832708768430569644906552884235562824979161795786742659273721548000927607714162503407927633179999301178640361879327038")),
                );
                let mut Z = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("572760840870178236853035772046735088944263883127408227743826525910146830870243892404319381470578310887213765508505455764803543911460257820884737481918507424123366311907049472911556841845632127052921477557404812452178270608325702011940")),
                    &Fp::decode_reduce(&bytes_from_str("234345867037187051610981381853782202706269938724656640077342440740797061625542324778369821023242484701745009558785082869683701592813823066867533105119081137278641881696726315579929247899922677257215353607153534417322160252557927055233")),
                );
                let Pa = PointX::new_xz(&X, &Z);

                X = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("240918786091575766840576186192357636557015078016470269736026936356922741088277795649683470521697315577783735936565273524638906792104597920177059922370304671071069165983929680649648434418196077112314159249136556207071450945526368441149")),
                    &Fp::decode_reduce(&bytes_from_str("538758818757349133002334203565947133564004939954432410267934456931526434527270870572394346051109484405705073114438783073697890177010198898338740113384887168556408028664889513891185249249962912110638216013619542334119366503201864418733")),
                );
                Z = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("1169575524930456511052422167045051820533219023082547967478869772682929446973890214841250759690721375644632220463740407305261311346976254039050568243335561380749084931602635687363313527871957929176745004143519314209571601982048697827207")),
                    &Fp::decode_reduce(&bytes_from_str("727944826454932951659000224915620402151366490329132155921526074300715688648246851218112947714875919996660710951266268356934298462332674623099042953085590030991142814710169476886351977037116530564392087825964419671160811151516418469213")),
                );
                let Qa = PointX::new_xz(&X, &Z);

                X = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("630576255341568418553692778103708953645275489050927473231058696076528648538755844196380112616730821819261821527773792327791767071905619304943613427900834120055267326015897098803345592293009449565330009764507857051778282758838436639333")),
                    &Fp::decode_reduce(&bytes_from_str("572692387501495165300465359133964248117792629963659040648955833588319400118204626515455147192613580589187207123673964110491328184007863670332627463948915442362224214373018410627024895127126781285809052280662624420632435167626977852644")),
                );
                Z = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("934528895432759183312079575240233098662340883181107972853542929751586703877430016790714857244407281607243861214466030319930641703666747625936825599016684951763351160072764395055810776103159553895722450143459640236417360808001692957459")),
                    &Fp::decode_reduce(&bytes_from_str("618491400089016280801250331471920178076947824186130910748776712282054053440359332042937904377257622451008426679933711959918774799830670634707033415329567326481809414407507931288819636182938880887462631358725785966646073477845194039682")),
                );
                let Pb = PointX::new_xz(&X, &Z);

                X = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("1241275950023586439013689116378720748876965969434280143905689586246085504432277909423009001535387784670214376517578777750967387508267287650936175160216476269222393215134005854840169391742818530966813407995990845951824045641644258265786")),
                    &Fp::decode_reduce(&bytes_from_str("1101836067764112893318898760255536598578741745542430318994644100635930427955200733485073786384797528174409015959712520818821819698289035672289507474595304828662744920366107776182474243030259998406783668391579868576605176009091160031096")),
                );
                Z = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("953162847760700614135716983545079449425897260816138491959529667060467709516658230763862569194684010481639043169670049747811826278040554674066204632376704705215149151318368876050200999934556596712125576456419739317465295641148490934137")),
                    &Fp::decode_reduce(&bytes_from_str("1283033347877881417316143767892577806674590848005883233315310958552843533960822989882679810681231221728350241609338522183118408443209615851670071439546639821816264337702538408299274885643855696867140155703750864991523815944370762899009")),
                );
                let Qb = PointX::new_xz(&X, &Z);

                X = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("472190568570193600451590234909149697753581192812203454700338316760710319453310411308211569122405471655764997766054731473888365635332171479648861858668437219774559857079771011481280551294299207956012067277720203479986035325368032472305")),
                    &Fp::decode_reduce(&bytes_from_str("212061112509645264568673213364926843870473123732873106269266603681497847878389656673733785638397605437505778987653010506437470058541901794627628431806506519175105245469937776321800146359662640116720955847430410967181950567515562279839")),
                );
                Z = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("1")),
                    &Fp::decode_reduce(&bytes_from_str("0")),
                );
                let PQb = PointX::new_xz(&X, &Z);

                X = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("681428482972751808421502760860138445549965666772715944635334458876283758730420344287506834395719537819680848817430471019869197617224247310992765624211495197606745854332486107960264387662719250795206481240660491271320650397338291904436")),
                    &Fp::decode_reduce(&bytes_from_str("686094387702535908251456749006509182617285743241612042069667305282006980567549501210799783333206185640773363882085008903752483156323040775491625814221735887801761527557050706356131140542996810814689862139735414910462883416070655471686")),
                );
                Z = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("549651389049981303959273744712410683422653015170651381258820367955385851183772617385738676402638842058150802333629643131133957187446599491582350896058866954391491278487022246231662375579308308507753660996138228372559530395518956429285")),
                    &Fp::decode_reduce(&bytes_from_str("594140730898402856982652133719012660418251534111898519533630719674400510698367172396414157516773853024840387555496274857687137449638381919922589366188690130776103136594462998404422670511036928782721015379588486453330645908904240398922")),
                );
                let Ra = PointX::new_xz(&X, &Z);

                X = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("583588841805073909298466372660327303892403696414982972312207811704082867847048015850482309802983335137923802583190428597570980862283822525259305550734207338147169326075876347168825742597777521727350543606408560469005974053228596528607")),
                    &Fp::decode_reduce(&bytes_from_str("59501724868485986852391135161052339555036198838525453741356605899596325968935860873597517022337401329003146359663891250722079439083601289061771246794503141711137613004208933521180985786715921019976041823909386473020146975847816887898")),
                );
                Z = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("1064260066063637903632038545888553490373190774852512190790118085534904041197000635007607003693104346854662827329109485068323447866204607933352124730047448850972391318710177965335256896922016330269833546087863996826469737636969355042979")),
                    &Fp::decode_reduce(&bytes_from_str("222116090344401317110181883770693901148552564276652322159640551495897146994175886234041859651320675316627817593856825724248535290586894022666441318558180874998194976220638958619868371462734328996832622297713976169391628276959880410209")),
                );
                let Pa1 = PointX::new_xz(&X, &Z);

                X = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("922041350990977902439740910427696803148400418164295124545744668819752239031089098744215222380270106443145277776224500713258511067863894848953676530735227248934395531307374685298605434851664664320462733720845143076469976626363672107950")),
                    &Fp::decode_reduce(&bytes_from_str("248262486030769717064579836874084781220487792355811110029078454000101673962230019724040407463823582718227674105397960105337298052773473040575232240173928815799778441114528946939385537809321887534615402977255685060704114656994158597448")),
                );
                Z = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("163468198729609337231370565767491520504294422490437778191866472036617878735048555470647663736952577509177877643891832856966068791789134215654467356025900134156042119712540558674531110613981803800163682385831216781588909966985860903513")),
                    &Fp::decode_reduce(&bytes_from_str("518116975018042521450289126740664762787745416169499180094027202510216521869659544726959448235622017884651155801708191137221801023842556987084134039114865842191362059071169604876473936930814369251195123424550563647621482985645801631875")),
                );
                let Qa1 = PointX::new_xz(&X, &Z);

                X = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("1169083748960531960282404934575646485148455266063335549895518806858388901793551593649790066393409849778342036989218552123730956634434049555737843342794304582717521332072755905707921316896624592183838146701749384442047709145469918269790")),
                    &Fp::decode_reduce(&bytes_from_str("959841901154456052665532892397224562587854277010811980159860431776881392833655031640128495985949889037258879131845676490817451104498309569050795582109822825714846837838471581078495892734145986464926641491547696911812387343867722648807")),
                );
                Z = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("141685170962052613720147216264376326782542335866625992645156528395298119573279611177286335475439268431862468006038753260929138688264609908065623788655973690713064510803397271771231912360552936924524880981066125581159927596884479885101")),
                    &Fp::decode_reduce(&bytes_from_str("950814153677230974898026235749747378246314546443534216497148874962503947550663370501826072422488636077872540760648865467376145080326794633190645781672175873388135999248656844681210560124206853482628808754281894086060675296703886390597")),
                );
                let Pb1 = PointX::new_xz(&X, &Z);

                X = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("1173624221962798595949466741895584279779224164419763451880508729030965606737307666684799791202143838296227374236097982366763676152073815722560105719421376368963600425447800897816857390282821669116717532227114917570229526929457471746765")),
                    &Fp::decode_reduce(&bytes_from_str("517153237664566289867241438802905180221547586832656689292825030094630508194793212185087472166426833392708310585634776427585470368277275966165584302820085745125710790781192440520900229495096378810894911800740278510944232888300757867200")),
                );
                Z = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("1029686899059005750081173825632177218496285814712328727759024224417934758728391292009209529184442450666015176273259529484310014745314638564482500534358328325646407366886256163959089393699786245141257295845522832039810327771038568325289")),
                    &Fp::decode_reduce(&bytes_from_str("434113451513889410347034797769081582096730498149096998084922694254840237792556144880650433988815015426150695719978506417986507464332954990397577029175678348038356350755455976071014915321409005061776465690949665396279837439483286238273")),
                );
                let Qb1 = PointX::new_xz(&X, &Z);

                X = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("675320109422614410017922021726129357343124136205342976845783996302770944909091828681989447344403810300069424082113546366152139636622996270584118191286256242594059874801904079665166250701538661494969395410514105769604458695452369079519")),
                    &Fp::decode_reduce(&bytes_from_str("313533281589695265556267747022713956127803568635445746210401178770898577264185037721879280879218347606474040625843461557762775323228472866027149922562419776190053029373432816221452702776372916423426869071304010148981028922050367362472")),
                );
                Z = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("991122466467826101337721802767080882579477425053081010190720596860136591151136679311692178324433236948406852703169229165653927356065379100592513920164515930624878713972457688131515088376167363141384414291556098720648727188224008381586")),
                    &Fp::decode_reduce(&bytes_from_str("399222155426583204172588466293783960681693228148523151499503578920785635105341853707906523811665994802623111132271239569056033462598692772128651467842253476238788067724433397719954542970799008404009920386261268585761845901366392106462")),
                );
                let PQb1 = PointX::new_xz(&X, &Z);

                X = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("1155393821577617518943596857514683084480775790668574661963955500833681080373903566830980332168688044618664145908607316402498956486742665722257567403945062201286148134557922834527715372954159020426064047107676892312664278601639533984188")),
                    &Fp::decode_reduce(&bytes_from_str("306872904021287076032110954866913899621609445045796878619671159182876670896437948040258070763443794817051177965506230007950177005907847059996110302820311463553521902123380697616394591309887362400196199618787607520918478195073606169128")),
                );
                Z = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str("947090807060019973010964395239624334008534941533364175653262895407398981403619063847817703050595043857231569887261350449928151677298696825648565457457058757104603824145231959932932517482278389360475526318157637694464998170906404964844")),
                    &Fp::decode_reduce(&bytes_from_str("367038986449595650266236756196535866599824873320257048955935382976782023986384669393292472866822595730161962237393970393086997575682658486581875620281629440188771988274316127034708118924890467063159055741495729472450083178836686944282")),
                );
                let Ra1 = PointX::new_xz(&X, &Z);

                let pubkey_points = PubKeyPoints::new(
                    Pa,
                    Qa,
                    Pb,
                    Qb,
                    PQb,
                    Ra
                );

                let pubkey_points1 = PubKeyPoints::new(
                    Pa1,
                    Qa1,
                    Pb1,
                    Qb1,
                    PQb1,
                    Ra1
                );

                let power_a = 3 * 128;
                let power_b = 162;
                let power_c = 56;

                (PubKey::new(
                    pubkey_points,
                    pubkey_points1,
                    power_a,
                    power_b,
                    power_c,
                ), (alice_secret0, alice_secret1, alpha))

            }

            pub fn encrypt(&self, pub_key: &PubKey, mu: Integer) -> Cipher {
                let start = Instant::now();

                let (points, points1) = (&pub_key.points, &pub_key.points1);
                let (mut Pa, mut Qa, Pb, Qb, PQb, mut R) = (points.Pa, points.Qa, points.Pb, points.Qb, points.PQb, points.R);
                let (mut Pa1, mut Qa1, Pb1, Qb1, PQb1, mut R1) = (points1.Pa, points1.Qa, points1.Pb, points1.Qb, points1.PQb, points1.R);

                let A24 = get_montgomery_A24(&Pb, &Qb, &PQb);
                let A24_1 = get_montgomery_A24(&Pb1, &Qb1, &PQb1);

                /*
                println!("");
                println!("");
                println!("A: {}", A24.0 / A24.1);
                println!("");
                */

                let s = self.l_b * generate_random_range(0.big(), self.l_b.big().pow(pub_key.power_b) - 1) +
                            generate_random_range(1.big(), 2.big()); // TODO: check if this is 1 or 2
                let s = "11".big(); // dbg

                let curve = Curve::new_fromA24(&A24.0, &A24.1);
                let curve1 = Curve::new_fromA24(&A24_1.0, &A24_1.1);

                let kernelx = curve.ladder_3pt(&Pb, &Qb, &PQb, s.clone());
                let kernel1x = curve1.ladder_3pt(&Pb1, &Qb1, &PQb1, s.clone());

                /*
                println!("");
                println!("");
                println!("");
                println!("Pb_to_be_mapped: {}", Pb.X / Pb.Z);
                println!("");
                println!("Qb_to_be_mapped: {}", Qb.X / Qb.Z);
                println!("");
                println!("PQb_to_be_mapped: {}", PQb.X / PQb.Z);

                println!("");
                println!("");
                println!("");
                println!("kernelx: {}", kernelx.X / kernelx.Z);
                println!("");
                */
                
                let n = self.strategy.len() + 1;

                let eval_points = [Pa, Qa, R];
                let (curve_new, image_points) = three_isogeny_chain(&curve, &kernelx, eval_points.to_vec(), n, &self.strategy);
                (Pa, Qa, R) = (image_points[0], image_points[1], image_points[2]);
 
                let eval_points = [Pa1, Qa1, R1];
                let (curve1_new, image_points) = three_isogeny_chain(&curve1, &kernel1x, eval_points.to_vec(), n, &self.strategy);
                (Pa1, Qa1, R1) = (image_points[0], image_points[1], image_points[2]);
 
                let beta = self.l_c * generate_random_range(0.big(), self.l_c.big().pow(pub_key.power_c - 1) - 1) +
                            generate_random_range(1.big(), 4.big()); // TODO: check if this can be 4
                let beta = 928347.big(); // dbg

                curve_new.xmul(&mut R, beta.clone());
                curve1_new.xmul(&mut R1, beta * mu);

                let bob_secret_1: Integer = self.l_a * generate_random_range(0.big(), self.l_a.big().pow(self.a + 2) - 1) + 1;
                let bob_secret_2: Integer = self.l_a * generate_random_range(0.big(), self.l_a.big().pow(self.a + 2) - 1) + 1;

                let bob_secret_1: Integer = 2 * 1324.big() + 1; // dbg
                let bob_secret_2: Integer = 2 * 8345.big() + 1; // dbg

                curve_new.xmul(&mut Pa, bob_secret_1.clone());
                curve_new.xmul(&mut Qa, bob_secret_2.clone());

                curve1_new.xmul(&mut Pa1, bob_secret_1);
                curve1_new.xmul(&mut Qa1, bob_secret_2);
 
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

            pub fn decrypt(&self, cipher: &Cipher, alice_secret: (Integer, Integer, Integer)) -> Integer {
                let start = Instant::now();

                let (points, points1) = (&cipher.points, &cipher.points1);
                let (mut Pa, mut Qa, mut R) = (points.Pa, points.Qa, points.R);
                let (Pa1, Qa1, R1) = (points1.Pa, points1.Qa, points1.R);

                let curve = Curve::new_fromA24(&cipher.A24_num, &cipher.A24_denom);
                let curve1 = Curve::new_fromA24(&cipher.A24_1_num, &cipher.A24_1_denom);

                let alpha = alice_secret.2;
                
                curve.xmul(&mut R, alpha);
                curve.xmul(&mut Pa, alice_secret.0);
                curve.xmul(&mut Qa, alice_secret.1);
 
                let (mut Pa_complete, _) = curve.complete_pointX(&Pa);
                let (Qa_complete, _) = curve.complete_pointX(&Qa);
                let (R_complete, _) = curve.complete_pointX(&R);

                let (Pa1_complete, _) = curve1.complete_pointX(&Pa1);
                let (Qa1_complete, _) = curve1.complete_pointX(&Qa1);
                let (R1_complete, _) = curve1.complete_pointX(&R1);


                let t = self.l_a.big().pow(self.a);
                let bytes = big_to_bytes(t.clone());
                let mut T1_0 = curve.mul(&Pa_complete, &bytes, bytes.len() * 8);
                let T1_1 = curve1.mul(&Pa1_complete, &bytes, bytes.len() * 8);

                let Qa_complete_mul = curve.mul(&Qa_complete, &bytes, bytes.len() * 8);
                let Qa1_complete_mul = curve1.mul(&Qa1_complete, &bytes, bytes.len() * 8);

                let t = 2;
                let (mut w1, ok) = curve.weil_pairing_2exp(t, &T1_0, &Qa_complete_mul);
                assert_eq!(ok, 0xFFFFFFFF);
                let bytes = big_to_bytes(3.big());
                w1.set_pow_simple(&bytes);

                let (w2, ok) = curve1.weil_pairing_2exp(t, &T1_1, &Qa1_complete_mul);
                assert_eq!(ok, 0xFFFFFFFF);

                if w1.equals(&w2) != 0xFFFFFFFF {
                    Pa_complete.Y.set_neg();
                    T1_0.Y.set_neg();
                }

                let R_shift = curve.add(&R_complete, &T1_0);
                let R_shiftX = PointX::new_xz(&R_shift.X, &R_shift.Z);

                let Pa_shift = curve.add(&Pa_complete, &T1_0);
                let Qa_shift = curve.add(&Qa_complete, &T1_0);

                let Pa1_shift = curve1.add(&Pa1_complete, &T1_1);
                let Qa1_shift = curve1.add(&Qa1_complete, &T1_1);

                let mut Pa_shiftX = PointX::new_xz(&Pa_shift.X, &Pa_shift.Z);
                let mut Qa_shiftX = PointX::new_xz(&Qa_shift.X, &Qa_shift.Z);

                let Pa1_shiftX = PointX::new_xz(&Pa1_shift.X, &Pa1_shift.Z);
                let Qa1_shiftX = PointX::new_xz(&Qa1_shift.X, &Qa1_shift.Z);

                let ell_product = EllipticProduct::new(&curve, &curve1);

                let b_points = [R];
                let b_shift_points = [R_shiftX];
                let points_R = compute_isogeny(
                        &ell_product,
                        &mut Pa,
                        &mut Qa,
                        &Pa1,
                        &Qa1,
                        &mut Pa_shiftX,
                        &mut Qa_shiftX,
                        &Pa1_shiftX,
                        &Qa1_shiftX,
                        b_points.to_vec(),
                        b_shift_points.to_vec(),
                        self.a as usize,
                        &self.strategy_2,
                    );

                let R1_2 = points_R[0];

                let (R1_2_complete, _) = curve1.complete_pointX(&R1_2);

                // let mu1 = dlog_5(&curve1, &R1_complete, &R1_2_complete, self.c.try_into().unwrap());
                let mu2 = dlog_5(&curve1, &R1_2_complete, &R1_complete, self.c.try_into().unwrap());

                /*
                println!("");
                println!("mu1: {:?}", mu1);
                println!("");
                println!("mu2: {:?}", mu2);
                println!("");
                */

                println!("decrypt: {:?}", start.elapsed());
                println!("");

                mu2
            }

            fn get_PQb_and_shift(&self, curve_1: &Curve, curve_2: &Curve, Pa: &Point, Qa: &mut Point,
                    Pa1: &Point, Qa1: &Point, torsion_b_minus: Integer, tau: Integer, dbg_index: usize) -> (Point, Point, Point, Point, Point, Point, Point, Point, Point, Point) {
                let mut Pb1;
                let mut Qb1;
                loop {
                    Pb1 = generate_random_fq(curve_1, torsion_b_minus.clone(), self.scalar_without_b.clone());
                    Qb1 = generate_random_fq(curve_1, torsion_b_minus.clone(), self.scalar_without_b.clone());

                    // Pb1 and Qb1 needs to be linear independent
                    let mut Pb1_check = PointX::new_xz(&Pb1.X, &Pb1.Z);
                    let mut Qb1_check = PointX::new_xz(&Qb1.X, &Qb1.Z);
                    curve_1.xmul(&mut Pb1_check, torsion_b_minus.clone());
                    curve_1.xmul(&mut Qb1_check, torsion_b_minus.clone());
                
                    if (Pb1_check.X * Qb1_check.Z).equals(&(Qb1_check.X * Pb1_check.Z)) != 0xFFFFFFFF {
                        break;
                    }
                } 

                // TODO: the following is only for debugging:

                // 128
                /*
                if dbg_index == 1 {
                    let px = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("589347475779467028048732638707745169074099572448781051271349481559473727012677648788277674070708563815384409822118544148404313695739142076641733342195789664410029747217412725418639602788946160431051681659016463369120955499530350884793")),
                        &Fp::decode_reduce(&bytes_from_str("295541731903830973966435177821009510545806235315633443160350494708236538905366964828599192114390194343094738330899036091771606752115327511842840580779568519473668508004731498100416200716785989577694849244176022777068709128659699958651")),
                    );

                    let pz = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("89542995344710952948107022790100922686334678202826420161095874206526834236077000668814682919417119967989382371825802002025554791683974908178374722677745710441261906652561600152816296300222128057559398717792476238492791294646212107584")),
                        &Fp::decode_reduce(&bytes_from_str("1028199781928077642554616804903240208550974242321502258355144963022928866609527418284098186948125250284199941214438305372352745234972686884241582592073140746296427566777630409750957270938545946624116302693122101933302505969963727803020")),
                    );
                    let PX = PointX::new_xz(&px, &pz);
                    (Pb1, _) = curve_1.complete_pointX(&PX);
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
                    (Qb1, _) = curve_1.complete_pointX(&QX);
                } else {
                    let px = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("284707190851302075926377063138643761676122518864665550138554350732573227293123587051907006783021460026585455871404018076970517581567031238032632831488051206693751540697662652780496985234073415652017084557602710706990588136599420102096")),
                        &Fp::decode_reduce(&bytes_from_str("4433206878219121884491950114955940050211700723058222307183932529106090176445532488192447704785491805751562283162365888008148621761674266348704859253404297466367931215163833005890081063719658319791338544511987560854732703312240947031")),
                    );

                    let pz = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("1093935580269457206149374076770843225706669129380470411936650859912240538078291466643328407199225712337866700580108763359406951491933814233135993306482146305732905615647772041901933910489239185311569121548211842583140047412454765613305")),
                        &Fp::decode_reduce(&bytes_from_str("1058746561392362583015971962998183602791552932464219152435984593356777921207689856807492624119867947182507749577424350815864146522253512104890649805341781640256245725471234537497460313528057066910946709776469323303225975507097522965590")),
                    );
                    let PX = PointX::new_xz(&px, &pz);
                    (Pb1, _) = curve_1.complete_pointX(&PX);
                    Pb1.Y.set_neg();

                    let qx = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("145906234163643449030383730881120164351395064575978090914954268500926786583728869199373317834271986476877315118658884030103672435546508158706589528861130614410530635777559247958672748047469744348747979052090700786837774220905750451334")),
                        &Fp::decode_reduce(&bytes_from_str("766927802926939947613940812232960547095107913377442979431559751914776504643902023279120677143016734429835058593038324957421113750100735428986654560453300619399286631099476511286144418104949927538098263242063404210826681074530878105377")),
                    );

                    let qz = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("17419782317540346937134817761607359768723700988741504660760338340476733698028437247332684816574260923297851007969869787984299743170956056497979887538897018574853409590185803871726852482940019308515889258533338726222427150787964900510")),
                        &Fp::decode_reduce(&bytes_from_str("891775174187706552397851193269110663997129891837664509985080668295184536400915054403320330055172222617408490414104838764654009673467101334621063315809437588065515559074427782962273905781376370930984854891614238782019433308352768046476")),
                    );

                    let QX = PointX::new_xz(&qx, &qz);
                    (Qb1, _) = curve_1.complete_pointX(&QX);
                }
                */

                /*
                if dbg_index == 1 {
                    let px = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("100753668991277171490355201963314361660277645075961895756320938950168994514482427293306239033660286623006285771602255294614047768873955453249009536038635480600609949967051972938597826938849228962334903205675411066846922940384643092092279790200426417290168896154871688513316825098806810230568949110978399701497397053374116351799534310886558731859695035")),
                        &Fp::decode_reduce(&bytes_from_str("92491204986494587062469574918792029832463703933911054081899559112869325256948727640481348400100974046002831075111319174084532620604355339075994447057770040806318499818104273395199086316466838953093133105852775798782177762798648120044926202821175782302660121264222153205705057638314529152060144033573462773770910925560526192345144512046921897361174455")),
                    );

                    let pz = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("24293388504956253077372164835795399721275131521456330722678329733580297003553177728701242059110522001002304721045756579397687672974167077731302249702631807053980403681828774847871792685839412925953940802345039701493033860963032605104836853290964917839269280826448654916072976331099812006292113423111931728348898201882959720154531942096808570830554846")),
                        &Fp::decode_reduce(&bytes_from_str("43193131590983153142628149447805190883166757103637037696639364321071101922462575285905642903061374594278575678085918592555482501700615245354855074621447418818071075599318865996257076907441678039615200455622152190660127030200760550190717864307233279429346277791029494764062475399613110935273631525431029296322950024647969375850340816648561647267325392")),
                    );
                    let PX = PointX::new_xz(&px, &pz);
                    (Pb1, _) = curve_1.complete_pointX(&PX);
                    // Pb1.Y.set_neg();

                    let qx = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("92800230721345441162567807132495527094605382065602766535910967876544586971842691712839102714824637428384554254754379959538413107014907462352067566929527779019185895766916159870537902829454602047029080398134950882584172060983785451580629223865090477781733518921514513602856580694516207736109612523928995149942659823699518566090343902142792932806359465")),
                        &Fp::decode_reduce(&bytes_from_str("17009880293362217888967415176967777605183176713298609826127896833072607034879406147555001951852622528855735101964893980079033771419931609675044074335255273634915213459776214947323530145315947765462918519568657464371409772689310797636928648506816499789238556591504127524784820075592657373174982784164226625062065860948910286098448181191730299744190572")),
                    );

                    let qz = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("67527619995723224728180604860988697475649633708037080082738323284235243594291683279084838163690353530491691421861704215504817711013715203968001095988641405790375897731286671338350633505408156887589812421863703757392526310199192819906051384545352358833970415961177037132038079246954736725532927296233309565083222857423759306842392386847884133310092827")),
                        &Fp::decode_reduce(&bytes_from_str("86157855271729660947311291427197838280989162217790018062054085227023835453186725992316786949148834193398279640580166117562123008802177447999933455172356327152343065666920033385043752744167092531197813461713116612836553217038357119946255007657328533780811581961708907868104889464415728150227200741976759222206481103363658542196037282203837877819517848")),
                    );

                    let QX = PointX::new_xz(&qx, &qz);
                    (Qb1, _) = curve_1.complete_pointX(&QX);
                    Qb1.Y.set_neg();
                } else {
                    let px = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("24123553726066587528282347391548844062398117349933605659150597882727115516644243785505737428080808079667896639775705254588405790832816990191780911594089407961117509498380673094842227466657610011310518180779497154085964531669703428981675171292181446413266424732304876912824483626937559661422274592653989312806472911262349932599817909950492252277419015")),
                        &Fp::decode_reduce(&bytes_from_str("86089270101478411532221744496882941055019147186485544032277556318986638540962961918928227824371721612938549288638310625282669014523991869763689639560965746712499817027091358084104324108310271865005399436786060817456995508628258164979516140461179174733998524173619183401994208437486785593265903471326134035167376968051813145337369344744781934453695240")),
                    );

                    let pz = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("7196527391637381021388815554108228805758427596218851809565786435006338144490652073568476690964235668974442272072706803458872944679364133302308925793592416048596473160942418533293416459394631927144737739314799548107677415164177473002984773049580306363885191542349898892244157374410426505004461991930354023868200240074227039839956975303050332194776378")),
                        &Fp::decode_reduce(&bytes_from_str("16876927515028497435210207082370409126907926471185051273726979891342039985456957544117071862615549559903128249880488248458024804850953653525682165041243622797143190803170275087931473602598195000144022898864990318827497017412135013221921082640427118060013393299669645231860362000610464727962449805772877286105805329339531630410497845128435767703832701")),
                    );
                    let PX = PointX::new_xz(&px, &pz);
                    (Pb1, _) = curve_1.complete_pointX(&PX);
                    // Pb1.Y.set_neg();

                    let qx = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("37079638735308167033772840183994506662677277205316204672463114378001211678933698119167329380356555301087704616996424541588148556389921385010934631173734669018707700682515435019999625692260558365936921128803907453783617552039960442428600564808718349921281176254878469006985309261650572994059409931043022721070701103293237134973294991374144086168748624")),
                        &Fp::decode_reduce(&bytes_from_str("54631547615660019645571509565663315774545005800320220680031247207053072129911562659602491911868546783931953782191565251149128611349805751051176509978919838495909622035604418583787196801349686848489316348365720281080234072077063838221139862690250796668206139826423061644509177295259570124641154342823730015365273641838639828717963785066417030247038695")),
                    );

                    let qz = Fq::new(
                        &Fp::decode_reduce(&bytes_from_str("59930351690211151079234231568884221288245475852798745202866104616940536070453537337773747732444050809023123050669952213414385119044810707523940976807488507703414366991904726761576974292407541361930693571980824267350037955069201714486060883800138167780523345795680878579169967037427405316015350719730952130846241209406029664093359491932504535179954893")),
                        &Fp::decode_reduce(&bytes_from_str("89733948583869801337456878283347821094749299555161383591850945170114318252579048065055284932598764652271867292532038398909072107626049786111894669979176320354892877250493343917471201283904584393336038857398851048670539241162672046372965129484645114342441514618179136524399886622659736422821305810441901001597577964501937323085631104758176578844077564")),
                    );

                    let QX = PointX::new_xz(&qx, &qz);
                    (Qb1, _) = curve_1.complete_pointX(&QX);
                }
                */

                // end of debugging
                /*
                println!("");
                println!("");
                println!("foo");
                println!("");
                println!("{}", Pb1.X / Pb1.Z);
                println!("");
                println!("{}", Pb1.Y / Pb1.Z);
                println!("");
                println!("{}", Qb1.X / Qb1.Z);
                println!("");
                println!("{}", Qb1.Y / Qb1.Z);
                println!("");
                println!("");
                */

                let PQb1 = curve_1.sub(&Pb1, &Qb1);
                // let PQb1 = curve_1.add(&Pb1, &Qb1); // TODO

                let t = (self.a+2).try_into().unwrap();
                let (mut w1, _) = curve_1.weil_pairing_2exp(t, Pa, Qa);
                let bytes = big_to_bytes(tau);
                w1.set_pow_simple(&bytes);
                let (w2, _) = curve_2.weil_pairing_2exp(t, Pa1, Qa1);
                if w1.equals(&w2) == 0xFFFFFFFF {
                    
                } else {
                    Qa.Y.set_neg();
                }

                // TODO: this is again computed using X coordinates in compute_isogeny, perhaps
                // use this result
                let shift_1 = curve_1.double_iter(Pa, self.a as usize);
                let shift_2 = curve_2.double_iter(Pa1, self.a as usize);
            
                /*
                println!("============ goo ==============");
                println!("");
                println!("Pa: {}", Pa.X / Pa.Z);
                println!("");
                println!("Pa: {}", Pa.Y / Pa.Z);
                println!("");

                println!("");
                println!("Pa1: {}", Pa1.X / Pa1.Z);
                println!("");
                println!("Pa1: {}", Pa1.Y / Pa1.Z);
                println!("");

                println!("");
                println!("shift 1: {}", shift_1.X / shift_1.Z);
                println!("");
                println!("shift 1: {}", shift_1.Y / shift_1.Z);
                println!("");

                println!("");
                println!("shift 2: {}", shift_2.X / shift_2.Z);
                println!("");
                println!("shift 2: {}", shift_2.Y / shift_2.Z);
                println!("");
                */

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
                    represent_integer(N.clone(), qa.clone(), order.clone(), bad_prod_primes.clone(), 1.big()).unwrap();
                
                /*
                println!("gamma:");
                println!("");
                println!("gamma[0] = {}", gamma.x);
                println!("");
                println!("gamma[1] = {}", gamma.y);
                println!("");
                println!("gamma[2] = {}", gamma.z);
                println!("");
                println!("gamma[3] = {}", gamma.t);
                println!("");
                println!("denom:");
                println!("{}", gamma.denom);

                let norm = gamma.reduced_norm(); 
                println!("");
                println!("???????");
                println!("");
                println!("norm:");
                println!("{:?}", norm.numer() / norm.denom());
                println!("");
                println!("N");
                println!("{}", N);
                */

                // assert!(norm.numer() == &(N * norm.denom()));

                gamma
            }
        }
    };
} // End of macro: define_litsigamal

pub(crate) use define_litsigamal;
