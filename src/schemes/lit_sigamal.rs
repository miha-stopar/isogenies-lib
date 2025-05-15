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
                    strategy  = &[65, 37, 23, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 9, 5, 4, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1];
                    strategy_2 = &[154, 93, 55, 33, 20, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 22, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 38, 22, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 65, 37, 21, 12, 7, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1];
                } else if param == 192 {
                    strategy = &[92, 65, 37, 21, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 37, 21, 13, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1];
                    strategy_2 = &[229, 146, 86, 49, 28, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 37, 21, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 60, 37, 21, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 23, 16, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 9, 5, 4, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 93, 55, 33, 20, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 22, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 38, 22, 13, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1];
                } else {
                    strategy = &[111, 89, 55, 34, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 34, 22, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1];
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
                let (Pa, _) = curve.complete_pointX(&PaX);
                let (Qa, _) = curve.complete_pointX(&QaX);
                let (PmQa, _) = curve.complete_pointX(&PmQaX);

                let (Pb, _) = curve.complete_pointX(&PbX);
                let (Qb, _) = curve.complete_pointX(&QbX);
                let (PmQb, _) = curve.complete_pointX(&PmQbX);

                let (Pc, _) = curve.complete_pointX(&PcX);
                let (Qc, _) = curve.complete_pointX(&QcX);
                let (PmQc, _) = curve.complete_pointX(&PmQcX);

                let l_a = 2;
                let l_b = 3;
                let l_c = 5;
 
                let mut a1 = generate_random_range(0.big(), l_a.big().pow(power_a) - 1);
                let mut a2 = generate_random_range(0.big(), l_a.big().pow(power_a - 1)) * l_a + 1;
                a1 = "67821213147103272121758535926151039005030912112029727951085328323311503610545609587044506281792717705600807906711
952".big();
                a2 = "10827571144677997620644424104650364372653000782849569633450213074630849473914074810111054720448115636003540948886
5063".big();

                let mut c1: Integer;
                let mut c2: Integer;

                let five = Integer::from(5);
                let zero = Integer::ZERO;
                loop {
                    c1 = generate_random_range(0.big(), l_c.big().pow(power_c) - 1);
                    c2 = generate_random_range(0.big(), l_c.big().pow(power_c) - 1);
                    if c1.clone().modulo(&five) != zero && c2.clone().modulo(&five) != zero {
                        break;
                    }
                }
                // Pc = c1 * Pc + c2 * Qc
                let mut bytes = big_to_bytes(c1.clone());
                let c1Pc = self.curve.mul(&Pc, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(c2.clone());
                let c2Qc = self.curve.mul(&Qc, &bytes, bytes.len() * 8);
                let Pc_rand = self.curve.add(&c1Pc, &c2Qc);
 
                let Rx = PointX::new_xz(&Pc_rand.X, &Pc_rand.Z);

                // N = (l_a**(2*a) - n**2) * l_b**b
                let tau = l_a.big().pow(self.a * 2) - self.n * self.n;
                let N = tau.clone() * l_b.big().pow(power_b);
                let gamma = self.generate_gamma(N.clone());

                let order = standard_maximal_extremal_order().order;
                let (coord, imprim) = gamma.factor_in_order(order.lattice.clone());

                // println!("1 (before apply endomorphism): {:?}", start.elapsed());
                // let second_part = Instant::now();

                let torsion_a = l_a.big().pow(self.a + 2);
                
                let (mut Pa_gamma, Qa_gamma, _) = apply_endomorphism_on_torsion_group(&curve, coord.clone(), imprim.clone(), torsion_a, mat2_2, mat2_3, mat2_4, &Pa, &Qa, &PmQa);

                let bytes1 = big_to_bytes(a1);
                let a1Pa_gamma = curve.mul(&Pa_gamma, &bytes1, bytes1.len() * 8);
                let bytes2 = big_to_bytes(a2);
                let a2Qa_gamma = curve.mul(&Qa_gamma, &bytes2, bytes2.len() * 8);
                let mut Qa_rand_gamma = curve.add(&a1Pa_gamma, &a2Qa_gamma);

                let torsion_c = l_c.big().pow(power_c);
                let (Pc_gamma, Qc_gamma, _) = apply_endomorphism_on_torsion_group(&curve, coord.clone(), imprim.clone(), torsion_c, mat5_2, mat5_3, mat5_4, &Pc, &Qc, &PmQc);

                let mut bytes = big_to_bytes(c1);
                let c1Pc = self.curve.mul(&Pc_gamma, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(c2);
                let c2Qc = self.curve.mul(&Qc_gamma, &bytes, bytes.len() * 8);
                let R = self.curve.add(&c1Pc, &c2Qc);
                let mut R_gamma = PointX::new_xz(&R.X, &R.Z);

                // Set [2^a]Pa_gamma = (1,*)
                let t = self.l_a.big().pow(self.a); // TODO: define once
                let bytes = big_to_bytes(t);
                let Pa_gamma_4 = curve.mul(&Pa_gamma, &bytes, bytes.len() * 8);
                if Pa_gamma_4.X.equals(&Pa_gamma_4.Z) != 0xFFFFFFFF {
                    Pa_gamma.X.set_neg();
                    Qa_rand_gamma.X.set_neg();
                    R_gamma.X.set_neg();
                    
                    let mut tmp = PointX::new_xz(&Pa_gamma.X, &Pa_gamma.Z);
                    (Pa_gamma, _) = curve.complete_pointX(&tmp);

                    tmp = PointX::new_xz(&Qa_rand_gamma.X, &Qa_rand_gamma.Z);
                    (Qa_rand_gamma, _) = curve.complete_pointX(&tmp);
                }
                let mut Pa1_to_be_mapped = PointX::new_xz(&Pa_gamma.X, &Pa_gamma.Z);
                let mut Ra1_to_be_mapped = PointX::new_xz(&R_gamma.X, &R_gamma.Z);

                let a1Pa = curve.mul(&Pa, &bytes1, bytes1.len() * 8);
                let a2Qa = curve.mul(&Qa, &bytes2, bytes2.len() * 8);
                let Qa_rand = curve.add(&a1Pa, &a2Qa);
                let Qa_rand_X = PointX::new_xz(&Qa_rand.X, &Qa_rand.Z);

                let mut Qa1_to_be_mapped = PointX::new_xz(&Qa_rand_gamma.X, &Qa_rand_gamma.Z);
            
                let torsion_b = l_b.big().pow(power_b);
                let torsion_b_minus = l_b.big().pow(power_b - 1);
                let (Pb_gamma, Qb_gamma, _) = apply_endomorphism_on_torsion_group(&curve, coord.clone(), imprim.clone(), torsion_b.clone(), mat3_2, mat3_3, mat3_4, &Pb, &Qb, &PmQb);

                // println!("2 (after apply endomorphism): {:?}", second_part.elapsed());
                // let third_part = Instant::now();

                let dlog1 = dlog_3(&curve, &Pb_gamma, &Qb_gamma, self.b.try_into().unwrap());

                // println!("3 (after dlog): {:?}", third_part.elapsed());
                // let fourth_part = Instant::now();

                let bytes = big_to_bytes(dlog1);
                let dlog_Pb = curve.mul(&Pb, &bytes, bytes.len() * 8);
                let kernel1 = curve.sub(&Qb, &dlog_Pb);

                let kernel1x = PointX::new_xz(&kernel1.X, &kernel1.Z);

                let n = self.strategy.len() + 1;

                // endomorphism gamma is of order (l_a**(2*a) - n**2) * l_b**b
                // let's take the isogeny gamma1 with kernel: ker(gamma) \cap E[l_b**b]

                let eval_points = [PaX, Qa_rand_X, Rx]; 
                let (mut codomain, mut image_points) = three_isogeny_chain(&curve, &kernel1x, eval_points.to_vec(), n, self.strategy);
                let (mut Pa_isog3X, mut Qa_rand_isog3X, R_isog3X) = (image_points[0], image_points[1], image_points[2]);

                let mut Pa_to_be_mapped = Pa_isog3X.clone(); // to be mapped by 3-isogeny chain
                let mut Qa_to_be_mapped = Qa_rand_isog3X.clone();
                let mut Ra_to_be_mapped = R_isog3X.clone();

                let (Pa, _) = codomain.complete_pointX(&Pa_isog3X);
                let (mut Qa_rand, _) = codomain.complete_pointX(&Qa_rand_isog3X);

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

                // start of debugging

                let Qb11 = generate_random_fq(&codomain, torsion_b_minus.clone(), self.scalar_without_b.clone());

                let mut bytes1 = big_to_bytes(Integer::from(self.n));
                let torsion8_Qb = self.curve.mul(&Qb11, &bytes1, bytes1.len() * 8);

                let factor = l_a.big().pow(self.a - 1);
                bytes1 = big_to_bytes(factor);
                let torsion8_Pa = self.curve.mul(&Pa_gamma, &bytes1, bytes1.len() * 8);
                let torsion8_Qa = self.curve.mul(&Qa_gamma, &bytes1, bytes1.len() * 8);

                let P1P2 = CouplePoint::new(&torsion8_Pa, &torsion8_Qa);
                let Q1Q2 = CouplePoint::new(&torsion8_Qa, &torsion8_Qb);


                // TODO:
                let image_points1 = vec![
                    P1P2,
                    P1P2,
                    P1P2
                    // CouplePoint::new(&self.two_dim.P, &self.two_dim.Q),
                    // CouplePoint::new(&self.two_dim.omegaP, &self.two_dim.omegaQ),
                ];

                // Precomputed with strategy.py
                // TODO: derive 383 from self.a
                // let strategy: [usize; 383] = [144, 89, 55, 34, 27, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 8, 6, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 55, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1];

                let chain2 = Instant::now();

                let (product1, points1) = product_isogeny(
                    &ell_product,
                    &P1P2,
                    &Q1Q2,
                    &image_points1,
                    // self.a as usize,
                    self.a as usize - 1,
                    &self.strategy_2,
                );

                println!("chain2: {:?}", chain2.elapsed());

                // end of debugging

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

                    let eval_points = [Pa_to_be_mapped, Qa_to_be_mapped, Ra_to_be_mapped, Qb_to_be_mapped];
                    (codomain, image_points) = three_isogeny_chain(&codomain, &kernel, eval_points.to_vec(), n, &self.strategy);
                    (Pa_to_be_mapped, Qa_to_be_mapped, Ra_to_be_mapped, backtracking_check) = (image_points[0], image_points[1], image_points[2], image_points[3]);

                    codomain.xmul(&mut backtracking_check, f_mul.clone());

                    let eval_points = [Pa1_to_be_mapped, Qa1_to_be_mapped, Ra1_to_be_mapped]; 
                    (curve, image_points) = three_isogeny_chain(&curve, &kernel1, eval_points.to_vec(), n, &self.strategy);
                    (Pa1_to_be_mapped, Qa1_to_be_mapped, Ra1_to_be_mapped) = (image_points[0], image_points[1], image_points[2]);

                    let (Pa, _) = codomain.complete_pointX(&Pa_to_be_mapped);
                    let (mut Qa, _) = codomain.complete_pointX(&Qa_to_be_mapped);
                    let (Pa1, _) = curve.complete_pointX(&Pa1_to_be_mapped);
                    let (Qa1, _) = curve.complete_pointX(&Qa1_to_be_mapped);

                    let (Pa_shift, Qa_shift, Pa1_shift, Qa1_shift, Pb, Qb, PQb, Pb_shift, Qb_shift, PQb_shift) =
                        self.get_PQb_and_shift(&codomain, &curve, &Pa, &mut Qa, &Pa1, &Qa1, torsion_b_minus.clone(), tau.clone(), 2);

                    Pb_to_be_mapped = PointX::new_xz(&Pb.X, &Pb.Z);
                    Qb_to_be_mapped = PointX::new_xz(&Qb.X, &Qb.Z);
                    PQb_to_be_mapped = PointX::new_xz(&PQb.X, &PQb.Z);

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

                let alpha: Integer = generate_random_range(0.big(), l_c.big().pow(power_c - 2)) * l_c + 
                    generate_random_range(1.big(), 5.big());
                curve.xmul(&mut Ra1_to_be_mapped, alpha.clone());

                let alice_secret0: Integer = generate_random_range(0.big(), l_a.big().pow(self.a + 2) - 2) * l_a + 1;
                let alice_secret1: Integer = generate_random_range(0.big(), l_a.big().pow(self.a + 2) - 2) * l_a + 1;

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

                println!("pub key gen: {:?}", start.elapsed());

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

                let s = self.l_b * generate_random_range(0.big(), self.l_b.big().pow(pub_key.power_b) - 2) +
                            generate_random_range(1.big(), 3.big());

                let curve = Curve::new_fromA24(&A24.0, &A24.1);
                let curve1 = Curve::new_fromA24(&A24_1.0, &A24_1.1);

                let kernelx = curve.ladder_3pt(&Pb, &Qb, &PQb, s.clone());
                let kernel1x = curve1.ladder_3pt(&Pb1, &Qb1, &PQb1, s.clone());

                let n = self.strategy.len() + 1;

                let eval_points = [Pa, Qa, R];
                let (curve_new, image_points) = three_isogeny_chain(&curve, &kernelx, eval_points.to_vec(), n, &self.strategy);
                (Pa, Qa, R) = (image_points[0], image_points[1], image_points[2]);
 
                let eval_points = [Pa1, Qa1, R1];
                let (curve1_new, image_points) = three_isogeny_chain(&curve1, &kernel1x, eval_points.to_vec(), n, &self.strategy);
                (Pa1, Qa1, R1) = (image_points[0], image_points[1], image_points[2]);
 
                let beta = self.l_c * generate_random_range(0.big(), self.l_c.big().pow(pub_key.power_c - 2) - 1) +
                            generate_random_range(1.big(), 5.big());

                curve_new.xmul(&mut R, beta.clone());
                curve1_new.xmul(&mut R1, beta * mu);

                let bob_secret_1: Integer = self.l_a * generate_random_range(0.big(), self.l_a.big().pow(self.a + 2) - 1) + 1;
                let bob_secret_2: Integer = self.l_a * generate_random_range(0.big(), self.l_a.big().pow(self.a + 2) - 1) + 1;

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
                let mu2 = dlog_5(&curve1, &R1_2_complete, &R1_complete, self.c.try_into().unwrap());

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

                let PQb1 = curve_1.sub(&Pb1, &Qb1);

                let t = (self.a+2).try_into().unwrap();
                let (mut w1, _) = curve_1.weil_pairing_2exp(t, Pa, Qa);
                let bytes = big_to_bytes(tau);
                w1.set_pow_simple(&bytes);
                let (w2, _) = curve_2.weil_pairing_2exp(t, Pa1, Qa1);
                if w1.equals(&w2) != 0xFFFFFFFF {
                    Qa.Y.set_neg();
                }

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
                    represent_integer(N.clone(), qa.clone(), order.clone(), bad_prod_primes.clone(), 1.big()).unwrap();
                
                // assert!(norm.numer() == &(N * norm.denom()));

                gamma
            }
        }
    };
} // End of macro: define_litsigamal

pub(crate) use define_litsigamal;
