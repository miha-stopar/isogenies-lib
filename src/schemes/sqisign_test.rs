
#[cfg(test)]
mod tests {
    use crate::{ec_lit, util::{big_to_bytes, Big}};

    #[test]
    fn dlog_3_test() {
        let curve = ec_lit::Curve::new(&ec_lit::Fq::ZERO);
        let fileName = "src/schemes/precomputed.json";
        let (P, _, _, _, _, _, power) = ec_lit::load_torsion_info(fileName, "lit-sigamal-128", 3);
        let (Pa, ok1) = curve.complete_pointX(&P);
        assert!(ok1 == 0xFFFFFFFF);

        let f: usize = power as usize;
        let test_data = vec![1.big(), 7.big(), 18.big(), 112123123.big(), "150094635296999120".big()];
        for i in test_data {
            let k_bytes = big_to_bytes(i.clone());
            let R = curve.mul(&Pa, &k_bytes, k_bytes.len() * 8);
            let dlog = ec_lit::dlog_3(&curve, &Pa, &R, f);
            assert!(dlog == i);
        }
    }

    #[test]
    fn dlog_5_test() {
        let curve = ec_lit::Curve::new(&ec_lit::Fq::ZERO);
        let fileName = "src/schemes/precomputed.json";
        let (_, P, _, _, _, _, power) = ec_lit::load_torsion_info(fileName, "lit-sigamal-128", 5);
        let (Pa, ok1) = curve.complete_pointX(&P);
        assert!(ok1 == 0xFFFFFFFF);

        let f: usize = power as usize;
        let test_data = vec![1.big(), 7.big(), 18.big(), 112123123.big(), "150094635296999120".big(), "1387778780781445675529539585113525390621".big()];
        for i in test_data {
            let k_bytes = big_to_bytes(i.clone());
            let R = curve.mul(&Pa, &k_bytes, k_bytes.len() * 8);
            let dlog = ec_lit::dlog_5(&curve, &Pa, &R, f);
            assert!(dlog == i);
        }
    }
    
}
