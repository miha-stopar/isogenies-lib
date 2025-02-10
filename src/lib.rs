#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

pub mod ec;
mod error;
pub mod fields;
pub mod finitefield;
pub mod linalg;
pub mod quaternion;
pub mod schemes;
mod theta;
pub mod util;

pub mod ec117 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp117::Fp;
    pub type Fq = crate::fields::Fp117Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec214 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp214::Fp;
    pub type Fq = crate::fields::Fp214Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec509 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp509::Fp;
    pub type Fq = crate::fields::Fp509Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec1757 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp1757::Fp;
    pub type Fq = crate::fields::Fp1757Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec5248 {
    pub type Fp = crate::fields::Fp5248::Fp;
    pub type Fq = crate::fields::Fp5248Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::ec::ec_helpers::define_ec_helpers! {}
    crate::theta::theta::define_theta_structure! {}
}

pub mod ec_lit {
    pub type Fp = crate::fields::FpLit128::Fp;
    pub type Fq = crate::fields::FpLit128Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::ec::ec_helpers::define_ec_helpers! {}
    
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::lit_sigamal::define_litsigamal! {}
}

pub mod ec_sqisign_lvl1 {
    pub type Fp = crate::fields::SQIsignLvl1::Fp;
    pub type Fq = crate::fields::SQIsignLvl1Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::ec::ec_helpers::define_ec_helpers! {}
}