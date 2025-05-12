use isogenies_lib::{ec_lit128::{self, Curve, LITSiGamal}, util::Big};
use criterion::{criterion_group, criterion_main, Criterion};

pub fn criterion_benchmark(cr: &mut Criterion) {
    let lit_sigamal = ec_lit128::LITSiGamal::new(128);

    cr.bench_function("generate_pub_key", |b| {
        b.iter(|| {
            lit_sigamal.generate_pub_key();
        });
    });
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10); // Adjust the sample size to limit iterations
    targets = criterion_benchmark
}
criterion_main!(benches);
