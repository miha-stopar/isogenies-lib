use isogenies_lib::{ec_lit::{self, Curve, LITSiGamal}, util::Big};
use criterion::{criterion_group, criterion_main, Criterion};

pub fn criterion_benchmark(cr: &mut Criterion) {
    let p = "1290217975993796939363993419446162388979006021159541007293712082644700121088673466685157498316158528176855539315411759315356741765308895915108991692829754882889263058278152142847999999999999999999999999999999999999999999999999999999999".big();
    let (a, b, c, f) = ec_lit::get_params(128);
    let A = ec_lit::Fq::ZERO;
    let curve = ec_lit::Curve::new(&A);
    let lit_sigamal = ec_lit::LITSiGamal::new(curve, p, 2, 3, 5, a, b, c, f, 3);

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
